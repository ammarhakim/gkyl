-- Gkyl ------------------------------------------------------------------------
--
-- Interpolation of a DG field, producing another DG field defined on a finer
-- grid.
--
-- 1) It's possible that in the future on wishes to pass a grid, or the nodes of
--    some arbitrary grid to interpolate to, instead of passing another field
--    to interpolate to.
-- 2) For now we'll restrict ourselves to the case in which the finer grid has
--    cell lengths that are smaller by an integer number.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local CartFldInterpDecl = require "Updater.interpolateCalcData.CartFieldInterpolateModDecl"
local Proto             = require "Lib.Proto"
local UpdaterBase       = require "Updater.Base"
local LinearDecomp      = require "Lib.LinearDecomp"
local Lin               = require "Lib.Linalg"
local ffi               = require "ffi"
local lume              = require "Lib.lume" 
local Range             = require "Lib.Range"

-- Interpolate updater object.
local CartFieldInterpolate = Proto(UpdaterBase)

function CartFieldInterpolate:init(tbl)
   CartFieldInterpolate.super.init(self, tbl) -- Setup base object.

   local toGrid    = assert(
      tbl.onGrid, "Updater.CartFieldInterpolate: Must provide grid we'll interpolate to using 'targetGrid'.")

   local fromGrid  = assert(
      tbl.fromGrid, "Updater.CartFieldInterpolate: Must provide grid we'll interpolate from using 'originGrid'.")

   local toBasis   = assert(
      tbl.onBasis, "Updater.CartFieldInterpolate: Must provide basis we'll interpolate to using 'targetBasis'.")

   local fromBasis = assert(
      tbl.fromBasis, "Updater.CartFieldInterpolate: Must provide basis we'll interpolate from using 'originBasis'.")

   -- Ensure dimensionality of two quantities is the same.
   assert(toBasis:ndim() == fromBasis:ndim(), "Updater.CartFieldInterpolate: the dimensionalities must be the same.")
   -- For now we'll restrict ourselves to equal bases as well.
   assert(toBasis:polyOrder() == fromBasis:polyOrder(), "Updater.CartFieldInterpolate: polynomial orders must be the same.")
   assert(toBasis:id() == fromBasis:id(), "Updater.CartFieldInterpolate: basis types must be the same.")

   self.dim        = fromBasis:ndim()        -- Dimension of space.
   local polyOrder = fromBasis:polyOrder()   -- Polynomial order.
   local basisID   = fromBasis:id()          -- Basis kind.

   -- Number of cells in coarse and fine grids (assuming fromGrid is the coarse-grid).
   self.numCellsC = {}
   self.numCellsF = {}
   for d = 1, self.dim do
      self.numCellsC[d] = fromGrid:numCells(d)
      self.numCellsF[d] = toGrid:numCells(d)
   end

   self.beta = {}   -- Ratio of cell-lengths in each direction.
   for d = 1, self.dim do self.beta[d] = fromGrid:dx(d)/toGrid:dx(d) end

   self.intStencilSize = {}   -- Largest interior (away from boundaries) stencil size, in each direction.
   for d = 1, self.dim do
      self.intStencilSize[d] = math.floor(self.beta[d]) + math.ceil(2*(self.beta[d] - math.floor(self.beta[d])))
   end

   -- This updater will loop through the coarse grid, and for each coarse grid it will pass
   -- the fine-grid cells to the kernel, one at a time, and the kernel will add
   -- the corresponding contributions to the fine (coarse) grid cells during prolongation (restriction).
   -- We will make a list of the cells that need to be passed in each region.
   self.stencils    = {}
   self.stencils[1] = {}
   self.stencils[1].stencilSize = self.intStencilSize
   -- For each cell in the stencil, the multidimensional index offset relative to th
   -- lower left corner index (of the stencil) will be given by the index of a range object.
   local rangeLimits = {{},{}}
   for d = 1, self.dim do
      rangeLimits[1][d] = 0
      rangeLimits[2][d] = self.stencils[1].stencilSize[d]-1
   end
   self.stencils[1].stencilRange = Range.Range(rangeLimits[1], rangeLimits[2])
   sI = 1
   for dI = 1, self.dim do
      for s = 1, #self.stencils do
         for mp = -1,1,2 do
            sI = sI+1
            self.stencils[sI] = {}
            self.stencils[sI].stencilSize     = lume.clone(self.stencils[s].stencilSize)
            self.stencils[sI].stencilSize[dI] = math.floor(self.beta[dI]) + math.ceil(self.beta[dI] - math.floor(self.beta[dI])) -- self.stencils[sI].stencilSize[dI]-math.ceil(self.beta[dI] - math.floor(self.beta[dI]))
--            self.stencils[sI].stencilSizeN = self.stencils[s].stencilSizeN*((self.intStencilSize[dI]-(math.ceil(self.beta[dI] - math.floor(self.beta[dI]))))/self.intStencilSize[dI])
            rangeLimits = {{},{}}
            for d = 1, self.dim do
               rangeLimits[1][d] = 0
               rangeLimits[2][d] = self.stencils[sI].stencilSize[d]-1
            end
            self.stencils[sI].stencilRange = Range.Range(rangeLimits[1], rangeLimits[2])
         end
      end
   end

   -- Select interpolation kernels.
   self._prolongation = CartFldInterpDecl.selectProlongation(basisID, self.dim, polyOrder)

   -- Cell lengths and cell centers in coarse and fine grids.
   self.dxC = Lin.Vec(self.dim)
   self.xcC = Lin.Vec(self.dim)
   self.dxF = Lin.Vec(self.dim)
   self.xcF = Lin.Vec(self.dim)

   self.fIdxLL = Lin.IntVec(self.dim)   -- Fine-grid index of the "lower-left" corner a coarse-grid cell contributes to.
   self.fIdx   = Lin.IntVec(self.dim)   -- Fine-grid index.

end

function CartFieldInterpolate:idx2stencil(idxIn, numCellsIn)
   -- Given a multi-dimensional index (idxIn) to a cell in a grid
   -- with numCellsIn cells, return the index of the stencil needed,
   -- within the table that holds relaxation/residue stencils.
   local stencilIdx = 1
   for d = 1, self.dim do
      if ((idxIn[d] == 1) or   -- First cell.
          (((idxIn[d]-1)*self.beta[d]-math.floor((idxIn[d]-1)*self.beta[d])) == 0)) then  -- Interior cell with a left-boundary-like stencil.
         stencilIdx = 2*stencilIdx + 3^(d-1) - 1
      elseif ((idxIn[d] == numCellsIn[d]) or   -- Last cell.
              ((idxIn[d]*self.beta[d]-math.floor(idxIn[d]*self.beta[d])) == 0)) then  -- Interior cell with a right-boundary-like stencil.
         stencilIdx = 2*stencilIdx + 3^(d-1)
      end
   end
   return stencilIdx
end

function CartFieldInterpolate:prolong(cFld,fFld)
   -- Prolongation of a coarse-grid field (cFld) to a fine-grid field (fFld). 

   local cGrid       = cFld:grid() 
   local fGrid       = fFld:grid() 

   localRangeDecomp  = LinearDecomp.LinearDecompRange {
      range = cFld:localRange(), numSplit = cGrid:numSharedProcs() }
   local tId         = cGrid:subGridSharedId()    -- Local thread ID.

   local cFldIndexer = cFld:genIndexer()
   local cFldItr     = cFld:get(1)

   local fFldIndexer = fFld:genIndexer()
   local fFldItr     = fFld:get(1)

   for cIdx in localRangeDecomp:rowMajorIter(tId) do

      cGrid:setIndex(cIdx)
      cGrid:getDx(self.dxC)
      cGrid:cellCenter(self.xcC)

      -- Compute the fine-grid index of the "lower-left" corner this coarse cell contributes to.
      for d = 1,self.dim do
         local eveOI    = self.beta[d]*(cIdx[d]-1)
         self.fIdxLL[d] = math.ceil(eveOI)+(math.ceil(eveOI-math.floor(eveOI))+1) % 2
      end

      cFld:fill(cFldIndexer(cIdx), cFldItr)   -- Coarse-grid field pointer.

      -- Loop over the fine-grid cells this coarse grid contributes to.
      for fIdxOff in self.stencils[self:idx2stencil(cIdx,self.numCellsC)].stencilRange:rowMajorIter() do
         for d=1,self.dim do self.fIdx[d] = self.fIdxLL[d]+fIdxOff[d] end

         fGrid:setIndex(self.fIdx)
         fGrid:getDx(self.dxF)
         fGrid:cellCenter(self.xcF)

         fFld:fill(fFldIndexer(self.fIdx), fFldItr)   -- Fine-grid field pointer.

         self._prolongation(self.xcC:data(), self.xcF:data(), self.dxC:data(), self.dxF:data(), cFldItr:data(), fFldItr:data())
      end
  
   end
end

-- Advance method.
function CartFieldInterpolate:_advance(tCurr, inFld, outFld)

   local coarseField = inFld[1]
   local fineField   = outFld[1]

   self:prolong(coarseField, fineField)

end

return CartFieldInterpolate
