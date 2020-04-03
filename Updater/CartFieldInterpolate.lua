-- Gkyl ------------------------------------------------------------------------
--
-- Interpolation of a DG field onto a finer or coarser grid.
-- 
-- Current limitations:
--    1) Only does whole grids with same domain (no grid subsets).
--    2) It may only work for (2^a)*(3^b) grids.
--    3) Does not always work in both directions, e.g. 48 -> 9 breaks.
--
-- Notes:
--    a] The code below makes references to a coarse and a fine grid. However,
--       the updater works in either direction.
--    b] It's possible that in the future on wishes to pass a grid, or the nodes of
--       some arbitrary grid to interpolate to, instead of passing another field
--       to interpolate to.
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
local PrimeFactor       = require "Lib.PrimeFactor"

-- Interpolate updater object.
local CartFieldInterpolate = Proto(UpdaterBase)

function CartFieldInterpolate:init(tbl)
   CartFieldInterpolate.super.init(self, tbl) -- Setup base object.

   local toGrid    = assert(
      tbl.onGrid, "Updater.CartFieldInterpolate: Must provide grid to interpolate to using 'onGrid'.")

   local fromGrid  = assert(
      tbl.fromGrid, "Updater.CartFieldInterpolate: Must provide grid to interpolate from using 'fromGrid'.")

   local toBasis   = assert(
      tbl.onBasis, "Updater.CartFieldInterpolate: Must provide basis to interpolate to using 'onBasis'.")

   local fromBasis = assert(
      tbl.fromBasis, "Updater.CartFieldInterpolate: Must provide basis to interpolate from using 'fromBasis'.")

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

      -- For now limit to (2^a)*(3^b) grids. Check:
      local pf = PrimeFactor.all(self.numCellsC[d])
      for i = 1, #pf do
         assert((pf[i]==2) or (pf[i]==3), string.format("Updater.CartFieldInterpolate: Only (2^a)*(3^b) grids supported. Prime factor: %d", pf[i]))
      end
      local pf = PrimeFactor.all(self.numCellsF[d])
      for i = 1, #pf do
         assert((pf[i]==2) or (pf[i]==3), string.format("Updater.CartFieldInterpolate: Only (2^a)*(3^b) grids supported. Prime factor: %d", pf[i]))
      end
   end

   self.beta = {}   -- Ratio of cell-lengths in each direction.
   for d = 1, self.dim do self.beta[d] = fromGrid:dx(d)/toGrid:dx(d) end

   self.intStencilSize = {}   -- Interior (away from boundaries) stencil size, in each direction.
   for d = 1, self.dim do
      if (self.beta[d] > 1) then   -- Mesh refinement.
--         self.intStencilSize[d] = math.ceil(self.beta[d]) + math.ceil(self.beta[d] - math.floor(self.beta[d]))
--         self.intStencilSize[d] = math.floor(self.beta[d]) + math.ceil(2*(self.beta[d] - math.floor(self.beta[d])))
--         self.intStencilSize[d] = math.floor(self.beta[d]) + math.ceil(self.beta[d] - math.floor(self.beta[d]))
--         if (math.ceil(self.beta[d] - math.floor(self.beta[d])) == 1) then
--            self.intStencilSize[d] = self.intStencilSize[d]+1
--         end
         -- Brute force search:
         -- Start with the size of the boundary stencil:
         local maxSizePossible = math.ceil(self.beta[d])+1
         local maxSize = math.floor(self.beta[d])+math.ceil(self.beta[d]-math.floor(self.beta[d]))
         for i = 2, self.numCellsC[d]-1 do
            local decimalL = 1-((i-1)*(self.beta[d])-math.floor((i-1)*(self.beta[d])))
            local decimalU = 1-(math.ceil(i*self.beta[d])-i*(self.beta[d]))
            local currSize = math.floor(self.beta[d]-decimalL-decimalU)+math.ceil(decimalL)+math.ceil(decimalU) 
            maxSize = math.max(maxSize, currSize)
         end
         self.intStencilSize[d] = maxSize
      else                         -- Mesh coarsening.
         self.intStencilSize[d] = 1+math.ceil(1*(1/self.beta[d]-math.floor(1/self.beta[d])))
      end
   end

   -- This updater will loop through the coarse grid, and for each coarse grid it will pass
   -- the fine-grid cells to the kernel, one at a time, and the kernel will add
   -- the corresponding contributions to the fine (coarse) grid cells during prolongation (restriction).
   -- We will make a list of the cells that need to be passed in each region.
   self.stencils    = {}
   self.stencils[1] = {}
   self.stencils[1].stencilSize = self.intStencilSize
   print(" stencil size = ", self.stencils[1].stencilSize[1])
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
            self.stencils[sI].stencilSize[dI] = math.floor(self.beta[dI]) + math.ceil(self.beta[dI] - math.floor(self.beta[dI]))
            print(" stencil size = ", self.stencils[sI].stencilSize[dI])
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
      local remDecL = (idxIn[d]-1)*(self.beta[d])-math.floor((idxIn[d]-1)*(self.beta[d]))
      local remDecU = math.ceil(idxIn[d]*self.beta[d])-idxIn[d]*(self.beta[d])
      if ((idxIn[d] == 1) or   -- First cell.
          (remDecL == 0) or   -- Interior cell with a left-boundary-like stencil.
          ((remDecL <= 0.5) and (remDecU <= 0.5))) then -- or
--          (((math.ceil(idxIn[d]*(self.beta[d]))-idxIn[d]*(self.beta[d]))+((idxIn[d]-1)*(self.beta[d])-math.floor((idxIn[d]-1)*(self.beta[d])))) <= 0.5) or
--          ((math.ceil((idxIn[d]-1)*self.beta[d]-math.floor(idxIn[d]*self.beta[d]))) == 1)) then  -- Interior cell with a left-boundary-like stencil.
         stencilIdx = 2*stencilIdx + 3^(d-1) - 1
      elseif ((idxIn[d] == numCellsIn[d]) or   -- Last cell.
              (remDecU == 0)) then   -- Interior cell with a right-boundary-like stencil.
--              ((idxIn[d]*self.beta[d]-math.floor(idxIn[d]*self.beta[d])) == 0)) then -- or
--              ((math.ceil((idxIn[d]-1)*self.beta[d]-math.floor(idxIn[d]*self.beta[d]))) == 1)) then  -- Interior cell with a right-boundary-like stencil.
         stencilIdx = 2*stencilIdx + 3^(d-1)
      end
   end
   return stencilIdx
end

-- Advance method.
function CartFieldInterpolate:_advance(tCurr, inFld, outFld)

   local cFld        = inFld[1]
   local fFld        = outFld[1]

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

      print(" cIdx = ",cIdx[1])

      -- Loop over the fine-grid cells this coarse grid contributes to.
      for fIdxOff in self.stencils[self:idx2stencil(cIdx,self.numCellsC)].stencilRange:rowMajorIter() do
         for d=1,self.dim do self.fIdx[d] = self.fIdxLL[d]+fIdxOff[d] end
         print("                 | fIdx = ",self.fIdx[1])

         fGrid:setIndex(self.fIdx)
         fGrid:getDx(self.dxF)
         fGrid:cellCenter(self.xcF)

         fFld:fill(fFldIndexer(self.fIdx), fFldItr)   -- Fine-grid field pointer.

         self._prolongation(self.xcC:data(), self.xcF:data(), self.dxC:data(), self.dxF:data(), cFldItr:data(), fFldItr:data())
      end
  
   end

end

return CartFieldInterpolate
