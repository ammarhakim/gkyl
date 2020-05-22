-- Gkyl ------------------------------------------------------------------------
--
-- Interpolation of a DG field defined on one grid, onto another field defined on
-- a Cartesian grid with a different resolution.
-- 
-- Current limitations:
--    1) Only does whole grids with same domain (no grid subsets).
--    2) It may only work for (2^a)*(3^b) grids.
--
-- Notes:
--    a] It's possible that in the future on wishes to pass a grid, or the nodes of
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
local lume              = require "Lib.lume" 
local Range             = require "Lib.Range"
local PrimeFactor       = require "Lib.PrimeFactor"

-- Interpolate updater object.
local CartFieldInterpolate = Proto(UpdaterBase)

local function IndexStencilMapRefine(dir, idx, nCells, beta, stencil)
   -- Given a index (idx) to a cell in the coarse grid
   -- with nCells cells, return the index of the refinement stencil needed,
   -- within the table that holds stencils (self.stencils).
   local remDecL = (idx-1)*(beta)-math.floor((idx-1)*(beta))
   local remDecU = math.ceil(idx*beta)-idx*(beta)
   if ((idx == 1) or   -- First cell.
       (remDecL == 0) or    -- Interior cell with a left-boundary-like stencil.
       ((remDecL <= 0.5) and (remDecU <= 0.5))) then -- or
      stencil = 2*stencil + 3^(dir-1) - 1
   elseif ((idx == nCells) or   -- Last cell.
           (remDecU == 0)) then             -- Interior cell with a right-boundary-like stencil.
      stencil = 2*stencil + 3^(dir-1)
   end
   return stencil
end

local function IndexStencilMapCoarsen(dir, idx, nCells, beta, stencil)
   -- Given a index (idx) to a cell in the fine grid
   -- with nCells cells, return the index of the coarsening stencil needed,
   -- within the table that holds stencils (self.stencils).

   local remDecL = (idx-1)*beta-math.floor(idx*beta)
   local remDecU = math.ceil(idx*beta)-idx*beta
   if ((idx == 1) or   -- First cell.
       (remDecL == 0) or    -- Interior cell with a left-boundary-like stencil.
       ((remDecL > 0) and (remDecU > 0))) then -- or
      stencil = 2*stencil + 3^(dir-1) - 1
   elseif ((idx == nCells) or   -- Last cell.
           (remDecU == 0)) then             -- Interior cell with a right-boundary-like stencil.
      stencil = 2*stencil + 3^(dir-1)
   end

   return stencil
end

function CartFieldInterpolate:init(tbl)
   CartFieldInterpolate.super.init(self, tbl) -- Setup base object.

   local outGrid  = assert(
      tbl.onGrid, "Updater.CartFieldInterpolate: Must provide grid to interpolate to using 'onGrid'.")

   local inGrid   = assert(
      tbl.fromGrid, "Updater.CartFieldInterpolate: Must provide grid to interpolate from using 'fromGrid'.")

   local outBasis = assert(
      tbl.onBasis, "Updater.CartFieldInterpolate: Must provide basis to interpolate to using 'onBasis'.")

   local inBasis  = assert(
      tbl.fromBasis, "Updater.CartFieldInterpolate: Must provide basis to interpolate from using 'fromBasis'.")

   assert(outBasis:ndim() == inBasis:ndim(), "Updater.CartFieldInterpolate: the dimensionalities must be the same.")
   -- For now we'll restrict ourselves to equal bases as well.
   assert(outBasis:polyOrder() == inBasis:polyOrder(), "Updater.CartFieldInterpolate: polynomial orders must be the same.")
   assert(outBasis:id() == inBasis:id(), "Updater.CartFieldInterpolate: basis types must be the same.")

   self.dim        = inBasis:ndim()
   local polyOrder = inBasis:polyOrder()
   local basisID   = inBasis:id()

   self.inNumCells  = {}
   self.outNumCells = {}
   for d = 1, self.dim do
      self.inNumCells[d]  = inGrid:numCells(d)
      self.outNumCells[d] = outGrid:numCells(d)

      -- For now limit to (2^a)*(3^b) grids. Check:
      local pf = PrimeFactor.all(self.inNumCells[d])
      for i = 1, #pf do
         assert((pf[i]==2) or (pf[i]==3), string.format("Updater.CartFieldInterpolate: Only (2^a)*(3^b) grids supported. Prime factor: %d", pf[i]))
      end
      local pf = PrimeFactor.all(self.outNumCells[d])
      for i = 1, #pf do
         assert((pf[i]==2) or (pf[i]==3), string.format("Updater.CartFieldInterpolate: Only (2^a)*(3^b) grids supported. Prime factor: %d", pf[i]))
      end
   end

   self.beta            = {}   -- Ratio of cell-lengths in each direction.
   self.indexStencilMap = {}
   for d = 1, self.dim do
      self.beta[d] = inGrid:dx(d)/outGrid:dx(d)
      if self.beta[d] > 1 then
         self.indexStencilMap[d] = IndexStencilMapRefine
      else
         self.indexStencilMap[d] = IndexStencilMapCoarsen
      end
   end

   self.intStencilSize = {}   -- Interior (away from boundaries) stencil size, in each direction.
   for d = 1, self.dim do
      if (self.beta[d] > 1) then   -- Mesh refinement.
         -- Brute force search:
         -- Start with the size of the boundary stencil:
         local maxSize = math.floor(self.beta[d])+math.ceil(self.beta[d]-math.floor(self.beta[d]))
         local maxSizePossible = math.ceil(self.beta[d])+1
         for i = 2, self.inNumCells[d]-1 do
            local decimalL = 1-((i-1)*(self.beta[d])-math.floor((i-1)*self.beta[d]))
            local decimalU = 1-(math.ceil(i*self.beta[d])-i*self.beta[d])
            local currSize = math.floor(self.beta[d]-decimalL-decimalU)+math.ceil(decimalL)+math.ceil(decimalU) 
            maxSize = math.max(maxSize, currSize)
         end
         self.intStencilSize[d] = maxSize
      else                         -- Mesh coarsening.
         self.intStencilSize[d] = 1+math.ceil(1*(1/self.beta[d]-math.floor(1/self.beta[d])))
      end
   end

   -- This updater will loop through the input-grid, and for each cell grid it will pass
   -- the output-grid cells to the kernel, one at a time, and the kernel will add
   -- the corresponding contributions to the output-grid cells.
   -- We will make a list of the cells that need to be passed in each region.
   self.stencils    = {}
   self.stencils[1] = {}
   self.stencils[1].stencilSize = self.intStencilSize
   -- For each cell in the stencil, the multidimensional index offset relative to the
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
   self._interpolation = CartFldInterpDecl.selectInterpolation(basisID, self.dim, polyOrder)

   -- Cell lengths and cell centers.
   self.inDx  = Lin.Vec(self.dim)
   self.inXc  = Lin.Vec(self.dim)
   self.outDx = Lin.Vec(self.dim)
   self.outXc = Lin.Vec(self.dim)

   self.outIdxLL = Lin.IntVec(self.dim)   -- Output-grid index of the "lower-left" corner a input-grid cell contributes to.
   self.outIdx   = Lin.IntVec(self.dim)   -- Output-grid index.

end

function CartFieldInterpolate:indexToStencil(idxIn, numCellsIn)
   -- Given a multi-dimensional index (idxIn) to a cell in a grid
   -- with numCellsIn cells, return the index of the stencil needed,
   -- within the table that holds relaxation/residue stencils.
   local stencilIdx = 1
   for d = 1, self.dim do
      stencilIdx = self.indexStencilMap[d](d,idxIn[d],numCellsIn[d],self.beta[d],stencilIdx)
   end
   return stencilIdx
end

-- Advance method.
function CartFieldInterpolate:_advance(tCurr, fldIn, fldOut)

   local inFld   = fldIn[1]
   local outFld  = fldOut[1]

   outFld:clear(0.0)

   local inGrid  = inFld:grid() 
   local outGrid = outFld:grid() 

   localRangeDecomp = LinearDecomp.LinearDecompRange {
      range = inFld:localRange(), numSplit = inGrid:numSharedProcs() }
   local tId        = inGrid:subGridSharedId()    -- Local thread ID.

   local inFldIndexer = inFld:genIndexer()
   local inFldItr     = inFld:get(1)

   local outFldIndexer = outFld:genIndexer()
   local outFldItr     = outFld:get(1)

   for inIdx in localRangeDecomp:rowMajorIter(tId) do

      inGrid:setIndex(inIdx)
      inGrid:getDx(self.inDx)
      inGrid:cellCenter(self.inXc)

      -- Compute the output-grid index of the "lower-left" corner this cell contributes to.
      for d = 1,self.dim do
         local eveOI    = self.beta[d]*(inIdx[d]-1)
         self.outIdxLL[d] = math.ceil(eveOI)+(math.ceil(eveOI-math.floor(eveOI))+1) % 2
      end

      inFld:fill(inFldIndexer(inIdx), inFldItr)   -- Pointer to input field.

      -- Loop over the output-grid cells this input-grid cell contributes to.
      for outIdxOff in self.stencils[self:indexToStencil(inIdx,self.inNumCells)].stencilRange:rowMajorIter() do
         for d=1,self.dim do self.outIdx[d] = self.outIdxLL[d]+outIdxOff[d] end

         outGrid:setIndex(self.outIdx)
         outGrid:getDx(self.outDx)
         outGrid:cellCenter(self.outXc)

         outFld:fill(outFldIndexer(self.outIdx), outFldItr)   -- Fine-grid field pointer.

         self._interpolation(self.inXc:data(), self.outXc:data(), self.inDx:data(), self.outDx:data(), inFldItr:data(), outFldItr:data())
      end
  
   end

end

return CartFieldInterpolate
