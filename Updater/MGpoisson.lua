-- Gkyl ------------------------------------------------------------------------
--
-- Multigrid solver for the Poisson equation
--     - L(phi) = rho
-- where L stands for Laplacian.
--
--
-- Question/development notes:
--   1) Do we need to create a grid object for each level? or create pseudo-grid object that is more lightweight and has the information needed?
--   2) Prolongation: so data doesn't get loaded multiple times, instead of iterating through each fine-grid cell, could we instead make
--                    the looping stencil-aware?
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local MGpoissonDecl    = require "Updater.mgPoissonCalcData.MGpoissonModDecl"
local Proto            = require "Lib.Proto"
local UpdaterBase      = require "Updater.Base"
local Grid             = require "Grid"
local DataStruct       = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local LinearDecomp     = require "Lib.LinearDecomp"
local Lin              = require "Lib.Linalg"

-- Multigrid updater object.
local MGpoisson = Proto(UpdaterBase)

function MGpoisson:init(tbl)
   MGpoisson.super.init(self, tbl) -- Setup base object.

   local topGrid = assert(
      tbl.onGrid, "Updater.MGpoisson: Must provide grid object using 'grid'.")

   local basis = assert(
      tbl.basis, "Updater.MGpoisson: Must provide the weak basis object using 'basis'.")

   self._dim       = basis:ndim()        -- Dimension of space.
   local polyOrder = basis:polyOrder()   -- Polynomial order.
   local basisID   = basis:id()          -- Basis kind.

   self._mgLevels  = 2   -- Number of multigrid levels (grids).

   -- Create a grid for each level.
   -- Not sure this is needed, but in general it probably is (e.g. unstructured, or even nonuniform meshes).
   self._mgGrids    = {}
   self._mgGrids[1] = topGrid
   periodicDirCount = 0
   self._rho        = {}   -- Right-side source field at each level.
   for i = 2, self._mgLevels do
      -- Determine the parameters of the coarse grid.
      local lowerC        = {}
      local upperC        = {}
      local cellsC        = {}
      local periodicDirsC = {}
      local decompCutsC   = {}
      for d = 1, self._dim do
         lowerC[d] = topGrid:lower(d)
         upperC[d] = topGrid:upper(d)
         if topGrid:isDirPeriodic(d) then
            periodicDirCount = periodicDirCount+1
            periodicDirsC[periodicDirCount] = d
         end

         -- The following two are more complicated. To be refined later.
         cellsC[d]      = (self._mgGrids[i-1]:numCells(d))/2
         decompCutsC[d] = self._mgGrids[i-1]:cuts(d) 
      end

      local isSharedC = topGrid:isShared()

      local decompC = DecompRegionCalc.CartProd {
         cuts      = decompCutsC,
         useShared = useSharedC,
      }

      self._mgGrids[i] = Grid.RectCart {
         lower         = lowerC, 
         upper         = upperC,
         cells         = cellsC,
         periodicDirs  = periodicDirsC,
         decomposition = decompC,
      }

      -- Allocate space for right-side source field, one for each level.
      self._rho[i] = DataStruct.Field {
         onGrid        = self._mgGrids[i],
         numComponents = basis:numBasis(),   -- NOTE: this will change if we do p-coarsening.
         ghost         = {1, 1},
      }
   end

   -- Select restriction and prolongation operator kernels.
   self._restrict = MGpoissonDecl.selectRestriction(basisID, self._dim, polyOrder)
   self._prolong  = MGpoissonDecl.selectProlongation(basisID, self._dim, polyOrder)

   -- Used to store fine and coarse grid indexes of neighboring cells needed.
   self.fIdx = {}
   self.cIdx = {}
   for i = 1, 3^(self._dim) do   -- The size 3^d is a worst case scenario for 2nd order BVPs.
      self.fIdx[i] = Lin.IntVec(self._dim)
      self.cIdx[i] = Lin.IntVec(self._dim)
   end

end

-- Restriction of a fine-grid field (fFld) to a coarse-grid field (cFld). 
function MGpoisson:restrict(fFld,cFld)

   local grid = cFld:grid() 

   localRangeDecomp = LinearDecomp.LinearDecompRange {
      range = cFld:localRange(), numSplit = grid:numSharedProcs() }

   local cFldIndexer = cFld:genIndexer()
   local cFldItr     = cFld:get(1)

   local fFldIndexer = fFld:genIndexer()
   -- Store pointers to all the neighbors needed. Ideally we would put
   -- all the data this points to into a single array or struct that
   -- we can pass to the kernel.
   local stencilSize = 2
   local fFldItr     = {}
   for i = 1, stencilSize do
      fFldItr[i] = fFld:get(1)
   end

   for cIdx in localRangeDecomp:rowMajorIter(tId) do

      grid:setIndex(cIdx)

      --if (cIdx[1] % 2 == 1) then
      -- Given the coarse-grid index, need to obtain the fine-grid index.
      for d = 1, self._dim do self.fIdx[1][d] = 2*(cIdx[d]-1)+1 end
      for d = 1, self._dim do self.fIdx[2][d] = self.fIdx[1][d]+1 end
  
  
      cFld:fill(cFldIndexer(cIdx), cFldItr)
  
      for i = 1, stencilSize do
         fFld:fill(fFldIndexer(self.fIdx[i]), fFldItr[i])
      end
  
      self._restrict(fFldItr[1]:data(), fFldItr[2]:data(), cFldItr:data())
      --end
   end
end

-- Function performing a single gamma-cycle.
--   gamma=1 V-cycle
--   gamma=2 W-cycle
-- What about F cycle? 
function MGpoisson:gammaCycle()
  return 0
end

-- Advance method.
function MGpoisson:_advance(tCurr, inFld, outFld)

   local srcFld  = inFld[1]

   -- Restrict the right-side source field.
   self._rho[1] = srcFld
   for i = 2, self._mgLevels do
      self:restrict(self._rho[i-1], self._rho[i])
      self._rho[i]:write(string.format("rho_%d.bp", i), tCurr, 0)
   end



   outFld[1]:copy(srcFld)   -- Temporary arbitrary output.

end

return MGpoisson
