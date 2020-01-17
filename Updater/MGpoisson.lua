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
local ffi              = require "ffi"

-- Multigrid updater object.
local MGpoisson = Proto(UpdaterBase)

function MGpoisson:init(tbl)
   MGpoisson.super.init(self, tbl) -- Setup base object.

   local topGrid = assert(
      tbl.onGrid, "Updater.MGpoisson: Must provide grid object using 'grid'.")

   local basis = assert(
      tbl.basis, "Updater.MGpoisson: Must provide the weak basis object using 'basis'.")

   self.dim        = basis:ndim()        -- Dimension of space.
   local polyOrder = basis:polyOrder()   -- Polynomial order.
   local basisID   = basis:id()          -- Basis kind.

   self.mgLevels   = 2   -- Number of multigrid levels (grids).

   -- Create a grid for each level.
   -- Not sure this is needed, but in general it probably is (e.g. unstructured, or even nonuniform meshes).
   self.mgGrids     = {}
   self.mgGrids[1]  = topGrid
   periodicDirCount = 0
   self.rho         = {}   -- Right-side source field at each level.
   for i = 2, self.mgLevels do
      -- Determine the parameters of the coarse grid.
      local lowerC        = {}
      local upperC        = {}
      local cellsC        = {}
      local periodicDirsC = {}
      local decompCutsC   = {}
      for d = 1, self.dim do
         lowerC[d] = topGrid:lower(d)
         upperC[d] = topGrid:upper(d)
         if topGrid:isDirPeriodic(d) then
            periodicDirCount = periodicDirCount+1
            periodicDirsC[periodicDirCount] = d
         end

         -- The following two are more complicated. To be refined later.
         cellsC[d]      = (self.mgGrids[i-1]:numCells(d))/2
         decompCutsC[d] = self.mgGrids[i-1]:cuts(d) 
      end

      local isSharedC = topGrid:isShared()

      local decompC = DecompRegionCalc.CartProd {
         cuts      = decompCutsC,
         useShared = useSharedC,
      }

      self.mgGrids[i] = Grid.RectCart {
         lower         = lowerC, 
         upper         = upperC,
         cells         = cellsC,
         periodicDirs  = periodicDirsC,
         decomposition = decompC,
      }

      -- Allocate space for right-side source field, one for each level.
      self.rho[i] = DataStruct.Field {
         onGrid        = self.mgGrids[i],
         numComponents = basis:numBasis(),   -- NOTE: this will change if we do p-coarsening.
         ghost         = {1, 1},
      }
   end

   -- Select restriction and prolongation operator kernels.
   self._restriction  = MGpoissonDecl.selectRestriction(basisID, self.dim, polyOrder)
   self._prolongation = MGpoissonDecl.selectProlongation(basisID, self.dim, polyOrder)
   -- Select relaction kernel.
   self._relaxation   = MGpoissonDecl.selectRelaxation(basisID, self.dim, polyOrder)

   -- Intergrid operator stencils: 
   self.igOpStencilWidth = 2
   self.igOpStencilSize  = self.igOpStencilWidth^self.dim
   -- Array of fine grid indexes of cells needed in inter-grid transfer.
   self.fineGridIdx = {}
   for i = 1, self.igOpStencilSize do
      self.fineGridIdx[i] = Lin.IntVec(self.dim)
   end
   -- List of pointers to the data in cells pointed to by the stencil.
   local DoublePtrVec = Lin.new_vec_ct(ffi.typeof("double*"))
   self.fineFldItr   = DoublePtrVec(self.igOpStencilSize)

   -- For now we'll only need one coarse-grid index.
   self.coarseGridIdx = Lin.IntVec(self.dim)

   -- Relxation stencil info: restrict ourselves to 'cross' relaxation stencils for now.
   self.relaxStencilWidth = 3
   self.relaxStencilSize  = (self.relaxStencilWidth-1)*self.dim+1
   self.relaxIdx = {}   -- List of cell indices pointed to by the stencil.
   for i = 1, self.relaxStencilSize do
      self.relaxIdx[i] = Lin.IntVec(self.dim)
   end
   -- List of pointers to the data in cells pointed to by the stencil.
   self.relaxItr      = DoublePtrVec(self.relaxStencilSize)
   self.relaxDxs      = DoublePtrVec(self.relaxStencilSize)


   -- Dimensions remaining when a dimension is removed.
   self.dimRemain = {}
   for d1 = 1, self.dim do
      self.dimRemain[d1] = {}
      for d2 = 1, self.dim do self.dimRemain[d1][d2] = d2 end
      table.remove(self.dimRemain[d1],d1)
   end

end

function MGpoisson:restrict(fFld,cFld)
   -- Restriction of a fine-grid field (fFld) to a coarse-grid field (cFld). 

   local grid = cFld:grid() 

   localRangeDecomp  = LinearDecomp.LinearDecompRange {
      range = cFld:localRange(), numSplit = grid:numSharedProcs() }
   local tId         = grid:subGridSharedId()    -- Local thread ID.

   local cFldIndexer = cFld:genIndexer()
   local cFldItr     = cFld:get(1)

   local fFldIndexer = fFld:genIndexer()
   local fFldItr     = fFld:get(1)

   for cIdx in localRangeDecomp:rowMajorIter(tId) do

      grid:setIndex(cIdx)

      -- Given the coarse-grid index, need to obtain the fine-grid index
      -- the of the cells needed by the coarsening stencil.
      -- For equal 2x coarsening in all directions, doubling the
      -- (zero-based) index puts you at the lower-corner
      -- of the hyper-cube needed in the fine grid. Each cell after that
      -- can be obtained by iterating through the dimensions, and stepping
      -- away from the cells already counted (starting with this corner one).
      for d = 1, self.dim do self.fineGridIdx[1][d] = 2*(cIdx[d]-1)+1 end
      local fIdxCount = 1
      for dir = 1, self.dim do
         for rI = 1, self.igOpStencilWidth^(dir-1) do
            fIdxCount = fIdxCount + 1
            for d = 1, self.dim do self.fineGridIdx[fIdxCount][d] = self.fineGridIdx[rI][d] end
            self.fineGridIdx[fIdxCount][dir] = self.fineGridIdx[rI][dir]+1
         end
      end
      
      cFld:fill(cFldIndexer(cIdx), cFldItr)  -- Coarse-grid field pointer.
  
      -- Array of pointers to fine-grid field data in cells pointed to by the coarsening stencil. 
      for i = 1, self.igOpStencilSize do
         fFld:fill(fFldIndexer(self.fineGridIdx[i]), fFldItr)
         self.fineFldItr[i] = fFldItr:data()
      end
  
      self._restriction(self.fineFldItr:data(), cFldItr:data())
   end
end

function MGpoisson:relaxStencilIndices(idxIn, stencilType)
   -- Given the index of the current cell (idxIn), return
   -- a table of indices of the cells in a stencil of
   -- type 'stencilType' used in the relaxtion operator.
   -- Stencil types are given by [t,w,l]:
   --   t: type of stencil
   --        0 : 'cross'.
   --        1 : partially filled cross (include nearest corner cells).
   --        2 : filled cross (include all corner cells).
   --   w: width of the stencil (along each dimension).
   --   l: location in grid along each dimension. For 1D
   --        [-1] = lower boundary cell.
   --        [ 0] = inner cell.
   --        [ 1] = upper boundary cell.
   --
   -- So stencilType=[0,[5,5],[0,0]] corresponds to 
   --                 8
   --                 6
   --         5   3   1   2   4
   --                 7
   --                 9

   -- First copy all indicies since (in higher dimensions)
   -- most of the stay the same for each cell.
   for i = 1, self.relaxStencilSize do
      idxIn:copyInto(self.relaxIdx[i])
   end

   if stencilType[1] == 0 then
      local sI = 1
      for d = 1, self.dim do
         for pm = 1,self.relaxStencilWidth-1 do
            sI = sI + 1
            self.relaxIdx[sI][d] = idxIn[d]-((-1)^(pm % 2))*((self.relaxStencilWidth-1)/2)
         end
      end
   end

end

function MGpoisson:relax(phiFld, rhoFld)
   -- Smoother/relaxation.

   local grid = fld:grid() 

   localRangeDecomp = LinearDecomp.LinearDecompRange {
      range = fld:localRange(), numSplit = grid:numSharedProcs() }
   local tId        = grid:subGridSharedId()    -- Local thread ID.

   local indexer = phiFld:genIndexer()

   local phiItr = phiFld:get(1)
   local rhoItr = rhoFld:get(1)

   for idx in localRangeDecomp:rowMajorIter(tId) do

      grid:setIndex(idx)

      -- Cell lengths and right-side source (rho) in this cell.
      grid:getDx(self.relaxDxs[1])
      rhoFld:fill(indexer(idx), rhoItr)   
  
      -- Get with indices of cells used by stencil.
      self:relaxStencilIndices(idx,{0,{3,3},{0,0}})

      -- Array of pointers to cell lengths and phi data in cells pointed to by the stencil. 
      for i = 1, sI do
         grid:setIndex(self.relaxIdx[i])
         grid:getDx(self.relaxDxs[i])
         phiFld:fill(indexer(self.relaxIdx[i]), phiItr)
         self.relaxItr[i] = phiItr:data()
      end

      self._relaxation(self.relaxDxs:data(), self.relaxItr:data(), rhoItr:data())
   end
end

function MGpoisson:areIndicesOdd(idxIn) 
   -- Identify if indices in all dimensions are odd.
   local areOdd = true
   for d = 1, self.dim do
      areOdd = areOdd and (idxIn[d] % 2 == 1)
   end
   return areOdd
end

function MGpoisson:prolong(cFld,fFld)
   -- Prolongation of a coarse-grid field (cFld) to a fine-grid field (fFld). 

   local grid = fFld:grid() 

   localRangeDecomp  = LinearDecomp.LinearDecompRange {
      range = fFld:localRange(), numSplit = grid:numSharedProcs() }
   local tId         = grid:subGridSharedId()    -- Local thread ID.

   local cFldIndexer = cFld:genIndexer()
   local cFldItr     = cFld:get(1)

   local fFldIndexer = fFld:genIndexer()
   local fFldItr     = fFld:get(1)

   for fIdx in localRangeDecomp:rowMajorIter(tId) do

      grid:setIndex(fIdx)

      if self:areIndicesOdd(fIdx) then

         -- Array of fine grid indexes to cells used by stencil.
         for d = 1, self.dim do self.fineGridIdx[1][d] = fIdx[d] end
         local fIdxCount = 1
         for dir = 1, self.dim do
            for rI = 1, self.igOpStencilWidth^(dir-1) do
               fIdxCount = fIdxCount + 1
               for d = 1, self.dim do self.fineGridIdx[fIdxCount][d] = self.fineGridIdx[rI][d] end
               self.fineGridIdx[fIdxCount][dir] = self.fineGridIdx[rI][dir]+1
            end
         end

         -- Given the fine-grid index, need to obtain the coarse-grid index.
         for d = 1, self.dim do self.coarseGridIdx[d] = (fIdx[d]-1)/2+1 end
  
  
         cFld:fill(cFldIndexer(self.coarseGridIdx), cFldItr)   -- Coarse field pointer.
  
         -- Array of pointers to fine-grid field data by stencil. 
         for i = 1, self.igOpStencilSize do
            fFld:fill(fFldIndexer(self.fineGridIdx[i]), fFldItr)
            self.fineFldItr[i] = fFldItr:data()
         end
  
         self._prolongation(cFldItr:data(), self.fineFldItr:data())
      end
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
   self.rho[1] = srcFld
   for i = 2, self.mgLevels do
      self:restrict(self.rho[i-1], self.rho[i])
      self.rho[i]:write(string.format("rho_%d.bp", i), tCurr, 0)
   end


   -- Prolong the right-side source field (just for testing).
   self:prolong(self.rho[2], outFld[1])

--   outFld[1]:copy(srcFld)   -- Temporary arbitrary output.

end

return MGpoisson
