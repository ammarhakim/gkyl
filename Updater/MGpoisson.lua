-- Gkyl ------------------------------------------------------------------------
--
-- Multigrid solver for the Poisson equation
--     - L(phi) = rho
-- where L stands for Laplacian.
--
--
-- Question/development notes:
--   0) I think it'd be better if the user does not have to indicate periodicDirs and
--      BCs in separate tables, but rather periodic is another type within the same table.
--   1) Do we need to create a grid object for each level? or create pseudo-grid object that is more lightweight and has the information needed?
--   2) Prolongation: so data doesn't get loaded multiple times, instead of iterating through each fine-grid cell, could we instead make
--                    the looping stencil-aware?
--   3) MF: I wish to change the stencil numbering to match kernel IDs.
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
local lume             = require "Lib.lume"

-- Boundary condition ID numbers.
local BVP_BC_PERIODIC  = 0
local BVP_BC_DIRICHLET = 1
local BVP_BC_NEUMANN   = 2

-- Multigrid updater object.
local MGpoisson = Proto(UpdaterBase)

local relaxationKinds = {
   "GaussSeidel", "DampedGaussSeidel"
}
local dampedRelaxationKinds = {
   "DampedGaussSeidel"
}
function MGpoisson:isRelaxKindGood(rlx)
   if lume.find(relaxationKinds, rlx) then
      return true
   end
   return false
end

function MGpoisson:isRelaxDamped(rlx)
   if lume.find(dampedRelaxationKinds, rlx) then
      return true
   end
   return false
end

function MGpoisson:init(tbl)
   MGpoisson.super.init(self, tbl) -- Setup base object.

   local topGrid = assert(
      tbl.onGrid, "Updater.MGpoisson: Must provide grid object using 'grid'.")

   local basis = assert(
      tbl.basis, "Updater.MGpoisson: Must provide the weak basis object using 'basis'.")


   -- Multigrid parameters.
   local relaxKind = assert(
      tbl.relaxType, "Updater.MGpoisson: Must provide the type of relaxations using 'relaxType'.")
   local nus = assert(
      tbl.relaxNum, "Updater.MGpoisson: Must provide the number of relaxations using 'relaxNum'.")
   self.gamma = assert(
      tbl.gamma, "Updater.MGpoisson: Must provide the MG gamma parameter using 'gamma'.")
   if self:isRelaxKindGood(relaxKind) then
      if self:isRelaxDamped(relaxKind) then
         self.omega = assert(
            tbl.relaxOmega, "Updater.MGpoisson: Must provide the relaxation damping using 'relaxOmega'.")
      else
         self.omega = 1.0
      end
   else
      assert(false, "Updater.MGpoisson: relaxType must be one of 'GaussSeidel', 'DampedGaussSeidel'")
   end

   -- Number of pre-, post- and coarsest-grid relaxations.
   self.nu1 = nus[1]
   self.nu2 = nus[2]
   self.nu3 = nus[3]

   self.dim        = basis:ndim()        -- Dimension of space.
   local polyOrder = basis:polyOrder()   -- Polynomial order.
   local basisID   = basis:id()          -- Basis kind.

   -- Establish periodicity of the domain.
   local periodicDirs    = {}
   local nonPeriodicDirs = {}
   for d = 1, self.dim do nonPeriodicDirs[d] = d end
   if tbl.periodicDirs then
      for i, d in ipairs(tbl.periodicDirs) do
         if d<1 or d>self.dim then
            assert(false, "Updater.MGpoisson: Directions in periodicDirs table should be 1 (for X), 2 (for Y), or 3 (for Z)")
         end
         periodicDirs[i]  = d
         lume.remove(nonPeriodicDirs, d)
      end
   end
   local isDirPeriodic    = {}
   for d = 1,self.dim do isDirPeriodic[d]=false end
   for _, d in ipairs(periodicDirs) do isDirPeriodic[d]=true end
   local isPeriodicDomain = lume.all(isDirPeriodic)

   -- Read non-periodic boundary conditions.
   -- Assume the input poissonBCs is a table, one entry for
   -- each non-periodic BC. Each of these entries is a table
   -- like {dir, {T = sL, V = vL}, {T = sU, V = vU}}, where these are
   --    dir:   integer indicating the direction.
   --    sL,sU: string indicating BC type ("P", "D", "N") on lower and upper boundaries.
   --    vL,vU: BC values on lower and upper boundaries.
   local function bcID(strIn)
      -- Given a string indicating BC type, return the corresponding
      -- BC ID number, defined a the top of this file.
      if strIn == "P" then return BVP_BC_PERIODIC
      elseif strIn == "D" then return BVP_BC_DIRICHLET
      elseif strIn == "N" then return BVP_BC_NEUMANN
      else
         assert(false, "Updater.MGpoisson: BC type (T) must be one of P, D, or N. Used " .. strIn .. " instead.")
      end
   end

   local topGridBCtypes  = {} 
   local topGridBCvalues = {} 
   if (tbl.poissonBCs) then
      for i = 1, #tbl.poissonBCs do
         topGridBCtypes[tbl.poissonBCs[i][1]]  = {bcID(tbl.poissonBCs[i][2]["T"]), bcID(tbl.poissonBCs[i][3]["T"])}
         topGridBCvalues[tbl.poissonBCs[i][1]] = {tbl.poissonBCs[i][2]["V"], tbl.poissonBCs[i][3]["V"]}
      end
      if #topGridBCtypes ~= #nonPeriodicDirs then
         assert(false, "Updater.MGpoisson: For a non-periodic domain must specify a BC for each non-peridic direction.")
      end
   elseif (not isPeriodicDomain) then
      assert(false, "Updater.MGpoisson: For a non-periodic domain must specify BCs in 'poissonBCs'.")
   end
   for _, d in ipairs(periodicDirs) do
      topGridBCtypes[d]  = {BVP_BC_PERIODIC, BVP_BC_PERIODIC}
      topGridBCvalues[d] = 0.0 
   end
   -- Translate BCvalues to a vector from which we can pass a pointer.
   self.topGridBCvals = Lin.Vec(self.dim*2)
   for d = 1,self.dim do 
     self.topGridBCvals[2*d-1] = topGridBCvalues[d][1]
     self.topGridBCvals[2*d]   = topGridBCvalues[d][2]
   end


   -- Create a grid for each level.
   -- Not sure this is needed, but in general it probably is (e.g. unstructured, or even nonuniform meshes).
   self.mgGrids     = {}
   self.mgGrids[1]  = topGrid
   periodicDirCount = 0
   -- Right-side source and residue fields at each level.
   self.phiAll     = {}
   self.rhoAll     = {}
   self.residueAll = {}

   self.mgLevels       = 1      -- Number of multigrid levels (grids).
   local notAtCoarsest = true
   while notAtCoarsest do
   
      self.mgLevels = self.mgLevels+1

      -- Determine parameters of the next coarse grid.
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
         cellsC[d]      = (self.mgGrids[self.mgLevels-1]:numCells(d))/2
         decompCutsC[d] = self.mgGrids[self.mgLevels-1]:cuts(d) 

         if cellsC[d] == 2 then notAtCoarsest = false end
      end

      local isSharedC = topGrid:isShared()

      local decompC = DecompRegionCalc.CartProd {
         cuts      = decompCutsC,
         useShared = useSharedC,
      }

      self.mgGrids[self.mgLevels] = Grid.RectCart {
         lower         = lowerC, 
         upper         = upperC,
         cells         = cellsC,
         periodicDirs  = periodicDirsC,
         decomposition = decompC,
      }
   end

   -- Allocate space for the iterate, right-side source field and
   -- the residue on each grid level. We don't need to allocate space for
   -- phi and rho on the top grid because we'll just use the fields
   -- given to the solver.
   self.residueAll[1] = DataStruct.Field {
      onGrid        = self.mgGrids[1],
      numComponents = basis:numBasis(),   -- NOTE: this will change if we do p-coarsening.
      ghost         = {1, 1},
   }
   for i = 2, self.mgLevels do
      -- Allocate space for the iterate, right-side source field and
      -- the residu on each grid level.
      self.phiAll[i] = DataStruct.Field {
         onGrid        = self.mgGrids[i],
         numComponents = basis:numBasis(),   -- NOTE: this will change if we do p-coarsening.
         ghost         = {1, 1},
      }
      self.rhoAll[i] = DataStruct.Field {
         onGrid        = self.mgGrids[i],
         numComponents = basis:numBasis(),   -- NOTE: this will change if we do p-coarsening.
         ghost         = {1, 1},
      }
      self.residueAll[i] = DataStruct.Field {
         onGrid        = self.mgGrids[i],
         numComponents = basis:numBasis(),   -- NOTE: this will change if we do p-coarsening.
         ghost         = {1, 1},
      }
   end

   -- Select restriction and prolongation operator kernels.
   self._restriction  = MGpoissonDecl.selectRestriction(basisID, self.dim, polyOrder)
   self._prolongation = MGpoissonDecl.selectProlongation(basisID, self.dim, polyOrder)
   -- Select kernels for relaxation and computing the residue.
   self._relaxation   = MGpoissonDecl.selectRelaxation(basisID, self.dim, polyOrder, relaxKind, topGridBCtypes)
   self._calcResidue  = MGpoissonDecl.selectResidueCalc(basisID, self.dim, polyOrder, topGridBCtypes)

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
   self.fineFldItr    = DoublePtrVec(self.igOpStencilSize)

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

   self.dxBuf = Lin.Vec(self.dim)    -- A buffer to store cell lengths.

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
   -- MF update: I wish to change the stencil numbering to match kernel IDs.
   -----------------
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

function MGpoisson:idx2stencil(idxIn, nCellsIn)
   -- Given a multi-dimensional index (idxIn) to a cell in a grid 
   -- with nCellsIn cells, return the index of the stencil needed, 
   -- within the table that holds relaxation/residue stencils.
   local stencilIdx = 1
   for d = 1, self.dim do
      if (idxIn[d] == 1) then                -- First cell.
         stencilIdx = 2*stencilIdx + 3^(d-1) - 1
      elseif (idxIn[d] == nCellsIn[d]) then  -- Last cell.
         stencilIdx = 2*stencilIdx + 3^(d-1)
      end
   end
   return stencilIdx
end

function MGpoisson:relax(numRelax, phiFld, rhoFld)
   -- Perform numRelax relaxations of the Poisson equation,

   local grid   = phiFld:grid() 
   local cellsN = {}
   for d = 1, self.dim do cellsN[d]=grid:numCells(d) end

   localRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phiFld:localRange(), numSplit = grid:numSharedProcs() }
   local tId        = grid:subGridSharedId()    -- Local thread ID.

   local indexer = phiFld:genIndexer()

   local phiItr = phiFld:get(1)
   local rhoItr = rhoFld:get(1)

   for nuI = 1, numRelax do    -- Relax numRelax times.
      for idx in localRangeDecomp:rowMajorIter(tId) do
   
         grid:setIndex(idx)
   
         -- Cell lengths and right-side source (rho) in this cell.
         grid:getDx(self.dxBuf)
         self.relaxDxs[1] = self.dxBuf:data()
         rhoFld:fill(indexer(idx), rhoItr)   
     
         -- Get with indices of cells used by stencil.
         self:relaxStencilIndices(idx,{0,{3,3},{0,0}})
   
         -- Array of pointers to cell lengths and phi data in cells pointed to by the stencil. 
         for i = 1, self.relaxStencilSize do
            grid:setIndex(self.relaxIdx[i])
            grid:getDx(self.dxBuf)
            self.relaxDxs[i] = self.dxBuf:data()

            phiFld:fill(indexer(self.relaxIdx[i]), phiItr)
            self.relaxItr[i] = phiItr:data()
         end

         self._relaxation[self:idx2stencil(idx,cellsN)](self.omega, self.relaxDxs:data(), self.topGridBCvals:data(), rhoItr:data(), self.relaxItr:data())
      end
   end
end

function MGpoisson:residue(phiFld, rhoFld, resFld)
   -- Compute the residue:
   --     r = rho + L(phi). 
   -- where L is the Laplacian.

   local grid = phiFld:grid() 
   local cellsN = {}
   for d = 1, self.dim do cellsN[d]=grid:numCells(d) end

   localRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phiFld:localRange(), numSplit = grid:numSharedProcs() }
   local tId        = grid:subGridSharedId()    -- Local thread ID.

   local indexer = phiFld:genIndexer()

   local phiItr = phiFld:get(1)
   local rhoItr = rhoFld:get(1)
   local resItr = rhoFld:get(1)

   for idx in localRangeDecomp:rowMajorIter(tId) do
   
      grid:setIndex(idx)
   
      -- Cell lengths, right-side source (rho) and residue in this cell.
      grid:getDx(self.dxBuf)
      self.relaxDxs[1] = self.dxBuf:data()
      rhoFld:fill(indexer(idx), rhoItr)   
      resFld:fill(indexer(idx), resItr)   
   
      -- Get with indices of cells used by stencil.
      self:relaxStencilIndices(idx,{0,{3,3},{0,0}})
   
      -- Array of pointers to cell lengths and phi data in cells pointed to by the stencil. 
      for i = 1, self.relaxStencilSize do
         grid:setIndex(self.relaxIdx[i])
         grid:getDx(self.dxBuf)
         self.relaxDxs[i] = self.dxBuf:data()
         phiFld:fill(indexer(self.relaxIdx[i]), phiItr)
         self.relaxItr[i] = phiItr:data()
      end
   
      self._calcResidue[self:idx2stencil(idx,cellsN)](self.relaxDxs:data(), self.topGridBCvals:data(), rhoItr:data(), self.relaxItr:data(), resItr:data())
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

function MGpoisson:gammaCycle(lCurr)
   -- Perform a single gamma-cycle at the lCurr grid level.
   --   gamma=1 V-cycle
   --   gamma=2 W-cycle
   -- What about F cycle? 

   if lCurr == self.mgLevels then

      -- Coarsest grid. Use a direct solver or many iterations.
      -- The latter is useful if the coarsest grid is large still.

      -- Relax nu3 times.
      self:relax(self.nu3, self.phiAll[lCurr], self.rhoAll[lCurr]) 

   else

      -- Relax nu1 times.
      self:relax(self.nu1, self.phiAll[lCurr], self.rhoAll[lCurr]) 

      -- Compute the residue.
      self:residue(self.phiAll[lCurr], self.rhoAll[lCurr], self.residueAll[lCurr]) 

      -- Restrict the residue to the next coarsest grid.
      self:restrict(self.residueAll[lCurr], self.residueAll[lCurr+1])

      -- Solve the problem on the coarser grid, gamma times, with an
      -- initial guess of zero.
      self.phiAll[lCurr+1]:clear(0.0)
      for gamI = 1, self.gamma do
         self:gammaCycle(lCurr+1) 
      end

      -- Prolong the error to the finer grid and correct the iterate.
      self:prolong(self.phiAll[lCurr+1], self.residueAll[lCurr])
      self.phiAll[lCurr]:accumulate(1.0,self.residueAll[lCurr])

      -- Relax nu2 times.
      self:relax(self.nu2, self.phiAll[lCurr], self.rhoAll[lCurr]) 

   end
   
end

-- Advance method.
function MGpoisson:_advance(tCurr, inFld, outFld)

   -- Phi iterate and right-side source fields on finest grid.
   -- Assume the given phi iterate is an initial guess, unless
   -- performing a full-multigrid (FMG).
   self.rhoAll[1] = inFld[1]
   self.phiAll[1] = outFld[1]

   -- Restrict the right-side source field.
   for i = 2, self.mgLevels do
      self:restrict(self.rhoAll[i-1], self.rhoAll[i])
   end

   -- Call a MG gamma-cycle starting at the finest grid.
   self:gammaCycle(1)


end

return MGpoisson
