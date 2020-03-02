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
--   3) MF: I wish to change the stencil numbering to match kernel IDs.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local MGpoissonDecl         = require "Updater.mgPoissonCalcData.MGpoissonModDecl"
local Proto                 = require "Lib.Proto"
local DirectDGPoissonSolver = require "Updater.DiscontPoisson"
local UpdaterBase           = require "Updater.Base"
local Grid                  = require "Grid"
local DataStruct            = require "DataStruct"
local DecompRegionCalc      = require "Lib.CartDecomp"
local LinearDecomp          = require "Lib.LinearDecomp"
local Lin                   = require "Lib.Linalg"
local ffi                   = require "ffi"
local lume                  = require "Lib.lume"
local IntQuantCalc          = require "Updater.CartFieldIntegratedQuantCalc"

-- Boundary condition ID numbers.
local BVP_BC_PERIODIC  = 0
local BVP_BC_DIRICHLET = 1
local BVP_BC_NEUMANN   = 2
local BVP_BC_ROBIN     = 3    -- Not supported yet.

-- Multigrid updater object.
local MGpoisson = Proto(UpdaterBase)

local relaxationKinds = {
   "Jacobi", "DampedJacobi", "GaussSeidel", "DampedGaussSeidel"
}
local dampedRelaxationKinds = {
   "DampedJacobi", "DampedGaussSeidel"
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

   local grid = assert(
      tbl.onGrid, "Updater.MGpoisson: Must provide grid object using 'grid'.")

   local basis = assert(
      tbl.basis, "Updater.MGpoisson: Must provide the weak basis object using 'basis'.")


   -- ~~............................ Multigrid parameters ..............................~~ --
   -- Relaxation method. 
   local relaxKind = ""
   if tbl.relaxType then
      relaxKind = tbl.relaxType
   else
      relaxKind = "DampedGaussSeidel"
   end
   if self:isRelaxKindGood(relaxKind) then
      if self:isRelaxDamped(relaxKind) then
         if tbl.relaxOmega then
            self.omega = tbl.relaxOmega
         else
            self.omega = 1.8   -- Arbitrary. Ideally the user inputs it, or we compute it later.
         end
      else
         self.omega = 1.0                    -- No damping. Set to 1.
         relaxKind  = "Damped" .. relaxKind  -- All kernels use "Damped".
      end
      if (relaxKind == "DampedJacobi") then
         -- Need to be Jacobi-aware as this involves an extra field copy.
         self.isJacobiRelax = true
      else
         self.isJacobiRelax = false
      end
   else
      assert(false, "Updater.MGpoisson: relaxType must be one of 'Jacobi', 'DampedJacobi', 'GaussSeidel', 'DampedGaussSeidel'")
   end

   -- Number of pre-, post- and coarsest-grid relaxations.
   if tbl.relaxNum then
      self.nu1 = tbl.relaxNum[1]
      self.nu2 = tbl.relaxNum[2]
      self.nu3 = tbl.relaxNum[3]
   else
      self.nu1 = 2
      self.nu2 = 2
      self.nu3 = 1000
   end

   -- The gamma parameter controls the type of MG cycle: =1 V-cycle, =2 W-cycle.
   if tbl.gamma then
      self.gamma = tbl.gamma
   else
      self.gamma = 1
   end

   if (tbl.tolerance) then
      -- User-defined tolerance (stopping point) if given.
      self.tol            = tbl.tolerance
      self.numGammaCycles = 1000
   else
      -- If not given perform the user-defined number of cycles
      -- or stop when the relative residual norm is 1e-12.
      self.numGammaCycles = assert(tbl.numCycles, "Updater.MGpoisson: if 'tolerance' is not specified, must provide the number of cycles with 'numCycles'.")
      self.tol = 1.0e-12
   end
   -- ~~.................... End of user-input multigrid parameters ......................~~ --

   -- Diagnostics flag indicates whether to write diagnostic info.
   if tbl.diagnostics then self.diagnostics=tbl.diagnostics else self.diagnostics=false end

   self.dim        = basis:ndim()        -- Dimension of space.
   local polyOrder = basis:polyOrder()   -- Polynomial order.
   local basisID   = basis:id()          -- Basis kind.

   -- Read boundary conditions.
   -- Assume the inputs are two tables, bcLower and bcUpper, with each table
   -- having a two element table for each dimension, For example, for a 2D sim
   -- use bcLower = {{T = sT1, V = val1}, {T = sT2, V = val2}} and analogously
   -- for bcUpper, where these are
   --    sT:  string indicating BC type ("P", "D", "N", "R").
   --    val: BC value.
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

   local bcLower = assert(
      tbl.bcLower, "Updater.MGpoisson: Must provide lower-boundary BCs along each direction in 'bcLower'.")
   local bcUpper = assert(
      tbl.bcUpper, "Updater.MGpoisson: Must provide upper-boundary BCs along each direction in 'bcUpper'.")
   assert(#bcLower == self.dim, "Updater.MGpoisson: number of lower BCs must equal number of dimensions.")
   assert(#bcUpper == self.dim, "Updater.MGpoisson: number of lower BCs must equal number of dimensions.")
   local bcTypes  = {} 
   local bcValues = {} 
   for i = 1, self.dim do
      bcTypes[i]  = {bcID(bcLower[i]["T"]), bcID(bcUpper[i]["T"])}
      bcValues[i] = {bcLower[i]["V"], bcUpper[i]["V"]}
   end
   -- Establish periodic directions.
   local periodicDirs  = {}
   local isDirPeriodic = {}
   for d = 1, self.dim do
      if (bcTypes[d] == 0) then
         lume.push(periodicDirs,d)
         isDirPeriodic[d] = true
         bcValues[d]      = {0.0, 0.0}   -- Not used, but a nil could cause problems.
      else
         isDirPeriodic[d] = false
      end
   end
   local isPeriodicDomain = lume.all(isDirPeriodic)

   -- Translate bcValues to a vector from which we can pass a pointer.
   self.bcValue = Lin.Vec(self.dim*2)
   for d = 1,self.dim do 
     self.bcValue[2*d-1] = bcValues[d][1]
     self.bcValue[2*d]   = bcValues[d][2]
   end


   -- Create a grid for each level.
   -- Not sure this is needed, but in general it probably is (e.g. unstructured, or even nonuniform meshes).
   self.mgGrids     = {}
   self.mgGrids[1]  = grid
   periodicDirCount = 0
   -- Iterate (phi), right-side source and residue fields at each level.
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
         lowerC[d] = grid:lower(d)
         upperC[d] = grid:upper(d)
         if grid:isDirPeriodic(d) then
            periodicDirCount = periodicDirCount+1
            periodicDirsC[periodicDirCount] = d
         end

         -- The following two are more complicated. To be refined later.
         cellsC[d]      = (self.mgGrids[self.mgLevels-1]:numCells(d))/2
         decompCutsC[d] = self.mgGrids[self.mgLevels-1]:cuts(d) 

         if cellsC[d] == 2 then notAtCoarsest = false end
      end

      local isSharedC = grid:isShared()

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

      if not notAtCoarsest then
         self.directSolver = DirectDGPoissonSolver {
            onGrid  = self.mgGrids[self.mgLevels],
            basis   = basis,
            bcLower = bcLower,
            bcUpper = bcUpper,
         }
      end
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
      -- the residue on each grid level.
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
   if self.isJacobiRelax then
      self.phiPrevAll = {}
      for i = 1, self.mgLevels do
         -- For Jacobi relaxation need an extrac copy of the field iterate.
         self.phiPrevAll[i] = DataStruct.Field {
            onGrid        = self.mgGrids[i],
            numComponents = basis:numBasis(),   -- NOTE: this will change if we do p-coarsening.
            ghost         = {1, 1},
         }
      end
   end

   -- Select restriction and prolongation operator kernels.
   self._restriction  = MGpoissonDecl.selectRestriction(basisID, self.dim, polyOrder)
   self._prolongation = MGpoissonDecl.selectProlongation(basisID, self.dim, polyOrder)
   -- Select kernels for relaxation and computing the residue.
   self._relaxation   = MGpoissonDecl.selectRelaxation(basisID, self.dim, polyOrder, relaxKind, bcTypes)
   self._calcResidue  = MGpoissonDecl.selectResidueCalc(basisID, self.dim, polyOrder, bcTypes)

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

   -- Relaxation stencil info: restrict ourselves to 'cross' relaxation stencils for now.
   self.opStencilWidth = 3
   self.opStencilSize  = (self.opStencilWidth-1)*self.dim+1
   self.opStencilIdx   = {}   -- List of cell indices pointed to by the stencil.
   for i = 1, self.opStencilSize do
      self.opStencilIdx[i] = Lin.IntVec(self.dim)
   end
   -- List of pointers to the data in cells pointed to by the stencil.
   self.opStencilItr   = DoublePtrVec(self.opStencilSize)
   self.opStencilDxs   = DoublePtrVec(self.opStencilSize)
   self.prevIterateItr = DoublePtrVec(self.opStencilSize)   -- Only used for Jacobi relaxations.

   self.dxBuf = Lin.Vec(self.dim)    -- A buffer to store cell lengths.

   -- Dimensions remaining when a dimension is removed.
   self.dimRemain = {}
   for d1 = 1, self.dim do
      self.dimRemain[d1] = {}
      for d2 = 1, self.dim do self.dimRemain[d1][d2] = d2 end
      table.remove(self.dimRemain[d1],d1)
   end

   -- Updater to compute the L2-norm of the residue.
   self.l2NormCalc = IntQuantCalc {
      onGrid   = grid, 
      basis    = basis,
      quantity = "RmsV",
   }
   self.relResNorm = DataStruct.DynVector {
      numComponents = 1,
   }
   self.residueNorm = DataStruct.DynVector {
      numComponents = 1,
   }
   self.rhoNorm = DataStruct.DynVector {
      numComponents = 1,
   }

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

function MGpoisson:opStencilIndices(idxIn, stencilType)
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
   for i = 1, self.opStencilSize do
      idxIn:copyInto(self.opStencilIdx[i])
   end

   if stencilType[1] == 0 then
      local sI = 1
      for d = 1, self.dim do
         for pm = 1,self.opStencilWidth-1 do
            sI = sI + 1
            self.opStencilIdx[sI][d] = idxIn[d]-((-1)^(pm % 2))*((self.opStencilWidth-1)/2)
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

function MGpoisson:jacobiCopyField(fldIn,fldOutAll)
   -- Need to copy the current iterate (phi). It would be easier
   -- to do so if the current level (lCurr) was available within
   -- the relax method, but currently it is not. So for now we will
   -- compare the grids until we find the right field to copy
   -- into and return the grid level so relax knows which one to use.
   local grid   = fldIn:grid()
   local cellsN = {}
   for d = 1, self.dim do cellsN[d]=grid:numCells(d) end

   local currLevel
   for l = 1, self.mgLevels do
      currLevel = l
      local otherGrid     = fldOutAll[l]:grid()
      local areGridsEqual = true
      for d = 1, self.dim do
         areGridsEqual = areGridsEqual and (cellsN[d] == otherGrid:numCells(d))
      end
      -- We will need to compare polyOrder as well when doing p-coarsening.
      if areGridsEqual then break end
   end

   fldOutAll[currLevel]:copy(fldIn)

   return currLevel
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

   local phiPrevItr, currLevel

   for nuI = 1, numRelax do    -- Relax numRelax times.

      if self.isJacobiRelax then
         currLevel  = self:jacobiCopyField(phiFld,self.phiPrevAll) 
         phiPrevItr = self.phiPrevAll[currLevel]:get(1)
      end

      for idx in localRangeDecomp:rowMajorIter(tId) do
   
         grid:setIndex(idx)
   
         -- Cell lengths and right-side source (rho) in this cell.
         grid:getDx(self.dxBuf)
         self.opStencilDxs[1] = self.dxBuf:data()
         rhoFld:fill(indexer(idx), rhoItr)   
     
         -- Get with indices of cells used by stencil. Store them in self.opStencilIdx.
         self:opStencilIndices(idx,{0,{3,3},{0,0}})
   
         -- Array of pointers to cell lengths and phi data in cells pointed to by the stencil. 
         for i = 1, self.opStencilSize do
            grid:setIndex(self.opStencilIdx[i])
            grid:getDx(self.dxBuf)
            self.opStencilDxs[i] = self.dxBuf:data()

            phiFld:fill(indexer(self.opStencilIdx[i]), phiItr)
            self.opStencilItr[i] = phiItr:data()

            if self.isJacobiRelax then
               self.phiPrevAll[currLevel]:fill(indexer(self.opStencilIdx[i]), phiPrevItr)
               self.prevIterateItr[i] = phiPrevItr:data()
            end
         end
         
         self._relaxation[self:idx2stencil(idx,cellsN)](self.omega, self.opStencilDxs:data(), self.bcValue:data(), rhoItr:data(), self.prevIterateItr:data(), self.opStencilItr:data())
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
   local resItr = resFld:get(1)

   for idx in localRangeDecomp:rowMajorIter(tId) do
   
      grid:setIndex(idx)
   
      -- Cell lengths, right-side source (rho) and residue in this cell.
      grid:getDx(self.dxBuf)
      self.opStencilDxs[1] = self.dxBuf:data()
      rhoFld:fill(indexer(idx), rhoItr)   
      resFld:fill(indexer(idx), resItr)   
   
      -- Get with indices of cells used by stencil.
      self:opStencilIndices(idx,{0,{3,3},{0,0}})
   
      -- Array of pointers to cell lengths and phi data in cells pointed to by the stencil. 
      for i = 1, self.opStencilSize do
         grid:setIndex(self.opStencilIdx[i])
         grid:getDx(self.dxBuf)
         self.opStencilDxs[i] = self.dxBuf:data()
         phiFld:fill(indexer(self.opStencilIdx[i]), phiItr)
         self.opStencilItr[i] = phiItr:data()
      end
   
      self._calcResidue[self:idx2stencil(idx,cellsN)](self.opStencilDxs:data(), self.bcValue:data(), rhoItr:data(), self.opStencilItr:data(), resItr:data())
   end
end

function MGpoisson:relResidueNorm(gamIdx)
   -- Compute the relative norm of the residue: ||rho + L(phi)||/||rho||.

   -- Compute the norm of the right-side source vector.
   self.l2NormCalc:advance(1,{self.rhoAll[1]},{self.rhoNorm})
   local _, rhsNorm = self.rhoNorm:lastData()
   -- Compute the norm of the residue.
   self:residue(self.phiAll[1], self.rhoAll[1], self.residueAll[1]) 
   self.l2NormCalc:advance(gamIdx,{self.residueAll[1]},{self.residueNorm})
   -- Compute the relative residue norm and store it in self.relResNorm.
   local _, resNorm    = self.residueNorm:lastData()
   local relResNormOut = resNorm[1]/rhsNorm[1]
   self.relResNorm:appendData(gamIdx, {relResNormOut})

   return relResNormOut
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

      if self.directSolver then
         -- Call the direct solver.
         self.directSolver:advance(0.0, {self.rhoAll[lCurr]}, {self.phiAll[lCurr]})
      else
         -- Relax nu3 times.
         self:relax(self.nu3, self.phiAll[lCurr], self.rhoAll[lCurr]) 
      end

   else

      -- Relax nu1 times.
      self:relax(self.nu1, self.phiAll[lCurr], self.rhoAll[lCurr]) 

      -- Compute the residue.
      self:residue(self.phiAll[lCurr], self.rhoAll[lCurr], self.residueAll[lCurr]) 

      -- Restrict the residue to the next coarsest grid.
      self:restrict(self.residueAll[lCurr], self.rhoAll[lCurr+1])

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
   -- Assume the given phi iterate is an initial guess if given
   -- within inFld. If not it will interpret that to mean that the
   -- user wishes to perform a full-multigrid (FMG) cycle.
   self.rhoAll[1]       = inFld[1]
   local initialGuess   = inFld[2]
   local relResNormCurr = 1.0e12    -- Current (relative) residue norm.
   if initialGuess then
      self.phiAll[1]  = initialGuess
   else
      -- No initial guess provided. Perform Full Multi-Grid (FMG).
      self.phiAll[1]  = outFld[1]

      -- FMG requires we restrict the right-side source field to all levels.
      for i = 2, self.mgLevels do
         self:restrict(self.rhoAll[i-1], self.rhoAll[i])
      end

      -- In FMG we start at the coarsest grid without an initial guess.
      self.phiAll[self.mgLevels]:clear(0.0)
      for i = self.mgLevels, 1, -1 do
         self:gammaCycle(i)
         if (i > 1) then self:prolong(self.phiAll[i], self.phiAll[i-1]) end
      end

   end

   local gI = 0         -- gamma cycle index.
   -- Compute the relative residue norm. It gets stored in self.relResNorm. 
   relResNormCurr = self:relResidueNorm(gI)

   -- Call MG gamma-cycles starting at the finest grid.
   while relResNormCurr > self.tol do
      gI = gI + 1
      self:gammaCycle(1)

      -- Compute the relative residue norm. It gets stored in self.relResNorm. 
      relResNormCurr = self:relResidueNorm(gI)

      if gI==self.numGammaCycles then break end
   end

   if self.diagnostics then
      -- Write out residue norm for each iteration.
      self.relResNorm:write(string.format("relResidue_RmsV.bp"), 0.0, 0)
   end

end

return MGpoisson
