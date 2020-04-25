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
local BVP_BC_ROBIN     = 3

-- Basis translation direction.
local DG_to_FEM = -1
local FEM_to_DG =  1

-- Multigrid updater object.
local MGpoisson = Proto(UpdaterBase)

local solverTypes           = { "DG", "FEM" }
local relaxationKinds       = {
   "Jacobi", "DampedJacobi", "GaussSeidel", "DampedGaussSeidel"
}
local dampedRelaxationKinds = {
   "DampedJacobi", "DampedGaussSeidel"
}
local function isSolverTypeGood(typeIn)
   idxFound = lume.find(solverTypes, typeIn)
   if idxFound then return true end
   return false
end

local function isRelaxKindGood(rlx)
   if lume.find(relaxationKinds, rlx) then
      return true
   end
   return false
end

local function isRelaxDamped(rlx)
   if lume.find(dampedRelaxationKinds, rlx) then
      return true
   end
   return false
end

local function createField(grid, basis, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {1, 1},
   }
   return fld
end

function MGpoisson:init(tbl)
   MGpoisson.super.init(self, tbl) -- Setup base object.

   local solverType = assert( 
      tbl.solverType, "Updater.MGpoisson: Must specify 'DG' or 'FEM' using 'solverType'.")
   assert(isSolverTypeGood(solverType), "Updater.MGpoisson: solverType must be one of 'DG' or 'FEM'.")
   if solverType == "DG" then
      self.isDG  = true
      self.isFEM = false
   elseif solverType == "FEM" then
      self.isDG  = false
      self.isFEM = true
   end

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
   if isRelaxKindGood(relaxKind) then
      if isRelaxDamped(relaxKind) then
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
      elseif strIn == "R" then return BVP_BC_ROBIN
      else
         assert(false, "Updater.MGpoisson: BC type (T) must be one of P, D, N or R. Used " .. strIn .. " instead.")
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
      for j = 1,2 do
         -- Ensure that for Robin BCs three values are given.
         if ( (bcTypes[i][j] == BVP_BC_ROBIN) and ((type(bcValues[i][j]) == "number") or 
              ((type(bcValues[i][j]) == "table") and (#bcValues[i][j] ~= 3))) ) then
            assert(false, "Updater.MGpoisson: for Robin BC please pass a table with 3 numbers in order to impose: bcValue1*f+bcValue2*df/dx=bcValue3.")
         end
         -- Turn BC value into a table if not Robin BC so they all are treated as tables below.
         if (type(bcValues[i][j]) == "number") then
            bcValues[i][j] = {bcValues[i][j]}
         end
      end
   end
   -- Establish periodic directions.
   local periodicDirs  = {}
   local isDirPeriodic = {}
   for d = 1, self.dim do
      if ((bcTypes[d][1] == 0) and (bcTypes[d][2] == 0))then
         lume.push(periodicDirs,d)
         isDirPeriodic[d] = true
         bcValues[d]      = {{0.0}, {0.0}}   -- Not used, but a nil could cause problems.
      elseif ( ((bcTypes[d][1] == 0) and (bcTypes[d][2] ~= 0)) or 
               ((bcTypes[d][1] ~= 0) and (bcTypes[d][2] == 0)) ) then
         assert(false, "Updater.MGpoisson: lower an upper BCs must both be T=\"P\" if periodic is desired.")
      else
         isDirPeriodic[d] = false
      end
   end
   local isPeriodicDomain = lume.all(isDirPeriodic)

   -- Translate bcValues to a vector from which we can pass a pointer.
   -- This vector has 3 values for each boundary in order to support a Robin BC like:
   --     bcValue1*f+bcValue2*df/dx=bcValue3.
   self.bcValue = Lin.Vec(self.dim*2*3)
   for d = 1,self.dim do 
      off = (d-1)*6
      for i = 1, 6 do
         self.bcValue[off+i] = 1.0
      end
      self.bcValue[off + 3 - (bcTypes[d][1] % 3)] = 0.0
      self.bcValue[off + 3]                       = bcValues[d][1][math.floor(1./3.)*2+1]
      if bcTypes[d][1] == BVP_BC_ROBIN then
         -- Robin BCs. First two values multiply boundary value and derivative, respectively.
         self.bcValue[off + 1] = bcValues[d][1][1]
         self.bcValue[off + 2] = bcValues[d][1][2]
      end
      
      self.bcValue[off + 6 - (bcTypes[d][2] % 3)] = 0.0
      self.bcValue[off + 6]                       = bcValues[d][2][math.floor(1./3.)*2+1]
      if bcTypes[d][2] == BVP_BC_ROBIN then
         -- Robin BCs. First two values multiply boundary value and derivative, respectively.
         self.bcValue[off + 4] = bcValues[d][2][1]
         self.bcValue[off + 5] = bcValues[d][2][2]
      end
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
   -- the residue on each grid level. For DG don't need to allocate space for
   -- phi and rho on the top grid because we'll just use the fields
   -- given to the solver.
   self.phiAll[1]     = nil 
   self.rhoAll[1]     = nil
   self.residueAll[1] = createField(self.mgGrids[1],basis)
   for i = 2, self.mgLevels do
      -- Allocate space for the iterate, right-side source field and
      -- the residue on each grid level.
      self.phiAll[i]     = createField(self.mgGrids[i],basis)
      self.rhoAll[i]     = createField(self.mgGrids[i],basis)
      self.residueAll[i] = createField(self.mgGrids[i],basis)
   end
   if self.isJacobiRelax then
      self.phiPrevAll = {}
      for i = 1, self.mgLevels do
         -- For Jacobi relaxation need an extrac copy of the field iterate.
         self.phiPrevAll[i] = createField(self.mgGrids[i],basis)
      end
   end

   if self.isDG then
-- temporary if-statement.
      -- Select restriction and prolongation operator kernels.
      self._restriction  = MGpoissonDecl.selectRestriction(solverType, basisID, self.dim, polyOrder)
      self._prolongation = MGpoissonDecl.selectProlongation(solverType, basisID, self.dim, polyOrder)
      -- Select kernels for relaxation and computing the residue.
      self._relaxation   = MGpoissonDecl.selectRelaxation(solverType, basisID, self.dim, polyOrder, relaxKind, bcTypes)
      self._calcResidue  = MGpoissonDecl.selectResidueCalc(solverType, basisID, self.dim, polyOrder, bcTypes)
   elseif self.isFEM then
      self._relaxation   = MGpoissonDecl.selectRelaxation(solverType, basisID, self.dim, polyOrder, relaxKind, bcTypes)
      self._prolongation = MGpoissonDecl.selectProlongation(solverType, basisID, self.dim, polyOrder)
   end

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

   -- We'll only need one coarse-grid index.
   self.coarseGridIdx = Lin.IntVec(self.dim)

   -- Stencil info: restrict ourselves to nearest neighbor stencils.
   self.phiStencilWidth = 3
   if self.isDG then
      -- 'Cross' stencils for DG (see opStencilIndices).
      self.phiStencilSize = (self.phiStencilWidth-1)*self.dim+1
      -- DG only needs the source in the current cell.
      self.rhoStencilSize = 1
   elseif self.isFEM then
      -- 'Filled' stencils for FEM (see opStencilIndices).
      self.phiStencilSize = self.phiStencilWidth^self.dim
      -- FEM uses the right-side source in neighboring cells as well.
      self.rhoStencilSize = self.phiStencilSize
   end
   -- List of cell indices pointed to by the stencils.
   self.phiStencilIdx = {}
   self.rhoStencilIdx = {}
   for i = 1, self.phiStencilSize do self.phiStencilIdx[i] = Lin.IntVec(self.dim) end
   for i = 1, self.rhoStencilSize do self.rhoStencilIdx[i] = Lin.IntVec(self.dim) end
   -- List of pointers to the data in cells pointed to by the stencil.
   self.phiStencil     = DoublePtrVec(self.phiStencilSize)
   self.prevPhiStencil = DoublePtrVec(self.phiStencilSize)   -- Only used for Jacobi relaxations.
   self.rhoStencil     = DoublePtrVec(self.rhoStencilSize)
   self.dxStencil      = DoublePtrVec(self.phiStencilSize)

   self.dxBuf = Lin.Vec(self.dim)    -- A buffer to store cell lengths.

   -- Dimensions remaining when a dimension is removed.
   self.dimRemain = {}
   for d1 = 1, self.dim do
      self.dimRemain[d1] = {}
      for d2 = 1, self.dim do self.dimRemain[d1][d2] = d2 end
      table.remove(self.dimRemain[d1],d1)
   end

   -- ..................... Things specific to the FEM solver ......................... --

   if self.isFEM then
      -- For FEM we do need fields at the finest grid because of the FEM-DG translations.
      self.phiAll[1] = createField(self.mgGrids[1],basis)
      self.rhoAll[1] = createField(self.mgGrids[1],basis)
   end

   -- For FEM solver, will need to translate (modal) DG coefficients into
   -- (nodal) FEM coefficients. Preselect the appropriate kernels here.
   self._dgToFEM = MGpoissonDecl.selectDGtoFEM(basisID, self.dim, polyOrder, bcTypes)

   self.dgToFEMstencilWidth = 2
   self.dgToFEMstencilSize  = self.dgToFEMstencilWidth^self.dim
   self.dgToFEMstencilIdx   = {}   -- List of cell indices pointed to by the stencil.
   for i = 1, self.dgToFEMstencilSize do
      self.dgToFEMstencilIdx[i] = Lin.IntVec(self.dim)
   end
   -- List of pointers to the data in cells pointed to by the stencil.
   self.dgToFEMstencilItr = DoublePtrVec(self.dgToFEMstencilSize)

   self.zeros  = {}
   self.threes = {}
   self.mOnes  = {}
   for d = 1, self.dim do
      self.zeros[d]  = 0
      self.threes[d] = 3
      self.mOnes[d]  = -1
   end

   -- ......................... End of FEM-specific things ............................ --

   -- Select MG components for FEM or DG solver.
   if self.isDG then
      self.relax   = function(numRelax, phiFld, rhoFld) MGpoisson['relaxDG'](self, numRelax, phiFld, rhoFld) end
      self.prolong = function(cFld,fFld) MGpoisson['prolongDG'](self,cFld,fFld) end
   else
      self.relax   = function(numRelax, phiFld, rhoFld) MGpoisson['relaxFEM'](self, numRelax, phiFld, rhoFld) end
      self.prolong = function(cFld,fFld) MGpoisson['prolongFEM'](self,cFld,fFld) end
   end

   -- Updater to compute the L2-norm of the residue.
   self.l2NormCalc = IntQuantCalc {
      onGrid   = grid, 
      basis    = basis,
      quantity = "RmsV",
   }
   self.relResNorm  = DataStruct.DynVector { numComponents = 1, }
   self.residueNorm = DataStruct.DynVector { numComponents = 1, }
   self.rhoNorm     = DataStruct.DynVector { numComponents = 1, }

end

function MGpoisson:opStencilIndices(idxIn, stencilType, stencilIdx)
   -----------------
   -- Given the index of the current cell (idxIn), return a table
   -- of indices of the cells in a stencil of type 'stencilType'
   -- Stencil types are given by [t,w,l]:
   --   t: type of stencil
   --        0 : 'cross'.
   --        1 : partially filled cross (include nearest corner cells).
   --        2 : filled (include all corner cells).
   --   w: width of the stencil (along each dimension).
   --   l: location in grid along each dimension. For 1D
   --        [-1] = lower boundary cell.
   --        [ 0] = inner cell.
   --        [ 1] = upper boundary cell.
   --      Note that depending on how the boundary kernels are codded
   --      it may be ok to pass an inner-cell stencil at the boundary.
   --
   -- Examples:
   -- a) stencilType=[0,[5,5],[0,0]] corresponds to 
   --                 9
   --                 7
   --         4   2   1   3   5
   --                 6
   --                 8
   -- b) stencilType=[1,[3,3],[0,0]] or [2,[3,3],[0,0]] correspond to 
   --             7   5   9
   --             2   1   3    
   --             6   4   8

   -- First copy all indicies since (in higher dimensions)
   -- most of the stay the same for each cell.
   for i = 1, #stencilIdx do
      idxIn:copyInto(stencilIdx[i])
   end

   if stencilType[1] == 0 then
      local sI = 1
      for d = 1, self.dim do
         for pm = 1,self.phiStencilWidth-1 do
            sI = sI + 1
            stencilIdx[sI][d] = idxIn[d]+((-1)^(pm % 2))*((self.phiStencilWidth-1)/2)
         end
      end
   elseif stencilType[1]==2 then
      local sI = 1
      for d = 1, self.dim do
         local prevDimCells = sI
         for pDC = 1, prevDimCells do
            for pm = 1-stencilType[3][d],self.phiStencilWidth-1 do
               sI = sI + 1
               for _, dr in ipairs(self.dimRemain[d]) do stencilIdx[sI][dr] = stencilIdx[pDC][dr] end
               stencilIdx[sI][d] = stencilIdx[pDC][d]+((-1)^(pm % 2))*((self.phiStencilWidth-1)/2)
            end
         end
      end
   end

end

function MGpoisson:idxToStencil(idxIn, nCellsIn)
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

-- ................................... Functions specific to FEM solver ................................... --

function MGpoisson:DG_FEM_coefTranslate(dgFld,femFld,dir)
  -- Translate the DG coefficients of a field into FEM expansion
  -- coefficients (dir=-1), and viceversa (dir=1).

   local grid   = dgFld:grid()
   local cellsN = {}
   for d = 1, self.dim do cellsN[d]=grid:numCells(d) end

   localRangeDecomp = LinearDecomp.LinearDecompRange {
      range = dgFld:localRange(), numSplit = grid:numSharedProcs() }
   local tId        = grid:subGridSharedId()    -- Local thread ID.

   local indexer    = dgFld:genIndexer()

   local dgFldItr   = dgFld:get(1)
   local femFldItr  = femFld:get(1)

   for idx in localRangeDecomp:rowMajorIter(tId) do

      grid:setIndex(idx)

      dgFld:fill(indexer(idx), dgFldItr)     -- DG field pointer.
      femFld:fill(indexer(idx), femFldItr)   -- FEM field pointer.
 
      -- Get with indices of cells used by stencil. Store them in self.phiStencilIdx.
      self:opStencilIndices(idx,{2,self.threes,self.mOnes},self.dgToFEMstencilIdx)
 
      -- Array of pointers to cell lengths and phi data in cells pointed to by the stencil.
      for i = 1, self.dgToFEMstencilSize do
         grid:setIndex(self.dgToFEMstencilIdx[i])
 
         femFld:fill(indexer(self.dgToFEMstencilIdx[i]), femFldItr)
         self.dgToFEMstencilItr[i] = femFldItr:data()
      end
 
      self._dgToFEM[self:idxToStencil(idx,cellsN)](dgFldItr:data(), self.dgToFEMstencilItr:data())
   end
end

-- ..................................... End of FEM solver functions ...................................... --

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

function MGpoisson:relaxDG(numRelax, phiFld, rhoFld)
   -- Perform numRelax DG relaxations of the Poisson equation.

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
         self.dxStencil[1] = self.dxBuf:data()
         rhoFld:fill(indexer(idx), rhoItr)   
     
         -- Get with indices of cells used by stencil. Store them in self.phiStencilIdx.
         self:opStencilIndices(idx,{0,self.threes,self.zeros},self.phiStencilIdx)
   
         -- Array of pointers to cell lengths and phi data in cells pointed to by the stencil. 
         for i = 1, self.phiStencilSize do
            grid:setIndex(self.phiStencilIdx[i])
            grid:getDx(self.dxBuf)
            self.dxStencil[i] = self.dxBuf:data()

            phiFld:fill(indexer(self.phiStencilIdx[i]), phiItr)
            self.phiStencil[i] = phiItr:data()

            if self.isJacobiRelax then
               self.phiPrevAll[currLevel]:fill(indexer(self.phiStencilIdx[i]), phiPrevItr)
               self.prevPhiStencil[i] = phiPrevItr:data()
            end
         end
         
         self._relaxation[self:idxToStencil(idx,cellsN)](self.omega, self.dxStencil:data(), self.bcValue:data(), rhoItr:data(), self.prevPhiStencil:data(), self.phiStencil:data())
      end
   end
end

function MGpoisson:relaxFEM(numRelax, phiFld, rhoFld)
   -- Perform numRelax FEM relaxations of the Poisson equation.

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
         self.dxStencil[1] = self.dxBuf:data()
         rhoFld:fill(indexer(idx), rhoItr)   
     
         -- Get with indices of cells used by stencil. Store them in self.phiStencilIdx.
         self:opStencilIndices(idx,{0,self.threes,self.zeros},self.phiStencilIdx)
   
         -- Array of pointers to cell lengths and phi data in cells pointed to by the stencil. 
         for i = 1, self.phiStencilSize do
            grid:setIndex(self.phiStencilIdx[i])
            grid:getDx(self.dxBuf)
            self.dxStencil[i] = self.dxBuf:data()

            phiFld:fill(indexer(self.phiStencilIdx[i]), phiItr)
            self.phiStencil[i] = phiItr:data()

            if self.isJacobiRelax then
               self.phiPrevAll[currLevel]:fill(indexer(self.phiStencilIdx[i]), phiPrevItr)
               self.prevPhiStencil[i] = phiPrevItr:data()
            end
         end
         
         -- Array of pointers to rho data in cells pointed to by the stencil. 
         for i = 1, self.rhoStencilSize do
            grid:setIndex(self.rhoStencilIdx[i])
            rhoFld:fill(indexer(self.rhoStencilIdx[i]), rhoItr)
            self.rhoStencil[i] = rhoItr:data()
         end
         
         self._relaxation[self:idxToStencil(idx,cellsN)](self.omega, self.dxStencil:data(), self.bcValue:data(), self.rhoStencil:data(), self.prevPhiStencil:data(), self.phiStencil:data())
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
      self.dxStencil[1] = self.dxBuf:data()
      rhoFld:fill(indexer(idx), rhoItr)   
      resFld:fill(indexer(idx), resItr)   
   
      -- Get with indices of cells used by stencil.
      self:opStencilIndices(idx,{0,self.threes,self.zeros},self.phiStencilIdx)
   
      -- Array of pointers to cell lengths and phi data in cells pointed to by the stencil. 
      for i = 1, self.phiStencilSize do
         grid:setIndex(self.phiStencilIdx[i])
         grid:getDx(self.dxBuf)
         self.dxStencil[i] = self.dxBuf:data()
         phiFld:fill(indexer(self.phiStencilIdx[i]), phiItr)
         self.phiStencil[i] = phiItr:data()
      end
   
      self._calcResidue[self:idxToStencil(idx,cellsN)](self.dxStencil:data(), self.bcValue:data(), rhoItr:data(), self.phiStencil:data(), resItr:data())
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

function MGpoisson:prolongDG(cFld,fFld)
   -- DG prolongation of a coarse-grid field (cFld) to a fine-grid field (fFld). 

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

function MGpoisson:prolongFEM(cFld,fFld)
   -- FEM prolongation of a coarse-grid field (cFld) to a fine-grid field (fFld). 

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

      -- Array of indexes to fine-grid cells this coarse-grid cell contributes to.
      for d = 1, self.dim do self.fineGridIdx[1][d] = 2*cIdx[d] end
      local parIdx     = 1
      local fCellCount = 1
      for dir = 1, self.dim do
         for rI = 1, self.igOpStencilWidth^(dir-1) do
            for k = 1, (self.igOpStencilWidth-1) do
               local newIdxInDir = self.fineGridIdx[parIdx][dir]-k
               if newIdxInDir<1 then break end
               fCellCount = fCellCount + 1
               for d = 1, self.dim do self.fineGridIdx[fCellCount][d] = self.fineGridIdx[parIdx][d] end
               self.fineGridIdx[fCellCount][dir] = newIdxInDir
            end
            parIdx = parIdx + self.igOpStencilWidth-1
         end
      end

      cFld:fill(cFldIndexer(cIdx), cFldItr)   -- Coarse field pointer.
  
      -- Array of pointers to fine-grid field data by stencil. 
      for i = 1, self.igOpStencilSize do
         fFld:fill(fFldIndexer(self.fineGridIdx[i]), fFldItr)
         self.fineFldItr[i] = fFldItr:data()
      end
  
      self._prolongation(cFldItr:data(), self.fineFldItr:data())
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
         self.relax(self.nu3, self.phiAll[lCurr], self.rhoAll[lCurr]) 
      end

   else

      -- Relax nu1 times.
      self.relax(self.nu1, self.phiAll[lCurr], self.rhoAll[lCurr]) 

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
      self.prolong(self.phiAll[lCurr+1], self.residueAll[lCurr])
      self.phiAll[lCurr]:accumulate(1.0,self.residueAll[lCurr])

      -- Relax nu2 times.
      self.relax(self.nu2, self.phiAll[lCurr], self.rhoAll[lCurr]) 

   end
   
end

-- Advance method.
function MGpoisson:_advance(tCurr, inFld, outFld)

   -- Phi iterate and right-side source fields on finest grid.
   -- Assume the given phi iterate is an initial guess if given
   -- within inFld. If not it will interpret that to mean that the
   -- user wishes to perform a full-multigrid (FMG) cycle.
   if self.isDG then
      self.rhoAll[1] = inFld[1]
   elseif self.isFEM then
      -- FEM solver. Translate RHS source DG coefficients to FEM.
      self.rhoAll[1]:clear(0.0)
      self:DG_FEM_coefTranslate(inFld[1], self.rhoAll[1],-1)
   end
   local initialGuess   = inFld[2]
   local relResNormCurr = 1.0e12    -- Current (relative) residue norm.
   if initialGuess then
      if self.isDG then
         self.phiAll[1] = initialGuess
      elseif self.isFEM then
         -- FEM solver. Translate initial guess DG coefficients to FEM.
         self.phiAll[1]:clear(0.0)
         self:DG_FEM_coefTranslate(initialGuess, self.phiAll[1],DG_to_FEM)
      end
   else
      -- No initial guess provided. Perform Full Multi-Grid (FMG).
      if self.isDG then self.phiAll[1] = outFld[1] end

      -- FMG requires we restrict the right-side source field to all levels.
      for i = 2, self.mgLevels do
         self:restrict(self.rhoAll[i-1], self.rhoAll[i])
      end

      -- In FMG we start at the coarsest grid without an initial guess.
      self.phiAll[self.mgLevels]:clear(0.0)
      for i = self.mgLevels, 1, -1 do
         self:gammaCycle(i)
         if (i > 1) then self.prolong(self.phiAll[i], self.phiAll[i-1]) end
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

--   if self.isFEM then
--      -- Translate final phi from FEM to DG.
--      outFld[1]:clear(0.0)
--      self:DG_FEM_coefTranslate(self.phiAll[1],outFld[1],FEM_to_DG)
--   end

end

return MGpoisson
