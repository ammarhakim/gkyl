-- Gkyl ------------------------------------------------------------------------
--
-- App support code: Maxwell field objects.
--
-- Contains functions needed by Maxwell equations' solvers, and the
-- funcField capability to impose external fields.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo  = require "Io.AdiosCartFieldIo"
local DataStruct        = require "DataStruct"
local FieldBase         = require "App.Field.FieldBase"
local LinearTrigger     = require "LinearTrigger"
local Mpi               = require "Comm.Mpi"
local PerfMaxwell       = require "Eq.PerfMaxwell"
local Proto             = require "Lib.Proto"
local Time              = require "Lib.Time"
local Updater           = require "Updater"
local BoundaryCondition = require "Updater.BoundaryCondition"
local xsys              = require "xsys"
local ffi               = require "ffi"
local lume              = require "Lib.lume"
local Lin               = require "Lib.Linalg"
local MGpoissonDecl     = require "Updater.mgPoissonCalcData.MGpoissonModDecl"
local LinearDecomp      = require "Lib.LinearDecomp"

-- MaxwellField ---------------------------------------------------------------------
--
-- Faraday's and Ampere's equations are evolved (electromagnetic), or the
-- Poisson equation is solved (electrostatic).
-------------------------------------------------------------------------------------

local MaxwellField = Proto(FieldBase.FieldBase)

-- Add constants to object indicate various supported boundary conditions.
local EM_BC_REFLECT = 1
local EM_BC_SYMMETRY = 2
local EM_BC_COPY = 3
-- AHH: This was 2 but seems that is unstable. So using plain copy.
local EM_BC_OPEN = EM_BC_COPY
local EM_BC_AXIS = 7

MaxwellField.bcOpen = EM_BC_OPEN    -- Zero gradient.
MaxwellField.bcCopy = EM_BC_COPY    -- Copy fields.
MaxwellField.bcReflect = EM_BC_REFLECT -- Perfect electric conductor.
MaxwellField.bcSymmetry = EM_BC_SYMMETRY
MaxwellField.bcAxis = EM_BC_AXIS

-- Function to check if BC type is good.
local function isBcGood(bcType)
   if bcType == EM_BC_OPEN or bcType == EM_BC_REFLECT or bcType == EM_BC_SYMMETRY then
      return true
   end
   if type(bcType) == "table" then
     return true
   end
   return false
end

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function MaxwellField:init(tbl)
   MaxwellField.super.init(self, tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function MaxwellField:fullInit(appTbl)
   local tbl = self.tbl -- Previously store table.
   
   self.epsilon0 = tbl.epsilon0
   self.mu0 = tbl.mu0
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- By default evolve field.
   -- By default, do not write data if evolve is false.
   self.forceWrite = xsys.pickBool(tbl.forceWrite, false)

   -- By default there is a magnetic field. Setting it to false runs Vlasov-Poisson.
   self.hasMagField = xsys.pickBool(tbl.hasMagneticField, true)

   if self.hasMagField then   -- Things not needed for Poisson solve.
      self.ce = tbl.elcErrorSpeedFactor and tbl.elcErrorSpeedFactor or 0.0
      self.cb = tbl.mgnErrorSpeedFactor and tbl.mgnErrorSpeedFactor or 1.0
   
      self.lightSpeed = 1/math.sqrt(self.epsilon0*self.mu0)
   
      -- tau parameter used for adding extra (less) diffusion to
      -- Ampere-Maxwell, while adding less (more) diffusion to Faraday
      -- equation if no tau parameter is specified, Eq object defaults to
      -- the speed of light.
      self.tau = tbl.tau

      self._hasSsBnd = tbl.hasSsBnd
      self._inOutFunc = tbl.inOutFunc

      -- No ghost current by default.
      self.useGhostCurrent = xsys.pickBool(tbl.useGhostCurrent, false)

      self.limiter = self.tbl.limiter and self.tbl.limiter or "monotonized-centered"

      -- numFlux used for selecting which type of numerical flux function to use
      -- defaults to "upwind" in Eq object, supported options: "central," "upwind"
      -- only used for DG Maxwell.
      self.numFlux = tbl.numFlux

      -- Store initial condition function (this is a wrapper around user
      -- supplied function as we need to add correction potential ICs here).
      self.initFunc = function (t, xn)
         local ex, ey, ez, bx, by, bz = tbl.init(t, xn)
         return ex, ey, ez, bx, by, bz, 0.0, 0.0
      end

   else

      -- Read in boundary conditions for the potential phi.
      -- Ideally this is not done through a separate infrastructue to bcx, bcy, bcz below.
      -- But changing those needs a little more work/careful thought. For now we reproduce
      -- the infrastructure in GkField.
      local ndim, periodicDirs, isDirPeriodic, periodicDomain = #appTbl.lower, nil, {}, true
      if appTbl.periodicDirs then periodicDirs = appTbl.periodicDirs else periodicDirs = {} end
      for d = 1, ndim do isDirPeriodic[d] = lume.find(periodicDirs,d) ~= nil end
      self.bcLowerPhi, self.bcUpperPhi = {}, {}
      for d=1,ndim do periodicDomain = periodicDomain and isDirPeriodic[d] end
      if tbl.bcLowerPhi==nil and tbl.bcUpperPhi==nil then
         if periodicDomain then
            for d=1,ndim do
               self.bcLowerPhi[d] = {T="P"}
               self.bcUpperPhi[d] = self.bcLowerPhi[d]
            end
         else
            assert(false, "App.Field.MaxwellField: must specify 'bcLower' and 'bcUpper.")
         end
      else
         assert(#tbl.bcLowerPhi==#tbl.bcUpperPhi, "App.Field.MaxwellField: number of entries in bcLower and bcUpper must be equal.")
         assert(#tbl.bcLowerPhi==ndim, "App.Field.MaxwellField: number of entries in bcLower/bcUpper must equal number of dimensions.")
         for d=1,ndim do
            assert((isDirPeriodic[d]==(tbl.bcLowerPhi[d].T=="P")) and (isDirPeriodic[d]==(tbl.bcUpperPhi[d].T=="P")),
                   string.format("App.Field.MaxwellField: direction %d is periodic. Must use {T='P'} in bcLower/bcUpper.",d))
            self.bcLowerPhi[d], self.bcUpperPhi[d] = tbl.bcLowerPhi[d], tbl.bcUpperPhi[d]
         end
      end

      self.esEnergyFirst = true

   end

   -- Create triggers to write fields.
   if tbl.nFrame then
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nFrame)
   else
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.ioFrame = 0 -- Frame number for IO.
   self.dynVecRestartFrame = 0 -- Frame number of restarts (for DynVectors only).

   self.hasNonPeriodicBc = false -- To indicate if we have non-periodic BCs.
   self.bcx, self.bcy, self.bcz = { }, { }, { }
   
   -- Read in boundary conditions.
   if tbl.bcx then
      self.bcx[1], self.bcx[2] = tbl.bcx[1], tbl.bcx[2]
      assert(isBcGood(self.bcx[1]) and isBcGood(self.bcx[2]), "MaxwellField: Incorrect X BC type specified!")
      self.hasNonPeriodicBc = true
   end
   if tbl.bcy then
      self.bcy[1], self.bcy[2] = tbl.bcy[1], tbl.bcy[2]
      assert(isBcGood(self.bcy[1]) and isBcGood(self.bcy[2]), "MaxwellField: Incorrect Y BC type specified!")
      self.hasNonPeriodicBc = true
   end
   if tbl.bcz then
      self.bcz[1], self.bcz[2] = tbl.bcz[1], tbl.bcz[2]
      assert(isBcGood(self.bcz[1]) and isBcGood(self.bcz[2]), "MaxwellField: Incorrect Z BC type specified!")
      self.hasNonPeriodicBc = true
   end

   self.ssBc = {}
   if tbl.ssBc then
      self.ssBc[1] = tbl.ssBc[1]
   end

   self.tmCurrentAccum = 0.0 -- Time spent in current accumulate.
   self.integratedEMTime = 0.0 -- Time spent integrating EM fields.

   -- Create trigger for how frequently to compute integrated EM fields.
   if appTbl.calcIntQuantEvery then
      self.calcIntEMQuantTrigger = LinearTrigger(0, appTbl.tEnd,  math.floor(1/appTbl.calcIntQuantEvery))
   else
      self.calcIntEMQuantTrigger = function(t) return true end
   end
end

-- Methods for EM field object.
function MaxwellField:hasEB() return true, self.hasMagField end
function MaxwellField:setCfl(cfl) self.cfl = cfl end
function MaxwellField:getCfl() return self.cfl end
function MaxwellField:setGrid(grid) self.grid = grid end

function MaxwellField:getEpsilon0() return self.epsilon0 end
function MaxwellField:getMu0() return self.mu0 end
function MaxwellField:getElcErrorSpeedFactor() return self.ce end
function MaxwellField:getMgnErrorSpeedFactor() return self.cb end

function MaxwellField:esEnergy(tCurr, fldIn, outDynV)
   -- Compute the electrostatic field energy given the potential. Here outDynV must
   -- be a DynVector with the same number of components as there are dimensions.
   local phiIn, esE = fldIn[1], outDynV[1]
   if GKYL_USE_GPU then
      phiIn:copyDeviceToHost()
   end

   local grid, dim = phiIn:grid(), phiIn:grid():ndim()

   if self.esEnergyFirst then
      -- Kernels and buffers used in computing electrostatic field energy.
      self.localEnergy, self.globalEnergy = Lin.Vec(dim), Lin.Vec(dim)
      self.dxBuf = Lin.Vec(dim)
      self._esEnergyCalc = MGpoissonDecl.selectESenergy(self.basis:id(), dim, self.basis:polyOrder(), nil)
   end

   local indexer, phiItr = phiIn:genIndexer(), phiIn:get(1)

   -- Clear local values.
   for d = 1, dim do self.localEnergy[d] = 0.0 end

   local phiRange = phiIn:localRange()

   for idx in phiRange:rowMajorIter() do
      grid:setIndex(idx)
      grid:getDx(self.dxBuf)

      phiIn:fill(indexer(idx), phiItr)

      self._esEnergyCalc(self.dxBuf:data(), phiItr:data(), self.localEnergy:data())
   end

   -- All-reduce across processors and push result into dyn-vector.
   Mpi.Allreduce(
      self.localEnergy:data(), self.globalEnergy:data(), dim, Mpi.DOUBLE, Mpi.SUM, grid:commSet().comm)

   esE:appendData(tCurr, self.globalEnergy)

   if self.esEnergyFirst then self.esEnergyFirst = false end
end

function MaxwellField:alloc(nRkDup)
   local nGhost = 2
   if self.basis:numBasis() > 1 then
      nGhost = 1
   end
   
   self.em = {}
   if self.hasMagField then   -- Maxwell's induction equations.
      -- Allocate fields needed in RK update.
      for i = 1, nRkDup do
         self.em[i] = DataStruct.Field {
            onGrid = self.grid,
            numComponents = 8*self.basis:numBasis(),
            ghost = {nGhost, nGhost}
         }
      end

      -- Array with one component per cell to store cflRate in each cell.
      self.cflRateByCell = DataStruct.Field {
         onGrid        = self.grid,
         numComponents = 1,
         ghost         = {1, 1},
      }
      self.cflRateByCell:clear(0.0)
      self.cflRatePtr = self.cflRateByCell:get(1)
      self.cflRateIdxr = self.cflRateByCell:genIndexer()
      self.dtGlobal = ffi.new("double[2]")
      
      -- For storing integrated energy components.
      self.emEnergy = DataStruct.DynVector { numComponents = 8 }

   else   -- Poisson equation.
      -- Electrostatic potential, phi.
      local phiFld = DataStruct.Field {
         onGrid = self.grid,
         numComponents = self.basis:numBasis(),
         ghost = {1, 1}
      }
      for i = 1, nRkDup do self.em[i] = phiFld end
      -- Keep copies of previous potentials so we can use three point
      -- extrapolation to form an initial guess for the MG solver.
      self.phiPrevNum = 3
      self.phiFldPrev = {}
      for i = 1, self.phiPrevNum do
         self.phiFldPrev[i] = {}
         self.phiFldPrev[i]["time"] = 0.0
         self.phiFldPrev[i]["fld"] = DataStruct.Field {
            onGrid = self.grid,
            numComponents = self.basis:numBasis(),
            ghost = {1, 1}
         }
      end

      -- For storing integrated energy components.
      self.emEnergy = DataStruct.DynVector { numComponents = self.grid:ndim() }
   end
end

function MaxwellField:createSolver(population)

   -- Create Adios object for field I/O.
   self.fieldIo = AdiosCartFieldIo {
      elemType = self.em[1]:elemType(),
      method   = self.ioMethod,
      metaData = {polyOrder = self.basis:polyOrder(),
                  basisType = self.basis:id(),
                  epsilon0  = self.epsilon0,
                  mu0       = self.mu0,
                  grid      = GKYL_OUT_PREFIX .. "_grid.bp"},
      writeRankInComm = {0, population:getComm(),},
   }

   if self.hasMagField then   -- Maxwell's induction equations.

      self.equation = PerfMaxwell {
         lightSpeed          = self.lightSpeed,  tau     = self.tau,
         elcErrorSpeedFactor = self.ce,          numFlux = self.numFlux,
         mgnErrorSpeedFactor = self.cb,          basis   = self.basis:numBasis() > 1 and self.basis or nil,
      }

      self.fieldSlvr, self.fieldHyperSlvr = nil, {}
      if self.basis:numBasis() > 1 then
         -- Using DG scheme.
         self.fieldSlvr = Updater.HyperDisCont {
            onGrid = self.grid,   cfl      = self.cfl,
            basis  = self.basis,  equation = self.equation,
         }
      else
         -- Using FV scheme.
         if self._hasSsBnd then
            self._inOut = DataStruct.Field {
               onGrid = self.grid,
               numComponents = 1,
               ghost = {2, 2}
            }
            local project = Updater.ProjectOnBasis {
               onGrid = self.grid,
               basis = self.basis,
               evaluate = self._inOutFunc,
               onGhosts = true,
            }
            project:advance(0.0, {}, {self._inOut})
            self.fieldIo:write(self._inOut, string.format("%s_inOut.bp", "field"), 0, 0)
         end

         local ndim = self.grid:ndim()
         for d = 1, ndim do
            self.fieldHyperSlvr[d] = Updater.WavePropagation {
               onGrid = self.grid,
               equation = self.equation,
               limiter = self.limiter,
               cfl = self.cfl,
               updateDirections = {d},
               hasSsBnd = self._hasSsBnd,
               inOut = self._inOut
            }
         end
      end

      self.emEnergyUpd = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.grid,   quantity      = "V2",
         basis  = self.basis,  numComponents = 8,
      }
      self.emEnergyCalc = function(tCurr, inFld, outDynV) self.emEnergyUpd:advance(tCurr, inFld, outDynV) end

   else   -- Poisson equation.

      self.isElliptic = true

      local isParallel = false
      for d=1,self.grid:ndim() do if self.grid:cuts(d)>1 then isParallel=true end end

--      MF 2022/08/15: disable multigrid for now.
--      if self.basis:polyOrder()>1 or isParallel or GKYL_USE_GPU then
         self.isSlvrMG = false
         self.fieldSlvr = Updater.FemPoisson {
            onGrid    = self.grid,   bcLower = self.bcLowerPhi,
            basis     = self.basis,  bcUpper = self.bcUpperPhi,
            epsilon_0 = self.epsilon0,
         }
         self.esEnergyUpd = Updater.CartFieldIntegratedQuantCalc {
            onGrid = self.grid,   quantity = "V2",
            basis  = self.basis,
         }
         self.emEnergyCalc = function(tCurr, inFld, outDynV) self:esEnergy(tCurr, inFld, outDynV) end
--      else
--         self.isSlvrMG = true
--         -- Multigrid parameters (hardcoded for now).
--         local gamma = 1                 -- V-cycles=1, W-cycles=2.
--         local relaxType = 'DampedJacobi'    -- DampedJacobi or DampedGaussSeidel
--         local relaxNum = {1,2,300}         -- Number of pre,post and coarsest-grid smoothings.
--         local relaxOmega                     -- Relaxation damping parameter.
--         local ndim = self.grid:ndim()
--         if ndim == 1 then
--            relaxOmega = 2./3.
--         elseif ndim == 2 then
--            relaxOmega = 4./5.
--         end
--         local tolerance = 1.e-6             -- Do cycles until reaching this relative residue norm.
--         -- After the 4th call to the advance method, we will start using a
--         -- 3-point extrapolation to obtain an initial guess for the MG solver. 
--         self.filledPhiPrev = false
--         self.phiPrevCount = 0
--         self.fieldSlvr = Updater.MGpoisson {                                   
--            solverType = 'FEM',                                                 
--            onGrid     = self.grid,
--            basis      = self.basis,
--            bcLower    = self.bcLowerPhi,
--            bcUpper    = self.bcUpperPhi,
--            gamma      = gamma,         -- V-cycles=1, W-cycles=2.
--            relaxType  = relaxType,     -- DampedJacobi or DampedGaussSeidel
--            relaxNum   = relaxNum,      -- Number of pre,post and coarsest-grid smoothings.
--            relaxOmega = relaxOmega,    -- Relaxation damping parameter.
--            tolerance  = tolerance,     -- Do cycles until reaching this relative residue norm.
--         }
--   
--         self.emEnergyCalc = function(tCurr, inFld, outDynV) self.fieldSlvr:esEnergy(tCurr, inFld, outDynV) end
--      end
   end

   -- Function to construct a BC updater.
   local function makeBcUpdater(dir, edge, bcList)
      return Updater.Bc {
         onGrid             = self.grid,
         boundaryConditions = bcList,
         dir                = dir,
         edge               = edge,
      }
   end

   -- Indices for tangent and normal components of E and B for dir.
   local idxEt = {{2, 3}, {1, 3}, {1, 2}}
   local idxEn = {1, 2, 3}
   local idxBt = {{5, 6}, {4, 6}, {4, 5}}
   local idxBn = {4, 5, 6}

   -- Various functions to apply BCs of different types.
   local function bcOpen(dir, tm, xc, fIn, fOut)
      local nb = self.basis:numBasis()
      local fInData, fOutData = fIn:data(), fOut:data()
      for i = 1, 8 do
	 self.basis:flipSign(dir, fInData+(i-1)*nb, fOutData+(i-1)*nb)
      end
   end
   local function bcCopy(dir, tm, xc, fIn, fOut)
      for i = 1, 8*self.basis:numBasis() do fOut[i] = fIn[i] end
   end 
   local function bcReflect(dir, tm, xc, fIn, fOut)
      local nb = self.basis:numBasis()
      local fInData, fOutData = fIn:data(), fOut:data()
      -- Zero gradient for all the components.
      for i = 1, 8 do
	 self.basis:flipSign(dir, fInData+(i-1)*nb, fOutData+(i-1)*nb)
      end
      for i = 1, self.basis:numBasis() do
	 -- Zero tangent for electric field.
	 fOutData[(idxEt[dir][1]-1)*nb + i - 1] = 
	    -1.0 * fOutData[(idxEt[dir][1]-1)*nb + i - 1]
	 fOutData[(idxEt[dir][2]-1)*nb + i - 1] = 
	    -1.0 * fOutData[(idxEt[dir][2]-1)*nb + i - 1]
	 -- Zero normal for magnetic field.
	 fOutData[(idxBn[dir]-1)*nb + i - 1] = 
	    -1.0 * fOutData[(idxBn[dir]-1)*nb + i - 1]
      end
   end
   local function bcSymmetry(dir, tm, xc, fIn, fOut)
      local nb = self.basis:numBasis()
      local fInData, fOutData = fIn:data(), fOut:data()
      -- Zero gradient for all the components.
      for i = 1, 8 do
	 self.basis:flipSign(dir, fInData+(i-1)*nb, fOutData+(i-1)*nb)
      end
      for i = 1, self.basis:numBasis() do
	 -- Zero normal for electric field.
	 fOutData[(idxEn[dir]-1)*nb + i - 1] = 
	    -1.0 * fOutData[(idxEn[dir]-1)*nb + i - 1]
	 -- Zero tangent for magnetic field.
	 fOutData[(idxBt[dir][1]-1)*nb + i - 1] = 
	    -1.0 * fOutData[(idxBt[dir][1]-1)*nb + i - 1]
	 fOutData[(idxBt[dir][2]-1)*nb + i - 1] = 
	    -1.0 * fOutData[(idxBt[dir][2]-1)*nb + i - 1]
      end
   end

   -- Functions to make life easier while reading in BCs to apply.
   self.boundaryConditions = { } -- List of Bcs to apply.
   local function appendBoundaryConditions(dir, edge, bcType)
      if bcType == EM_BC_OPEN then
	 table.insert(self.boundaryConditions,
		      makeBcUpdater(dir, edge, { bcCopy }))
      elseif bcType == EM_BC_COPY then
	 table.insert(self.boundaryConditions,
		      makeBcUpdater(dir, edge, { bcCopy }))
      elseif bcType == EM_BC_REFLECT then
	 table.insert(self.boundaryConditions,
		      makeBcUpdater(dir, edge, { bcReflect }))
      elseif bcType == EM_BC_SYMMETRY then
	 table.insert(self.boundaryConditions,
		      makeBcUpdater(dir, edge, { bcSymmetry }))
      elseif bcType == EM_BC_AXIS then
	 table.insert(self.boundaryConditions,
		      makeBcUpdater(dir, edge,  PerfMaxwell.bcAxis ))
      elseif type(bcType) == "table" then
	 table.insert(self.boundaryConditions,
		      makeBcUpdater(dir, edge, bcType))
      else
	 assert(false, "MaxwellField: Unsupported BC type!")
      end
   end

   local function handleBc(dir, bc)
      if bc[1] then
	 appendBoundaryConditions(dir, "lower", bc[1])
      end
      if bc[2] then
	 appendBoundaryConditions(dir, "upper", bc[2])
      end
   end
   
   -- Add various BCs to list of BCs to apply.
   handleBc(1, self.bcx)
   handleBc(2, self.bcy)
   handleBc(3, self.bcz)

   self.ssBoundaryConditions = { }
   if self._hasSsBnd then
      function makeSsBcUpdater(dir, inOut, bcList)
         return Updater.StairSteppedBc {
            onGrid = self.grid,
            inOut = inOut,
            boundaryConditions = bcList,
            dir = dir,
         }
      end
      local function appendSsBoundaryConditions(dir, inOut, bcType)
         if bcType == EM_BC_OPEN then
            table.insert(self.ssBoundaryConditions,
                         makeSsBcUpdater(dir, inOut, { bcCopy }))
         elseif bcType == EM_BC_COPY then
            table.insert(self.ssBoundaryConditions,
                         makeSsBcUpdater(dir, inOut, { bcCopy }))
         elseif bcType == EM_BC_REFLECT then
            table.insert(self.ssBoundaryConditions,
                         makeSsBcUpdater(dir, inOut, { bcReflect }))
         elseif bcType == EM_BC_SYMMETRY then
            table.insert(self.ssBoundaryConditions,
                         makeSsBcUpdater(dir, inOut, { bcSymmetry }))
         elseif bcType == EM_BC_AXIS then
            table.insert(self.ssBoundaryConditions,
                         makeSsBcUpdater(dir, inOut,  PerfMaxwell.bcAxis ))
         elseif type(bcType) == "table" then
            table.insert(self.ssBoundaryConditions,
                         makeSsBcUpdater(dir, inOut, bcType))
         else
            assert(false, "MaxwellField.ssBoundaryConditions: Unsupported BC type!")
         end
      end
      for dir = 1, self.grid:ndim() do
         appendSsBoundaryConditions(dir, self._inOut, self.ssBc[1])
      end
   end

   self.bcTime = 0.0 -- Timer for BCs.
end

function MaxwellField:createDiagnostics() end

function MaxwellField:initField(population)
   if self.hasMagField then   -- Maxwell's induction equations.
      -- Create field for total current density. Need to do this
      -- here because field object does not know about vdim.
      self.currentDensLocal, self.currentDens = nil, nil
      local vdim = nil
      for _, s in population.iterLocal() do
         vdim = vdim and vdim or s.vdim
         self.currentDens = self.currentDens and self.currentDens or s:allocVectorMoment(vdim)
         self.currentDensLocal = self.currentDensLocal and self.currentDensLocal or s:allocVectorMoment(vdim)
         assert(vdim == s.vdim, "MaxwellField: currently don't support species with different vdim.")
      end

      local project = Updater.ProjectOnBasis {
         onGrid = self.grid,   evaluate = self.initFunc,
         basis  = self.basis,
      }
      project:advance(0.0, {}, {self.em[1]})
      self:applyBc(0.0, self.em[1])
   else   -- Poisson equation. Solve for initial phi.
      -- Create field for total charge density.
      self.chargeDensLocal, self.chargeDens = nil, nil
      for _, s in population.iterLocal() do
         self.chargeDens = self.chargeDens and self.chargeDens or s:allocMoment()
         self.chargeDensLocal = self.chargeDensLocal and self.chargeDensLocal or s:allocMoment()
         break
      end

      self:advance(0.0, population, 1, 1)
      local emStart = self:rkStepperFields()[1]
      for i = 1, self.phiPrevNum do self.phiFldPrev[i]["fld"]:copy(emStart) end
   end
end

function MaxwellField:write(tm, force)
   if self.evolve or self.forceWrite then
      local tmStart = Time.clock()
      -- Compute EM energy integrated over domain.
      if self.calcIntEMQuantTrigger(tm) then
         self.emEnergyCalc(tm, { self.em[1] }, { self.emEnergy })
      end
      -- Time computation of integrated moments.
      self.integratedEMTime = self.integratedEMTime + Time.clock() - tmStart

      if self.ioTrigger(tm) or force then
	 self.fieldIo:write(self.em[1], string.format("field_%d.bp", self.ioFrame), tm, self.ioFrame)
	 self.emEnergy:write(string.format("fieldEnergy.bp"), tm, self.ioFrame)
	 
	 self.ioFrame = self.ioFrame+1
      end
   else
      -- If not evolving species, don't write anything except initial conditions.
      if self.ioFrame == 0 then
	 self.fieldIo:write(self.em[1], string.format("field_%d.bp", self.ioFrame), tm, self.ioFrame)
      end
      self.ioFrame = self.ioFrame+1
   end
end

function MaxwellField:writeRestart(tm)
   -- (the final "false" prevents writing of ghost cells)
   self.fieldIo:write(self.em[1], "field_restart.bp", tm, self.ioFrame, false)

   -- (the first "false" prevents flushing of data after write, the second "false" prevents appending)
   self.emEnergy:write("fieldEnergy_restart.bp", tm, self.dynVecRestartFrame, false, false)
   self.dynVecRestartFrame = self.dynVecRestartFrame + 1
end

function MaxwellField:readRestart()
   local tm, fr = self.fieldIo:read(self.em[1], "field_restart.bp")
   self:applyBc(tm, self.em[1])
   self.em[1]:sync() -- Must get all ghost-cell data correct.
     
   self.emEnergy:read("fieldEnergy_restart.bp", tm)
   self.ioFrame = fr
   -- Iterate triggers.
   self.ioTrigger(tm)

   return tm
end

function MaxwellField:rkStepperFields()
   return self.em
end

-- For RK timestepping.
function MaxwellField:copyRk(outIdx, aIdx)
   if self:rkStepperFields()[aIdx] then self:rkStepperFields()[outIdx]:copy(self:rkStepperFields()[aIdx]) end
end
-- For RK timestepping
function MaxwellField:combineRk(outIdx, a, aIdx, ...)
   if self:rkStepperFields()[aIdx] and self.hasMagField then 
      local args = {...} -- Package up rest of args as table.
      local nFlds = #args/2
      self:rkStepperFields()[outIdx]:combine(a, self:rkStepperFields()[aIdx])
      for i = 1, nFlds do -- Accumulate rest of the fields.
         self:rkStepperFields()[outIdx]:accumulate(args[2*i-1], self:rkStepperFields()[args[2*i]])
      end	 
   end
end

function MaxwellField:suggestDt()
   if self.hasMagField then 
      return math.min(self.cfl/self.cflRateByCell:reduce('max')[1], GKYL_MAX_DOUBLE)
   else
      return GKYL_MAX_DOUBLE
   end
end

function MaxwellField:clearCFL()
   -- Clear cflRateByCell for next cfl calculation.
   if self.hasMagField then self.cflRateByCell:clear(0.0) end
end

function MaxwellField:accumulateCurrent(current, emRhs)
   if current == nil then return end

   local tmStart = Time.clock()

--   local cItr, eItr = current:get(1), emRhs:get(1)
--   local cIdxr, eIdxr = current:genIndexer(), emRhs:genIndexer()

   -- If we are to use ghost currents, compute mean current first.
   local ghostCurrent = 0.0
   if self.useGhostCurrent then
      assert(false, "ghost current not yet implemented in g0 merge")
      --local nx = self.grid:numCells(1)
      --local localMeanCurrent = ffi.new("double[2]")
      --for idx in emRhs:localRangeIter() do
      --   current:fill(cIdxr(idx), cItr)
      --   localMeanCurrent[0] = localMeanCurrent[0]+cItr[1]
      --end
      --local globalMeanCurrent = ffi.new("double[2]")
      --Mpi.Allreduce(localMeanCurrent, globalMeanCurrent, 1, Mpi.DOUBLE, Mpi.SUM, self.grid:commSet().comm)
      --ghostCurrent = globalMeanCurrent[0]/nx
   end

   --for idx in emRhs:localRangeIter() do
   --   current:fill(cIdxr(idx), cItr)
   --   emRhs:fill(eIdxr(idx), eItr)
   --   --eItr[1] = eItr[1]-1.0/self.epsilon0*(cItr[1]-ghostCurrent)
   --   for i = 1, current:numComponents() do
   --      eItr[i] = eItr[i]-1.0/self.epsilon0*cItr[i]
   --   end
   --end
   emRhs:accumulateRange(-1.0/self.epsilon0, current, emRhs:localRange())
   self.tmCurrentAccum = self.tmCurrentAccum + Time.clock()-tmStart
end

function MaxwellField:advance(tCurr, population, inIdx, outIdx)
   local emIn = self:rkStepperFields()[inIdx]
   local emRhsOut = self:rkStepperFields()[outIdx]

   if self.hasMagField then   -- Maxwell's induction equations.
      if self.evolve then
         self.fieldSlvr:advance(tCurr, {emIn}, {emRhsOut, self.cflRateByCell})
         if self.currentDens then -- No currents for source-free Maxwell.
            self.currentDensLocal:clear(0.0)
            -- Accumulate currents of local species.
            for _, s in population.iterLocal() do
               self.currentDensLocal:accumulate(s:getCharge(), s:getMomDensity())
            end
            -- Reduce currents across species communicator.
            self.currentDens:clear(0.0)
            self.currentDens:reduceByCell('sum', population:getComm(), self.currentDensLocal)
            -- Add species current to curl{B} in Amepere's equation.
            self:accumulateCurrent(self.currentDens, emRhsOut)
         end
      else
         emRhsOut:clear(0.0)   -- No RHS.
      end
   else   -- Poisson equation. Solve for phi.
      -- Accumulate the charge density of local species.
      self.chargeDensLocal:clear(0.0)
      for _, s in population.iterLocal() do
         self.chargeDensLocal:accumulate(s:getCharge(), s:getNumDensity())
      end
      -- Reduce charge density across species communicator.
      self.chargeDens:clear(0.0)
      self.chargeDens:reduceByCell('sum', population:getComm(), self.chargeDensLocal)

--      if self.isSlvrMG then
--         self.chargeDens:scale(1.0/self.epsilon0)
--         if inIdx == 1 then
--            -- In the first RK stage shuffle the storage of previous potentials.
--            for i = 1, self.phiPrevNum-1 do
--               self.phiFldPrev[i]["fld"]:copy(self.phiFldPrev[i+1]["fld"])
--               self.phiFldPrev[i]["time"] = self.phiFldPrev[i+1]["time"]
--            end
--            self.phiFldPrev[self.phiPrevNum]["fld"]:copy(emIn)
--            self.phiFldPrev[self.phiPrevNum]["time"] = tCurr
--            -- Count until phiPrevNum time steps have been taken.
--            if not self.filledPhiPrev then
--               self.phiPrevCount = self.phiPrevCount+1
--               if self.phiPrevCount > self.phiPrevNum then self.filledPhiPrev = true end
--            end
--         end
--         if self.filledPhiPrev then
--            -- Form an initial guess with 3-point Lagrange extrapolation.
--            local tMt1 = tCurr-self.phiFldPrev[1]["time"]
--            local tMt2 = tCurr-self.phiFldPrev[2]["time"]
--            local tMt3 = tCurr-self.phiFldPrev[3]["time"]
--            local t1Mt2, t1Mt3 = tMt2-tMt1, tMt3-tMt1
--            local t2Mt1, t2Mt3 = -t1Mt2, tMt3-tMt2
--            local t3Mt1, t3Mt2 = -t1Mt3, -t2Mt3
--            local f1 = tMt2*tMt3/(t1Mt2*t1Mt3)
--            local f2 = tMt1*tMt3/(t2Mt1*t2Mt3)
--            local f3 = tMt1*tMt2/(t3Mt1*t3Mt2)
--            emIn:combine(f1,self.phiFldPrev[1]["fld"],
--                         f2,self.phiFldPrev[2]["fld"],
--                         f3,self.phiFldPrev[3]["fld"]) 
--         end
--      end

      -- Solve for the potential.
      self.fieldSlvr:advance(tCurr, {self.chargeDens}, {emIn})
   end
end

function MaxwellField:updateInDirection(dir, tCurr, dt, fIn, fOut)
   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   if self.evolve then
      self:applyBc(tCurr, fIn, dir)
      self.fieldHyperSlvr[dir]:setDtAndCflRate(dt, nil)
      status, dtSuggested = self.fieldHyperSlvr[dir]:advance(tCurr, {fIn}, {fOut})
   else
      fOut:copy(fIn)
   end
   return status, dtSuggested   
end

function MaxwellField:applyBcIdx(tCurr, idx)
   self:applyBc(tCurr, self:rkStepperFields()[idx])
end 

function MaxwellField:applyBc(tCurr, emIn, dir)
   local tmStart = Time.clock()
   if self.hasNonPeriodicBc then
      for _, bc in ipairs(self.boundaryConditions) do
         if (not dir) or dir == bc:getDir() then
	    bc:advance(tCurr, {}, {emIn})
         end
      end
   end   

   if self._hasSsBnd then
      emIn:sync()
      for _, bc in ipairs(self.ssBoundaryConditions) do
         if (not dir) or dir == bc:getDir() then
            bc:advance(tCurr, {}, {emIn})
         end
      end
   end

   emIn:sync()
   self.bcTime = self.bcTime + Time.clock()-tmStart
end
 
MaxwellField.bcConst = function(Ex, Ey, Ez, Bx, By, Bz, phiE, phiB)
   local bc = BoundaryCondition.Const {
      components = {1, 2, 3, 4, 5, 6, 7, 8},
      values = {Ex, Ey, Ez, Bx, By, Bz, phiE, phiB}
   }
   return { bc }
end
  
function MaxwellField:totalSolverTime()
   local ftm = 0.0
   if self.fieldSlvr then
      ftm = self.fieldSlvr.totalTime
   else
      for d = 1, self.grid:ndim() do
	 ftm = ftm+self.fieldHyperSlvr[d].totalTime
      end
   end
   return ftm + self.tmCurrentAccum
end

function MaxwellField:totalBcTime()
   return self.bcTime
end

function MaxwellField:energyCalcTime()
   return self.integratedEMTime
end

-- ExternalMaxwellField ---------------------------------------------------------------------
--
-- A field object with external EM fields specified as a time-dependent function.
-----------------------------------------------------------------------------------------

local ExternalMaxwellField = Proto(FieldBase.ExternalFieldBase)

-- Methods for no field object.
function ExternalMaxwellField:init(tbl)
   ExternalMaxwellField.super.init(self, tbl)
   self.tbl = tbl
end

function ExternalMaxwellField:fullInit(appTbl)
   local tbl = self.tbl -- Previously store table.

   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- By default evolve field.

   self.hasMagField = xsys.pickBool(tbl.hasMagneticField, true) -- By default there is a magnetic field.

   -- Create triggers to write fields.
   if tbl.nFrame then
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nFrame)
   else
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.ioFrame = 0 -- Frame number for IO.
   
   -- Store function to compute EM field.
   if self.hasMagField then
      self.emFunc = function (t, xn)
         local ex, ey, ez, bx, by, bz = tbl.emFunc(t, xn)
         return ex, ey, ez, bx, by, bz, 0.0, 0.0
      end
   else
      self.emFunc = function (t, xn)
         local ex, ey, ez = tbl.emFunc(t, xn)
         return ex, ey, ez
      end
   end
end

function ExternalMaxwellField:hasEB() return true, self.hasMagField end
function ExternalMaxwellField:setGrid(grid) self.grid = grid end

function ExternalMaxwellField:alloc(nField)
   local nGhost = 2
   if self.basis:numBasis() > 1 then
      nGhost = 1
   end

   -- Allocate fields needed in RK update.
   local emVecComp = 8
   if not self.hasMagField then emVecComp = 1 end  -- Electric field only.
   self.em = DataStruct.Field {
      onGrid = self.grid,
      numComponents = emVecComp*self.basis:numBasis(),
      ghost         = {nGhost, nGhost},
   }
end

function ExternalMaxwellField:createSolver(population)
   self.fieldSlvr = Updater.ProjectOnBasis {
      onGrid = self.grid,   evaluate = self.emFunc,
      basis  = self.basis,  onGhosts = true,
   }
      
   -- Create Adios object for field I/O.
   self.fieldIo = AdiosCartFieldIo {
      elemType = self.em:elemType(),
      method   = self.ioMethod,
      metaData = {polyOrder = self.basis:polyOrder(),
                  basisType = self.basis:id(),
                  epsilon0  = self.epsilon0,
                  mu0       = self.mu0,
                  grid      = GKYL_OUT_PREFIX .. "_grid.bp"},
      writeRankInComm = {0, population:getComm(),},
   }   
end

function ExternalMaxwellField:createDiagnostics()
end

function ExternalMaxwellField:initField()
   self.fieldSlvr:advance(0.0, {}, {self.em})
   self:applyBc(0.0, self.em)
end

function ExternalMaxwellField:write(tm, force)
   if self.evolve or self.forceWrite then
      if self.ioTrigger(tm) or force then
	 self.fieldIo:write(self.em, string.format("externalField_%d.bp", self.ioFrame), tm, self.ioFrame)
         self.ioFrame = self.ioFrame+1
      end
   else
      -- If not evolving species, don't write anything except initial conditions.
      if self.ioFrame == 0 then
	 self.fieldIo:write(self.em, string.format("externalField_%d.bp", self.ioFrame), tm, self.ioFrame)
      end
      self.ioFrame = self.ioFrame+1
   end
end

function ExternalMaxwellField:writeRestart(tm)
   self.fieldIo:write(self.em, "externalField_restart.bp", tm, self.ioFrame)
end

function ExternalMaxwellField:readRestart()
   local tm, fr = self.fieldIo:read(self.em, "externalField_restart.bp")
   self.em:sync() -- Must get all ghost-cell data correct.
   
   self.ioFrame = fr
end

function ExternalMaxwellField:rkStepperFields()
   return { self.em, self.em, self.em, self.em }
end

function ExternalMaxwellField:advance(tCurr)
   local emOut = self:rkStepperFields()[1]
   if self.evolve then
      self.fieldSlvr:advance(tCurr, {}, {emOut})
      self:applyBc(tCurr, emOut)
   end
end

function ExternalMaxwellField:applyBcIdx(tCurr, idx)
   self:applyBc(tCurr, self:rkStepperFields()[1])
end

function ExternalMaxwellField:applyBc(tCurr, emIn)
   emIn:sync()
end

function ExternalMaxwellField:totalSolverTime()
   return self.fieldSlvr.totalTime
end

function ExternalMaxwellField:totalBcTime() return 0.0 end
function ExternalMaxwellField:energyCalcTime() return 0.0 end

return {
   MaxwellField = MaxwellField,
   ExternalMaxwellField = ExternalMaxwellField,
   FuncMaxwellField = ExternalMaxwellField,   -- For backwards compatibility.
}
