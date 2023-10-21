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
local Range             = require "Lib.Range"
local FieldBase         = require "App.Field.FieldBase"
local Grid              = require "Grid"
local DecompRegionCalc  = require "Lib.CartDecomp"
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

MaxwellField.bcOpen     = EM_BC_OPEN    -- Zero gradient.
MaxwellField.bcCopy     = EM_BC_COPY    -- Copy fields.
MaxwellField.bcReflect  = EM_BC_REFLECT -- Perfect electric conductor.
MaxwellField.bcSymmetry = EM_BC_SYMMETRY
MaxwellField.bcAxis     = EM_BC_AXIS

-- Function to check if BC type is good.
local function isBcGood(bcType)
   if bcType == EM_BC_OPEN or bcType == EM_BC_REFLECT or bcType == EM_BC_SYMMETRY then
      return true
   end
   if type(bcType) == "table" then return true end
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
function MaxwellField:fullInit(appTbl, plasmaField)
   local tbl = self.tbl -- Previously store table.
   
   self.epsilon0 = tbl.epsilon0
   self.mu0      = tbl.mu0
   self.evolve   = xsys.pickBool(tbl.evolve, true) -- By default evolve field.

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

      self._inOutFunc = tbl.inOutFunc

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

   -- Create trigger for how frequently to compute integrated EM fields.
   if appTbl.calcIntQuantEvery then
      self.calcIntEMQuantTrigger = LinearTrigger(0, appTbl.tEnd,  math.floor(1/appTbl.calcIntQuantEvery))
   else
      self.calcIntEMQuantTrigger = function(t) return true end
   end

   self.timers = {advance = 0.,   bc = 0.}
end

-- Methods for EM field object.
function MaxwellField:hasEB() return true, self.hasMagField end
function MaxwellField:setCfl(cfl) self.cfl = cfl end
function MaxwellField:getCfl() return self.cfl end
function MaxwellField:setGrid(grid, myadios)
   self.grid    = grid
   self.myadios = myadios

   self.ndim = self.grid:ndim()

   local keepDims = {};  for i = 1, self.ndim do keepDims[i] = i end
   local gridInfo = grid:childGrid(keepDims)
   local GridConstructor = grid.mapc2p and Grid.MappedCart or Grid.RectCart

   -- Create global grid for Poisson solver.
   -- Need a communicator with just ourselves.
   local confComm = self.grid:getMessenger():getConfComm_host()
   local confGroup = Mpi.Comm_group(confComm)
   local noDecompGroupRanks = Lin.IntVec(1);  noDecompGroupRanks[1] = self.grid:getMessenger():getConfRank()
   local noDecompGroup = Mpi.Group_incl(confGroup, 1, noDecompGroupRanks:data());
   local tag = self.grid:getMessenger():getConfRank()
   self.noDecompComm = Mpi.Comm_create_group(confComm, noDecompGroup, tag);
   Mpi.Group_free(confGroup);  Mpi.Group_free(noDecompGroup)
   local noCuts = {};  for d = 1,self.ndim do noCuts[d] = 1 end
   local noDecomp = DecompRegionCalc.CartProd {  -- Decomp in x-y but not z.
      cuts = noCuts,  comm = self.noDecompComm,
   }
   self.gridGlobal = GridConstructor {
      lower = gridInfo.lower,  world    = gridInfo.world,
      upper = gridInfo.upper,  mappings = gridInfo.coordinateMap,
      cells = gridInfo.cells,  mapc2p   = gridInfo.mapc2p,
      periodicDirs  = gridInfo.periodicDirs,
      messenger     = gridInfo.messenger,
      decomposition = noDecomp,
   }
end

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

local function createField(grid, basis, ghostCells, ncomp, periodicSync, useDevice)
   local metadata = basis and {polyOrder = basis:polyOrder(), basisType = basis:id()} or nil
   local fld = DataStruct.Field {
      onGrid   = grid,        numComponents    = ncomp,
      metaData = metadata,    syncPeriodicDirs = periodicSync,
      ghost    = ghostCells,  useDevice        = useDevice,
   }
   fld:clear(0.0)
   return fld
end

function MaxwellField:alloc(nRkDup)
   
   local ghostNum = {1,1}

   self.em = {}
   if self.hasMagField then   -- Maxwell's induction equations.
      -- Allocate fields needed in RK update.
      for i = 1, nRkDup do
         self.em[i] = createField(self.grid, self.basis, ghostNum, 8*self.basis:numBasis())
      end

      -- Array with one component per cell to store cflRate in each cell.
      self.cflRateByCell = createField(self.grid, self.basis, ghostNum, 1)
      self.cflRateByCell:clear(0.0)
      self.cflRatePtr  = self.cflRateByCell:get(1)
      self.cflRateIdxr = self.cflRateByCell:genIndexer()
      self.dtGlobal    = ffi.new("double[2]")
      
      -- For storing integrated energy components.
      self.emEnergy = DataStruct.DynVector { numComponents = 8, adiosSystem = self.myadios, }

   else   -- Poisson equation.
      -- Electrostatic potential, phi, and external magnetic potential A_ext (computed by ExternalField).
      local potentials = createField(self.grid, self.basis, ghostNum, 4*self.basis:numBasis())
      for i = 1, nRkDup do self.em[i] = potentials end

      -- Extra fields needed for the distributed parallel field solve.
      self.localBuffer  = createField(self.grid, self.basis, {0,0}, self.basis:numBasis())
      self.globalBuffer = createField(self.gridGlobal, self.basis, {0,0}, self.basis:numBasis())
      self.globalSol    = createField(self.gridGlobal, self.basis, {1,1}, self.basis:numBasis())
   
      -- For storing integrated energy components.
      self.emEnergy = DataStruct.DynVector { numComponents = self.grid:ndim(), adiosSystem = self.myadios, }
   end
end

function MaxwellField:createSolver(population)

   -- Create Adios object for field I/O.
   self.fieldIo = AdiosCartFieldIo {
      elemType = self.em[1]:elemType(),
      metaData = {polyOrder = self.basis:polyOrder(),
                  basisType = self.basis:id(),
                  epsilon0  = self.epsilon0,
                  mu0       = self.mu0,
                  grid      = GKYL_OUT_PREFIX .. "_grid.bp"},
      writeRankInComm = {0, population:getComm_host(),},
   }

   if self.hasMagField then   -- Maxwell's induction equations.

      self.equation = PerfMaxwell {
         lightSpeed          = self.lightSpeed,  tau     = self.tau,
         elcErrorSpeedFactor = self.ce,          numFlux = self.numFlux,
         mgnErrorSpeedFactor = self.cb,          basis   = self.basis,
      }

      self.fieldSlvr = Updater.HyperDisCont {
         onGrid = self.grid,   cfl      = self.cfl,
         basis  = self.basis,  equation = self.equation,
      }

      self.emEnergyUpd = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.grid,   quantity      = "V2",
         basis  = self.basis,  numComponents = 8,
      }
      self.emEnergyCalc = function(tCurr, inFld, outDynV) self.emEnergyUpd:advance(tCurr, inFld, outDynV) end

   else   -- Poisson equation.

      self.isElliptic = true

      local decompRange = self.grid:decomposedRange()
      local subdomIdx   = {};  decompRange:cutsInvIndexer()(self.grid:subGridId(), subdomIdx)

      -- Range of MPI processes participating in field solve.
      local ndim = self.grid:ndim()
      self.cutsRange = decompRange:cutsRange()  -- Range of MPI processes sharing a field solve.
      -- Ranges in a global field into which we will copy a buffer obtained from other processes.
      self.destRange, self.bufferOffset = {}, {}
      self.cutsRangeIdxr = self.cutsRange:indexer(Range.colMajor)
      for idx in self.cutsRange:colMajorIter() do
         local linIdx    = decompRange:cutsIndexer()(idx)
         local destRange = decompRange:subDomain(linIdx)
         local lv, uv    = destRange:lowerAsVec(), destRange:upperAsVec()
         self.destRange[linIdx]    = self.globalSol:localExtRange():subRange(lv,uv)
         self.bufferOffset[linIdx] = (linIdx-1)*self.localBuffer:size()
      end
      -- Range in a global field we'll copy the local solution from.
      local lv, uv = self.em[1]:localRange():lowerAsVec(), self.em[1]:localRange():upperAsVec()
      self.srcRange = self.globalSol:localExtRange():subRange(lv,uv)

      self.fieldSlvr = Updater.FemPoisson {
         onGrid  = self.gridGlobal,  bcLower = self.bcLowerPhi,
         basis   = self.basis,       bcUpper = self.bcUpperPhi,
         epsilon = self.epsilon0,
      }
      self.esEnergyUpd = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.grid,   quantity = "V2",
         basis  = self.basis,
      }
      self.emEnergyCalc = function(tCurr, inFld, outDynV) self:esEnergy(tCurr, inFld, outDynV) end
   end

   -- Function to construct a BC updater.
   local function makeBcUpdater(dir, edge, bcList)
      return Updater.Bc {
         onGrid             = self.grid,  dir  = dir,
         boundaryConditions = bcList,     edge = edge,
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

   self.bcTime = 0.0 -- Timer for BCs.
end

function MaxwellField:AllgatherField(fldLocal, fldGlobal)
   -- Gather a CartField distributed in configuration space onto a global field.
   fldLocal:copyRangeToBuffer(fldLocal:localRange(), self.localBuffer:data())

   -- Gather flat buffers from other processes.
   local mess = self.grid:getMessenger()
   local comm = mess:getComms()["conf"]
   mess:Allgather(self.localBuffer, self.globalBuffer, comm)

   -- Rearrange into a global field, i.e. reconstruct the multi-D array
   -- using the collection of flat buffers we gathered from each process.
   for idx in self.cutsRange:colMajorIter() do
      local linIdx = self.cutsRangeIdxr(idx[1],idx[2],idx[3])
      fldGlobal:copyRangeFromBuffer(self.destRange[linIdx],
         self.globalBuffer:data()+self.bufferOffset[linIdx])
   end
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
   end
end

function MaxwellField:write(tm, force)
   if self.evolve or self.forceWrite then
      -- Compute EM energy integrated over domain.
      if self.calcIntEMQuantTrigger(tm) then
         self.emEnergyCalc(tm, { self.em[1] }, { self.emEnergy })
      end

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

function MaxwellField:rkStepperFields() return self.em end

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
   emRhs:accumulateRange(-1.0/self.epsilon0, current, emRhs:localRange())
end

function MaxwellField:advance(tCurr, population, inIdx, outIdx)
   local tmStart = Time.clock()
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
	    population:AllreduceByCell(self.currentDensLocal, self.currentDens, 'sum')
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
      population:AllreduceByCell(self.chargeDensLocal, self.chargeDens, 'sum')

      -- Solve for the potential.
      self:AllgatherField(self.chargeDens, self.globalSol)  -- Gather charge density in global field.

      self.fieldSlvr:advance(tCurr, {self.globalSol}, {self.globalSol})

      -- Copy portion of global field belonging to this process.
      emIn:copyRangeToRange(self.globalSol, emIn:localRange(), self.srcRange)
   end
   self.timers.advance = self.timers.advance + Time.clock() - tmStart
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

   emIn:sync()
   self.timers.bc = self.timers.bc + Time.clock() - tmStart
end

function MaxwellField:clearTimers()
   for nm, _ in pairs(self.timers) do self.timers[nm] = 0. end
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

function ExternalMaxwellField:fullInit(appTbl, plasmaField)
   local tbl = self.tbl -- Previously store table.

   self.evolve = xsys.pickBool(tbl.evolve, true) -- By default evolve field.

   -- By default there is a plasma-generated magnetic field, unless running Vlasov-Poisson.
   if plasmaField then
      self.hasPlasmaMagField = xsys.pickBool(plasmaField.hasMagField, true)
   else
      self.hasPlasmaMagField = true -- Assumption
   end

   self.hasMagField = xsys.pickBool(tbl.hasMagneticField, true) -- By default there is a magnetic field.

   -- Create triggers to write fields.
   if tbl.nFrame then
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nFrame)
   else
      self.ioTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.ioFrame = 0 -- Frame number for IO.
   
   -- Store function to compute the external field.
   if self.hasPlasmaMagField then
      if self.hasMagField then
         self.emFunc = function (t, xn)
            local ex, ey, ez, bx, by, bz = tbl.emFunc(t, xn)
            return ex, ey, ez, bx, by, bz, 0., 0.
         end
      else
         self.emFunc = function (t, xn)
            local ex, ey, ez = tbl.emFunc(t, xn)
            return ex, ey, ez, 0., 0., 0., 0., 0.
         end
      end
   else
      if self.hasMagField then
         self.emFunc = function (t, xn)
            local phi, ax, ay, az = tbl.emFunc(t, xn)
            return phi, ax, ay, az
         end
      else
         self.emFunc = function(t, xn)
            local phi = tbl.emFunc(t, xn)
            return phi, 0., 0., 0.
         end
      end
   end

   self.timers = {advance = 0.,   bc = 0.}
end

function ExternalMaxwellField:hasEB() return true, self.hasMagField end
function ExternalMaxwellField:setGrid(grid, myadios)
   self.grid    = grid
   self.myadios = myadios
end

function ExternalMaxwellField:alloc(nField)
   local ghostNum = {1,1}

   -- Allocate fields needed in RK update.
   local emVecComp = self.hasPlasmaMagField and 8 or 4
   self.em = createField(self.grid, self.basis, ghostNum, emVecComp*self.basis:numBasis())
end

function ExternalMaxwellField:createSolver(population)
   self.fieldSlvr = Updater.EvalOnNodes {
      onGrid = self.grid,   evaluate = self.emFunc,
      basis  = self.basis,  onGhosts = false,
   }
      
   -- Create Adios object for field I/O.
   self.fieldIo = AdiosCartFieldIo {
      elemType = self.em:elemType(),
      metaData = {polyOrder = self.basis:polyOrder(),
                  basisType = self.basis:id(),
                  epsilon0  = self.epsilon0,
                  mu0       = self.mu0,
                  grid      = GKYL_OUT_PREFIX .. "_grid.bp"},
      writeRankInComm = {0, population:getComm(),},
   }   
end

function ExternalMaxwellField:createDiagnostics() end

function ExternalMaxwellField:initField(population)
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
   local tmStart = Time.clock()
   local emOut = self:rkStepperFields()[1]
   if self.evolve then
      self.fieldSlvr:advance(tCurr, {}, {emOut})
      self:applyBc(tCurr, emOut)
   end
   self.timers.advance = self.timers.advance + Time.clock() - tmStart
end

function ExternalMaxwellField:applyBcIdx(tCurr, idx)
   self:applyBc(tCurr, self:rkStepperFields()[1])
end

function ExternalMaxwellField:applyBc(tCurr, emIn)
   local tmStart = Time.clock()
   emIn:sync()
   self.timers.bc = self.timers.bc + Time.clock() - tmStart
end

function ExternalMaxwellField:clearTimers()
   for nm, _ in pairs(self.timers) do self.timers[nm] = 0. end
end

return {
   MaxwellField         = MaxwellField,
   ExternalMaxwellField = ExternalMaxwellField,
   FuncMaxwellField     = ExternalMaxwellField,   -- For backwards compatibility.
}
