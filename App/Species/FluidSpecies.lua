-- Gkyl ------------------------------------------------------------------------
--
-- App support code: FluidSpecies object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Basis = require "Basis"
local Collisions = require "App.Collisions"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local ffi = require "ffi"
local Grid = require "Grid"
local LinearTrigger = require "Lib.LinearTrigger"
local Proto = require "Lib.Proto"
local SpeciesBase = require "App.Species.SpeciesBase"
local Time = require "Lib.Time"
local Updater = require "Updater"
local xsys = require "xsys"

-- function to create basis functions
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   end
end

-- base class for kinetic species
local FluidSpecies = Proto(SpeciesBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function FluidSpecies:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function FluidSpecies:fullInit(appTbl)
   local tbl = self.tbl -- previously store table

   self.cfl =  0.1
   self.charge = tbl.charge and tbl.charge or 1.0
   self.mass = tbl.mass and tbl.mass or 1.0
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default, evolve species
   self.evolveCollisionless = xsys.pickBool(tbl.evolveCollisionless,
                                            self.evolve)
   self.evolveCollisions = xsys.pickBool(tbl.evolveCollisions, self.evolve)
   self.confBasis = nil -- Will be set later

   -- create triggers to write diagnostics
   if tbl.nDiagnosticFrame then
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nDiagnosticFrame)
   else
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.diagIoFrame = 0 -- frame number for diagnostics

   -- for storing integrated moments
   self.integratedMoments = nil -- allocated in alloc() method

   -- store initial condition function
   self.initFunc = tbl.init

   -- default to a single moment
   self.nMoments = 1
   self.nGhost = 1 -- default is 1 ghost-cell in each direction

   self.hasNonPeriodicBc = false -- to indicate if we have non-periodic BCs
   self.bcx, self.bcy, self.bcz = { }, { }, { }

   -- read in boundary conditions
   -- check to see if bc type is good is now done in createBc
   if tbl.bcx then
      self.bcx[1], self.bcx[2] = tbl.bcx[1], tbl.bcx[2]
      self.hasNonPeriodicBc = true
   end
   if tbl.bcy then
      self.bcy[1], self.bcy[2] = tbl.bcy[1], tbl.bcy[2]
      self.hasNonPeriodicBc = true
   end
   if tbl.bcz then
      self.bcz[1], self.bcz[2] = tbl.bcz[1], tbl.bcz[2]
      self.hasNonPeriodicBc = true
   end

   self.boundaryConditions = { } -- list of Bcs to apply
   self.zeroFluxDirections = {}

   self.bcTime = 0.0 -- timer for BCs

   -- Collisions: currently used for a diffusion term.
   self.collisions = {}
   for nm, val in pairs(tbl) do
      if Collisions.CollisionsBase.is(val) then
         self.collisions[nm] = val
         self.collisions[nm]:setName(nm)
         val:setSpeciesName(self.name)
         val:fullInit(tbl)    -- Initialize collisions (diffusion).
      end
   end

   self.useShared = xsys.pickBool(appTbl.useShared, false)
   self.positivity = xsys.pickBool(tbl.applyPositivity, false)
   self.positivityDiffuse = xsys.pickBool(tbl.positivityDiffuse, self.positivity)
   self.positivityRescale = xsys.pickBool(tbl.positivityRescale, false)
   self.deltaF = xsys.pickBool(appTbl.deltaF, false)

   self.tCurr = 0.0
end

function FluidSpecies:getCharge() return self.charge end
function FluidSpecies:getMass() return self.mass end
function FluidSpecies:getEvolve() return self.evolve end

function FluidSpecies:getNdim()
   return self.ndim
end
function FluidSpecies:vdim()
   return 0
end
function FluidSpecies:setName(nm)
   self.name = nm
end
function FluidSpecies:setCfl(cfl)
   self.cfl = cfl
   for _, c in pairs(self.collisions) do
      c:setCfl(cfl)
   end
end
function FluidSpecies:setIoMethod(ioMethod)
   self.ioMethod = ioMethod
end
function FluidSpecies:setConfBasis(basis)
   self.confBasis = basis
   for _, c in pairs(self.collisions) do
      c:setConfBasis(basis)
   end
end
function FluidSpecies:setConfGrid(cgrid)
   self.confGrid = cgrid
   for _, c in pairs(self.collisions) do
      c:setConfGrid(cgrid)
   end
end

function FluidSpecies:createGrid(cLo, cUp, cCells, cDecompCuts, cPeriodicDirs)
   self.cdim = #cCells
   self.ndim = self.cdim

   -- create decomposition
   local decompCuts = {}
   for d = 1, self.cdim do table.insert(decompCuts, cDecompCuts[d]) end
   self.decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      useShared = self.useShared,
   }

   -- create computational domain
   local lower, upper, cells = {}, {}, {}
   for d = 1, self.cdim do
      table.insert(lower, cLo[d])
      table.insert(upper, cUp[d])
      table.insert(cells, cCells[d])
   end
   self.grid = Grid.RectCart {
      lower = lower,
      upper = upper,
      cells = cells,
      periodicDirs = cPeriodicDirs,
      decomposition = self.decomp,
   }
end

function FluidSpecies:createBasis(nm, polyOrder)
   self.basis = createBasis(nm, self.ndim, polyOrder)
end

function FluidSpecies:allocMoment()
   local m = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost = {self.nGhost, self.nGhost}
   }
   return m
end
function FluidSpecies:allocVectorMoment(dim)
   local m = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.confBasis:numBasis()*dim,
      ghost = {self.nGhost, self.nGhost}
   }
   return m
end

function FluidSpecies:allocMomCouplingFields()
   return {self:allocVectorMoment(self.nMoments)}
end

function FluidSpecies:bcAbsorbFunc(dir, tm, idxIn, fIn, fOut)
   -- note that for bcAbsorb there is no operation on fIn,
   -- so skinLoop (which determines indexing of fIn) does not matter 
   for i = 1, self.nMoments*self.basis:numBasis() do
      fOut[i] = 0.0
   end
end

function FluidSpecies:bcCopyFunc(dir, tm, idxIn, fIn, fOut)
   for i = 1, self.nMoments*self.basis:numBasis() do
      fOut[i] = fIn[i]
   end
end

-- function to construct a BC updater
function FluidSpecies:makeBcUpdater(dir, edge, bcList)
   return Updater.Bc {
      onGrid = self.grid,
      boundaryConditions = bcList,
      dir = dir,
      edge = edge,
   }
end

function FluidSpecies:createBCs()
   -- functions to make life easier while reading in BCs to apply
   -- note: appendBoundaryConditions defined in sub-classes
   local function handleBc(dir, bc)
      if bc[1] then
	 self:appendBoundaryConditions(dir, "lower", bc[1])
      end
      if bc[2] then
	 self:appendBoundaryConditions(dir, "upper", bc[2])
      end
   end

   -- add various BCs to list of BCs to apply
   handleBc(1, self.bcx)
   handleBc(2, self.bcy)
   handleBc(3, self.bcz)
end

function FluidSpecies:createSolver(funcField)
   -- Create solvers for collisions (diffusion).
   for _, c in pairs(self.collisions) do
      c:createSolver(funcField)
   end
end

function FluidSpecies:alloc(nRkDup)
   -- allocate fields needed in RK update
   self.moments = {}
   for i = 1, nRkDup do
      self.moments[i] = self:allocVectorMoment(self.nMoments)
   end
   -- create Adios object for field I/O
   self.momIo = AdiosCartFieldIo {
      elemType = self.moments[1]:elemType(),
      method = self.ioMethod,
   }
   self.couplingMoments = self:allocVectorMoment(self.nMoments)
   self.integratedMoments = DataStruct.DynVector { numComponents = self.nMoments }

   if self.positivity then
      self.fPos = self:allocVectorMoment(self.nMoments)
   end

   -- array with one component per cell to store cflRate in each cell
   self.cflRateByCell = DataStruct.Field {
	onGrid = self.grid,
	numComponents = 1,
	ghost = {1, 1},
   }
   self.cflRateByCell:clear(0.0)
   self.cflRatePtr = self.cflRateByCell:get(1)
   self.cflRateIdxr = self.cflRateByCell:genIndexer()
   self.dt = ffi.new("double[2]")
   self.dtGlobal = ffi.new("double[2]")

   self:createBCs()
end

function FluidSpecies:initDist()
   local project = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      evaluate = self.initFunc,
      projectOnGhosts = true,
   }
   project:advance(0.0, {}, {self.moments[1]})

   if self.positivityRescale or self.positivityDiffuse then
     self.posRescaler:advance(0.0, {self.moments[1]}, {self.moments[1]})
   end
end

function FluidSpecies:rkStepperFields()
   return self.moments
end

-- for RK timestepping 
function FluidSpecies:copyRk(outIdx, aIdx)
   self:rkStepperFields()[outIdx]:copy(self:rkStepperFields()[aIdx])
end
-- for RK timestepping 
function FluidSpecies:combineRk(outIdx, a, aIdx, ...)
   local args = {...} -- package up rest of args as table
   local nFlds = #args/2
   self:rkStepperFields()[outIdx]:combine(a, self:rkStepperFields()[aIdx])
   for i = 1, nFlds do -- accumulate rest of the fields
      self:rkStepperFields()[outIdx]:accumulate(args[2*i-1], self:rkStepperFields()[args[2*i]])
   end	 
   self:applyBc(nil, self:rkStepperFields()[outIdx])
   if self.positivityDiffuse and a<=self.dtGlobal[0] then -- only diffuse when this combine is a forwardEuler 
      self.posRescaler:advance(self.tCurr, {self:rkStepperFields()[outIdx]}, {self:rkStepperFields()[outIdx]})
   end
end

function FluidSpecies:suggestDt()
   return GKYL_MAX_DOUBLE
end

function FluidSpecies:clearCFL()
end

function FluidSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   self.tCurr = tCurr
   local fIn = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      local em = emIn[1]:rkStepperFields()[inIdx]
      if self.positivityRescale then
         self.posRescaler:advance(tCurr, {fIn}, {self.fPos})
         self.solver:advance(tCurr, {self.fPos, em}, {fRhsOut})
      else
         self.solver:advance(tCurr, {fIn, em}, {fRhsOut})
      end
   else
      fRhsOut:clear(0.0) -- no RHS
   end

   -- Perform the collision (diffusion) update.
   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
         c.diffusionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
         c:advance(tCurr, fIn, species, fRhsOut)
         -- the full 'species' list is needed for the cross-species
         -- collisions
      end
   end
end

function FluidSpecies:applyBc(tCurr, fIn)
   local tmStart = Time.clock()
   if self.evolve then
      if self.hasNonPeriodicBc then
         for _, bc in ipairs(self.boundaryConditions) do
            bc:advance(tCurr, {}, {fIn})
         end
      end
      fIn:sync()
   end
   self.bcTime = self.bcTime + Time.clock()-tmStart
end

function FluidSpecies:createDiagnostics()
   -- create updater to compute volume-integrated moments
   self.intMom2Calc = Updater.CartFieldIntegratedQuantCalc {
      onGrid = self.grid,
      basis = self.basis,
      numComponents = self.nMoments,
      quantity = "V"
   }
end

function FluidSpecies:write(tm)
   if self.evolve then
      -- compute integrated diagnostics
      self.intMom2Calc:advance(tm, { self.moments[1] }, { self.integratedMoments })
      
      -- only write stuff if triggered
      if self.diagIoTrigger(tm) then
	 self.momIo:write(
	    self.moments[1], string.format("%s_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame)
         self.integratedMoments:write(
            string.format("%s_intMom_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame)
         self.diagIoFrame = self.diagIoFrame+1
      end
   else
      -- if not evolving species, don't write anything except initial conditions
      if self.diagIoFrame == 0 then
         self.momIo:write(self.moments[1], string.format("%s_%d.bp", self.name, 0), tm, 0)
      end
      self.diagIoFrame = self.diagIoFrame+1
   end
end

function FluidSpecies:writeRestart(tm)
   self.momIo:write(
      self.moments[1], string.format("%s_restart.bp", self.name), tm, self.diagIoFrame)
   self.integratedMoments:write(
      string.format("%s_intMom_restart.bp", self.name), tm, self.diagIoFrame, false)
end

function FluidSpecies:readRestart()
   local tm, fr = self.momIo:read(self.moments[1], string.format("%s_restart.bp", self.name))
   self.diagIoFrame = fr -- reset internal frame counter
   self.integratedMoments:read(string.format("%s_intMom_restart.bp", self.name))   
   
   self:applyBc(tm, 0.0, self.moments[1])
   self.moments[1]:sync() -- must get all ghost-cell data correct

   return tm
end

-- timers
function FluidSpecies:totalSolverTime()
   if self.solver then
      return self.solver.totalTime
   end
   return 0
end
function FluidSpecies:totalBcTime()
   return self.bcTime
end
function FluidSpecies:momCalcTime()
   return 0
end
function FluidSpecies:intMomCalcTime()
   return self.intMom2Calc.totalTime
end

return FluidSpecies
