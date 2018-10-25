-- Gkyl ------------------------------------------------------------------------
--
-- App support code: FluidSpecies object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Basis = require "Basis"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
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

   self.name = "name"
   self.cfl =  0.1
   self.charge = tbl.charge and tbl.charge or 1.0
   self.mass = tbl.mass and tbl.mass or 1.0
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default, evolve species
   self.basis = nil -- Will be set later

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

   self.bcTime = 0.0 -- timer for BCs

   self.useShared = xsys.pickBool(appTbl.useShared, false)
   self.positivity = xsys.pickBool(tbl.applyPositivity, false)
   self.deltaF = xsys.pickBool(appTbl.deltaF, false)
end

function FluidSpecies:getCharge() return self.charge end
function FluidSpecies:getMass() return self.mass end

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
end
function FluidSpecies:setIoMethod(ioMethod)
   self.ioMethod = ioMethod
end
function FluidSpecies:setConfBasis(basis)
   self.basis = basis
end
function FluidSpecies:setConfGrid(cgrid)
   self.grid = cgrid
   self.ndim = self.grid:ndim()
end

function FluidSpecies:createGrid(cgrid)
   self.cdim = #cCells
   self.ndim = self.cdim

   -- create decomposition
   local decompCuts = {}
   for d = 1, self.cdim do table.insert(decompCuts, cgrid:cuts(d)) end
   self.decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      useShared = self.useShared,
   }

   -- create computational domain
   local lower, upper, cells = {}, {}, {}
   for d = 1, self.cdim do
      table.insert(lower, cgrid:lower(d))
      table.insert(upper, cgrid:upper(d))
      table.insert(cells, cgrid:numCells(d))
   end
   self.grid = Grid.RectCart {
      lower = lower,
      upper = upper,
      cells = cells,
      periodicDirs = cgrid:getPeriodicDirs(),
      decomposition = self.decomp,
   }
end

function FluidSpecies:createBasis(nm, polyOrder)
end

function FluidSpecies:allocMoment()
   local m = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {self.nGhost, self.nGhost}
   }
   return m
end
function FluidSpecies:allocVectorMoment(dim)
   local m = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis()*dim,
      ghost = {self.nGhost, self.nGhost}
   }
   return m
end

function FluidSpecies:allocMomCouplingFields()
   return {self:allocVectorMoment(self.nMoments)}
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
   self:createBCs()
end

function FluidSpecies:initDist()
   local project = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      evaluate = self.initFunc,
      projectOnGhosts = true,
   }
   project:advance(0.0, 0.0, {}, {self.moments[1]})
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
end

function FluidSpecies:forwardEuler(tCurr, dt, species, emIn, inIdx, outIdx)
   local fIn = self:rkStepperFields()[inIdx]
   local fOut = self:rkStepperFields()[outIdx]

   if self.evolve then
      local em = emIn[1]:rkStepperFields()[inIdx]
      local myStatus, myDt = self.solver:advance(tCurr, dt, {fIn, em}, {fOut})

      if self.positivity then self.positivityRescale:advance(tCurr, dt, {fOut}, {fOut}) end

      -- apply BCs
      self:applyBc(tCurr, dt, fOut)
      return myStatus, myDt
   else
      fOut:copy(fIn) -- just copy stuff over
      return true, GKYL_MAX_DOUBLE
   end
end

function FluidSpecies:applyBc(tCurr, dt, fIn)
   local tmStart = Time.clock()
   if self.evolve then
      if self.hasNonPeriodicBc then
         for _, bc in ipairs(self.boundaryConditions) do
            bc:advance(tCurr, dt, {}, {fIn})
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
      self.intMom2Calc:advance(tm, 0.0, { self.moments[1] }, { self.integratedMoments })
      
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
