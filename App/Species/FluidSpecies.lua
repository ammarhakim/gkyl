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
local Updater = require "Updater"
local xsys = require "xsys"
local SpeciesBase = require "App.Species.SpeciesBase"

-- function to create basis functions
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   end
end

-- add constants to object indicate various supported boundary conditions
local SP_BC_ABSORB = 1
local SP_BC_OPEN = 2
local SP_BC_REFLECT = 3

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
   self.confBasis = nil -- Will be set later

   -- create triggers to write diagnostics
   if tbl.nDiagnosticFrame then
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nDiagnosticFrame)
   else
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.diagIoFrame = 0 -- frame number for diagnostics

   -- for storing integrated moments
   self.integratedMoments = DataStruct.DynVector { numComponents = 1 }

   -- store initial condition function
   self.initFunc = tbl.init

   -- default to a single moment
   self.nMoments = 1

   self.hasNonPeriodicBc = false -- to indicate if we have non-periodic BCs
   self.bcx, self.bcy, self.bcz = { }, { }, { }

   -- function to check if BC type is good
   local function isBcGood(bcType)
      if bcType == SP_BC_ABSORB or bcType == SP_BC_OPEN or bcType == SP_BC_REFLECT then
         return true
      end
      return false
   end
   
   -- read in boundary conditions
   if tbl.bcx then
      self.bcx[1], self.bcx[2] = tbl.bcx[1], tbl.bcx[2]
      assert(isBcGood(self.bcx[1]) and isBcGood(self.bcx[2]), "FluidSpecies: Incorrect X BC type specified!")
      self.hasNonPeriodicBc = true
   end
   if tbl.bcy then
      self.bcy[1], self.bcy[2] = tbl.bcy[1], tbl.bcy[2]
      assert(isBcGood(self.bcy[1]) and isBcGood(self.bcy[2]), "FluidSpecies: Incorrect Y BC type specified!")
      self.hasNonPeriodicBc = true
   end
   if tbl.bcz then
      self.bcz[1], self.bcz[2] = tbl.bcz[1], tbl.bcz[2]
      assert(isBcGood(self.bcz[1]) and isBcGood(self.bcz[2]), "FluidSpecies: Incorrect Z BC type specified!")
      self.hasNonPeriodicBc = true
   end
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
   self.confBasis = basis
end
function FluidSpecies:setConfGrid(cgrid)
   self.confGrid = cgrid
end

function FluidSpecies:createGrid(cLo, cUp, cCells, cDecompCuts, cPeriodicDirs)
   self.cdim = #cCells
   self.ndim = self.cdim

   -- create decomposition
   local decompCuts = {}
   for d = 1, self.cdim do table.insert(decompCuts, cDecompCuts[d]) end
   self.decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      shared = false,
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
	ghost = {1, 1}
   }
   return m
end
function FluidSpecies:allocVectorMoment(dim)
   local m = DataStruct.Field {
	onGrid = self.confGrid,
	numComponents = self.confBasis:numBasis()*dim,
	ghost = {1, 1}
   }
   return m
end

function FluidSpecies:createBCs()
   -- function to construct a BC updater
   local function makeBcUpdater(dir, edge, bcList)
      return Updater.Bc {
	 onGrid = self.grid,
	 boundaryConditions = bcList,
	 dir = dir,
	 edge = edge,
      }
   end

   -- various functions to apply BCs of different types
   local function bcAbsorb(dir, tm, xc, fIn, fOut)
      for i = 1, self.basis:numBasis() do
	 fOut[i] = 0.0
      end
   end
   local function bcOpen(dir, tm, xc, fIn, fOut)
      self.basis:flipSign(dir, fIn, fOut)
   end

   -- functions to make life easier while reading in BCs to apply
   self.boundaryConditions = { } -- list of Bcs to apply
   local function appendBoundaryConditions(dir, edge, bcType)
      if bcType == SP_BC_ABSORB then
	 table.insert(self.boundaryConditions, makeBcUpdater(dir, edge, { bcAbsorb }))
      elseif bcType == SP_BC_OPEN then
	 table.insert(self.boundaryConditions, makeBcUpdater(dir, edge, { bcOpen }))
      else
	 assert(false, "FluidSpecies: Unsupported BC type!")
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

   self:createBCs()
end

function FluidSpecies:initDist()
   local project = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      evaluate = self.initFunc
   }
   project:advance(0.0, 0.0, {}, {self.moments[1]})
   self:applyBc(0.0, 0.0, self.moments[1])
end

function FluidSpecies:rkStepperFields()
   return self.moments
end

function FluidSpecies:forwardEuler(tCurr, dt, fIn, emIn, fOut)
   if self.evolve then
      return self.solver:advance(tCurr, dt, {fIn, emIn}, {fOut})
   else
      fOut:copy(fIn) -- just copy stuff over
      return true, GKYL_MAX_DOUBLE
   end
end

function FluidSpecies:applyBc(tCurr, dt, fIn)
   if self.hasNonPeriodicBc then
      for _, bc in ipairs(self.boundaryConditions) do
	 bc:advance(tCurr, dt, {}, {fIn})
      end
   end
   fIn:sync()
end

function FluidSpecies:createDiagnostics()
   -- create updater to compute volume-integrated moments
   self.intMom2Calc = Updater.CartFieldIntegratedQuantCalc {
      onGrid = self.grid,
      basis = self.basis,
      numComponents = self.nMoments,
      quantity = "V2"
   }
end

function FluidSpecies:write(tm)
   if self.evolve then
      -- compute integrated diagnostics
      self.intMom2Calc:advance(tm, 0.0, { self.moments[1] }, { self.integratedMoments })
      
      -- only write stuff if triggered
      if self.diagIoTrigger(tm) then
	 self.momIo:write(self.moments[1], string.format("%s_%d.bp", self.name, self.diagIoFrame), tm)
         self.integratedMoments:write(
            string.format("%s_intMom_%d.bp", self.name, self.diagIoFrame), tm)
         self.diagIoFrame = self.diagIoFrame+1
      end
   else
      -- if not evolving species, don't write anything except initial conditions
      if self.diagIoFrame == 0 then
         self.momIo:write(self.moments[1], string.format("%s_%d.bp", self.name, 0), tm)
      end
      self.diagIoFrame = self.diagIoFrame+1
   end
end

-- timers
function FluidSpecies:totalSolverTime()
   return self.solver.totalTime
end
function FluidSpecies:totalBcTime()
   local tm = 0.0
   for _, bc in ipairs(self.boundaryConditions) do
      tm = tm + bc.totalTime
   end
   return tm
end
function FluidSpecies:momCalcTime()
   return 0
end
function FluidSpecies:intMomCalcTime()
   return self.intMom2Calc.totalTime
end

return FluidSpecies
