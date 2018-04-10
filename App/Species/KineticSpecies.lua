-- Gkyl ------------------------------------------------------------------------
--
-- App support code: KineticSpecies object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Basis = require "Basis"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local GaussQuadRules = require "Lib.GaussQuadRules"
local Grid = require "Grid"
local LinearTrigger = require "Lib.LinearTrigger"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
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
local SP_BC_SEE = 4

-- base class for kinetic species
local KineticSpecies = Proto(SpeciesBase)

KineticSpecies.bcAbsorb = SP_BC_ABSORB -- absorb all particles
KineticSpecies.bcOpen = SP_BC_OPEN -- zero gradient
KineticSpecies.bcReflect = SP_BC_REFLECT -- specular reflection
KineticSpecies.bcSEE = SP_BC_SEE -- specular reflection

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function KineticSpecies:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function KineticSpecies:fullInit(appTbl)
   local tbl = self.tbl -- previously store table

   self.name = "name"
   self.cfl =  0.1
   self.charge = tbl.charge and tbl.charge or 1.0
   self.mass = tbl.mass and tbl.mass or 1.0
   self.lower, self.upper = tbl.lower, tbl.upper
   self.cells = tbl.cells
   self.vdim = #self.cells -- velocity dimensions
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default, evolve species
   self.confBasis = nil -- Will be set later

   assert(#self.lower == self.vdim, "'lower' must have " .. self.vdim .. " entries")
   assert(#self.upper == self.vdim, "'upper' must have " .. self.vdim .. " entries")
   self.coordinateMap = tbl.coordinateMap

   self.decompCuts = {}
   -- parallel decomposition stuff
   if tbl.decompCuts then
      assert(self.vdim == #tbl.decompCuts, "decompCuts should have exactly " .. self.vdim .. " entries")
      self.decompCuts = tbl.decompCuts
   else
      -- if not specified, use 1 processor
      for d = 1, self.vdim do self.decompCuts[d] = 1 end
   end

   -- create triggers to write distribution functions and moments
   if tbl.nDistFuncFrame then
      self.distIoTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nDistFuncFrame)
   else
      self.distIoTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   if tbl.nDiagnosticFrame then
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nDiagnosticFrame)
   else
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.distIoFrame = 0 -- frame number for distrubution fucntion
   self.diagIoFrame = 0 -- frame number for diagnostcia

   -- read in which diagnostic moments to compute on output
   self.diagnosticMoments = { }
   if tbl.diagnosticMoments then
      for i, nm in ipairs(tbl.diagnosticMoments) do
	 self.diagnosticMoments[i] = nm
      end
   end

   -- for storing integrated moments
   self.integratedMoments = DataStruct.DynVector { numComponents = 5 }

   -- store initial condition function
   self.initBackgroundFunc = tbl.initBackground
   if self.initBackgroundFunc then 
      self.initFunc = function (t, xn)
         return tbl.init(t, xn, tbl)
      end
   else 
      self.initFunc = tbl.init
   end

   self.fluctuationBCs = xsys.pickBool(tbl.fluctuationBCs, false)
   if self.fluctuationBCs then 
      assert(self.initBackgroundFunc, [[KineticSpecies: must specify an initial
        background distribution with 'initBackground' in order to use fluctuation-only BCs]]) 
   end

   self.hasNonPeriodicBc = false -- to indicate if we have non-periodic BCs
   self.bcx, self.bcy, self.bcz = { }, { }, { }

   -- function to check if BC type is good
   local function isBcGood(bcType)
      if bcType == SP_BC_ABSORB or bcType == SP_BC_OPEN or bcType == SP_BC_REFLECT or bcType == SP_BC_SEE then
         return true
      end
      return false
   end

   -- read in boundary conditions
   if tbl.bcx then
      self.bcx[1], self.bcx[2] = tbl.bcx[1], tbl.bcx[2]
      assert(isBcGood(self.bcx[1]) and isBcGood(self.bcx[2]), "KineticSpecies: Incorrect X BC type specified!")
      self.hasNonPeriodicBc = true
   end
   if tbl.bcy then
      self.bcy[1], self.bcy[2] = tbl.bcy[1], tbl.bcy[2]
      assert(isBcGood(self.bcy[1]) and isBcGood(self.bcy[2]), "KineticSpecies: Incorrect Y BC type specified!")
      self.hasNonPeriodicBc = true
   end
   if tbl.bcz then
      self.bcz[1], self.bcz[2] = tbl.bcz[1], tbl.bcz[2]
      assert(isBcGood(self.bcz[1]) and isBcGood(self.bcz[2]), "KineticSpecies: Incorrect Z BC type specified!")
      self.hasNonPeriodicBc = true
   end
end

function KineticSpecies:createQuadratureData()
   -- Gauss-Legendre quadrature data
   self._N = self.basis:polyOrder() + 1
   local ordinates1D = GaussQuadRules.ordinates[self._N]
   local weights1D = GaussQuadRules.weights[self._N]

   local numDims = self.basis:ndim()
   local l, d = {}, {}
   for d = 1, numDims do l[d], u[d] = 1, self._N end
   quadRange = Range.Range(l, u) -- for looping over quadrature nodes
   self._numOrdinates = quadRange:volume() -- number of ordinates
   -- construct weights and ordinates for integration in multiple dimensions
   self._ordinates = Lin.Mat(self._numOrdinates, numDims)
   self._weights = Lin.Vec(self._numOrdinates)
   local nodeNum = 1
   for idx in quadRange:colMajorIter() do
      self._weights[nodeNum] = 1.0
      for d = 1, numDims do
	 self._weights[nodeNum] =
	    self._weights[nodeNum] * weights1D[idx[d]]
	 self._ordinates[nodeNum][d] = ordinates1D[idx[d]]
      end
      nodeNum = nodeNum + 1
   end
   numBasis = self.basis:numBasis()
   self._basisAtOrdinates = Lin.Mat(self._numOrdinates, numBasis)
   -- pre-compute values of basis functions at quadrature nodes
   for n = 1, self._numOrdinates do
      self.basis:evalBasis(self._ordinates[n],
			   self._basisAtOrdinates[n])
   end
end

function KineticSpecies:getCharge() return self.charge end
function KineticSpecies:getMass() return self.mass end

function KineticSpecies:getNdim()
   return self.ndim
end
function KineticSpecies:getVdim()
   return self.vdim
end
function KineticSpecies:setName(nm)
   self.name = nm
end
function KineticSpecies:setCfl(cfl)
   self.cfl = cfl
end
function KineticSpecies:setIoMethod(ioMethod)
   self.ioMethod = ioMethod
end
function KineticSpecies:setConfBasis(basis)
   self.confBasis = basis
end
function KineticSpecies:setConfGrid(cgrid)
   self.confGrid = cgrid
end

function KineticSpecies:createGrid(cLo, cUp, cCells, cDecompCuts, cPeriodicDirs, cMap)
   self.cdim = #cCells
   self.ndim = self.cdim+self.vdim

   -- create decomposition
   local decompCuts = {}
   for d = 1, self.cdim do table.insert(decompCuts, cDecompCuts[d]) end
   for d = 1, self.vdim do table.insert(decompCuts, self.decompCuts[d]) end
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
   for d = 1, self.vdim do
      table.insert(lower, self.lower[d])
      table.insert(upper, self.upper[d])
      table.insert(cells, self.cells[d])
   end

   local GridConstructor = Grid.RectCart
   local coordinateMap = {} -- table of functions
   -- construct comp -> phys mappings if they exist
   if self.coordinateMap or cMap then
      if cMap and self.coordinateMap then
         for d = 1, self.cdim do
            table.insert(coordinateMap, cMap[d])
         end
         for d = 1, self.vdim do
            table.insert(coordinateMap, self.coordinateMap[d])
         end
      elseif cMap then
         for d = 1, self.cdim do
            table.insert(coordinateMap, cMap[d])
         end
         for d = 1, self.vdim do
            table.insert(coordinateMap, function (z) return z end)
         end
      else
         for d = 1, self.cdim do
            table.insert(coordinateMap, function (z) return z end)
         end
         for d = 1, self.vdim do
            table.insert(coordinateMap, self.coordinateMap[d])
         end
      end
      GridConstructor = Grid.NonUniformRectCart
   end
   self.grid = GridConstructor {
      lower = lower,
      upper = upper,
      cells = cells,
      periodicDirs = cPeriodicDirs,
      decomposition = self.decomp,
      mappings = coordinateMap,
   }
end

function KineticSpecies:createBasis(nm, polyOrder)
   self.basis = createBasis(nm, self.ndim, polyOrder)
end

function KineticSpecies:allocDistf()
   local f = DataStruct.Field {
	onGrid = self.grid,
	numComponents = self.basis:numBasis(),
	ghost = {1, 1},
   }
   return f
end
function KineticSpecies:allocMoment()
   local m = DataStruct.Field {
	onGrid = self.confGrid,
	numComponents = self.confBasis:numBasis(),
	ghost = {1, 1},
   }
   return m
end
function KineticSpecies:allocVectorMoment(dim)
   local m = DataStruct.Field {
	onGrid = self.confGrid,
	numComponents = self.confBasis:numBasis()*dim,
	ghost = {1, 1},
   }
   return m
end

function KineticSpecies:createBCs()
   -- function to construct a BC updater
   local function makeBcUpdater(dir, edge, bcList, skinLoop)
      return Updater.Bc {
	 onGrid = self.grid,
	 boundaryConditions = bcList,
	 dir = dir,
	 edge = edge,
	 skinLoop = skinLoop,
	 cdim = self.cdim,
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
   local function bcReflect(dir, tm, xc, fIn, fOut)
      -- get handle to function to compute basis functions at specified coordinates
      local bName = ""
      if self.basis:id() == "serendipity" then
	 bName = "Serendip"
      else
	 bName = "MaxOrder"
      end
      local fName = "Basis._data.Modal" .. bName .. "BasisReflect" .. self.cdim .. "x" .. self.vdim .. "v"
      local _m = require(fName)
      local polyOrder = self.basis:polyOrder()
      _m[polyOrder](dir, fIn, fOut) -- function to flip sign of both configuration and velocity component
   end

   local function bcSEE(dir, tm, xcIn, xcOut, fIn, fOut)

   end

   -- functions to make life easier while reading in BCs to apply
   self.boundaryConditions = { } -- list of Bcs to apply
   local function appendBoundaryConditions(dir, edge, bcType)
      if bcType == SP_BC_ABSORB then
	 table.insert(self.boundaryConditions, makeBcUpdater(dir, edge, { bcAbsorb }, "pointwise"))
      elseif bcType == SP_BC_OPEN then
	 table.insert(self.boundaryConditions, makeBcUpdater(dir, edge, { bcOpen }, "pointwise"))
      elseif bcType == SP_BC_REFLECT then
	 table.insert(self.boundaryConditions, makeBcUpdater(dir, edge, { bcReflect }, "flip"))
      elseif bcType == SP_BC_SEE then
	 table.insert(self.boundaryConditions, makeBcUpdater(dir, edge, { bcSEE }, "integrate"))
      else
	 assert(false, "KineticSpecies: Unsupported BC type!")
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

function KineticSpecies:alloc(nRkDup)
   -- allocate fields needed in RK update
   self.distf = {}
   for i = 1, nRkDup do
      self.distf[i] = self:allocDistf()
   end
   -- create Adios object for field I/O
   self.distIo = AdiosCartFieldIo {
      elemType = self.distf[1]:elemType(),
      method = self.ioMethod,
   }

   -- background (or initial) distribution
   self.f0 = self:allocDistf()

   self:createBCs()
end

function KineticSpecies:initDist()
   local syncPeriodicDirs = true
   if self.fluctuationBCs then syncPeriodicDirs = false end

   local project = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      evaluate = self.initFunc,
      projectOnGhosts = true
   }
   project:advance(0.0, 0.0, {}, {self.distf[1]})

   if self.initBackgroundFunc then
      local projectBackground = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = self.initBackgroundFunc,
         projectOnGhosts = true
      }
      projectBackground:advance(0.0, 0.0, {}, {self.f0})
      self.f0:sync(syncPeriodicDirs)
   else
      -- if no specified background distribution, use initial
      -- distribution function as background
      self.f0:copy(self.distf[1])
      self:applyBc(0.0, 0.0, self.f0)
   end
end

function KineticSpecies:rkStepperFields()
   return self.distf
end

function KineticSpecies:applyBc(tCurr, dt, fIn)
   -- fIn is total distribution function

   local syncPeriodicDirsTrue = true

   if self.fluctuationBCs then
     -- if fluctuation-only BCs, subtract off background before applying BCs
     fIn:accumulate(-1.0, self.f0)
   end

   -- apply non-periodic BCs (to only fluctuations if fluctuation BCs)
   if self.hasNonPeriodicBc then
      for _, bc in ipairs(self.boundaryConditions) do
	 bc:advance(tCurr, dt, {}, {fIn})
      end
   end

   -- apply periodic BCs (to only fluctuations if fluctuation BCs)
   fIn:sync(syncPeriodicDirsTrue)

   if self.fluctuationBCs then
     -- put back together total distribution
     fIn:accumulate(1.0, self.f0)

     -- update ghosts in total distribution, without enforcing periodicity
     fIn:sync(not syncPeriodicDirsTrue)
   end

end

function KineticSpecies:createDiagnostics()
   -- create updater to compute volume-integrated moments
   self.intMomentCalc = Updater.DistFuncIntegratedMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
   }
end

function KineticSpecies:calcDiagnosticMoments()
   local numMoms = #self.diagnosticMoments
   for i = 1, numMoms do
      self.diagnosticMomentUpdaters[i]:advance(
	 0.0, 0.0, {self.distf[1]}, {self.diagnosticMomentFields[i]})
   end
end

function KineticSpecies:write(tm)
   if self.evolve then
      -- compute integrated diagnostics
      self.intMomentCalc:advance(tm, 0.0, { self.distf[1] }, { self.integratedMoments })

      -- only write stuff if triggered
      if self.distIoTrigger(tm) then
	 self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, self.distIoFrame), tm)
         if self.initBackgroundFunc then
            self.distIo:write(self.f0, string.format("%s_f0_%d.bp", self.name, self.distIoFrame), tm)
            self.distf[1]:accumulate(-1, self.f0)
            self.distIo:write(self.distf[1], string.format("%s_f1_%d.bp", self.name, self.distIoFrame), tm)
            self.distf[1]:accumulate(1, self.f0)
            self.distIo:write(self.distf[1], string.format("%s_ftot_%d.bp", self.name, self.distIoFrame), tm)
         end
	 self.distIoFrame = self.distIoFrame+1
      end


      if self.diagIoTrigger(tm) then
         -- compute moments and write them out
         self:calcDiagnosticMoments()
         for i, mom in ipairs(self.diagnosticMoments) do
            -- should one use AdiosIo object for this?
            self.diagnosticMomentFields[i]:write(
               string.format("%s_%s_%d.bp", self.name, mom, self.diagIoFrame), tm)
         end
         self.integratedMoments:write(
            string.format("%s_intMom_%d.bp", self.name, self.diagIoFrame), tm)
         self.diagIoFrame = self.diagIoFrame+1
      end
   else
      -- if not evolving species, don't write anything except initial conditions
      if self.distIoFrame == 0 then
	 self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, 0), tm)

	 -- compute moments and write them out
	 self:calcDiagnosticMoments()
	 for i, mom in ipairs(self.diagnosticMoments) do
	    -- should one use AdiosIo object for this?
	    self.diagnosticMomentFields[i]:write(string.format("%s_%s_%d.bp", self.name, mom, 0), tm)
	 end
      end
      self.distIoFrame = self.distIoFrame+1
   end
end

-- timers
function KineticSpecies:totalSolverTime()
   return self.solver.totalTime
end
function KineticSpecies:totalBcTime()
   local tm = 0.0
   for _, bc in ipairs(self.boundaryConditions) do
      tm = tm + bc.totalTime
   end
   return tm
end
function KineticSpecies:intMomCalcTime()
   return self.intMomentCalc.totalTime
end

return KineticSpecies
