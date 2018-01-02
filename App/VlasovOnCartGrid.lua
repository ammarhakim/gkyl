-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov solver on a Cartesian grid. Works in arbitrary CDIM/VDIM
-- (VDIM>=CDIM) with either Maxwell, Poisson or specified EM fields.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Basis = require "Basis"
local BoundaryCondition = require "Updater.BoundaryCondition"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local Time = require "Lib.Time"
local Updater = require "Updater"
local date = require "Lib.date"
local xsys = require "xsys"

-- function to create basis functions
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   end   
end

-- For returning module table
local M = {}

-- Class to store species-specific info
local Species = {}
function Species:new(tbl)
   local self = setmetatable({}, Species)
   self.type = "species" -- to identify objects of this (Species) type

   self.name = "name"
   self.cfl =  0.1
   self.charge, self.mass = tbl.charge, tbl.mass
   self.qbym = self.charge/self.mass
   self.lower, self.upper = tbl.lower, tbl.upper
   self.cells = tbl.cells
   self.vdim = #self.cells -- velocity dimensions
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default, evolve species
   self.confBasis = nil -- Will be set later

   assert(#self.lower == self.vdim, "'lower' must have " .. self.vdim .. " entries")
   assert(#self.upper == self.vdim, "'upper' must have " .. self.vdim .. " entries")

   self.decompCuts = {}
   -- parallel decomposition stuff
   if tbl.decompCuts then
      assert(self.vdim == #tbl.decompCuts, "decompCuts should have exactly " .. self.vdim .. " entries")
      self.decompCuts = tbl.decompCuts
   else
      -- if not specified, use 1 processor
      for d = 1, self.vdim do self.decompCuts[d] = 1 end
   end

   -- store initial condition function
   self.initFunc = tbl.init
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(Species, { __call = function (self, o) return self.new(self, o) end })

-- methods for species object
Species.__index = {
   ndim = function(self)
      return self.vdim
   end,
   setName = function(self, nm)
      self.name = nm
   end,
   setCfl = function(self, cfl)
      self.cfl = cfl
   end,
   setIoMethod = function(self, ioMethod)
      self.ioMethod = ioMethod
   end,
   setConfBasis = function (self, basis)
      self.confBasis = basis
   end,   
   createGrid = function(self, cLo, cUp, cCells, cDecompCuts, cPeriodicDirs)
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
      self.grid = Grid.RectCart {
	 lower = lower,
	 upper = upper,
	 cells = cells,
	 periodicDirs = cPeriodicDirs,
	 decomposition = self.decomp,
      }
   end,
   createBasis = function(self, nm, polyOrder)
      self.basis = createBasis(nm, self.ndim, polyOrder)
   end,
   alloc = function(self)
      -- allocate fields needed in RK update
      self.distf = {}
      for i = 1, 3 do
	 self.distf[i] = DataStruct.Field {
	    onGrid = self.grid,
	    numComponents = self.basis:numBasis(),
	    ghost = {1, 1}
	 }
      end
      -- create Adios object for field I/O
      self.distIo = AdiosCartFieldIo {
	 elemType = self.distf[1]:elemType(),
	 method = self.ioMethod,
      }
   end,
   createSolver = function (self, hasE, hasB)
      -- create updater to advance solution by one time-step
      self.vlasovSlvr = Updater.VlasovDisCont {
	 onGrid = self.grid,
	 phaseBasis = self.basis,
	 confBasis = self.confBasis,
	 charge = self.charge,
	 mass = self.mass,
	 cfl = self.cfl,
	 hasElectricField = hasE,
	 hasMagneticField = hasB,
      }
   end,
   initDist = function(self)
      local project = Updater.ProjectOnBasis {
	 onGrid = self.grid,
	 basis = self.basis,
	 evaluate = self.initFunc
      }
      project:advance(0.0, 0.0, {}, {self.distf[1]})
   end,
   write = function(self, frame, tm)
      if self.evolve then
	 self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, frame), tm)
      else
	 -- if not evolving species, don't write anything except initial conditions
	 if frame == 0 then
	    self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, frame), tm)
	 end
      end
   end,
   rkStepperFields = function(self)
      return self.distf[1], self.distf[2], self.distf[3]
   end,
   forwardEuler = function(self, tCurr, dt, fIn, fOut)
      if self.evolve then
	 return self.vlasovSlvr:advance(tCurr, dt, {fIn}, {fOut})
      else
	 fOut:copy(fIn) -- just copy stuff over
	 return true, GKYL_MAX_DOUBLE
      end
   end,
   applyBc = function(self, tCurr, dt, fIn)
      fIn:sync()
   end,
   totalSolverTime = function(self)
      return self.vlasovSlvr.totalTime
   end,
   volTime = function(self)
      return self.vlasovSlvr:volTime()
   end,
   surfTime = function(self)
      return self.vlasovSlvr:surfTime()
   end,
}

-- add to module table
M.Species = Species

-- Class to store Maxwell field info
local EmField = {}
function EmField:new(tbl)
   local self = setmetatable({}, EmField)
   self.type = "field" -- to identify objects of this type

   self.epsilon0 = tbl.epsilon0
   self.mu0 = tbl.mu0
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default evolve field

   -- store initial condition function
   self.initFunc = tbl.init

   self.lightSpeed = 1/math.sqrt(self.epsilon0*self.mu0)

   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(EmField, { __call = function (self, o) return self.new(self, o) end })
-- methods for EM field object
EmField.__index = {
   hasEB = function (self)
      return true, false
   end,
   setCfl = function(self, cfl)
      self.cfl = cfl
   end,
   setIoMethod = function(self, ioMethod)
      self.ioMethod = ioMethod
   end,
   setBasis = function (self, basis)
      self.basis = basis
   end,
   setGrid = function(self, grid)
      self.grid = grid
   end,
   alloc = function(self)
      -- allocate fields needed in RK update
      self.em = {}
      for i = 1, 3 do
	 self.em[i] = DataStruct.Field {
	    onGrid = self.grid,
	    numComponents = self.basis:numBasis(),
	    ghost = {1, 1}
	 }
      end
      
      -- create Adios object for field I/O
      self.fieldIo = AdiosCartFieldIo {
	 elemType = self.em[1]:elemType(),
	 method = self.ioMethod,
      }
   end,
   initField = function(self)
      local project = Updater.ProjectOnBasis {
	 onGrid = self.grid,
	 basis = self.basis,
	 evaluate = self.initFunc
      }
      project:advance(0.0, 0.0, {}, {self.em[1]})
   end,
   write = function(self, frame, tm)
      if self.evolve then
	 self.fieldIo:write(self.em[1], string.format("field_%d.bp", frame), tm)
      else
	 -- if not evolving species, don't write anything except initial conditions
	 if frame == 0 then
	    self.fieldIo:write(self.em[1], string.format("field_%d.bp", frame), tm)
	 end
      end
   end,
   rkStepperFields = function(self)
      return self.em[1], self.em[2], self.em[3]
   end,
   forwardEuler = function(self, tCurr, dt, emIn, emOut)
      if self.evolve then
	 assert(false, "EmField:forwardEuler. NYI!")
      else
	 emOut:copy(emIn) -- just copy stuff over
	 return true, GKYL_MAX_DOUBLE
      end
   end,
   applyBc = function(self, tCurr, dt, emIn)
      emIn:sync()
   end,
   totalSolverTime = function(self)
      return 0.0
   end,
   volTime = function(self)
      return 0.0
   end,
   surfTime = function(self)
      return 0.0
   end,   
}

-- add to module table
M.EmField = EmField

-- function to check "type" of obj.
local function isType(obj, typeString)
   if type(obj) == typeString then
      return true
   elseif type(obj) == "table" then
      if obj.type == typeString then
	 return true
      end
   end
   return false
end

-- top-level method to build application "run" method
local function buildApplication(self, tbl)
   -- create logger
   local log = Logger {
      logToFile = xsys.pickBool(tbl.logToFile, true)
   }

   log(date(false):fmt()) -- time-stamp for sim start

   -- function to warn user about default values
   local function warnDefault(varVal, varNm, default)
      if varVal then return varVal end
      log(string.format(" ** WARNING: %s not specified, assuming %s", varNm, tostring(default)))
      return default
   end

   log("Initializing VlasovOnCartGrid simulation ...")
   local tmStart = Time.clock()

   local cdim = #tbl.lower -- configuration space dimensions
   assert(cdim == #tbl.upper, "upper should have exactly " .. cdim .. " entries")
   assert(cdim == #tbl.cells, "cells should have exactly " .. cdim .. " entries")

   -- basis function name
   local basisNm = warnDefault(tbl.basis, "basis", "serendipity")
   if basisNm ~= "serendipity" and basisNm ~= "maximal-order" then
      assert(false, "Incorrect basis type " .. basisNm .. " specified")
   end

   -- polynomial order
   local polyOrder = tbl.polyOrder

   -- create basis function for configuration space
   local confBasis = createBasis(basisNm, cdim, polyOrder)

   -- I/O method
   local ioMethod = tbl.ioMethod and tbl.ioMethod or "MPI"
   if ioMethod ~= "POSIX" and ioMethod ~= "MPI" then
      assert(false, "ioMethod must be one of 'MPI' or 'POSIX'. Provided '" .. ioMethod .. "' instead")
   end

   -- time-stepper
   local timeStepperNm = warnDefault(tbl.timeStepper, "timeStepper", "rk3")
   if timeStepperNm ~= "rk2" and timeStepperNm ~= "rk3" and timeStepperNm ~= "rk3s4" then
      assert(false, "Incorrect timeStepper type " .. timeStepperNm .. " specified")
   end

   -- CFL fractions for various steppers
   local stepperCFLFracs = { rk2 = 1.0, rk3 = 1.0, rk3s4 = 2.0 }

   local cflFrac = tbl.cflFrac
   -- Compute CFL fraction if not specified
   if  cflFrac == nil then
      cflFrac = stepperCFLFracs[timeStepperNm]
   end

   -- parallel decomposition stuff
   local decompCuts = tbl.decompCuts
   if tbl.decompCuts then
      assert(cdim == #tbl.decompCuts, "decompCuts should have exactly " .. cdim .. " entries")
   else
      -- if not specified, use 1 processor
      decompCuts = {}
      for d = 1, cdim do decompCuts[d] = 1 end
   end
   local useShared = xsys.pickBool(tbl.useShared, false)

   -- extract periodic directions
   local periodicDirs = {}
   if tbl.periodicDirs then
      for i, d in ipairs(tbl.periodicDirs) do
	 if d<1 and d>3 then
	    assert(false, "Directions in periodicDirs table should be 1 (for X), 2 (for Y) or 3 (for Z)")
	 end
	 periodicDirs[i] = d
      end
   end

   -- read in information about each species
   local species = {}
   for nm, val in pairs(tbl) do
      if isType(val, "species") then
	 species[nm] = val
	 species[nm]:setName(nm)
	 species[nm]:setIoMethod(ioMethod)
      end
   end

   local cflMin = GKYL_MAX_DOUBLE
   -- compute CFL numbers
   for _, s in pairs(species) do
      local ndim = cdim+s:ndim()
      local myCfl = tbl.cfl and tbl.cfl or cflFrac/(ndim*(2*polyOrder+1))
      cflMin = math.min(cflMin, myCfl)
      s:setCfl(cflMin)
   end
   log(string.format("Using CFL number %g", cflMin))

   -- setup each species
   for _, s in pairs(species) do
      s:createGrid(tbl.lower, tbl.upper, tbl.cells, decompCuts, periodicDirs)
      s:setConfBasis(confBasis)
      s:createBasis(basisNm, polyOrder)
      s:alloc()
   end

   -- setup configuration space grid
   local grid = Grid.RectCart {
      lower = tbl.lower,
      upper = tbl.upper,
      cells = tbl.cells,
      periodicDirs = periodicDirs,
      -- NEED TO CREATE DECOMPOSITION
   }

   -- setup information about fields
   local field = tbl.field -- must always be called "field"
   field:setIoMethod(ioMethod)
   field:setBasis(confBasis)
   field:setGrid(grid)
   do
      local myCfl = tbl.cfl and tbl.cfl or cflFrac/(cdim*(2*polyOrder+1))
      field:setCfl(myCfl)
   end
   -- allocate field data
   field:alloc()

   -- initialize species solvers
   for _, s in pairs(species) do
      local hasE, hasB = field:hasEB()
      s:createSolver(hasE, hasB)
   end   

   -- initialize species distributions and write them out
   for _, s in pairs(species) do
      s:initDist()
   end
   -- initialize fields
   field:initField()

   -- function to write data to file
   local function writeData(frame, tCurr)
      for _, s in pairs(species) do s:write(frame, tCurr) end
      field:write(frame, tCurr)
   end

   writeData(0, 0.0) -- write initial conditions
   
   -- store fields used in RK time-stepping for each species
   local speciesRkFields = { }
   for nm, s in pairs(species) do
      speciesRkFields[nm] = { s:rkStepperFields() } -- package them up as a table
   end

   -- function to take a single forward-euler time-step
   local function fowardEuler(tCurr, dt, inIdx, outIdx)
      local status, dtSuggested = true, GKYL_MAX_DOUBLE
      for nm, s in pairs(species) do
	 local myStatus, myDtSuggested = s:forwardEuler(tCurr, dt, speciesRkFields[nm][inIdx], speciesRkFields[nm][outIdx])
	 status = status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
	 s:applyBc(tCurr, dt, speciesRkFields[nm][outIdx])
      end
      return status, dtSuggested
   end

   -- function to increment fields
   local function increment(a, aIdx, b, bIdx, outIdx)
      for nm, s in pairs(species) do
	 speciesRkFields[nm][outIdx]:combine(a, speciesRkFields[nm][aIdx], b, speciesRkFields[nm][bIdx])
      end
   end

   local timeSteppers = {}

   -- function to advance solution using SSP-RK2 scheme
   function timeSteppers.rk2(tCurr, dt)
      local status, dtSuggested
      -- RK stage 1
      status, dtSuggested = fowardEuler(tCurr, dt, 1, 2)
      if status == false then return status, dtSuggested end

      -- RK stage 2
      status, dtSuggested = fowardEuler(tCurr, dt, 2, 3)
      if status == false then return status, dtSuggested end

      -- final update
      increment(0.5, 1, 0.5, 3, 2)
      for nm, s in pairs(species) do
	 speciesRkFields[nm][1]:copy(speciesRkFields[nm][2])
      end
      return status, dtSuggested
   end

   -- function to advance solution using SSP-RK3 scheme
   function timeSteppers.rk3(tCurr, dt)

   end

   -- function to advance solution using 4-stage SSP-RK3 scheme
   function timeSteppers.rk3s4(tCurr, dt)
      assert(false, "timeSteppers.rk3s4: NYI!")
   end   

   local tmEnd = Time.clock()
   log(string.format("Initializing completed in %g sec\n", tmEnd-tmStart))

   -- return function that runs main simulation loop   
   return function(self)
      log("Starting main loop of HyperEqnOnCartGrid simulation ...")
      local tStart, tEnd, nFrame = 0, tbl.tEnd, tbl.nFrame
      local initDt =  tbl.suggestedDt and tbl.suggestedDt or tEnd-tStart -- initial time-step
      local frame = 1
      local tFrame = (tEnd-tStart)/nFrame
      local nextIOt = tFrame
      local step = 1
      local tCurr = tStart
      local myDt = initDt

      local tmSimStart = Time.clock()
      -- main simulation loop
      while true do
	 -- if needed adjust dt to hit tEnd exactly
	 if tCurr+myDt > tEnd then myDt = tEnd-tCurr end
	 log (string.format(" Taking step %5d at time %6g with dt %g", step, tCurr, myDt))
	 local status, dtSuggested = timeSteppers[timeStepperNm](tCurr, myDt) -- take a time-step

	 if not status then
	    -- updater failed, time-step too large
	    log (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	    myDt = dtSuggested
	 else
	    -- write out data if needed
	    if tCurr+myDt > nextIOt or tCurr+myDt >= tEnd then
	       log (string.format(" Writing data at time %g (frame %d) ...\n", tCurr+myDt, frame))
	       writeData(frame, tCurr+myDt)
	       frame = frame + 1
	       nextIOt = nextIOt + tFrame
	       step = 0
	    end
	    tCurr = tCurr + myDt
	    myDt = dtSuggested
	    step = step + 1
	    if (tCurr >= tEnd) then
	       break
	    end
	 end 
      end -- end of time-step loop
      local tmSimEnd = Time.clock()

      local tmVlasovSlvr = 0.0 -- total time in Vlasov solver
      for _, s in pairs(species) do
	 tmVlasovSlvr = tmVlasovSlvr+s:totalSolverTime()
      end

      local tmVlasovVol, tmVlasovSurf = 0.0, 0.0
      for _, s in pairs(species) do
	 tmVlasovVol = tmVlasovVol + s:volTime()
	 tmVlasovSurf = tmVlasovSurf + s:surfTime()
      end

      log(string.format("Vlasov solver took %g sec", tmVlasovSlvr))
      log(string.format("  [Vol updates took %g sec. Surf updates took %g sec]", tmVlasovVol, tmVlasovSurf))
      log(string.format("Main loop completed in %g sec", tmSimEnd-tmSimStart))
      log(date(false):fmt()) -- time-stamp for sim end
   end
end

-- VlasovOnCartGrid application object
local App = {}
-- constructor
function App:new(tbl)
   local self = setmetatable({}, App)
   self._runApplication = buildApplication(self, tbl)
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(App, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
App.__index = {
   run = function (self)
      return self:_runApplication()
   end
}

-- add to table
M.App = App

return M
