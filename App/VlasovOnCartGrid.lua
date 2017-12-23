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

local function showTbl(nm, tbl)
   io.write(string.format("--- %s ---\n", nm))
   for _, v in ipairs(tbl) do
      io.write(v); io.write(" ")
   end
   io.write("\n---\n")
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
      if nm == "serendipity" then
	 self.basis = Basis.CartModalSerendipity { ndim = self.ndim, polyOrder = polyOrder }
	 self.confBasis = Basis.CartModalSerendipity { ndim = self.cdim, polyOrder = polyOrder }
      elseif nm == "maximal-order" then
	 self.basis = Basis.CartModalMaxOrder { ndim = self.ndim, polyOrder = polyOrder }
	 self.confBasis = Basis.CartModalMaxOrder { ndim = self.cdim, polyOrder = polyOrder }
      end
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
      -- create updater to advance solution by one time-step
      self.vlasovSlvr = Updater.VlasovDisCont {
	 onGrid = self.grid,
	 phaseBasis = self.basis,
	 confBasis = self.confBasis,
	 charge = self.charge,
	 mass = self.mass,
	 cfl = self.cfl
      }
      
      -- create Adios object for field I/O
      self.distIo = AdiosCartFieldIo {
	 elemType = self.distf[1]:elemType()
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
      self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, frame), tm)
   end,
   rkStepperFields = function(self)
      return self.distf[1], self.distf[2], self.distf[3]
   end,
   forwardEuler = function(self, tCurr, dt, fIn, fOut)
      return self.vlasovSlvr:advance(tCurr, dt, {fIn}, {fOut})
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
      logToFile = tbl.logToFile and tbl.logToFile or false
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

   -- time-stepper
   local timeStepperNm = warnDefault(tbl.timeStepper, "timeStepper", "rk3")
   if timeStepperNm ~= "rk1" and timeStepperNm ~= "rk2" and timeStepperNm ~= "rk3" then
      assert(false, "Incorrect timeStepper type " .. timeStepperNm .. " specified")
   end

   -- CFL fractions for various steppers
   local stepperCFLFracs = { rk2 = 1.0, rk3 = 1.0, rk4 = 2.0 }

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
   local useShared = tbl.useShared and tbl.useShared or false

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
      end
   end

   local cflMin = GKYL_MAX_DOUBLE
   -- compute CFL numbers
   for _, s in pairs(species) do
      local ndim = cdim+s:ndim()
      local myCfl = tbl.cfl and tbl.cfl or 1/(ndim*(2*polyOrder+1))
      cflMin = math.min(cflMin, myCfl)
      s:setCfl(myCfl)
   end
   log(string.format("Using CFL number %g", cflMin))

   -- setup each species
   for _, s in pairs(species) do
      s:createGrid(tbl.lower, tbl.upper, tbl.cells, decompCuts, periodicDirs)
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

   -- initialize species distributions and write them out
   for _, s in pairs(species) do
      s:initDist()
   end

   -- function to write data to file
   local function writeData(frame, tCurr)
      for _, s in pairs(species) do s:write(frame, tCurr) end
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
   function timeSteppers.rk1(tCurr, dt)
      local status, dtSuggested
      -- RK stage 1
      status, dtSuggested = fowardEuler(tCurr, dt, 1, 2)
      if status == false then return status, dtSuggested end

      for nm, s in pairs(species) do
	 speciesRkFields[nm][1]:copy(speciesRkFields[nm][2])
      end
      return status, dtSuggested
   end   

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
