-- Gkyl ------------------------------------------------------------------------
--
-- Plasma solver on a Cartesian grid. Works in arbitrary CDIM/VDIM
-- (VDIM>=CDIM) with either Vlasov, gyrokinetic, fuilds and Maxwell,
-- Poisson or specified EM fields.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Basis = require "Basis"
local Collisions = require "App.Collisions"
local DecompRegionCalc = require "Lib.CartDecomp"
local Field = require "App.Field"
local Grid = require "Grid"
local LinearTrigger = require "Lib.LinearTrigger"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Sources = require "App.Sources"
local Species = require "App.Species"
local Time = require "Lib.Time"
local date = require "xsys.date"
local lume = require "Lib.lume"
local xsys = require "xsys"
local Projection = require "App.Projection"
local lfs = require "lfs"

-- function to create basis functions
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "tensor" then
      return Basis.CartModalTensor { ndim = ndim, polyOrder = polyOrder }
   end
end

-- function to check if file exists
local function file_exists(name)
   if lfs.attributes(name) then return true else return false end
end

-- top-level method to build application "run" method
local function buildApplication(self, tbl)
   local log = Logger {
      logToFile = xsys.pickBool(tbl.logToFile, true)
   }
      
   log(date(false):fmt()); log("\n") -- time-stamp for sim start
   if GKYL_HG_CHANGESET then
      log(string.format("Gkyl built with %s\n", GKYL_HG_CHANGESET))
   end
   if GKYL_BUILD_DATE then
      log(string.format("Gkyl built on %s\n", GKYL_BUILD_DATE))
   end

   -- function to warn user about default values
   local function warnDefault(varVal, varNm, default)
      if varVal then return varVal end
      log(string.format(" ** WARNING: %s not specified, assuming %s\n", varNm, tostring(default)))
      return default
   end

   log("Initializing PlasmaOnCartGrid simulation ...\n")
   local tmStart = Time.clock()

   local cdim = #tbl.lower -- configuration space dimensions
   assert(cdim == #tbl.upper, "upper should have exactly " .. cdim .. " entries")
   assert(cdim == #tbl.cells, "cells should have exactly " .. cdim .. " entries")

   -- basis function name
   local basisNm = tbl.basis and tbl.basis or "serendipity"
   if basisNm ~= "serendipity" and basisNm ~= "maximal-order" and basisNm ~= "tensor" then
      assert(false, "Incorrect basis type " .. basisNm .. " specified")
   end

   -- FV methods don't need to specify polyOrder
   local polyOrder = tbl.polyOrder and tbl.polyOrder or 0 -- polynomial order

   -- create basis function for configuration space
   local confBasis = createBasis(basisNm, cdim, polyOrder)

   -- I/O method
   local ioMethod = tbl.ioMethod and tbl.ioMethod or "MPI"
   if ioMethod ~= "POSIX" and ioMethod ~= "MPI" then
      assert(false, "ioMethod must be one of 'MPI' or 'POSIX'. Provided '" .. ioMethod .. "' instead")
   end

   local goodStepperNames = { "rk1", "rk2", "rk3", "rk3s4", "fvDimSplit" }
   -- time-stepper
   local timeStepperNm = warnDefault(tbl.timeStepper, "timeStepper", "rk3")
   if not lume.find(goodStepperNames, timeStepperNm) then
      assert(false, "Incorrect timeStepper type " .. timeStepperNm .. " specified")
   end

   -- CFL fractions for various steppers
   local stepperCFLFracs = { rk1 = 1.0, rk2 = 1.0, rk3 = 1.0, rk3s4 = 2.0, fvDimSplit = 1.0 }

   local cflFrac = tbl.cflFrac
   -- Compute CFL fraction if not specified
   if  cflFrac == nil then
      cflFrac = stepperCFLFracs[timeStepperNm]
   end

   -- Number of fields needed for each stepper type
   local stepperNumFields = { rk1 = 3, rk2 = 3, rk3 = 3, rk3s4 = 4, fvDimSplit = 3 }

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
	 if d<1 or d>cdim then
	    assert(false, "Directions in periodicDirs table should be 1 (for X), 2 (for Y) or 3 (for Z)")
	 end
	 periodicDirs[i] = d
      end
   end

   -- configuration space decomp object (eventually, this will be
   -- slaved to the phase-space decomp)
   local decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      useShared = useShared,
   }

   -- pick grid ctor based on uniform/non-uniform grid
   local GridConstructor = Grid.RectCart
   if tbl.coordinateMap then
      GridConstructor = Grid.NonUniformRectCart
   elseif tbl.mapc2p then 
      GridConstructor = Grid.MappedCart
   end
   -- setup configuration space grid
   local confGrid = GridConstructor {
      lower = tbl.lower,
      upper = tbl.upper,
      cells = tbl.cells,
      periodicDirs = periodicDirs,
      decomposition = decomp,
      mappings = tbl.coordinateMap,
      mapc2p = tbl.mapc2p,
   }
   confGrid:write("grid.bp")

   -- read in information about each species
   local species = {}
   for nm, val in pairs(tbl) do
      if Species.SpeciesBase.is(val) then
	 val:fullInit(tbl) -- initialize species
	 species[nm] = val
	 species[nm]:setName(nm)
	 species[nm]:setIoMethod(ioMethod)
      end
   end

   -- set up each species
   for _, s in pairs(species) do
      -- set up conf grid and basis
      s:setConfGrid(confGrid)
      s:setConfBasis(confBasis)
      -- set up phase grid and basis
      s:createGrid(confGrid)
      s:createBasis(basisNm, polyOrder)
   end

   -- read in information about each species
   local sources = {}
   for nm, val in pairs(tbl) do
      if Sources.SourceBase.is(val) then
	 val:fullInit(tbl) -- initialize sources
	 sources[nm] = val
	 sources[nm]:setName(nm)
      end
   end

   -- configuration space decomp object (eventually, this will be
   -- slaved to the phase-space decomp)
   local decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      useShared = useShared,
   }

   -- pick grid ctor based on uniform/non-uniform grid
   local GridConstructor = Grid.RectCart
   if tbl.coordinateMap then
      GridConstructor = Grid.NonUniformRectCart
   end
   -- setup configuration space grid
   local grid = GridConstructor {
      lower = tbl.lower,
      upper = tbl.upper,
      cells = tbl.cells,
      periodicDirs = periodicDirs,
      decomposition = decomp,
      mappings = tbl.coordinateMap,
   }

   -- add grid to app object
   self._confGrid = grid

   -- set conf grid for each species
   for _, s in pairs(species) do
      s:setConfGrid(grid)
      s:alloc(stepperNumFields[timeStepperNm])
   end

   -- set conf grid for each source
   for _, s in pairs(sources) do
      s:setConfGrid(grid)
   end   

   local cflMin = GKYL_MAX_DOUBLE
   -- compute CFL numbers
   for _, s in pairs(species) do
      local ndim = s:getNdim()
      local myCfl = tbl.cfl and tbl.cfl or cflFrac/(2*polyOrder+1)
      cflMin = math.min(cflMin, myCfl)
      s:setCfl(cflMin)
   end

   local function completeFieldSetup(fld)
      fld:fullInit(tbl) -- complete initialization
      fld:setIoMethod(ioMethod)
      fld:setBasis(confBasis)
      fld:setGrid(confGrid)
      do
	 local myCfl = tbl.cfl and tbl.cfl or cflFrac/(2*polyOrder+1)
	 if fld.isElliptic then
	    local myCfl = tbl.cfl and tbl.cfl or cflFrac/(cdim*(2*polyOrder+1))
	 end
	 cflMin = math.min(cflMin, myCfl)
	 fld:setCfl(myCfl)
      end
      log(string.format("Using CFL number %g\n", cflMin))
      
      -- allocate field data
      fld:alloc(stepperNumFields[timeStepperNm])

      -- initialize field solvers and diagnostics
      fld:createDiagnostics()
   end

   -- setup information about fields: if this is not specified, it is
   -- assumed there are no force terms (neutral particles)
   local field = nil
   local nfields = 0
   for _, val in pairs(tbl) do
      if Field.FieldBase.is(val) then
        field = val
        completeFieldSetup(field)
        nfields = nfields + 1
      end
   end
   assert(nfields<=1, "PlasmaOnCartGrid: can only specify one Field object!")
   if field == nil then field = Field.NoField {} end

   -- initialize funcField, which is sometimes needed to initialize species
   local funcField = nil
   nfields = 0
   for _, val in pairs(tbl) do
      if Field.FuncFieldBase.is(val) then
        funcField = val
        completeFieldSetup(funcField)
        nfields = nfields + 1
      end
   end
   assert(nfields<=1, "PlasmaOnCartGrid: can only specify one FuncField object!")
   if funcField == nil then funcField = Field.NoField {} end
   funcField:createSolver()
   funcField:initField()
   
   -- initialize species solvers and diagnostics
   for nm, s in pairs(species) do
      local hasE, hasB = field:hasEB()
      local funcHasE, funcHasB = funcField:hasEB()
      s:createSolver(hasE or funcHasE, hasB or funcHasB, funcField)
      s:initDist()
      s:createDiagnostics()
   end

   -- initialize source solvers
   for nm, s in pairs(sources) do
      s:createSolver(species, field)
   end   

   -- initialize field (sometimes requires species to have been initialized)
   for nm, s in pairs(species) do
      s:calcCouplingMoments(0, 0, 1)
   end
   field:createSolver(species, funcField)
   field:initField(species)

   -- apply species BCs 
   for nm, s in pairs(species) do
      -- this is a dummy forwardEuler call because some BCs require 
      -- auxFields to be set, which is controlled by species solver
      s:forwardEuler(0, 0, species, {field, funcField}, 1, 2)
      s:applyBc(0, 0, s:rkStepperFields()[1])
   end

   -- function to write data to file
   local function writeData(tCurr, force)
      for _, s in pairs(species) do s:write(tCurr, force) end
      field:write(tCurr, force)
      funcField:write(tCurr, force)
   end

   -- function to write restart frames to file
   local function writeRestart(tCurr)
      for _, s in pairs(species) do s:writeRestart(tCurr) end
      field:writeRestart(tCurr)
      funcField:writeRestart(tCurr)
   end

   -- function to read from restart frame
   local function readRestart() --> time at which restart was written
      local rTime = 0.0
      -- read fields first, in case needed for species init or BCs
      field:readRestart()
      funcField:readRestart()
      for _, s in pairs(species) do
	 rTime = s:readRestart()
      end
      return rTime
   end

   local tStart = 0.0 -- by default start at t=0
   if GKYL_COMMANDS[1] == "restart" then
      -- give everyone a chance to adjust ICs based on restart frame
      -- and adjust tStart accordingly
      tStart = readRestart()
   else
      writeData(0.0) -- write initial conditions
   end

   -- determine whether we need two steps in forwardEuler
   local step2 = false
   if field.step2 ~= nil then step2 = field.step2 end

   -- function to take a single forward-euler time-step
   local function forwardEuler(tCurr, dt, inIdx, outIdx)
      local status, dtSuggested = true, GKYL_MAX_DOUBLE

      -- update EM field
      for nm, s in pairs(species) do
	 -- compute moments needed in coupling with fields and
	 -- collisions (the species should update internal datastructures). 
         s:calcCouplingMoments(tCurr, dt, inIdx)
      end
      -- note that this can be either an elliptic solve, which updates inIdx
      -- or a hyperbolic solve, which updates outIdx, or a combination of both
      local myStatus, myDtSuggested = field:forwardEuler(tCurr, dt, species, inIdx, outIdx)
      status = status and myStatus
      dtSuggested = math.min(dtSuggested, myDtSuggested)

      -- compute functional field (if any)
      funcField:forwardEuler(tCurr, dt)
      
      -- update species
      for nm, s in pairs(species) do
         local myStatus, myDtSuggested = s:forwardEuler(tCurr, dt, species, {field, funcField}, inIdx, outIdx)

	 status = status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end

      -- some systems (e.g. EM GK) require a second step to complete the forward Euler
      if step2 then      
         -- update EM field.. step 2 (if necessary). 
         -- note: no calcCouplingMoments call because field:forwardEulerStep2 either reuses already calculated moments, 
         --       or other moments are calculated in field:forwardEulerStep2
         local myStatus, myDtSuggested = field:forwardEulerStep2(tCurr, dt, species, inIdx, outIdx)
         status = status and myStatus
         dtSuggested = math.min(dtSuggested, myDtSuggested)

         -- update species.. step 2 (if necessary)
         for nm, s in pairs(species) do
            local myStatus, myDtSuggested = s:forwardEulerStep2(tCurr, dt, species, {field, funcField}, inIdx, outIdx)

            status = status and myStatus
            dtSuggested = math.min(dtSuggested, myDtSuggested)
         end
      end

      return status, dtSuggested
   end

   -- various functions to copy/increment fields
   local function copy(outIdx, aIdx)
      for nm, s in pairs(species) do
         s:copyRk(outIdx, aIdx)
      end
      field:copyRk(outIdx, aIdx)
   end
   local function combine(outIdx, a, aIdx, ...)
      for nm, s in pairs(species) do
         s:combineRk(outIdx, a, aIdx, ...)
      end
      field:combineRk(outIdx, a, aIdx, ...)
   end

   -- various time-steppers. See gkyl docs for formulas for various
   -- SSP-RK schemes:
   -- http://gkyl.readthedocs.io/en/latest/dev/ssp-rk.html
   local timeSteppers = {}
   local stepperTime = 0.0

   -- function to advance solution using RK1 scheme (UNSTABLE! Only for testing)
   function timeSteppers.rk1(tCurr, dt)
      local status, dtSuggested
      status, dtSuggested = forwardEuler(tCurr, dt, 1, 2)
      if status == false then return status, dtSuggested end
      local tm = Time.clock()
      copy(1, 2)
      stepperTime = stepperTime + (Time.clock() - tm)

      return status, dtSuggested 
   end

   -- function to advance solution using SSP-RK2 scheme (mildly
   -- unstable and in general should not be used)
   function timeSteppers.rk2(tCurr, dt)
      local status, dtSuggested
      -- RK stage 1
      status, dtSuggested = forwardEuler(tCurr, dt, 1, 2)
      if status == false then return status, dtSuggested end

      -- RK stage 2
      status, dtSuggested = forwardEuler(tCurr+dt, dt, 2, 3)
      if status == false then return status, dtSuggested end
      local tm = Time.clock()
      combine(2, 1.0/2.0, 1, 1.0/2.0, 3)
      copy(1, 2)
      stepperTime = stepperTime + (Time.clock() - tm)

      return status, dtSuggested
   end

   -- function to advance solution using SSP-RK3 scheme
   function timeSteppers.rk3(tCurr, dt)
      local status, dtSuggested, tm
      -- RK stage 1
      status, dtSuggested = forwardEuler(tCurr, dt, 1, 2)
      if status == false then return status, dtSuggested end

      -- RK stage 2
      status, dtSuggested = forwardEuler(tCurr+dt, dt, 2, 3)
      if status == false then return status, dtSuggested end
      tm = Time.clock()
      combine(2, 3.0/4.0, 1, 1.0/4.0, 3)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 3
      status, dtSuggested = forwardEuler(tCurr+dt/2, dt, 2, 3)
      if status == false then return status, dtSuggested end
      tm = Time.clock()
      combine(2, 1.0/3.0, 1, 2.0/3.0, 3)
      copy(1, 2)
      stepperTime = stepperTime + (Time.clock() - tm)

      return status, dtSuggested
   end

   -- function to advance solution using 4-stage SSP-RK3 scheme
   function timeSteppers.rk3s4(tCurr, dt)
      local status, dtSuggested, tm
      -- RK stage 1
      status, dtSuggested = forwardEuler(tCurr, dt, 1, 2)
      if status == false then return status, dtSuggested end
      tm = Time.clock()
      combine(3, 1.0/2.0, 1, 1.0/2.0, 2)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 2
      status, dtSuggested = forwardEuler(tCurr+dt/2, dt, 3, 4)
      if status == false then return status, dtSuggested end
      tm = Time.clock()
      combine(2, 1.0/2.0, 3, 1.0/2.0, 4)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 3
      status, dtSuggested = forwardEuler(tCurr+dt, dt, 2, 3)
      if status == false then return status, dtSuggested end
      tm = Time.clock()
      combine(4, 2.0/3.0, 1, 1.0/6.0, 2, 1.0/6.0, 3)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 4
      status, dtSuggested = forwardEuler(tCurr+dt/2, dt, 4, 3)
      if status == false then return status, dtSuggested end
      tm = Time.clock()
      combine(1, 1.0/2.0, 4, 1.0/2.0, 3)
      stepperTime = stepperTime + (Time.clock() - tm)

      return status, dtSuggested
   end

   -- update solution in specified direction
   local function updateInDirection(dir, tCurr, dt)
      local status, dtSuggested = true, GKYL_MAX_DOUBLE
      local fIdx = { {1,2}, {2,1}, {1,2} } -- for indexing inp/out fields

      -- update species
      for nm, s in pairs(species) do
	 local vars = s:rkStepperFields()
	 local inp, out = vars[fIdx[dir][1]], vars[fIdx[dir][2]]
	 local myStatus, myDtSuggested = s:updateInDirection(dir, tCurr, dt, inp, out)
	 status =  status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end
      do
	 -- update field
	 local vars = field:rkStepperFields()
	 local inp, out = vars[fIdx[dir][1]], vars[fIdx[dir][2]]
	 local myStatus, myDtSuggested = field:updateInDirection(dir, tCurr, dt, inp, out)
	 status =  status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end

      return status, dtSuggested
   end

   -- update sources
   local function updateSource(dataIdx, tCurr, dt)
      -- make list of species data to operate on
      local speciesVar = {}
      for nm, s in pairs(species) do
	 speciesVar[nm] = s:rkStepperFields()[dataIdx]
      end
      -- field data to operate on
      local fieldVar = field:rkStepperFields()[dataIdx]

      local status, dtSuggested = true, GKYL_MAX_DOUBLE
      -- update sources
      for nm, s in pairs(sources) do
	 local myStatus, myDtSuggested = s:updateSource(tCurr, dt, speciesVar, fieldVar)
	 status =  status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end

      return status, dtSuggested
   end

   -- function to advance solution using FV dimensionally split scheme
   function timeSteppers.fvDimSplit(tCurr, dt)
      local status, dtSuggested = true, GKYL_MAX_DOUBLE
      local fIdx = { {1,2}, {2,1}, {1,2} } -- for indexing inp/out fields      

      -- copy in case we need to take this step again
      copy(3, 1)

      -- update source by half time-step
      do
	 local myStatus, myDtSuggested = updateSource(1, tCurr, dt/2)
	 status = status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end

      -- update solution in each direction
      for d = 1, cdim do
	 local myStatus, myDtSuggested = updateInDirection(d, tCurr, dt)
	 status =  status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end

      -- update source by half time-step
      if status then
	 local myStatus, myDtSuggested
	 if fIdx[cdim][2] == 2 then
	    myStatus, myDtSuggested = updateSource(2, tCurr, dt/2)
	 else
	    myStatus, myDtSuggested = updateSource(1, tCurr, dt/2)
	 end
	 status = status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end

      if not status then
	 copy(1, 3) -- restore old solution in case of failure
      else
	 -- if solution not already in field[1], copy for use in next
	 -- time-step
	 if fIdx[cdim][2] == 2 then copy(1, 2) end
      end
      

      return status, dtSuggested
   end

   local tmEnd = Time.clock()
   log(string.format("Initializing completed in %g sec\n\n", tmEnd-tmStart))

   -- read some info about restarts (default is to write restarts 1/5 of sim)
   local restartFrameEvery = tbl.restartFrameEvery and tbl.restartFrameEvery or 0.2
   local restartFrameAfter = tbl.restartFrameAfter and tbl.restartFrameAfter or GKYL_MAX_DOUBLE

   -- return function that runs main simulation loop
   return function(self)
      log("Starting main loop of PlasmaOnCartGrid simulation ...\n\n")
      local tStart, tEnd = tStart, tbl.tEnd

      -- sanity check: don't run if not needed
      if tStart >= tEnd then return end

      local maxDt = tbl.maximumDt and tbl.maximumDt or tEnd-tStart -- max time-step
      local initDt =  tbl.suggestedDt and tbl.suggestedDt or maxDt -- initial time-step
      local step = 1
      local tCurr = tStart
      local myDt = initDt

      -- triggers for 10% and 1% loggers
      local logTrigger = LinearTrigger(0, tEnd, 10)
      local logTrigger1p = LinearTrigger(0, tEnd, 100)
      local tenth = 0
      if tStart > 0 then
	 tenth = math.floor(tStart/tEnd*10)
      end

      -- triggers for restarts
      local restartTrigger = LinearTrigger(0.0, tEnd, math.floor(1/restartFrameEvery))
      local nRestart = 0
      -- function to check if restart frame should happen
      local function checkWriteRestart(t)
	 if restartTrigger(t) then
	    if nRestart > 1 then return true end
	 end
	 nRestart = nRestart+1
	 return false
      end

      local p1c = 0
      if tStart > 0 then
	 p1c = math.floor(tStart/tEnd*100) % 10
      end

      local logCount = 0 -- this is needed to avoid initial log message
      -- for writing out log messages
      local function writeLogMessage(tCurr)
	 if logTrigger(tCurr) then
	    if logCount > 0 then
	       log (string.format(
		       " Step %5d at time %g. Time step %g. Completed %g%s\n", step, tCurr, myDt, tenth*10, "%"))
	    else
	       logCount = logCount+1
	    end
	    tenth = tenth+1
	 end
	 if logTrigger1p(tCurr) then
	    log(string.format("%d", p1c))
	    p1c = (p1c+1)%10
	 end
      end

      local tmSimStart = Time.clock()
      local first = true
      local failcount = 0
      local stopfile = GKYL_OUT_PREFIX .. ".stop"

      -- main simulation loop
      while true do
	 -- call time-stepper
	 local status, dtSuggested = timeSteppers[timeStepperNm](tCurr, myDt)
   
         -- if stopfile exists, break
         if (file_exists(stopfile)) then
            writeData(tCurr+myDt, true)
            writeRestart(tCurr+myDt)
            break
         end

	 -- check status and determine what to do next
	 if status then
            if first then initDt = dtSuggested; first = false end
	    writeLogMessage(tCurr+myDt)
	    -- we must write data first before calling writeRestart in
	    -- order not to mess up numbering of frames on a restart
	    writeData(tCurr+myDt)
	    if checkWriteRestart(tCurr+myDt) then
	       writeRestart(tCurr+myDt)
	    end	    
	    
	    tCurr = tCurr + myDt
	    myDt = math.min(dtSuggested, maxDt)
	    step = step + 1
	    if (tCurr >= tEnd) then
	       break
	    end
	 else
	    log (string.format(" ** Time step %g too large! Will retake with dt %g\n", myDt, dtSuggested))
	    myDt = dtSuggested
            if (myDt < 1e-3*initDt) then 
               failcount = failcount + 1
               if failcount > 20 then
                  writeData(tCurr+myDt, true)
                  log(string.format("ERROR: Timestep below 1e-3*initDt for 20 consecutive steps. Exiting...\n"))
                  break
               end
            end
	 end
      end
      local tmSimEnd = Time.clock()

      -- compute time spent in various parts of code
      local tmSlvr = 0.0
      for _, s in pairs(species) do
	 tmSlvr = tmSlvr+s:totalSolverTime()
      end

      local tmMom, tmIntMom, tmBc, tmColl = 0.0, 0.0, 0.0, 0.0
      for _, s in pairs(species) do
         tmMom = tmMom + s:momCalcTime()
         tmIntMom = tmIntMom + s:intMomCalcTime()
         tmBc = tmBc + s:totalBcTime()
         if s.collisions then
	    for _, c in pairs(s.collisions) do
	       tmColl = tmColl + c:totalTime()
	    end
         end
      end

      local tmSrc = 0.0
      for _, s in pairs(sources) do
         tmSrc = tmSrc + s:totalTime()
      end

      local tmTotal = tmSimEnd-tmSimStart
      log(string.format("\nTotal number of time-steps %s\n\n", step))
      log(string.format(
	     "Solver took				%9.5f sec   (%7.6f s/step)   (%6.3f%%)\n",
	     tmSlvr, tmSlvr/step, 100*tmSlvr/tmTotal))
      log(string.format(
	     "Solver BCs took 			%9.5f sec   (%7.6f s/step)   (%6.3f%%)\n",
	     tmBc, tmBc/step, 100*tmBc/tmTotal))
      log(string.format(
	     "Field solver took 			%9.5f sec   (%7.6f s/step)   (%6.3f%%)\n",
	     field:totalSolverTime(), field:totalSolverTime()/step, 100*field:totalSolverTime()/tmTotal))
      log(string.format(
	     "Field solver BCs took			%9.5f sec   (%7.6f s/step)   (%6.3f%%)\n",
	     field:totalBcTime(), field:totalBcTime()/step, 100*field:totalBcTime()/tmTotal))
      log(string.format(
	     "Function field solver took		%9.5f sec   (%7.6f s/step)   (%6.3f%%)\n",
	     funcField:totalSolverTime(), funcField:totalSolverTime()/step, 100*funcField:totalSolverTime()/tmTotal))
      log(string.format(
	     "Moment calculations took		%9.5f sec   (%7.6f s/step)   (%6.3f%%)\n",
	     tmMom, tmMom/step, 100*tmMom/tmTotal))
      log(string.format(
	     "Integrated moment calculations took	%9.5f sec   (%7.6f s/step)   (%6.3f%%)\n",
	     tmIntMom, tmIntMom/step, 100*tmIntMom/tmTotal))
      log(string.format(
	     "Field energy calculations took		%9.5f sec   (%7.6f s/step)   (%6.3f%%)\n",
	     field:energyCalcTime(), field:energyCalcTime()/step, 100*field:energyCalcTime()/tmTotal))
      log(string.format(
	     "Collision solver(s) took		%9.5f sec   (%7.6f s/step)   (%6.3f%%)\n",
	     tmColl, tmColl/step, 100*tmColl/tmTotal))
      log(string.format(
	     "Source updaters took 			%9.5f sec   (%7.6f s/step)   (%6.3f%%)\n",
	     tmSrc, tmSrc/step, 100*tmSrc/tmTotal))
      log(string.format(
	     "Stepper combine/copy took		%9.5f sec   (%7.6f s/step)   (%6.3f%%)\n",
	     stepperTime, stepperTime/step, 100*stepperTime/tmTotal))
      log(string.format(
	     "Main loop completed in			%9.5f sec   (%7.6f s/step)   (%6.f%%)\n\n",
	     tmTotal, tmTotal/step, 100*tmTotal/tmTotal))
      log(date(false):fmt()); log("\n") -- time-stamp for sim end

      if file_exists(stopfile) then os.remove(stopfile) end -- clean up
   end
end

-- PlasmaOnCartGrid application object
local App = Proto()

function App:init(tbl)
   self._runApplication = buildApplication(self, tbl)
end

function App:getConfGrid()
   return self._confGrid
end

function App:run()
   -- by default command is "run"
   if #GKYL_COMMANDS == 0 then GKYL_COMMANDS[1] = "run" end

   -- take action
   if GKYL_COMMANDS[1] == "run" or GKYL_COMMANDS[1] == "restart" then
      return self:_runApplication()
   elseif GKYL_COMMANDS[1] == "init" then
      return function (...) end
   end
end

return {
   AdiabaticSpecies = Species.AdiabaticSpecies,
   App = App,
   BgkCollisions = Collisions.BgkCollisions,   
   FuncMaxwellField = Field.FuncMaxwellField,
   FuncVlasovSpecies = Species.FuncVlasovSpecies,
   GkField = Field.GkField,
   GkGeometry = Field.GkGeometry,
   GkSpecies = Species.GkSpecies,
   HamilVlasovSpecies = Species.HamilVlasovSpecies,
   IncompEulerSpecies = Species.IncompEulerSpecies,
   MaxwellField = Field.MaxwellField,
   MomentSpecies = Species.MomentSpecies,
   NoField = Field.NoField,
   VlasovSpecies = Species.VlasovSpecies,
   VmLBOCollisions = Collisions.VmLBOCollisions,
   GkLBOCollisions = Collisions.GkLBOCollisions,
   VoronovIonization = Collisions.VoronovIonization,
   Projection = Projection,

   -- valid pre-packaged species-field systems
   Gyrokinetic = {
      App = App, Species = Species.GkSpecies, Field = Field.GkField, Geometry = Field.GkGeometry,
      FunctionProjection = Projection.GkProjection.FunctionProjection, 
      MaxwellianProjection = Projection.GkProjection.MaxwellianProjection,
      BgkCollisions = Collisions.BgkCollisions,
      LBOCollisions = Collisions.GkLBOCollisions,
      AdiabaticSpecies = Species.AdiabaticSpecies,
   },
   IncompEuler = {
      App = App, Species = Species.IncompEulerSpecies, Field = Field.GkField
   },
   VlasovMaxwell = {
      App = App, Species = Species.VlasovSpecies, FuncSpecies = Species.FuncVlasovSpecies,
      Field = Field.MaxwellField,
      FunctionProjection = Projection.VlasovProjection.FunctionProjection, 
      MaxwellianProjection = Projection.VlasovProjection.MaxwellianProjection,
      BgkCollisions = Collisions.BgkCollisions,
      LBOCollisions = Collisions.VmLBOCollisions,
   },
   Moments = {
      App = App, Species = Species.MomentSpecies, Field = Field.MaxwellField,
      CollisionlessEmSource = Sources.CollisionlessEmSource
   } 
}
