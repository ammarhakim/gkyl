-- Gkyl ------------------------------------------------------------------------
--
-- Plasma solver on a Cartesian grid. Works in arbitrary CDIM/VDIM
-- (VDIM>=CDIM) with either Vlasov, gyrokinetic, fluids and Maxwell,
-- Poisson or specified EM fields.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Infrastructure loads
local Alloc            = require "Alloc"
local Basis            = require "Basis"
local DataStruct       = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid             = require "Grid"
local Lin              = require "Lib.Linalg"
local LinearTrigger    = require "Lib.LinearTrigger"
local Messenger        = require "Comm.Messenger"
local Logger           = require "Lib.Logger"
local Mpi              = require "Comm.Mpi"
local Proto            = require "Lib.Proto"
local Time             = require "Lib.Time"
local date             = require "xsys.date"
local lfs              = require "lfs"
local lume             = require "Lib.lume"
local xsys             = require "xsys"
local ffi              = require "ffi"

math = require("sci.math").generic -- this is global so that it affects input file

-- App loads (do not load specific app objects here, but only things
-- needed to run the App itself. Specific objects should be loaded in
-- the  methods defined at the bottom of this file)
local SpeciesBase = require "App.Species.SpeciesBase"
local FieldBase = require ("App.Field.FieldBase").FieldBase
local ExternalFieldBase = require ("App.Field.FieldBase").ExternalFieldBase
local NoField = require ("App.Field.FieldBase").NoField
local PopApp = require "App.Species.Population"

-- Function to create basis functions.
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "tensor" then
      return Basis.CartModalTensor { ndim = ndim, polyOrder = polyOrder }
   end
end

-- Function to check if file exists.
local function file_exists(name)
   if lfs.attributes(name) then return true else return false end
end

-- Function to check if runtime is less than walltime.
local function out_of_walltime(startTime, wallTime, boolC)
   boolC[0] = (Time.clock()-startTime) > (wallTime-300.)  -- Exit 5 min early.
   -- Need to reduce tsofar across MPI ranks, because they keep different clocks.
   Mpi.Bcast(boolC, 1, Mpi.C_BOOL, 0, Mpi.COMM_WORLD)
   return boolC[0]
end

-- Top-level method to build application "run" method.
local function buildApplication(self, tbl)
   local log = Logger {
      logToFile = xsys.pickBool(tbl.logToFile, true)
   }
      
   log(date(false):fmt()); log("\n") -- Time-stamp for sim start.
   if GKYL_GIT_CHANGESET then
      log(string.format("Gkyl built with %s\n", GKYL_GIT_CHANGESET))
   end
   if GKYL_BUILD_DATE then
      log(string.format("Gkyl built on %s\n", GKYL_BUILD_DATE))
   end

   -- Function to warn user about default values.
   local function warnDefault(varVal, varNm, default)
      if varVal then return varVal end
      log(string.format(" ** WARNING: %s not specified, assuming %s\n", varNm, tostring(default)))
      return default
   end

   local timers = {
      total = 0.,   init = 0.,   timeloop = 0.,   unaccounted = 0.,
      diagnostics = 0.,   restart = 0.,   dydt = 0.,   fwdEuler = 0.,
      combine = 0.,   copy = 0.,   combineFwdEuler = 0.,   dtCalc = 0.,
      bc = 0.,   timeloop_unaccounted = 0.,   dydt_unaccounted = 0.,   fwdEuler_unaccounted = 0.,
      dydt_extfield = 0.,   dydt_mom = 0.,   dydt_momcross = 0.,   dydt_field = 0.,
      dydt_speciesadv = 0.,   dydt_speciesadvcross = 0.,   dydt_collisionless = 0.,   dydt_collisions = 0.,
      dydt_boundflux = 0.,   dydt_sources = 0.,   dydt_speciesadv_unaccounted = 0.,
   }
   timers.total = Time.clock()
   timers.init  = timers.total
   log(string.format("Initializing %s simulation ...\n", self.label))

   local cdim = #tbl.lower -- Configuration space dimensions.
   assert(tbl.cells, "Must specify the number of cells in 'cells'.")
   assert(cdim == #tbl.upper, "upper should have exactly " .. cdim .. " entries")
   assert(cdim == #tbl.cells, "cells should have exactly " .. cdim .. " entries")

   -- Basis function name.
   local basisNm = tbl.basis and tbl.basis or "serendipity"
   if basisNm ~= "serendipity" and basisNm ~= "maximal-order" and basisNm ~= "tensor" then
      assert(false, "Incorrect basis type " .. basisNm .. " specified")
   end

   local polyOrder = assert(tbl.polyOrder, "Must specify polynomial order in 'polyOrder'.")

   -- Create basis function for configuration space.
   local confBasis = createBasis(basisNm, cdim, polyOrder)

   -- I/O method
   local ioMethod = tbl.ioMethod and tbl.ioMethod or "MPI"
   if ioMethod ~= "POSIX" and ioMethod ~= "MPI" then
      assert(false, "ioMethod must be one of 'MPI' or 'POSIX'. Provided '" .. ioMethod .. "' instead")
   end

   -- Optional wallclock time allotted for this simulation.
   local maxWallTime = tbl.maxWallTime and tbl.maxWallTime or GKYL_MAX_DOUBLE

   -- Time-stepper.
   local goodStepperNames = { "rk1", "rk2", "rk3", "rk3s4" }
   local timeStepperNm = warnDefault(tbl.timeStepper, "timeStepper", "rk3")
   local timeIntegrator
   if lume.find(goodStepperNames, timeStepperNm) then
      if timeStepperNm == "rk1" then
         timeIntegrator = require "App.TimeSteppers.RK1"
      elseif timeStepperNm == "rk2" then
         timeIntegrator = require "App.TimeSteppers.SSP_RK2"
      elseif timeStepperNm == "rk3" then
         timeIntegrator = require "App.TimeSteppers.SSP_RK3"
      elseif timeStepperNm == "rk3s4" then
         timeIntegrator = require "App.TimeSteppers.SSP_RK3s4"
      else
         assert(false, "Time stepper not implemented.")
      end
   else
      assert(false, "Incorrect timeStepper type " .. timeStepperNm .. " specified.")
   end
   local timeStepper = timeIntegrator{}

   local maxDt = tbl.maximumDt or GKYL_MAX_DOUBLE

   local cflFrac = tbl.cflFrac or timeStepper.cflFrac   -- CFL fraction.

   -- Tracker for timestep
   local dtTracker = DataStruct.DynVector { numComponents = 1, }
   local dtPtr = Lin.Vec(1)

   -- Used in reducing time step across species communicator.
   local dtMinLocal, dtMinGlobal = Lin.Vec(1), Lin.Vec(1)

   -- Count the number species, needed for decomposition.
   local numSpecies = 0
   for _, val in pairs(tbl) do
      if SpeciesBase.is(val) then numSpecies = numSpecies+1 end
   end

   -- Object that handles MPI/NCCL communication (some comms are still elsewhere).
   local commManager = Messenger{
      cells      = tbl.cells,   decompCutsConf     = tbl.decompCuts,
      numSpecies = numSpecies,  parallelizeSpecies = tbl.parallelizeSpecies,
   }

   -- Extract periodic directions.
   local periodicDirs = {}
   if tbl.periodicDirs then
      for i, d in ipairs(tbl.periodicDirs) do
	 if d<1 or d>cdim then
	    assert(false, "Directions in periodicDirs table should be 1 (for X), 2 (for Y) or 3 (for Z)")
	 end
	 periodicDirs[i] = d
      end
   end

   -- Pick grid ctor based on uniform/non-uniform grid.
   local GridConstructor = Grid.RectCart
   if tbl.coordinateMap then
      GridConstructor = Grid.NonUniformRectCart
   elseif tbl.mapc2p then 
      GridConstructor = Grid.MappedCart
   end
   -- Setup configuration space grid.
   local confGrid = GridConstructor {
      lower = tbl.lower,  decomposition = commManager:getConfDecomp(),
      upper = tbl.upper,  mappings = tbl.coordinateMap,
      cells = tbl.cells,  mapc2p = tbl.mapc2p,
      periodicDirs = periodicDirs,  world = tbl.world, 
      messenger = commManager,
   }
   -- Read in information about each species.
   local population = PopApp{ messenger = commManager }
   for nm, val in pairs(tbl) do
      if SpeciesBase.is(val) then
	 population.species[nm] = val
	 population.species[nm]:setName(nm)
	 population.species[nm]:setIoMethod(ioMethod)
      end
   end
   lume.setOrder(population:getSpecies())  -- Save order in metatable to loop in the same order (w/ orderedIter, better for I/O).

   -- Distribute species across MPI ranks in species communicator.
   population:decompose()

   -- Create other subcommunicators.
   commManager:createSubComms(confGrid)

   if tbl.coordinateMap or tbl.mapc2p then 
      -- Placing this writing step after population:decompose so that only one species rank writes.
      confGrid:write(confBasis, {0, population:getComm_host(),})
   end

   for _, s in population.iterGlobal() do
      s:fullInit(tbl) -- Initialize species.
   end

   -- Setup each species.
   for _, s in population.iterGlobal() do
      -- Set up conf grid and basis.
      s:setConfGrid(confGrid)
      s:setConfBasis(confBasis)
      -- Set up phase grid and basis.
      s:createGrid(confGrid)
      s:createBasis(basisNm, polyOrder)
   end
   for _, s in population.iterLocal() do -- Only allocate fields for species in this rank.
      s:alloc(timeStepper.numFields)
   end

   -- Add grid to app object.
   self._confGrid = confGrid

   local cflMin = GKYL_MAX_DOUBLE
   -- Compute CFL numbers.
   for _, s in population.iterLocal() do
      local myCfl = tbl.cfl and tbl.cfl or cflFrac
      cflMin = math.min(cflMin, myCfl)
      s:setCfl(cflMin)
   end

   local function completeFieldSetup(fld)
      fld:fullInit(tbl) -- Complete initialization.
      fld:setIoMethod(ioMethod)
      fld:setBasis(confBasis)
      fld:setGrid(confGrid)
      do
	 local myCfl = tbl.cfl and tbl.cfl or cflFrac
	 if fld.isElliptic then
	    myCfl = tbl.cfl and tbl.cfl or cflFrac
	 end
	 cflMin = math.min(cflMin, myCfl)
	 fld:setCfl(myCfl)
      end
      
      -- Allocate field data.
      fld:alloc(timeStepper.numFields)

      -- Initialize field solvers and diagnostics.
      fld:createDiagnostics()
   end

   -- Setup information about fields: if this is not specified, it is
   -- assumed there are no force terms (neutral particles).
   local field = nil
   local nfields = 0
   for _, val in pairs(tbl) do
      if FieldBase.is(val) then
         field = val
         completeFieldSetup(field)
         nfields = nfields + 1
      end
   end
   assert(nfields<=1, "PlasmaOnCartGrid: can only specify one Field object!")
   if field == nil then field = NoField {} end

   -- Initialize externalField, which is sometimes needed to initialize species.
   local externalField = nil
   nfields = 0
   for _, val in pairs(tbl) do
      if ExternalFieldBase.is(val) then
         externalField = val
         completeFieldSetup(externalField)
         nfields = nfields + 1
      end
   end
   assert(nfields<=1, "PlasmaOnCartGrid: can only specify one ExternalField object!")
   if externalField == nil then externalField = NoField {} end
   externalField:createSolver(population)
   externalField:initField()
   
   -- Initialize species solvers and diagnostics.
   for _, s in population.iterGlobal() do
      s:initCrossSpeciesCoupling(population)    -- Call this before createSolver if updaters are all created in createSolver.
   end
   for _, s in population.iterLocal() do
      s:createSolver(field, externalField)
      s:initDist(externalField, population:getSpecies())
   end
   -- Create field solver (sometimes requires species solver to have been created).
   field:createSolver(population, externalField)

   -- Some objects require an additional createSolver step after self-species
   -- solvers are created due to cross-species interactions.
   for _, s in population.iterGlobal() do
      s:createCouplingSolver(population, field, externalField)
   end

   for _, s in population.iterLocal() do
      -- Compute moments needed in coupling with fields and collisions.
      s:calcCouplingMoments(0.0, 1, population:getSpecies())
   end

   -- Initialize field (sometimes requires species to have been initialized).
   field:initField(population)

   commManager:startCommGroup()  -- Needed by NCCL.
   for _, s in population.iterGlobal() do
      -- Compute/comunicate moments needed for cross-species interactions.
      -- Needs to happen after field:advance because we'll initiate (multi-GPU)
      -- communication here that needs to finish before starting other comms.
      s:calcCrossCouplingMoments(0., 1, population)
   end
   commManager:endCommGroup()  -- Needed by NCCL.

   -- Initialize diagnostic objects.
   for _, s in population.iterLocal() do s:createDiagnostics(field) end

   -- Apply species BCs.
   for _, s in population.iterLocal() do
      -- This is a dummy forwardEuler call because some BCs require 
      -- auxFields to be set, which is controlled by species solver.
      if s.charge == 0.0 then
      	 s:advance(0, population, {NoField {}, NoField {}}, 1, 2)
      else
	 s:advance(0, population, {field, externalField}, 1, 2)
      end
      s:applyBcInitial(0, field, externalField, 1, 1)
   end
   for _, s in population.iterGlobal() do
      s:advanceCrossSpeciesCoupling(0, population, {field, externalField}, 1, 2)
   end

   -- Function to write data to file.
   local function writeData(tCurr, force)
      for _, s in population.iterLocal() do s:write(tCurr, field, force) end
      field:write(tCurr, force)
      externalField:write(tCurr, force)
   end

   -- Function to write restart frames to file.
   local function writeRestart(tCurr)
      for _, s in population.iterLocal() do s:writeRestart(tCurr) end
      field:writeRestart(tCurr)
      externalField:writeRestart(tCurr)
   end
   
   -- Function to read from restart frame.
   local function readRestart() --> Time at which restart was written.
      local rTime = 0.0
      dtTracker:read(string.format("dt.bp"))
      local _, dtLast = dtTracker:lastData()
      -- Read fields first, in case needed for species init or BCs.
      field:readRestart()
      externalField:readRestart()
      for _, s in population.iterLocal() do
         -- This is a dummy forwardEuler call because some BCs require 
         -- auxFields to be set, which is controlled by species solver.
         if s.charge == 0 then
            s:advance(0, population, {NoField {}, NoField {}}, 1, 2)
         else
            s:advance(0, population, {field, externalField}, 1, 2)
         end
         s:setDtGlobal(dtLast[1])
         rTime = s:readRestart(field, externalField)
      end
      for _, s in population.iterGlobal() do
         s:advanceCrossSpeciesCoupling(0, population, {field, externalField}, 1, 2)
      end
      return rTime
   end

   local tStart = 0.0 -- By default start at t=0.
   if GKYL_COMMANDS[1] == "restart" then
      -- Give everyone a chance to adjust ICs based on restart frame
      -- and adjust tStart accordingly.
      tStart = readRestart()
   else
      writeData(0.0)      -- Write initial conditions.
      writeRestart(0.0)   -- Write initial conditions as a restart file.
   end

   -- Determine whether we need two steps in forwardEuler.
   local nstep = 1
   if field.nstep ~= nil then nstep = field.nstep end

   -- Various functions to copy/increment fields.
   local function copy(outIdx, aIdx)
      local tmCopy = Time.clock()
      for _, s in population.iterLocal() do s:copyRk(outIdx, aIdx) end
      field:copyRk(outIdx, aIdx)
      timers.copy = timers.copy + Time.clock() - tmCopy
   end
   local function combineFunc(outIdx, a, aIdx, ...)
      for _, s in population.iterLocal() do s:combineRk(outIdx, a, aIdx, ...) end
      field:combineRk(outIdx, a, aIdx, ...)
   end
   local function combine(outIdx, a, aIdx, ...)
      local tmCombine = Time.clock()
      combineFunc(outIdx, a, aIdx, ...)
      timers.combine = timers.combine + Time.clock() - tmCombine
   end
   local function combineFwdEuler(outIdx, a, aIdx, ...)
      local tmCombineFwdEuler = Time.clock()
      combineFunc(outIdx, a, aIdx, ...)
      timers.combineFwdEuler = timers.combineFwdEuler + Time.clock() - tmCombineFwdEuler
   end
   local function applyBc(tCurr, inIdx, outIdx, ...)
      local tmBc = Time.clock()
      for _, s in population.iterLocal() do s:applyBcIdx(tCurr, field, externalField, inIdx, outIdx, ...) end
      field:applyBcIdx(tCurr, outIdx)
      timers.bc = timers.bc + Time.clock() - tmBc
   end

   local appStatus = {
      -- Store some flags and info about the state and work done by the app.
      success = true,
      step = 0,
      -- For diagnostics:
      nFwdEuler = 0,
      -- Below: an entry per stage (max 4 stages + 2 for operator splitting).
      nFail = {0, 0, 0, 0, 0, 0},
      dtDiff = {{GKYL_MAX_DOUBLE, 0.}, {GKYL_MAX_DOUBLE, 0.},
                {GKYL_MAX_DOUBLE, 0.}, {GKYL_MAX_DOUBLE, 0.},
                {GKYL_MAX_DOUBLE, 0.}, {GKYL_MAX_DOUBLE, 0.}},
   }

   local function dydt(tCurr, inIdx, outIdx)
      -- Compute the time rate of change (dy/dt).
      local tmDydt = Time.clock()

      field:clearCFL()
      for _, s in population.iterLocal() do s:clearCFL() end

      -- Compute functional field (if any).
      externalField:advance(tCurr)
      
      for _, s in population.iterLocal() do
         -- Compute moments needed in coupling with fields and collisions.
         s:calcCouplingMoments(tCurr, inIdx, population:getSpecies())
      end

      -- Update EM field.
      -- Note that this can be either an elliptic solve, which updates inIdx
      -- or a hyperbolic solve, which updates outIdx = RHS, or a combination of both.
      field:advance(tCurr, population, inIdx, outIdx)

      commManager:startCommGroup()  -- Needed by NCCL.
      for _, s in population.iterGlobal() do
         -- Compute/comunicate moments needed for cross-species interactions.
	 -- Needs to happen after field:advance because we'll initiate (multi-GPU)
	 -- communication here that needs to finish before starting other comms.
         s:calcCrossCouplingMoments(tCurr, inIdx, population)
      end
      commManager:endCommGroup()  -- Needed by NCCL.

      -- Update species.
      for _, s in population.iterLocal() do
         if s.charge == 0 then
            s:advance(tCurr, population, {NoField {}, NoField {}}, inIdx, outIdx)
         else
            s:advance(tCurr, population, {field, externalField}, inIdx, outIdx)
         end
      end
      for _, s in population.iterGlobal() do
         s:advanceCrossSpeciesCoupling(tCurr, population, {field, externalField}, inIdx, outIdx)
      end

      -- Some systems (e.g. EM GK) require additional step(s) to complete the forward Euler.
      for istep = 2, nstep do      
         -- Update EM field.. step 2 (if necessary). 
         -- Note: no calcCouplingMoments call because field:forwardEulerStep2
         -- either reuses already calculated moments, or other moments are
         -- calculated in field:forwardEulerStep2.
         local advanceString = "advanceStep" .. istep
         field[advanceString](field, tCurr, population, inIdx, outIdx)

         -- Update species.. step 2 (if necessary).
         for _, s in population.iterLocal() do
            s[advanceString](s, tCurr, population, {field, externalField}, inIdx, outIdx)
         end
      end
      timers.dydt = timers.dydt + Time.clock() - tmDydt
   end

   -- Function to take a single forward-euler time-step.
   local function forwardEuler(tCurr, dt, inIdx, outIdx, stat)
      local tmFwdEuler = Time.clock()

      appStatus.nFwdEuler = appStatus.nFwdEuler + 1

      local tmDtCalc = Time.clock()
      -- Get suggested dt from each field and species.
      -- Species dt:
      dtMinLocal[1], dtMinGlobal[1] = GKYL_MAX_DOUBLE, GKYL_MAX_DOUBLE
      for _, s in population.iterLocal() do dtMinLocal[1] = math.min(dtMinLocal[1], s:suggestDt()) end
      -- Reduce dtMin across species communicator.
      Mpi.Allreduce(dtMinLocal:data(), dtMinGlobal:data(), 1, Mpi.DOUBLE, Mpi.MIN, commManager:getSpeciesComm_host())

      -- Field dt:
      local dtMin = math.min(dtMinGlobal[1], field:suggestDt())

      -- MF 2021/08/04: We will disable this criteria for now, so that the dt is
      --                as it's been in g2 and not quite like it is in g0 now.
      --local dt_maxRelDiff = 0.01
      ---- Check if dtMin is slightly smaller than dt. Use dt if it is
      ---- (avoids retaking steps if dt changes are very small).
      --local dt_relDiff = (dt-dtMin)/dt
      --if (dt_relDiff > 0 and dt_relDiff < dt_maxRelDiff) then dtMin = dt end

      -- Don't take a time-step larger that input dt.
      stat.dt_actual = dt < dtMin and dt or dtMin
      stat.dt_suggested = dtMin

      stat.dt_actual = math.min(math.min(stat.dt_actual, tbl.tEnd - tCurr + 1e-20), maxDt)
      
      -- After deciding global dt, tell species.
      for _, s in population.iterLocal() do s:setDtGlobal(stat.dt_actual) end
      timers.dtCalc = timers.dtCalc + Time.clock() - tmDtCalc

      -- Take forward Euler step in fields and species
      -- NOTE: order of these arguments matters... outIdx must come before inIdx.
      combineFwdEuler(outIdx, stat.dt_actual, outIdx, 1.0, inIdx)
      applyBc(tCurr, inIdx, outIdx, calcCflFlag)

      timers.fwdEuler = timers.fwdEuler + Time.clock() - tmFwdEuler
   end

   -- Set functions in time stepper object.
   timeStepper:createSolver(appStatus, {combine, copy, dydt, forwardEuler},
                            {population:getSpecies(), field, externalField})

   local devDiagnose = function()
      -- Perform performance/numerics-related diagnostics.
      field:printDevDiagnostics()
   end

   timers.init = Time.clock()-timers.init
   log(string.format("Initialization completed in %g sec\n\n", timers.init))

   -- Read some info about restarts (default is to write restarts 1/20 (5%) of sim, 
   -- but no need to write restarts more frequently than regular diagnostic output).
   local restartFrameEvery = tbl.restartFrameEvery and tbl.restartFrameEvery or math.max(1/20.0, 1/tbl.nFrame)

   -- Zero out species, field and extField timers (so they don't contain init times).
   for _, s in population.iterLocal() do s:clearTimers() end
   field:clearTimers()
   externalField:clearTimers()

   -- Return function that runs main simulation loop.
   return function(self)
      log(string.format("Starting main loop of %s simulation ...\n\n", self.label))
      local tStart, tEnd = tStart, tbl.tEnd

      -- Sanity check: don't run if not needed.
      if tStart >= tEnd then return end

      local dt_max = tbl.maximumDt and tbl.maximumDt or tEnd-tStart -- max time-step
      local dt_init = tbl.suggestedDt and tbl.suggestedDt or dt_max -- initial time-step
      local tCurr = tStart
      local dt_next = dt_init
      appStatus.step = 1

      -- Triggers for 10% and 1% loggers.
      local logTrigger = LinearTrigger(0, tEnd, 10)
      local logTrigger1p = LinearTrigger(0, tEnd, 100)
      local tenth = 0
      if tStart > 0 then
	 tenth = math.floor(tStart/tEnd*10)
      end

      -- Triggers for restarts.
      local restartTrigger = LinearTrigger(0.0, tEnd, math.floor(1/restartFrameEvery))
      local nRestart = 0
      -- Function to check if restart frame should happen.
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

      local logCount = 0 -- This is needed to avoid initial log message.
      -- For writing out log messages.
      local function writeLogMessage(tCurr)
	 if logTrigger(tCurr) then
	    if logCount > 0 then
	       log (string.format(
		       " Step %6d at time  %#11.8g.  Time step  %.6e.  Completed %g%s\n", 
                       appStatus.step, tCurr, dt_next, tenth*10, "%"))
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

      local first = true
      local failcount = 0
      local irestart = 0
      local stopfile = GKYL_OUT_PREFIX .. ".stop"
      local timesUp = ffi.new("bool [?]", 1)

      -- Main simulation loop.
      timers.timeloop = Time.clock()
      while true do
	 -- Call time-stepper.
	 local stepStatus = timeStepper:advance(tCurr, dt_next)
    
         -- If stopfile exists, break.
         if file_exists(stopfile) or out_of_walltime(timers.total, maxWallTime, timesUp) then
            local tlatest = tCurr+stepStatus.dt_actual
            local tmDiags = Time.clock()
            writeData(tlatest, true)
            timers.diagnostics = timers.diagnostics + Time.clock() - tmDiags
            local tmRestart = Time.clock()
            writeRestart(tlatest)
            timers.restart = timers.restart + Time.clock() - tmRestart
            dtPtr:data()[0] = stepStatus.dt_actual
            dtTracker:appendData(tlatest, dtPtr)
            dtTracker:write(string.format("dt.bp"), tlatest, irestart)
            break
         end

         -- Abort simulation if the suggested timestep is 0, which means there are likely NaNs.
         -- Don't write anything.
         if (stepStatus.dt_suggested == 0.0) then
            log(string.format(" ERROR: dt is zero, there are likely NaNs. Terminating without writing files."))
            break
         end

	 if appStatus.success then
            if first then 
               log(string.format(" Step 0 at time %g. Time step %g. Completed 0%%\n", tCurr, stepStatus.dt_actual))
               dt_init = math.min(dt_max, stepStatus.dt_actual); first = false
            end
	    tCurr = tCurr + stepStatus.dt_actual
            -- Track dt.
            dtPtr:data()[0] = stepStatus.dt_actual
            dtTracker:appendData(tCurr, dtPtr)
            -- Write log
	    writeLogMessage(tCurr)
	    -- We must write data first before calling writeRestart in
	    -- order not to mess up numbering of frames on a restart.
            local tmDiags = Time.clock()
	    writeData(tCurr)
            timers.diagnostics = timers.diagnostics + Time.clock() - tmDiags
	    if checkWriteRestart(tCurr) then
               local tmRestart = Time.clock()
	       writeRestart(tCurr)
               dtTracker:write(string.format("dt.bp"), tCurr, irestart)
               irestart = irestart + 1
               timers.restart = timers.restart + Time.clock() - tmRestart
	    end	    
	    
	    dt_next = math.min(stepStatus.dt_suggested, dt_max)
	    appStatus.step = appStatus.step + 1
	    if (tCurr >= tEnd) then break end
	 else
	    log(string.format(" ** Step failed with dt=%g! Will retake with dt=%g\n", dt_next, stepStatus.dt_suggested))
	    dt_next = stepStatus.dt_suggested
	 end

         if (dt_next < 1e-4*dt_init) then 
            failcount = failcount + 1
            log(string.format("WARNING: Timestep dt = %g is below 1e-4*dt_init. Fail counter = %d...\n", dt_next, failcount))
            if failcount > 20 then
               local tmDiags = Time.clock()
               writeData(tCurr+stepStatus.dt_actual, true)
               dtTracker:write(string.format("dt.bp"), tCurr+stepStatus.dt_actual)
               timers.diagnostics = timers.diagnostics + Time.clock() - tmDiags
               log(string.format("ERROR: Timestep below 1e-4*dt_init for 20 consecutive steps. Exiting...\n"))
               break
            end
         else
            failcount = 0
         end
      end
      local tmEnd = Time.clock()
      timers.timeloop = tmEnd - timers.timeloop
      timers.total = tmEnd - timers.total

      -- Compute time spent in various parts of code.
      local tmSlvr = 0.0
      local tmMom, tmIntMom, tmBc, tmColl = 0.0, 0.0, 0.0, 0.0
      local tmSrc, tmCollNonSlvr = 0.0, 0.0

      log(string.format("\n\nTotal number of time-steps %s\n", appStatus.step))
      log(string.format("   Number of forward-Euler calls %s\n", appStatus.nFwdEuler))
      for stI = 2, 3 do
         log(string.format("   Number of RK stage-"..stI.." failures %s\n", appStatus.nFail[stI]))
         if appStatus.nFail[stI] > 0 then
            log(string.format("     Min rel dt diff for RK stage-"..stI.." failures %s\n", appStatus.dtDiff[stI][1]))
            log(string.format("     Max rel dt diff for RK stage-"..stI.." failures %s\n", appStatus.dtDiff[stI][2]))
         end
      end

      -- Compute the time spent inside of dydt components.
      timers.dydt_extfield = externalField:getTimer('advance')
      timers.dydt_field = field:getTimer('advance')
      for _, s in population.iterLocal() do
         timers.dydt_mom             = timers.dydt_mom             + s:getTimer('mom')
         timers.dydt_momcross        = timers.dydt_momcross        + s:getTimer('momcross')
         timers.dydt_speciesadv      = timers.dydt_speciesadv      + s:getTimer('advance')
         timers.dydt_speciesadvcross = timers.dydt_speciesadvcross + s:getTimer('advancecross')
         timers.dydt_collisionless   = timers.dydt_collisionless   + s:getTimer('collisionless')
         timers.dydt_collisions      = timers.dydt_collisions      + s:getTimer('collisions')
         timers.dydt_boundflux       = timers.dydt_boundflux       + s:getTimer('boundflux')
         timers.dydt_sources         = timers.dydt_sources         + s:getTimer('sources')
      end
      log("\n")
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Initialization took", timers.init, 0., 100*timers.init/timers.total))
      timers.unaccounted = timers.unaccounted + timers.init
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Time loop took", timers.timeloop, timers.timeloop/appStatus.step, 100*timers.timeloop/timers.total))
      timers.unaccounted = timers.unaccounted + timers.timeloop
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "  - dydt took", timers.dydt, timers.dydt/appStatus.step, 100*timers.dydt/timers.total))
      timers.timeloop_unaccounted = timers.timeloop_unaccounted + timers.dydt
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "    * moments:", timers.dydt_mom, timers.dydt_mom/appStatus.step, 100*timers.dydt_mom/timers.total))
      timers.dydt_unaccounted = timers.dydt_unaccounted + timers.dydt_mom
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "    * cross moments:", timers.dydt_momcross, timers.dydt_momcross/appStatus.step, 100*timers.dydt_momcross/timers.total))
      timers.dydt_unaccounted = timers.dydt_unaccounted + timers.dydt_momcross
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "    * field solver:", timers.dydt_field, timers.dydt_field/appStatus.step, 100*timers.dydt_field/timers.total))
      timers.dydt_unaccounted = timers.dydt_unaccounted + timers.dydt_field
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "    * external field solver:", timers.dydt_extfield, timers.dydt_extfield/appStatus.step, 100*timers.dydt_extfield/timers.total))
      timers.dydt_unaccounted = timers.dydt_unaccounted + timers.dydt_extfield
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "    * species advance:", timers.dydt_speciesadv, timers.dydt_speciesadv/appStatus.step, 100*timers.dydt_speciesadv/timers.total))
      timers.dydt_unaccounted = timers.dydt_unaccounted + timers.dydt_speciesadv
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "      > collisionless:", timers.dydt_collisionless, timers.dydt_collisionless/appStatus.step, 100*timers.dydt_collisionless/timers.total))
      timers.dydt_speciesadv_unaccounted = timers.dydt_speciesadv_unaccounted + timers.dydt_collisionless
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "      > collisions:", timers.dydt_collisions, timers.dydt_collisions/appStatus.step, 100*timers.dydt_collisions/timers.total))
      timers.dydt_speciesadv_unaccounted = timers.dydt_speciesadv_unaccounted + timers.dydt_collisions
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "      > boundary flux:", timers.dydt_boundflux, timers.dydt_boundflux/appStatus.step, 100*timers.dydt_boundflux/timers.total))
      timers.dydt_speciesadv_unaccounted = timers.dydt_speciesadv_unaccounted + timers.dydt_boundflux
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "      > sources:", timers.dydt_sources, timers.dydt_sources/appStatus.step, 100*timers.dydt_sources/timers.total))
      timers.dydt_speciesadv_unaccounted = timers.dydt_speciesadv_unaccounted + timers.dydt_sources
      timers.dydt_speciesadv_unaccounted = timers.dydt_speciesadv - timers.dydt_speciesadv_unaccounted
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "      > [Unaccounted]", timers.dydt_speciesadv_unaccounted, timers.dydt_speciesadv_unaccounted/appStatus.step, 100*timers.dydt_speciesadv_unaccounted/timers.dydt_speciesadv))
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "    * species advanceCross:", timers.dydt_speciesadvcross, timers.dydt_speciesadvcross/appStatus.step, 100*timers.dydt_speciesadvcross/timers.total))
      timers.dydt_unaccounted = timers.dydt_unaccounted + timers.dydt_speciesadvcross
      timers.dydt_unaccounted = timers.dydt - timers.dydt_unaccounted
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "    * [Unaccounted]", timers.dydt_unaccounted, timers.dydt_unaccounted/appStatus.step, 100*timers.dydt_unaccounted/timers.dydt))
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "  - forward Euler took", timers.fwdEuler, timers.fwdEuler/appStatus.step, 100*timers.fwdEuler/timers.total))
      timers.timeloop_unaccounted = timers.timeloop_unaccounted + timers.fwdEuler
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "    * dt calc took", timers.dtCalc, timers.dtCalc/appStatus.step, 100*timers.dtCalc/timers.total))
      timers.fwdEuler_unaccounted = timers.fwdEuler_unaccounted + timers.dtCalc
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "    * combine took", timers.combineFwdEuler, timers.combineFwdEuler/appStatus.step, 100*timers.combineFwdEuler/timers.total))
      timers.fwdEuler_unaccounted = timers.fwdEuler_unaccounted + timers.combineFwdEuler
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "    * BCs took", timers.bc, timers.bc/appStatus.step, 100*timers.bc/timers.total))
      timers.fwdEuler_unaccounted = timers.fwdEuler_unaccounted + timers.bc
      timers.fwdEuler_unaccounted = timers.fwdEuler - timers.fwdEuler_unaccounted
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "    * [Unaccounted]", timers.fwdEuler_unaccounted, timers.fwdEuler_unaccounted/appStatus.step, 100*timers.fwdEuler_unaccounted/timers.fwdEuler))
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "  - copy took", timers.copy, timers.copy/appStatus.step, 100*timers.copy/timers.total))
      timers.timeloop_unaccounted = timers.timeloop_unaccounted + timers.copy
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "  - combine took", timers.combine, timers.combine/appStatus.step, 100*timers.combine/timers.total))
      timers.timeloop_unaccounted = timers.timeloop_unaccounted + timers.combine
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "  - diagnostic write took", timers.diagnostics, timers.diagnostics/appStatus.step, 100*timers.diagnostics/timers.total))
      timers.timeloop_unaccounted = timers.timeloop_unaccounted + timers.diagnostics
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "  - restart write took", timers.restart, timers.restart/appStatus.step, 100*timers.restart/timers.total))
      timers.timeloop_unaccounted = timers.timeloop_unaccounted + timers.restart
      timers.timeloop_unaccounted = timers.timeloop - timers.timeloop_unaccounted
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "  - [Unaccounted]", timers.timeloop_unaccounted, timers.timeloop_unaccounted/appStatus.step, 100*timers.timeloop_unaccounted/timers.timeloop))

      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Whole simulation took", timers.total, timers.total/appStatus.step, 100*timers.total/timers.total))
      timers.unaccounted = timers.total - timers.unaccounted
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "[Unaccounted]", timers.unaccounted, timers.unaccounted/appStatus.step, 100*timers.unaccounted/timers.total))

      log("\n")
      log(date(false):fmt()); log("\n") -- Time-stamp for sim end.

      -- Perform other numerical/performance diagnostics.
      devDiagnose()

      if file_exists(stopfile) then os.remove(stopfile) end -- Clean up.
   end
end

-- PlasmaOnCartGrid application object.
local App = Proto()

function App:init(tbl) self._runApplication = buildApplication(self, tbl) end

function App:getConfGrid() return self._confGrid end

function App:run()
   -- By default command is "run".
   if #GKYL_COMMANDS == 0 then GKYL_COMMANDS[1] = "run" end

   -- Take action.
   if GKYL_COMMANDS[1] == "run" or GKYL_COMMANDS[1] == "restart" then
      return self:_runApplication()
   elseif GKYL_COMMANDS[1] == "init" then
      return function (...) end
   end
end

return {
   Gyrofluid = function ()
      App.label = "Gyrofluid"
      return  {
	 AdiabaticSpecies = require ("App.Species.AdiabaticSpecies"),
         AdiabaticBasicBC = require "App.BCs.AdiabaticBasic",
	 App = App,
         BasicBC = require ("App.BCs.GyrofluidBasic").GyrofluidBasic,
         AbsorbBC = require ("App.BCs.GyrofluidBasic").GyrofluidAbsorb,
         CopyBC = require ("App.BCs.GyrofluidBasic").GyrofluidCopy,
         SheathBC = require ("App.BCs.GyrofluidBasic").GyrofluidSheath,
         ZeroFluxBC = require ("App.BCs.GyrofluidBasic").GyrofluidZeroFlux,
	 Field = require ("App.Field.GkField").GkField,
	 FunctionProjection = require ("App.Projection.GyrofluidProjection").FunctionProjection, 
	 Geometry = require ("App.Field.GkField").GkGeometry,
	 GyrofluidProjection = require ("App.Projection.GyrofluidProjection").GyrofluidProjection, 
         HeatFlux = require "App.Collisions.GfHeatFlux",
         PASCollisions = require "App.Collisions.GfPitchAngleScattering",
         Source = require "App.Sources.GyrofluidSource",
	 Species = require "App.Species.GyrofluidSpecies",
      }
   end,

   Gyrokinetic = function ()
      App.label = "Gyrokinetic"
      return  {
	 AdiabaticSpecies = require ("App.Species.AdiabaticSpecies"),
         AdiabaticBasicBC = require "App.BCs.AdiabaticBasic",
	 App = App,
         BasicBC = require ("App.BCs.GkBasic").GkBasic,
         MaxwellianBC = require("App.BCs.GkMaxwellianBC"),
         AbsorbBC = require ("App.BCs.GkBasic").GkAbsorb,
         CopyBC = require ("App.BCs.GkBasic").GkCopy,
         NeutralRecyclingBC = require "App.BCs.NeutralRecycling",
         OpenBC = require ("App.BCs.GkBasic").GkOpen,
         ReflectBC = require ("App.BCs.GkBasic").GkReflect,
         SheathBC = require "App.BCs.GkSheath",
         ZeroFluxBC = require ("App.BCs.GkBasic").GkZeroFlux,
         TwistShiftBC = require "App.BCs.TwistShift",
	 VmAbsorbBC = require ("App.BCs.VlasovBasic").VlasovAbsorb,
	 VmReflectBC = require ("App.BCs.VlasovBasic").VlasovReflect,
	 BGKCollisions = require "App.Collisions.GkBGKCollisions",
	 BgkCollisions = require "App.Collisions.GkBGKCollisions",
	 Diffusion = require ("App.Collisions.Diffusion").DiffusionGyrokinetic,
	 ChargeExchange = require "App.Collisions.GkChargeExchange",
	 Field = require ("App.Field.GkField").GkField,
	 AmbipolarSheathField = require "App.Field.AmbipolarSheathField",
	 FunctionProjection = require ("App.Projection.GkProjection").FunctionProjection, 
	 Geometry = require ("App.Field.GkField").GkGeometry,
	 Ionization = require "App.Collisions.GkIonization",
	 LBOCollisions = require "App.Collisions.GkLBOCollisions",
	 LboCollisions = require "App.Collisions.GkLBOCollisions",
	 MaxwellianProjection = require ("App.Projection.GkProjection").MaxwellianProjection,
	 Species = require "App.Species.GkSpecies",
	 Source = require "App.Sources.GkSource",
	 Vlasov = require ("App.Species.VlasovSpecies"),
	 VmMaxwellianProjection = require ("App.Projection.VlasovProjection").MaxwellianProjection,
	 VmSource = require "App.Sources.VmSource",
      }
   end,

   IncompEuler = function ()
      App.label = "Incompressible Euler"
      return {
	 App = App,
         BasicBC = require ("App.BCs.FluidBasic").FluidBasic,
         AbsorbBC = require ("App.BCs.FluidBasic").FluidAbsorb,
         CopyBC = require ("App.BCs.FluidBasic").FluidCopy,
         ZeroFluxBC = require ("App.BCs.FluidBasic").FluidZeroFlux,
	 Diffusion = require "App.Collisions.Diffusion",
	 Field = require ("App.Field.GkField").GkField,
         Source = require "App.Sources.FluidSource",
	 Species = require "App.Species.IncompEulerSpecies",
      }
   end,

   HasegawaWakatani = function ()
      App.label = "Hasegawa-Wakatani"
      return {
	 App = App,
	 Field = require ("App.Field.GkField").GkField,
	 Species = require "App.Species.HasegawaWakataniSpecies",
	 Diffusion = require "App.Collisions.Diffusion",
      }
   end,

   VlasovMaxwell = function ()
      App.label = "Vlasov-Maxwell"
      GenSpecies = {
         -- Non-relativistic Vlasov flavors
         Vlasov = require ("App.Species.VlasovSpecies").VlasovMaxwell,
         VlasovPoisson = require ("App.Species.VlasovSpecies").VlasovPoisson,
         VlasovPoissonA = require ("App.Species.VlasovSpecies").VlasovPoissonA,
         VlasovNeutral = require ("App.Species.VlasovSpecies").VlasovNeutral,
         VlasovGenGeoNeutral = require ("App.Species.VlasovSpecies").VlasovGenGeoNeutral,
         -- Relativistic Vlasov flavors
         VlasovSR = require ("App.Species.VlasovSpecies").VlasovSRMaxwell,
         VlasovSRNeutral = require ("App.Species.VlasovSpecies").VlasovSRNeutral,
      }
      GenField = {
         Maxwell = require ("App.Field.MaxwellField").MaxwellField,
      }
      return {
         App = App,
         BasicBC = require ("App.BCs.VlasovBasic").VlasovBasic,
         AbsorbBC = require ("App.BCs.VlasovBasic").VlasovAbsorb,
         CopyBC = require ("App.BCs.VlasovBasic").VlasovCopy,
         NeutralRecyclingBC = require "App.BCs.NeutralRecycling",
         OpenBC = require ("App.BCs.VlasovBasic").VlasovOpen,
         ReflectBC = require ("App.BCs.VlasovBasic").VlasovReflect,
         ZeroFluxBC = require ("App.BCs.VlasovBasic").VlasovZeroFlux,
         BronoldFehskeBC = require "App.BCs.BronoldFehskeReflection",
         -- Backwards compatible species and field objects
         Species = require ("App.Species.VlasovSpecies").VlasovMaxwell,
         FuncSpecies = require "App.Species.FuncVlasovSpecies",
         Field = require ("App.Field.MaxwellField").MaxwellField,
         ExternalField = require ("App.Field.MaxwellField").ExternalMaxwellField,
         -- Flexible species and field objects
         GenSpecies = GenSpecies,
         GenField = GenField,
         FunctionProjection = require ("App.Projection.VlasovProjection").FunctionProjection,
         MaxwellianProjection = require ("App.Projection.VlasovProjection").MaxwellianProjection,
         BGKCollisions = require "App.Collisions.VmBGKCollisions",
         LBOCollisions = require "App.Collisions.VmLBOCollisions",
         BgkCollisions = require "App.Collisions.VmBGKCollisions",
         LboCollisions = require "App.Collisions.VmLBOCollisions",
         ChargeExchange = require "App.Collisions.VmChargeExchange",
         Ionization = require "App.Collisions.VmIonization",
	 Diffusion = require ("App.Collisions.Diffusion").DiffusionVlasov,
         SteadySource = require "App.Sources.VmSteadyStateSource",
         Source = require "App.Sources.VmSource",
      }
   end,
   
   PassiveAdvection = function ()
      App.label = "Passively-advected scalar fluid"
      return {
	 App = App,
         BasicBC = require ("App.BCs.FluidBasic").FluidBasic,
         AbsorbBC = require ("App.BCs.FluidBasic").FluidAbsorb,
         CopyBC = require ("App.BCs.FluidBasic").FluidCopy,
         ZeroFluxBC = require ("App.BCs.FluidBasic").FluidZeroFlux,
         TwistShiftBC = require "App.BCs.TwistShift",
	 Diffusion = require "App.Collisions.Diffusion",
         Source = require "App.Sources.FluidSource",
	 Species = require "App.Species.PassiveAdvectionSpecies",
      }
   end,
   
   EulerIso = function ()
      App.label = "Isothermal Euler"
      return {
        App = App,
        BasicBC = require ("App.BCs.FluidBasic").FluidBasic,
        AbsorbBC = require ("App.BCs.FluidBasic").FluidAbsorb,
        CopyBC = require ("App.BCs.FluidBasic").FluidCopy,
        ZeroFluxBC = require ("App.BCs.FluidBasic").FluidZeroFlux,
        Species = require "App.Species.EulerSpecies",
      }
   end,

   Moments = function ()
      -- this is a mere redirection to the stand-alone Moments App
      -- that wraps the G0 Moments App
      local App = require("App.Moments.Moments")
      return App
   end,
}
