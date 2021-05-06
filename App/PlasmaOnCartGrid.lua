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
local Alloc = require "Alloc"
local AllocShared = require "AllocShared"
local Basis = require "Basis"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local LinearTrigger = require "Lib.LinearTrigger"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local date = require "xsys.date"
local lfs = require "lfs"
local lume = require "Lib.lume"
local xsys = require "xsys"
math = require("sci.math").generic -- this is global so that it affects input file

-- App loads (do not load specific app objects here, but only things
-- needed to run the App itself. Specific objects should be loaded in
-- the  methods defined at the bottom of this file)
local SpeciesBase       = require "App.Species.SpeciesBase"
local FluidSourceBase   = require "App.FluidSources.FluidSourceBase"
local FieldBase         = require ("App.Field.FieldBase").FieldBase
local ExternalFieldBase = require ("App.Field.FieldBase").ExternalFieldBase
local NoField           = require ("App.Field.FieldBase").NoField

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

   log(string.format("Initializing %s simulation ...\n", self.label))
   local tmStart = Time.clock()

   local cdim = #tbl.lower -- Configuration space dimensions.
   assert(cdim == #tbl.upper, "upper should have exactly " .. cdim .. " entries")
   assert(cdim == #tbl.cells, "cells should have exactly " .. cdim .. " entries")

   -- Basis function name.
   local basisNm = tbl.basis and tbl.basis or "serendipity"
   if basisNm ~= "serendipity" and basisNm ~= "maximal-order" and basisNm ~= "tensor" then
      assert(false, "Incorrect basis type " .. basisNm .. " specified")
   end

   -- FV methods don't need to specify polyOrder.
   local polyOrder = tbl.polyOrder and tbl.polyOrder or 0 -- Polynomial order.

   -- Create basis function for configuration space.
   local confBasis = createBasis(basisNm, cdim, polyOrder)

   -- I/O method
   local ioMethod = tbl.ioMethod and tbl.ioMethod or "MPI"
   if ioMethod ~= "POSIX" and ioMethod ~= "MPI" then
      assert(false, "ioMethod must be one of 'MPI' or 'POSIX'. Provided '" .. ioMethod .. "' instead")
   end

   local goodStepperNames = { "rk1", "rk2", "rk3", "rk3s4", "fvDimSplit" }
   -- Time-stepper.
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

   -- Tracker for timestep
   local dtTracker = DataStruct.DynVector {
      numComponents = 1,
   }
   local dtPtr = Lin.Vec(1)

   -- Parallel decomposition stuff.
   local useShared  = xsys.pickBool(tbl.useShared, false)   
   local decompCuts = tbl.decompCuts
   if tbl.decompCuts then
      assert(cdim == #tbl.decompCuts, "decompCuts should have exactly " .. cdim .. " entries")
   else
      if not useShared then
	 -- if not specified and not using shared, construct a
	 -- decompCuts automatically
	 local numRanks = Mpi.Comm_size(Mpi.COMM_WORLD)
	 decompCuts = DecompRegionCalc.makeCuts(cdim, numRanks, tbl.cells)
      else
	 assert(false, "Must specify decompCuts when useShared = true")
      end
   end

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

   -- Configuration space decomp object.
   local decomp = DecompRegionCalc.CartProd {
      cuts      = decompCuts,
      useShared = useShared,
   }

   -- Some timers.
   local stepperTime         = 0.
   local fwdEulerCombineTime = 0.
   local writeDataTime       = 0.
   local writeRestartTime    = 0.

   -- Pick grid ctor based on uniform/non-uniform grid.
   local GridConstructor = Grid.RectCart
   if tbl.coordinateMap then
      GridConstructor = Grid.NonUniformRectCart
   elseif tbl.mapc2p then 
      GridConstructor = Grid.MappedCart
   end
   -- Setup configuration space grid.
   local confGrid = GridConstructor {
      lower = tbl.lower,
      upper = tbl.upper,
      cells = tbl.cells,
      periodicDirs  = periodicDirs,
      decomposition = decomp,
      mappings      = tbl.coordinateMap,
      mapc2p        = tbl.mapc2p,
      world         = tbl.world,
   }
   if tbl.coordinateMap or tbl.mapc2p then 
      local metaData = {polyOrder = confBasis:polyOrder(),
                        basisType = confBasis:id(),
                        grid      = GKYL_OUT_PREFIX .. "_grid.bp"}
      confGrid:write("grid.bp", 0.0, metaData)
   end

   -- Read in information about each species.
   local species = {}
   for nm, val in pairs(tbl) do
      if SpeciesBase.is(val) then
	 species[nm] = val
	 species[nm]:setName(nm)
	 species[nm]:setIoMethod(ioMethod)
	 val:fullInit(tbl) -- Initialize species.
      end
   end
   -- Create a keys entry in species so we always loop in the same order.
   local species_keys = {}
   for k in pairs(species) do table.insert(species_keys, k) end
   table.sort(species_keys)
   setmetatable(species, species_keys)

   -- Setup each species.
   for _, s in lume.orderedIter(species) do
      -- Set up conf grid and basis.
      s:setConfGrid(confGrid)
      s:setConfBasis(confBasis)
      -- Set up phase grid and basis.
      s:createGrid(confGrid)
      s:createBasis(basisNm, polyOrder)
      s:alloc(stepperNumFields[timeStepperNm])
   end

   -- Read in information about each fluid source
   local fluidSources = {}
   for nm, val in pairs(tbl) do
      if FluidSourceBase.is(val) then
	 fluidSources[nm] = val
	 fluidSources[nm]:setName(nm)
	 val:fullInit(tbl) -- Initialize fluid sources.
      end
   end

   -- Create a keys entry in sources so we always loop in the same order.
   local fluidSources_keys = {}
   for k in pairs(fluidSources) do table.insert(fluidSources_keys, k) end
   table.sort(fluidSources_keys)
   setmetatable(fluidSources, fluidSources_keys)

   -- Add grid to app object.
   self._confGrid = confGrid

   -- Set conf grid for each fluid source.
   for _, s in lume.orderedIter(fluidSources) do
      s:setConfGrid(confGrid)
   end  

   local cflMin = GKYL_MAX_DOUBLE
   -- Compute CFL numbers.
   for _, s in lume.orderedIter(species) do
      local ndim = s:getNdim()
      local myCfl = tbl.cfl and tbl.cfl or cflFrac/(2*polyOrder+1)
      cflMin = math.min(cflMin, myCfl)
      s:setCfl(cflMin)
   end

   local function completeFieldSetup(fld)
      fld:fullInit(tbl) -- Complete initialization.
      fld:setIoMethod(ioMethod)
      fld:setBasis(confBasis)
      fld:setGrid(confGrid)
      do
	 local myCfl = tbl.cfl and tbl.cfl or cflFrac/(2*polyOrder+1)
	 if fld.isElliptic then
	    myCfl = tbl.cfl and tbl.cfl or cflFrac/(cdim*(2*polyOrder+1))
	 end
	 cflMin = math.min(cflMin, myCfl)
	 fld:setCfl(myCfl)
      end
      
      -- Allocate field data.
      fld:alloc(stepperNumFields[timeStepperNm])

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
   externalField:createSolver()
   externalField:initField()
   
   -- Initialize species solvers and diagnostics.
   for _, s in lume.orderedIter(species) do
      local hasE, hasB = field:hasEB()
      local extHasE, extHasB = externalField:hasEB()
      s:initCrossSpeciesCoupling(species)    -- Call this before createSolver if updaters are all created in createSolver.
      s:createSolver(hasE or extHasE, hasB or extHasB, externalField, hasB)
      s:initDist(externalField)
      s:createDiagnostics()
   end

   -- Initialize fluid source solvers.
   for _, s in lume.orderedIter(fluidSources) do
      s:createSolver(species, field)
   end   

   -- Compute the coupling moments.
   for _, s in lume.orderedIter(species) do
      s:clearMomentFlags(species)
      s:calcCouplingMoments(0.0, 1, species)
   end

   -- Initialize field (sometimes requires species to have been initialized).
   field:createSolver(species, externalField)
   field:initField(species)

   -- Apply species BCs.
   for _, s in lume.orderedIter(species) do
      -- This is a dummy forwardEuler call because some BCs require 
      -- auxFields to be set, which is controlled by species solver.
      if s.charge == 0.0 then
      	 s:advance(0, species, {NoField {}, NoField {}}, 1, 2)
      else
	 s:advance(0, species, {field, externalField}, 1, 2)
      end
      s:applyBc(0, s:rkStepperFields()[1])
   end

   -- Function to write data to file.
   local function writeData(tCurr, force)
      for _, s in lume.orderedIter(species) do s:write(tCurr, force) end
      for _, src in lume.orderedIter(fluidSources) do src:write(tCurr) end 
      field:write(tCurr, force)
      externalField:write(tCurr, force)
   end

   -- Function to write restart frames to file.
   local function writeRestart(tCurr)
      for _, s in lume.orderedIter(species) do s:writeRestart(tCurr) end
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
      for _, s in lume.orderedIter(species) do
         -- This is a dummy forwardEuler call because some BCs require 
         -- auxFields to be set, which is controlled by species solver.
	 if s.charge == 0 then
	    s:advance(0, species, {NoField {}, NoField {}}, 1, 2)
	 else
	    s:advance(0, species, {field, externalField}, 1, 2)
	 end
         s:setDtGlobal(dtLast[1])
	 rTime = s:readRestart()
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
      for _, s in lume.orderedIter(species) do
         s:copyRk(outIdx, aIdx)
      end
      field:copyRk(outIdx, aIdx)
   end
   local function combine(outIdx, a, aIdx, ...)
      for _, s in lume.orderedIter(species) do
         s:combineRk(outIdx, a, aIdx, ...)
      end
      field:combineRk(outIdx, a, aIdx, ...)
   end
   local function applyBc(tCurr, idx, ...)
      for _, s in lume.orderedIter(species) do
         s:applyBcIdx(tCurr, idx, ...)
      end
      field:applyBcIdx(tCurr, idx)
   end

   -- Function to take a single forward-euler time-step.
   local function forwardEuler(tCurr, dt, inIdx, outIdx)
      local calcCflFlag = false
      local dtSuggested
      if dt == nil then calcCflFlag = true end
      field:clearCFL()
      for _, s in lume.orderedIter(species) do
         s:clearCFL()
         s:clearMomentFlags(species)
      end
      -- Compute functional field (if any).
      externalField:advance(tCurr)
      
      for _, s in lume.orderedIter(species) do
         -- Compute moments needed in coupling with fields and
         -- collisions (the species should update internal datastructures). 
         s:calcCouplingMoments(tCurr, inIdx, species)
      end

      -- Update EM field.
      -- Note that this can be either an elliptic solve, which updates inIdx
      -- or a hyperbolic solve, which updates outIdx = RHS, or a combination of both.
      field:advance(tCurr, species, inIdx, outIdx)

      -- Update species.
      for _, s in lume.orderedIter(species) do
         if s.charge == 0 then
            s:advance(tCurr, species, {NoField {}, NoField {}}, inIdx, outIdx)
         else
            s:advance(tCurr, species, {field, externalField}, inIdx, outIdx)
         end
      end

      -- Some systems (e.g. EM GK) require additional step(s) to complete the forward Euler.
      for istep = 2, nstep do      
         -- Update EM field.. step 2 (if necessary). 
         -- Note: no calcCouplingMoments call because field:forwardEulerStep2
         -- either reuses already calculated moments, or other moments are
         -- calculated in field:forwardEulerStep2.
         local advanceString = "advanceStep" .. istep
         field[advanceString](field, tCurr, species, inIdx, outIdx)

         -- Update species.. step 2 (if necessary).
         for _, s in lume.orderedIter(species) do
            s[advanceString](s, tCurr, species, {field, externalField}, inIdx, outIdx)
         end
      end

      if calcCflFlag then
         dtSuggested = tbl.tEnd - tCurr + 1e-20
         if tbl.maximumDt then dtSuggested = math.min(dtSuggested, tbl.maximumDt) end
         
         -- Get suggested dt from each field and species.
         dtSuggested = math.min(dtSuggested, field:suggestDt())
         for _, s in lume.orderedIter(species) do
            dtSuggested = math.min(dtSuggested, s:suggestDt())
         end
         
         -- After deciding global dt, tell species.
         for _, s in lume.orderedIter(species) do
            s:setDtGlobal(dtSuggested)
         end
      else 
         dtSuggested = dt -- From argument list.
         -- If calcCflFlag not being used, need to barrier before doing the RK combine.
         -- When running with calcCflFlag, an all-reduce is done on the time-step to find
         -- the smallest time step, giving us an implicit barrier before we combine RK steps.
         Mpi.Barrier(self._confGrid:commSet().sharedComm)
      end
      -- Take forward Euler step in fields and species
      -- NOTE: order of these arguments matters... outIdx must come before inIdx.
      local tm = Time.clock()
      combine(outIdx, dtSuggested, outIdx, 1.0, inIdx)
      fwdEulerCombineTime = fwdEulerCombineTime + Time.clock() - tm
      applyBc(tCurr, outIdx, calcCflFlag)

      return dtSuggested
   end

   -- Various time-steppers. See gkyl docs for formulas for various
   -- SSP-RK schemes:
   -- http://gkyl.readthedocs.io/en/latest/dev/ssp-rk.html
   local timeSteppers = {}

   -- Function to advance solution using RK1 scheme (UNSTABLE! Only for testing).
   function timeSteppers.rk1(tCurr)
      local dt = forwardEuler(tCurr, nil, 1, 2)
      local tm = Time.clock()
      copy(1, 2)
      stepperTime = stepperTime + (Time.clock() - tm)

      return true, dt
   end

   -- Function to advance solution using SSP-RK2 scheme (mildly
   -- unstable and in general should not be used).
   function timeSteppers.rk2(tCurr)
      -- RK stage 1.
      local dt = forwardEuler(tCurr, nil, 1, 2)

      -- RK stage 2.
      forwardEuler(tCurr+dt, dt, 2, 3)
      local tm = Time.clock()
      combine(2, 1.0/2.0, 1, 1.0/2.0, 3)
      copy(1, 2)
      stepperTime = stepperTime + (Time.clock() - tm)

      return true, dt
   end

   -- Function to advance solution using SSP-RK3 scheme.
   function timeSteppers.rk3(tCurr)
      -- RK stage 1.
      local dt = forwardEuler(tCurr, nil, 1, 2)

      -- RK stage 2.
      forwardEuler(tCurr+dt, dt, 2, 3)
      local tm = Time.clock()
      combine(2, 3.0/4.0, 1, 1.0/4.0, 3)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 3.
      forwardEuler(tCurr+dt/2, dt, 2, 3)
      tm = Time.clock()
      combine(2, 1.0/3.0, 1, 2.0/3.0, 3)
      copy(1, 2)
      stepperTime = stepperTime + (Time.clock() - tm)

      return true, dt
   end

   -- Function to advance solution using 4-stage SSP-RK3 scheme.
   function timeSteppers.rk3s4(tCurr)
      -- RK stage 1.
      local dt = forwardEuler(tCurr, nil, 1, 2)
      local tm = Time.clock()
      combine(3, 1.0/2.0, 1, 1.0/2.0, 2)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 2.
      forwardEuler(tCurr+dt/2, dt, 3, 4)
      tm = Time.clock()
      combine(2, 1.0/2.0, 3, 1.0/2.0, 4)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 3.
      forwardEuler(tCurr+dt, dt, 2, 3)
      tm = Time.clock()
      combine(4, 2.0/3.0, 1, 1.0/6.0, 2, 1.0/6.0, 3)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 4.
      forwardEuler(tCurr+dt/2, dt, 4, 3)
      tm = Time.clock()
      combine(1, 1.0/2.0, 4, 1.0/2.0, 3)
      stepperTime = stepperTime + (Time.clock() - tm)

      return true, dt
   end

   -- Update solution in specified direction.
   local function updateInDirection(dir, tCurr, dt, tryInv)
      local status, dtSuggested = true, GKYL_MAX_DOUBLE
      local fIdx = { {1,2}, {2,1}, {1,2} } -- For indexing inp/out fields.

      local tryInv_next = {}
      -- Update species.
      for _, s in lume.orderedIter(species) do
	 local vars = s:rkStepperFields()
	 local inp, out = vars[fIdx[dir][1]], vars[fIdx[dir][2]]
	 local myStatus, myDtSuggested, myTryInv = s:updateInDirection(
	    dir, tCurr, dt, inp, out, tryInv[s])
	 tryInv_next[s] = myTryInv
	 status =  status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end
      do
	 -- Update field.
	 local vars = field:rkStepperFields()
	 local inp, out = vars[fIdx[dir][1]], vars[fIdx[dir][2]]
	 local myStatus, myDtSuggested = field:updateInDirection(dir, tCurr, dt, inp, out)
	 status =  status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end

      return status, dtSuggested, tryInv_next
   end

   -- Update fluid sources.
   local function updateFluidSource(dataIdx, tCurr, dt)
      -- Make list of species data to operate on.
      local speciesVar = {}
      for nm, s in lume.orderedIter(species) do
	 speciesVar[nm] = s:rkStepperFields()[dataIdx]
      end
      -- Field data to operate on.
      local fieldVar = field:rkStepperFields()[dataIdx]

      -- Expose freely-available array space space. Useful for storing
      -- intermediate quantities like spatial gradient.
      local bufIdx = (dataIdx==1) and 2 or 1
      local speciesBuf = {}
      for nm, s in lume.orderedIter(species) do
         speciesBuf[nm] = s:rkStepperFields()[bufIdx]
      end
      -- Field data to operate on.
      local fieldBuf = field:rkStepperFields()[bufIdx]

      local status, dtSuggested = true, GKYL_MAX_DOUBLE
      -- Update fluid source terms.
      for _, s in lume.orderedIter(fluidSources) do
	 local myStatus, myDtSuggested = s:updateFluidSource(
            tCurr, dt, speciesVar, fieldVar, speciesBuf, fieldBuf, species,
            field, externalField.em)
	 status =  status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end

      return status, dtSuggested
   end

   -- Function to advance solution using FV dimensionally split scheme.
   function timeSteppers.fvDimSplit(tCurr, dt, tryInv)
      local status, dtSuggested = true, GKYL_MAX_DOUBLE
      local fIdx = { {1,2}, {2,1}, {1,2} } -- For indexing inp/out fields.

      -- Copy in case we need to take this step again.
      copy(3, 1)

      -- Update fluid source by half time-step.
      do
	 local myStatus, myDtSuggested = updateFluidSource(1, tCurr, dt/2)
	 status = status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end

      -- Update solution in each direction.
      local isInv = true
      for d = 1, cdim do
	 local myStatus, myDtSuggested, myTryInv = updateInDirection(d, tCurr, dt, tryInv)
	 status =  status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
	 if not status then
	    log(" ** Time step too large! Aborting this step!")
	    break
	 else
	    -- If an updated species is invalid, plan to use lax flux for THIS
	    -- species in the re-taken step.
            for nm, s in lume.orderedIter(species) do
	       if myTryInv[s] then
		  isInv = false
		  tryInv[s] = true
		  log(string.format(
			 "\n ** Invalid values in %s; Will re-update using Lax flux!\n", nm))
	       end
	    end
	    -- Break the loop if any species is invalid.
	    if not isInv then
	       break
	    end
	 end
      end
      -- Is all species is valid, do not use lax in the next step.
      if isInv then
         for _, s in lume.orderedIter(species) do tryInv[s] = false end
      end

      -- Update fluid source by half time-step.
      if status and isInv then
	 local myStatus, myDtSuggested
	 if fIdx[cdim][2] == 2 then
	    myStatus, myDtSuggested = updateFluidSource(2, tCurr, dt/2)
	 else
	    myStatus, myDtSuggested = updateFluidSource(1, tCurr, dt/2)
	 end
	 status = status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
      end

      if not (status and isInv) then
	 copy(1, 3) -- Restore old solution in case of failure.
      else
	 -- If solution not already in field[1], copy for use in next
	 -- time-step.
	 if fIdx[cdim][2] == 2 then copy(1, 2) end
      end
      
      return status, dtSuggested, isInv
   end

   local devDiagnose = function()
      -- Perform performance/numerics-related diagnostics.
      field:printDevDiagnostics()
   end

   local tmEnd = Time.clock()
   log(string.format("Initialization completed in %g sec\n\n", tmEnd-tmStart))

   -- Read some info about restarts (default is to write restarts 1/20 (5%) of sim, 
   -- but no need to write restarts more frequently than regular diagnostic output).
   local restartFrameEvery = tbl.restartFrameEvery and tbl.restartFrameEvery or math.max(1/20.0, 1/tbl.nFrame)

   -- Return function that runs main simulation loop.
   return function(self)
      log(string.format("Starting main loop of %s simulation ...\n\n", self.label))
      local tStart, tEnd = tStart, tbl.tEnd

      -- Sanity check: don't run if not needed.
      if tStart >= tEnd then return end

      local maxDt = tbl.maximumDt and tbl.maximumDt or tEnd-tStart -- max time-step
      local initDt =  tbl.suggestedDt and tbl.suggestedDt or maxDt -- initial time-step
      local step = 1
      local tCurr = tStart
      local myDt = initDt

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
		       " Step %6d at time  %#11.8g.  Time step  %.6e.  Completed %g%s\n", step, tCurr, myDt, tenth*10, "%"))
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
      local irestart = 0
      local stopfile = GKYL_OUT_PREFIX .. ".stop"

      -- For the fvDimSplit updater, tryInv contains for indicators for each
      -- species whether the domain-invariant equation should be used in the
      -- next step; they might be changed during fvDimSplit calls.
      local tryInv = {}
      for _, s in lume.orderedIter(species) do tryInv[s] = false end
      local isInv = true
      -- Main simulation loop.
      while true do
	 -- Call time-stepper.
         local status, dtSuggested
         if timeStepperNm == "fvDimSplit" then
	    status, dtSuggested, isInv = timeSteppers[timeStepperNm](tCurr, myDt, tryInv)
         else
            status, myDt = timeSteppers[timeStepperNm](tCurr)
            dtSuggested = myDt
         end
    
         -- If stopfile exists, break.
         if (file_exists(stopfile)) then
            writeData(tCurr+myDt, true)
            writeRestart(tCurr+myDt)
            break
         end

         -- Abort simulation if the suggested timestep is 0, which means there are likely NaNs.
         -- Don't write anything.
         if (myDt == 0.0) then
            log(string.format(" ERROR: dt is zero, there are likely NaNs. Terminating without writing files."))
            break
         end

	 -- Check status and determine what to do next.
	 if status and isInv then
            if first then 
               log(string.format(" Step 0 at time %g. Time step %g. Completed 0%%\n", tCurr, myDt))
               initDt = math.min(maxDt, dtSuggested); first = false
            end
            -- Track dt.
            dtPtr:data()[0] = myDt
            dtTracker:appendData(tCurr+myDt, dtPtr)
            -- Write log
	    writeLogMessage(tCurr+myDt)
	    -- We must write data first before calling writeRestart in
	    -- order not to mess up numbering of frames on a restart.
            local tmWrite = Time.clock()
	    writeData(tCurr+myDt)
            writeDataTime = writeDataTime + Time.clock() - tmWrite
	    if checkWriteRestart(tCurr+myDt) then
               local tmRestart = Time.clock()
	       writeRestart(tCurr+myDt)
               dtTracker:write(string.format("dt.bp"), tCurr+myDt, irestart)
               irestart = irestart + 1
               writeRestartTime = writeRestartTime + Time.clock() - tmRestart
	    end	    
	    
	    tCurr = tCurr + myDt
	    myDt = math.min(dtSuggested, maxDt)
	    step = step + 1
	    if (tCurr >= tEnd) then break end
	 elseif not status then
	    log (string.format(" ** Time step %g too large! Will retake with dt %g\n", myDt, dtSuggested))
	    myDt = dtSuggested
	 elseif not isInv then
	    log (string.format(" ** Invalid values detected! Will retake with dt %g\n", dtSuggested))
	    myDt = dtSuggested
	 end

         if (myDt < 1e-4*initDt) then 
            failcount = failcount + 1
            log(string.format("WARNING: Timestep dt = %g is below 1e-4*initDt. Fail counter = %d...\n", myDt, failcount))
            if failcount > 20 then
               writeData(tCurr+myDt, true)
               dtTracker:write(string.format("dt.bp"), tCurr+myDt)
               log(string.format("ERROR: Timestep below 1e-4*initDt for 20 consecutive steps. Exiting...\n"))
               break
            end
         else
            failcount = 0
         end
      end
      local tmSimEnd = Time.clock()

      -- Compute time spent in various parts of code.
      local tmSlvr = 0.0
      for _, s in lume.orderedIter(species) do
	 tmSlvr = tmSlvr+s:totalSolverTime()
      end

      local tmMom, tmIntMom, tmBc, tmColl = 0.0, 0.0, 0.0, 0.0
      local tmSrc, tmCollNonSlvr = 0.0, 0.0
      for _, s in lume.orderedIter(species) do
         tmMom = tmMom + s:momCalcTime()
         tmIntMom = tmIntMom + s:intMomCalcTime()
         tmBc = tmBc + s:totalBcTime()
         if s.collisions then
	    for _, c in pairs(s.collisions) do
	       tmColl = tmColl + c:slvrTime()
               tmCollNonSlvr = tmCollNonSlvr + c:nonSlvrTime()
	    end
         end
         if s.sources then
	    for _, src in pairs(s.sources) do
               tmSrc = tmSrc + src:srcTime()
	    end
	 end
      end

      for _, s in lume.orderedIter(fluidSources) do tmSrc = tmSrc + s:totalTime() end

      local tmTotal = tmSimEnd-tmSimStart
      local tmAccounted = 0.0
      log(string.format("\nTotal number of time-steps %s\n", step))
      --log(string.format(
	--     "Number of barriers %d barriers (%g barriers/step)\n\n",
	--     Mpi.getNumBarriers(), Mpi.getNumBarriers()/step))
      
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Solver took", tmSlvr, tmSlvr/step, 100*tmSlvr/tmTotal))
      tmAccounted = tmAccounted + tmSlvr
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Solver BCs took", tmBc, tmBc/step, 100*tmBc/tmTotal))
      tmAccounted = tmAccounted + tmBc
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Field solver took",
	     field:totalSolverTime(), field:totalSolverTime()/step, 100*field:totalSolverTime()/tmTotal))
      tmAccounted = tmAccounted + field:totalSolverTime()
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Field solver BCs",
	     field:totalBcTime(), field:totalBcTime()/step, 100*field:totalBcTime()/tmTotal))
      tmAccounted = tmAccounted + field:totalBcTime()
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Function field solver took",
	     externalField:totalSolverTime(), externalField:totalSolverTime()/step, 100*externalField:totalSolverTime()/tmTotal))
      tmAccounted = tmAccounted + externalField:totalSolverTime()
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Moment calculations took",
	     tmMom, tmMom/step, 100*tmMom/tmTotal))
      tmAccounted = tmAccounted + tmMom
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Integrated moment calculations took",
	     tmIntMom, tmIntMom/step, 100*tmIntMom/tmTotal))
      tmAccounted = tmAccounted + tmIntMom
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Field energy calculations took",
	     field:energyCalcTime(), field:energyCalcTime()/step, 100*field:energyCalcTime()/tmTotal))
      tmAccounted = tmAccounted + field:energyCalcTime()
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Collision solver(s) took",
	     tmColl, tmColl/step, 100*tmColl/tmTotal))
      tmAccounted = tmAccounted + tmColl
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Collision (other) took",
	     tmCollNonSlvr, tmCollNonSlvr/step, 100*tmCollNonSlvr/tmTotal))
      tmAccounted = tmAccounted + tmCollNonSlvr
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Source updaters took",
	     tmSrc, tmSrc/step, 100*tmSrc/tmTotal))
      tmAccounted = tmAccounted + tmSrc
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Stepper combine/copy took",
	     stepperTime, stepperTime/step, 100*stepperTime/tmTotal))
      tmAccounted = tmAccounted + stepperTime
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Forward Euler combine took",
	     fwdEulerCombineTime, fwdEulerCombineTime/step, 100*fwdEulerCombineTime/tmTotal))
      tmAccounted = tmAccounted + fwdEulerCombineTime
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
      	     "Time spent in barrier function",
      	     Mpi.getTimeBarriers(), Mpi.getTimeBarriers()/step, 100*Mpi.getTimeBarriers()/tmTotal))      
      tmUnaccounted = tmTotal - tmAccounted
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Data write took",
	     writeDataTime, writeDataTime/step, 100*writeDataTime/tmTotal))
      tmAccounted = tmAccounted + writeDataTime
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n",
	     "Write restart took",
	     writeRestartTime, writeRestartTime/step, 100*writeRestartTime/tmTotal))
      tmAccounted = tmAccounted + writeRestartTime
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.3f%%)\n\n",
	     "[Unaccounted for]",
	     tmUnaccounted, tmUnaccounted/step, 100*tmUnaccounted/tmTotal))
      
      log(string.format(
	     "%-40s %13.5f s   (%9.6f s/step)   (%6.f%%)\n\n",
	     "Main loop completed in",
	     tmTotal, tmTotal/step, 100*tmTotal/tmTotal))      
      log(date(false):fmt()); log("\n") -- Time-stamp for sim end.

      -- Perform other numerical/performance diagnostics.
      devDiagnose()

      if file_exists(stopfile) then os.remove(stopfile) end -- Clean up.
   end
end

-- PlasmaOnCartGrid application object.
local App = Proto()

function App:init(tbl)
   self._runApplication = buildApplication(self, tbl)
end

function App:getConfGrid()
   return self._confGrid
end

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
	 AdiabaticSpecies    = require ("App.Species.AdiabaticSpecies"),
	 App                 = App,
	 Field               = require ("App.Field.GkField").GkField,
	 FunctionProjection  = require ("App.Projection.GyrofluidProjection").FunctionProjection, 
	 Geometry            = require ("App.Field.GkField").GkGeometry,
	 GyrofluidProjection = require ("App.Projection.GyrofluidProjection").GyrofluidProjection, 
         PASCollisions       = require "App.Collisions.GfPitchAngleScattering",
         Source              = require "App.Sources.GyrofluidSource",
	 Species             = require "App.Species.GyrofluidSpecies",
      }
   end,

   Gyrokinetic = function ()
      App.label = "Gyrokinetic"
      return  {
	 AdiabaticSpecies       = require ("App.Species.AdiabaticSpecies"),
	 App                    = App,
	 BGKCollisions          = require "App.Collisions.GkBGKCollisions",
	 BgkCollisions          = require "App.Collisions.GkBGKCollisions",
	 ChargeExchange         = require "App.Collisions.GkChargeExchange",
	 Field                  = require ("App.Field.GkField").GkField,
	 FunctionProjection     = require ("App.Projection.GkProjection").FunctionProjection, 
	 Geometry               = require ("App.Field.GkField").GkGeometry,
	 Ionization             = require "App.Collisions.GkIonization",
	 LBOCollisions          = require "App.Collisions.GkLBOCollisions",
	 LboCollisions          = require "App.Collisions.GkLBOCollisions",
	 MaxwellianProjection   = require ("App.Projection.GkProjection").MaxwellianProjection,
	 Species                = require "App.Species.GkSpecies",
	 Source                 = require "App.Sources.GkSource",
	 Vlasov                 = require ("App.Species.VlasovSpecies"),
	 VmMaxwellianProjection = require ("App.Projection.VlasovProjection").MaxwellianProjection,
      }
   end,

   IncompEuler = function ()
      App.label = "Incompressible Euler"
      return {
	 App       = App,
	 Diffusion = require "App.Collisions.Diffusion",
	 Field     = require ("App.Field.GkField").GkField,
         Source    = require "App.Sources.FluidSource",
	 Species   = require "App.Species.IncompEulerSpecies",
      }
   end,
   
   Moments = function ()
      App.label = "Multi-fluid"
      return {
         App = App,
         Species = require "App.Species.MomentSpecies",
         Field = require ("App.Field.MaxwellField").MaxwellField,
         CollisionlessEmSource = require "App.FluidSources.CollisionlessEmSource",
         TenMomentRelaxSource  = require "App.FluidSources.TenMomentRelaxSource",
         AxisymmetricMomentSource = require "App.FluidSources.AxisymmetricMomentSource",
         AxisymmetricPhMaxwellSource = require "App.FluidSources.AxisymmetricPhMaxwellSource",
         BraginskiiHeatConductionSource = require "App.FluidSources.BraginskiiHeatConductionSource",
         BraginskiiViscosityDiffusionSource = require "App.FluidSources.BraginskiiViscosityDiffusionSource",
      }
   end,

   VlasovMaxwell = function ()
      App.label = "Vlasov-Maxwell"
      return {
	 App = App,
	 Species = require "App.Species.VlasovSpecies",
	 FuncSpecies = require "App.Species.FuncVlasovSpecies",
	 Field = require ("App.Field.MaxwellField").MaxwellField,
	 ExternalField = require ("App.Field.MaxwellField").ExternalMaxwellField,
	 FuncField = require ("App.Field.MaxwellField").ExternalMaxwellField, -- for backwards compat
	 FunctionProjection = require ("App.Projection.VlasovProjection").FunctionProjection,
	 MaxwellianProjection = require ("App.Projection.VlasovProjection").MaxwellianProjection,
	 BGKCollisions = require "App.Collisions.VmBGKCollisions",
	 LBOCollisions = require "App.Collisions.VmLBOCollisions",
	 BgkCollisions = require "App.Collisions.VmBGKCollisions",
	 LboCollisions = require "App.Collisions.VmLBOCollisions",
	 ChargeExchange = require "App.Collisions.VmChargeExchange",
	 Ionization = require "App.Collisions.VmIonization",
	 Diffusion = require "App.Collisions.Diffusion",
	 SteadySource = require "App.Sources.VmSteadyStateSource",
	 Source = require "App.Sources.VmSource",
      }
   end,
   
   Moments = function ()
      App.label = "Multi-fluid"
      return {
         App = App,
         Species = require "App.Species.MomentSpecies",
         Field = require ("App.Field.MaxwellField").MaxwellField,
         ExternalField = require ("App.Field.MaxwellField").ExternalMaxwellField,
         CollisionlessEmSource = require "App.FluidSources.CollisionlessEmSource",
         TenMomentRelaxSource  = require "App.FluidSources.TenMomentRelaxSource",
         AxisymmetricMomentSource = require "App.FluidSources.AxisymmetricMomentSource",
         AxisymmetricPhMaxwellSource = require "App.FluidSources.AxisymmetricPhMaxwellSource",
         BraginskiiHeatConductionSource = require "App.FluidSources.BraginskiiHeatConductionSource",
         BraginskiiViscosityDiffusionSource = require "App.FluidSources.BraginskiiViscosityDiffusionSource",
      }
   end
}
