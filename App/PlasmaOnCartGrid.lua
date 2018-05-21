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
local Collisions = require "App.Collisions"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Field = require "App.Field"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local LinearTrigger = require "Lib.LinearTrigger"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Species = require "App.Species"
local Time = require "Lib.Time"
local Updater = require "Updater"
local date = require "xsys.date"
local xsys = require "xsys"

-- function to create basis functions
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   end
end

-- top-level method to build application "run" method
local function buildApplication(self, tbl)
   -- create logger
   local log = Logger {
      logToFile = xsys.pickBool(tbl.logToFile, true)
   }

   log(date(false):fmt()); log("\n") -- time-stamp for sim start

   -- function to warn user about default values
   local function warnDefault(varVal, varNm, default)
      if varVal then return varVal end
      log(string.format(" ** WARNING: %s not specified, assuming %s", varNm, tostring(default)))
      return default
   end

   log("Initializing PlasmaOnCartGrid simulation ...\n")
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
   if timeStepperNm ~= "rk1" and timeStepperNm ~= "rk2" and timeStepperNm ~= "rk3" and timeStepperNm ~= "rk3s4" then
      assert(false, "Incorrect timeStepper type " .. timeStepperNm .. " specified")
   end

   -- CFL fractions for various steppers
   local stepperCFLFracs = { rk1 = 1.0, rk2 = 1.0, rk3 = 1.0, rk3s4 = 2.0 }

   local cflFrac = tbl.cflFrac
   -- Compute CFL fraction if not specified
   if  cflFrac == nil then
      cflFrac = stepperCFLFracs[timeStepperNm]
   end

   -- Number of fields needed for each stepper type
   local stepperNumFields = { rk1 = 3, rk2 = 3, rk3 = 3, rk3s4 = 4 }

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
      if Species.SpeciesBase.is(val) then
	 val:fullInit(tbl) -- initialize species
	 species[nm] = val
	 species[nm]:setName(nm)
	 species[nm]:setIoMethod(ioMethod)
      end
   end

   -- setup each species
   for _, s in pairs(species) do
      s:createGrid(tbl.lower, tbl.upper, tbl.cells, decompCuts, periodicDirs, tbl.coordinateMap)
      s:setConfBasis(confBasis)
      s:createBasis(basisNm, polyOrder)
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

   -- set conf grid for each species
   for _, s in pairs(species) do
      s:setConfGrid(grid)
      s:alloc(stepperNumFields[timeStepperNm])
   end

   local cflMin = GKYL_MAX_DOUBLE
   -- compute CFL numbers
   for _, s in pairs(species) do
      local ndim = s:getNdim()
      local myCfl = tbl.cfl and tbl.cfl or cflFrac/(2*polyOrder+1)
      cflMin = math.min(cflMin, myCfl)
      s:setCfl(cflMin)
   end


   -- read in information about collisions
   local collisions = {}
   for nm, val in pairs(tbl) do
      if Collisions.CollisionsBase.is(val) then
	 val:fullInit(tbl) -- initialize species
	 collisions[nm] = val
	 collisions[nm]:setName(nm)
	 collisions[nm]:setConfGrid(grid)
	 collisions[nm]:setConfBasis(confBasis)
	 collisions[nm]:setPhaseGrid(species)
	 collisions[nm]:setPhaseBasis(species)
	 collisions[nm]:createSolver(species)
      end
   end

   local function completeFieldSetup(fld)
      fld:fullInit(tbl) -- complete initialization
      fld:setIoMethod(ioMethod)
      fld:setBasis(confBasis)
      fld:setGrid(grid)
      do
	 local myCfl = tbl.cfl and tbl.cfl or cflFrac/(cdim*(2*polyOrder+1))
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
   -- store fields from EM field for RK time-stepper
   local emRkFields = field:rkStepperFields()

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
   -- store fields from funcField for RK time-stepper
   local emRkFuncFields = funcField:rkStepperFields()
   funcField:createSolver()
   funcField:initField()
   funcField:applyBc(0, 0, emRkFuncFields[1])
   
   -- initialize species solvers and diagnostics
   local speciesRkFields = { }
   for nm, s in pairs(species) do
      local hasE, hasB = field:hasEB()
      local funcHasE, funcHasB = funcField:hasEB()
      speciesRkFields[nm] = s:rkStepperFields()
      s:createSolver(hasE or funcHasE, hasB or funcHasB, funcField)
      s:initDist()
      s:createDiagnostics()
   end

   -- initialize field (sometimes requires species to have been initialized)
   for nm, s in pairs(species) do
      s:calcCouplingMoments(0, 0, speciesRkFields[nm][1])
   end
   field:createSolver(species)
   field:initField(species)
   field:applyBc(0, 0, emRkFields[1])

   -- apply species BCs 
   for nm, s in pairs(species) do
      -- this is a dummy forwardEuler call because some BCs require 
      -- auxFields to be set, which is controlled by species solver
      s:forwardEuler(0, 0, speciesRkFields[nm][1], {emRkFields[1], emRkFuncFields[1]}, speciesRkFields[nm][2])
      -- restore initial condition
      s:initDist()
      -- apply BCs
      s:applyBc(0, 0, speciesRkFields[nm][1])
   end

   -- function to write data to file
   local function writeData(tCurr)
      for _, s in pairs(species) do s:write(tCurr) end
      field:write(tCurr)
      funcField:write(tCurr)
   end

   writeData(0.0) -- write initial conditions

   -- determine if field equations are elliptic 
   local ellipticFieldEqn = false
   if field.isElliptic ~= nil then ellipticFieldEqn = field.isElliptic end

   -- function to take a single forward-euler time-step
   local function forwardEuler(tCurr, dt, inIdx, outIdx)
      local status, dtSuggested = true, GKYL_MAX_DOUBLE

      -- update EM field
      for nm, s in pairs(species) do
	 -- compute moments needed in coupling with fields and
	 -- collisions (the species should update internal datastructures). 
         s:calcCouplingMoments(tCurr, dt, speciesRkFields[nm][inIdx])
      end
      if ellipticFieldEqn then
        -- if field equation is elliptic, calculate field 
        -- that is self-consistent with speciesRkFields[inIdx] 
        -- to use in species update
        local myStatus, myDtSuggested = field:forwardEuler(
           tCurr, dt, emRkFields[inIdx], species, emRkFields[inIdx])
        field:applyBc(tCurr, dt, emRkFields[inIdx])
        status = status and myStatus
        dtSuggested = math.min(dtSuggested, myDtSuggested)
      else
        -- otherwise evolve field inIdx -> outIdx
        local myStatus, myDtSuggested = field:forwardEuler(
           tCurr, dt, emRkFields[inIdx], species, emRkFields[outIdx])
        field:applyBc(tCurr, dt, emRkFields[outIdx])
        status = status and myStatus
        dtSuggested = math.min(dtSuggested, myDtSuggested)
      end

      -- compute functional field (if any)
      funcField:forwardEuler(tCurr, dt, nil, nil, emRkFuncFields[1])
      
      -- update species
      for nm, s in pairs(species) do
	 local myStatus, myDtSuggested = s:forwardEuler(
	    tCurr, dt, speciesRkFields[nm][inIdx], {emRkFields[inIdx], emRkFuncFields[1]}, speciesRkFields[nm][outIdx])

	 status = status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
	 -- s:applyBc(tCurr, dt, speciesRkFields[nm][outIdx]) -- see comment below
      end
      --update species with collisions
      for _, c in pairs(collisions) do
	 local myStatus, myDtSuggested = c:forwardEuler(
	    tCurr, dt, inIdx, outIdx, species)
	 status = status and myStatus
	 dtSuggested = math.min(dtSuggested, myDtSuggested)
	 -- apply BC
         -- NRM: @Petr, does collision updater require that species BCs have been set prior to call? 
         -- NRM: if not, can remove applyBc from update species loop above and this loop below, and just have
         -- NRM: one loop over all species that sets BCs before returning (as commented below)
	 -- for _, nm in ipairs(c.speciesList) do
	 --    species[nm]:applyBc(tCurr, dt, speciesRkFields[nm][outIdx])
	 -- end
      end

      -- if collisions doesn't require BCs to be set before update, then just do BCs down here
      for nm, s in pairs(species) do
        s:applyBc(tCurr, dt, speciesRkFields[nm][outIdx])
      end

      return status, dtSuggested
   end

   -- various functions to copy/increment fields
   local function copy1(aIdx, outIdx)
      for nm, s in pairs(species) do
	 speciesRkFields[nm][outIdx]:copy(speciesRkFields[nm][aIdx])
      end
      if emRkFields[aIdx] and not ellipticFieldEqn then -- only increment EM fields if there are any
	 emRkFields[outIdx]:copy(emRkFields[aIdx])
      end
   end
   local function combine2(a, aIdx, b, bIdx, outIdx)
      for nm, s in pairs(species) do
	 speciesRkFields[nm][outIdx]:combine(a, speciesRkFields[nm][aIdx], b, speciesRkFields[nm][bIdx])
      end
      if emRkFields[aIdx] and not ellipticFieldEqn then -- only increment EM fields if there are any
	 emRkFields[outIdx]:combine(a, emRkFields[aIdx], b, emRkFields[bIdx])
      end
   end
   local function combine3(a, aIdx, b, bIdx, c, cIdx, outIdx)
      for nm, s in pairs(species) do
	 speciesRkFields[nm][outIdx]:combine(
	    a, speciesRkFields[nm][aIdx], b, speciesRkFields[nm][bIdx], c, speciesRkFields[nm][cIdx])
      end
      if emRkFields[aIdx] and not ellipticFieldEqn then -- only increment EM fields if there are any
	 emRkFields[outIdx]:combine(
	    a, emRkFields[aIdx], b, emRkFields[bIdx], c, emRkFields[cIdx])
      end
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
      copy1(2, 1)
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
      combine2(1.0/2.0, 1, 1.0/2.0, 3, 2)
      copy1(2, 1)
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
      combine2(3.0/4.0, 1, 1.0/4.0, 3, 2)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 3
      status, dtSuggested = forwardEuler(tCurr+dt/2, dt, 2, 3)
      if status == false then return status, dtSuggested end
      tm = Time.clock()
      combine2(1.0/3.0, 1, 2.0/3.0, 3, 2)
      copy1(2, 1)
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
      combine2(1.0/2.0, 1, 1.0/2.0, 2, 3)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 2
      status, dtSuggested = forwardEuler(tCurr+dt/2, dt, 3, 4)
      if status == false then return status, dtSuggested end
      tm = Time.clock()
      combine2(1.0/2.0, 3, 1.0/2.0, 4, 2)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 3
      status, dtSuggested = forwardEuler(tCurr+dt, dt, 2, 3)
      if status == false then return status, dtSuggested end
      tm = Time.clock()
      combine3(2.0/3.0, 1, 1.0/6.0, 2, 1.0/6.0, 3, 4)
      stepperTime = stepperTime + (Time.clock() - tm)

      -- RK stage 4
      status, dtSuggested = forwardEuler(tCurr+dt/2, dt, 4, 3)
      if status == false then return status, dtSuggested end
      tm = Time.clock()
      combine2(1.0/2.0, 4, 1.0/2.0, 3, 1)
      stepperTime = stepperTime + (Time.clock() - tm)

      return status, dtSuggested
   end

   local tmEnd = Time.clock()
   log(string.format("Initializing completed in %g sec\n\n", tmEnd-tmStart))

   -- return function that runs main simulation loop
   return function(self)
      log("Starting main loop of PlasmaOnCartGrid simulation ...\n\n")
      local tStart, tEnd = 0, tbl.tEnd
      local initDt =  tbl.suggestedDt and tbl.suggestedDt or tEnd-tStart -- initial time-step
      local frame = 1
      local step = 1
      local tCurr = tStart
      local myDt = initDt

      -- triggers for 10% and 1% loggers
      local logTrigger = LinearTrigger(tStart, tEnd, 10)
      local logTrigger1p = LinearTrigger(tStart, tEnd, 100)
      local tenth = 0

      local p1c = 0
      -- for writing out log messages
      local function writeLogMessage(tCurr, myDt)
	 if logTrigger(tCurr) then
	    log (string.format(
		    " Step %5d at time %g. Time step %g. Completed %g%s\n", step, tCurr, myDt, tenth*10, "%"))
	    tenth = tenth+1
	 end
	 if logTrigger1p(tCurr) then
	    log(string.format("%d", p1c))
	    p1c = (p1c+1)%10
	 end
      end

      local tmSimStart = Time.clock()
      -- main simulation loop
      while true do
	 -- if needed adjust dt to hit tEnd exactly
	 if tCurr+myDt > tEnd then myDt = tEnd-tCurr end
	 -- take a time-step
	 local status, dtSuggested = timeSteppers[timeStepperNm](tCurr, myDt)
	 -- check if step was successful
	 if status then
	    writeLogMessage(tCurr, myDt)
	    writeData(tCurr+myDt) -- give chance to everyone to write data
	    tCurr = tCurr + myDt
	    myDt = dtSuggested
	    step = step + 1
	    if (tCurr >= tEnd) then
	       break
	    end
	 else
	    log (string.format(" ** Time step %g too large! Will retake with dt %g\n", myDt, dtSuggested))
	    myDt = dtSuggested
	 end
      end -- end of time-step loop
      writeLogMessage(tCurr, myDt)
      local tmSimEnd = Time.clock()

      -- compute time spent in various parts of code
      local tmSlvr = 0.0 -- total time in ptcl solver
      local tmVol = 0.0 -- total time in ptcl solver vol terms
      local tmSurf = 0.0 -- total time in ptcl solver surf terms
      for _, s in pairs(species) do
	 tmSlvr = tmSlvr+s:totalSolverTime()
         if s.solverVolTime then tmVol = tmVol+s:solverVolTime() end
         if s.solverSurfTime then tmSurf = tmSurf+s:solverSurfTime() end
      end

      local tmMom, tmIntMom, tmBc = 0.0, 0.0, 0.0
      for _, s in pairs(species) do
         tmMom = tmMom + s:momCalcTime()
         tmIntMom = tmIntMom + s:intMomCalcTime()
         tmBc = tmBc + s:totalBcTime()
      end
      local tmColl, tmCollEvalMom, tmCollProjectMaxwell = 0.0, 0.0, 0.0
      for _, c in pairs(collisions) do
         tmColl = tmColl + c:totalSolverTime()
         tmCollEvalMom = tmCollEvalMom + c:evalMomTime()
         tmCollProjectMaxwell = tmCollProjectMaxwell + c:projectMaxwellTime()
      end

      log(string.format("\nTotal number of time-steps %s\n", step))
      log(string.format("\nSolver took %g sec\n", tmSlvr))
      log(string.format("  [Vol kernels %g sec. Surf kernels %g sec. Surf-to-Vol ratio = %g]\n", tmVol, tmSurf, tmSurf/tmVol))
      log(string.format("Solver BCs took %g sec\n", tmBc))
      log(string.format("Field solver took %g sec\n", field:totalSolverTime()))
      log(string.format("Field solver BCs took %g sec\n", field:totalBcTime()))
      log(string.format("Function field solver took %g sec\n", funcField:totalSolverTime()))
      log(string.format("Moment calculations took %g sec\n", tmMom))
      log(string.format("Integrated moment calculations took %g sec\n", tmIntMom))
      log(string.format("Field energy calculations took %g sec\n", field:energyCalcTime()))
      log(string.format("Collision solver took %g sec\n", tmColl))
      log(string.format(
	     "  [Moment evaluation %g sec. Maxwellian projection %g sec]\n",
	     tmCollEvalMom, tmCollProjectMaxwell))
      log(string.format("Stepper combine/copy took %g sec\n", stepperTime))
      log(string.format("Main loop completed in %g sec\n\n", tmSimEnd-tmSimStart))
      log(date(false):fmt()); log("\n") -- time-stamp for sim end
   end
end

-- PlasmaOnCartGrid application object
local App = Proto()

function App:init(tbl)
   self._runApplication = buildApplication(self, tbl)
end

function App:run()
   if GKYL_COMMAND == "run" then
      return self:_runApplication()
   elseif GKYL_COMMAND == "init" then
      return function (...) end
   end
end

return {
   AdiabaticSpecies = Species.AdiabaticSpecies,
   App = App,
   BgkCollisions = Collisions.BgkCollisions,   
   FuncMaxwellField = Field.FuncMaxwellField,
   GkField = Field.GkField,
   GkGeometry = Field.GkGeometry,
   GkSpecies = Species.GkSpecies,
   HamilVlasovSpecies = Species.HamilVlasovSpecies,
   IncompEulerSpecies = Species.IncompEulerSpecies,
   MaxwellField = Field.MaxwellField,
   NoField = Field.NoField,
   VlasovSpecies = Species.VlasovSpecies,
   VoronovIonization = Collisions.VoronovIonization,

   -- valid pre-packaged species-field systems
   Gyrokinetic = {Species = Species.GkSpecies, Field = Field.GkField, Geometry = Field.GkGeometry},
   IncompEuler = {Species = Species.IncompEulerSpecies, Field = Field.GkField},
   VlasovMaxwell = {Species = Species.VlasovSpecies, Field = Field.MaxwellField},
}
