-- Gkyl ------------------------------------------------------------------------
--
-- Hyperbolic solver on a Cartesian grid. Works in 1D, 2D and 3D.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local BoundaryCondition = require "Updater.BoundaryCondition"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local Proto = require "Proto"
local Time = require "Lib.Time"
local Updater = require "Updater"
local date = require "xsys.date"
local xsys = require "xsys"

-- For returning module table
local M = {}

-- Boundary condition types
local bcCopyId = 1
local bcConstId = 2
local bcCustomId = 3

-- Boundary condition objects
M.bcCopy = { id = bcCopyId }
M.bcConst = function (tbl)
   local bc = { id = bcConstId }
   bc.components = tbl.components
   bc.values = tbl.values
   assert(#bc.components == #bc.values, "In HyperEqn.bcConst 'components' and 'values' must have same size")
   return bc
end
M.bcCustom = function(tbl)
   local bc = { id = bcCustomId, bcList = tbl }
   return bc
end

-- top-level method to build application "run" method
local function buildApplication(self, tbl)
   -- create logger
   local log = Logger {
      logToFile = xsys.pickBool(tbl.logToFile, false)
   }

   log(date(false):fmt()); log("\n") -- time-stamp for sim start
   
   -- function to warn user about default values
   local function warnDefault(varVal, varNm, default)
      if varVal then return varVal end
      log(string.format(" ** WARNING: %s not specified, assuming %s\n", varNm, tostring(default)))
      return default
   end

   log("Initializing HyperEqnOnCartGrid simulation ...\n")
   local tmStart = Time.clock()

   -- basic parameters and checks
   local ndim = #tbl.lower -- simulation dimension
   assert(ndim == #tbl.upper, "upper should have exactly " .. ndim .. " entries")
   assert(ndim == #tbl.cells, "cells should have exactly " .. ndim .. " entries")

   -- CFL number
   local cfl = warnDefault(tbl.cfl, "cfl", 0.9)
   -- limiter
   local limiter = warnDefault(tbl.limiter, "limiter", "monotonized-centered")

   -- parallel decomposition stuff
   local decompCuts = tbl.decompCuts
   if tbl.decompCuts then
      assert(ndim == #tbl.decompCuts, "decompCuts should have exactly " .. ndim .. " entries")
   else
      -- if not specified, use 1 processor
      decompCuts = { }
      for d = 1, ndim do decompCuts[d] = 1 end
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

   -- create decomposition and grid
   local decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      useShared = useShared,
   }
   -- computational domain
   local grid = Grid.RectCart {
      lower = tbl.lower,
      upper = tbl.upper,
      cells = tbl.cells,
      periodicDirs = periodicDirs,
      decomposition = decomp,
   }

   local hyperEqn
   -- equation object: provides Reimann solver
   if tbl.equation then
      hyperEqn = tbl.equation and tbl.equation
   else
      assert("HyperEqn: Must specify equation to solve using 'equation'!")
   end
   local meqn = hyperEqn.numEquations()

   -- allocate space for use in dimensional sweeps
   local field = {}
   for d = 1, 2 do -- we only two fields for dimensional split algorithm
      field[d] = DataStruct.Field {
	 onGrid = grid,
	 numComponents = meqn,
	 ghost = {2, 2},
      }
   end
   -- in case we need to take time-step again
   local fieldDup = DataStruct.Field {
      onGrid = grid,
      numComponents = meqn,
      ghost = {2, 2},
   }

   -- create solvers for updates in each direction
   local hyperSlvr = {}
   for d = 1, ndim do
      hyperSlvr[d] = Updater.WavePropagation {
	 onGrid = grid,
	 equation = hyperEqn,
	 limiter = limiter,
	 cfl = cfl,
	 updateDirections = {d}
      }
   end

   -- set flags to indicate which directions are periodic
   local isDirPeriodic = {false, false, false}
   for _, d in ipairs(periodicDirs) do isDirPeriodic[d] = true end

   -- initialize component list
   local cList = {}
   for m = 1, meqn do cList[m] = m end
   -- for copy BC
   local bcCopyAll = BoundaryCondition.Copy { components = cList }

   -- function to construct a BC updater
   local function makeBcUpdater(dir, edge, bcList)
      return Updater.Bc {
	 onGrid = grid,
	 boundaryConditions = bcList,
	 dir = dir,
	 edge = edge,
      }
   end

   local boundaryConditions = { } -- list of Bcs to apply
   -- function to determine what BC to apply and insert into list
   local function appendBoundaryConditions(dir, edge, bcType)
      if bcType.id == bcCopyId then
	 table.insert(boundaryConditions, makeBcUpdater(dir, edge, {bcCopyAll}))
      elseif bcType.id == bcConstId then
	 local bc = BoundaryCondition.Const { components = bcType.components, values = bcType.values  }
	 table.insert(boundaryConditions, makeBcUpdater(dir, edge, { bc }))
      elseif bcType.id == bcCustomId then
	 -- equation specific BCs
	 table.insert(boundaryConditions, makeBcUpdater(dir, edge, bcType.bcList))
      end
   end

   -- determine BCs
   if not isDirPeriodic[1] then -- X
      appendBoundaryConditions(1, "lower", tbl.bcx[1])
      appendBoundaryConditions(1, "upper", tbl.bcx[2])
   end
   if ndim > 1 then -- Y
      if not isDirPeriodic[2] then
	 appendBoundaryConditions(2, "lower", tbl.bcy[1])
	 appendBoundaryConditions(2, "upper", tbl.bcy[2])
      end
   end
   if ndim > 2 then -- Z
      if not isDirPeriodic[3] then
	 appendBoundaryConditions(3, "lower", tbl.bcz[1])
	 appendBoundaryConditions(3, "upper", tbl.bcz[2])
      end
   end

   -- function to apply boundary conditions
   local function applyBc(fld, tCurr, t)
      for _, bc in ipairs(boundaryConditions) do
	 bc:advance(tCurr, t-tCurr, {}, {fld})
      end
      fld:sync()
   end

   -- construct list of diagnostics to compute
   local diagnostics = { }
   if tbl.diagnostics then
      for _, d in ipairs(tbl.diagnostics) do
	 table.insert(diagnostics, {
			 name = d.name, dynFunc = d.diagnostic,
			 dynVec = DataStruct.DynVector { numComponents = 1 }
	 })
      end
   end

   local tmDiag = 0.0
   -- function to compute diagnostics
   local function calcDiagnostics(fld, tCurr, t)
      if #diagnostics == 0 then return end -- no diagnostics to compute

      local ndv = #diagnostics
      local dv = Lin.Vec(ndv)
      local tms = Time.clock()
      for i = 1, ndv do dv[i] = 0.0 end

      local fItr = fld:get(1)
      local indexer = fld:genIndexer()
      -- loop over field, updating diagnostic
      for idx in fld:localRangeIter() do
	 grid:setIndex(idx)
	 local vol = grid:cellVolume()
	 fld:fill(indexer(idx), fItr) -- pointer to data in cell
	 for i = 1, ndv do
	    dv[i] = dv[i] + vol*diagnostics[i].dynFunc(t, fItr)
	 end
      end

      local dvGlobal = Lin.Vec(#diagnostics)
      -- all-reduce across global communicator. NOTE: If DynVector
      -- storage changes to use shared memory on shared nodes, code
      -- below will need modification
      Mpi.Allreduce(dv:data(), dvGlobal:data(), ndv, Mpi.DOUBLE, Mpi.SUM, grid:commSet().comm)
      for i = 1, ndv do
	 diagnostics[i].dynVec:appendData(t, { dvGlobal[i] })
      end
      tmDiag = tmDiag+(Time.clock()-tms)
   end

   -- function to write diagnostics
   local function writeDiagnostics(frame, t)
      for _, d in ipairs(diagnostics) do
	 d.dynVec:write(string.format("%s_%d.bp", d.name, frame), t)
      end
   end

   -- function to update hyperbolic equations
   local function updateHyper(tCurr, t)
      local status = true
      local suggestedDt = GKYL_MAX_DOUBLE
      local fIdx = { {1,2}, {2,1}, {1,2} } -- for indexing inp/out fields

      -- update solution in each direction using dimensional splitting
      for d = 1, ndim do
	 local inpField, outField = field[fIdx[d][1]], field[fIdx[d][2]]
	 local myStatus, mySuggestedDt = hyperSlvr[d]:advance(tCurr, t-tCurr, {inpField}, {outField})

	 status =  status and myStatus
	 suggestedDt = math.min(suggestedDt, mySuggestedDt)
	 if not status then
	    return status, suggestedDt
	 end
	 applyBc(outField, tCurr, t)
      end
       -- if solution is not already in field[1], copy for use in next time-step
      if fIdx[ndim][2] == 2 then
	 field[1]:copy(field[2])
      end
      return status, suggestedDt
   end

   local ioMethod = tbl.ioMethod and tbl.ioMethod or "MPI"
   if ioMethod ~= "POSIX" and ioMethod ~= "MPI" then
      assert(false, "ioMethod must be one of 'MPI' or 'POSIX'. Provided '" .. ioMethod .. "' instead")
   end
   
   -- create Adios object for field I/O
   local fieldIo = AdiosCartFieldIo {
      elemType = field[1]:elemType(),
      method = ioMethod,
   }

   -- function to apply initial conditions
   local function init(fld)
      local xn = Lin.Vec(ndim)
      local fItr = fld:get(1)      
      local indexer = fld:genIndexer()
      for idx in fld:localExtRangeIter() do
	 grid:setIndex(idx)
	 grid:cellCenter(xn)
	 -- evaluate supplied IC function
	 local v = { tbl.init(0.0, xn) } -- braces around function put return values in table
	 -- set values in cell
	 fld:fill(indexer(idx), fItr) -- pointer to data in cell
	 for c = 1, meqn do fItr[c] = v[c] end
      end
   end

   init(field[1]) -- initialize field
   applyBc(field[1], 0.0, 0.0) -- apply Bc to initial conditions

   -- write out ICs
   fieldIo:write(field[1], "field_0.bp", 0.0)

   -- compute diagnostics from ICs and write them out
   calcDiagnostics(field[1], 0.0, 0.0)
   writeDiagnostics(0, 0.0)

   local tmEnd = Time.clock()
   log(string.format("Initializing completed in %g sec\n\n", tmEnd-tmStart))

   -- return function that runs main simulation loop
   return function (self)
      log("Starting main loop of HyperEqnOnCartGrid simulation ...\n")
      local tStart, tEnd, nFrame = 0, tbl.tEnd, tbl.nFrame
      local initDt =  tbl.suggestedDt and tbl.suggestedDt or tEnd-tStart -- initial time-step
      local frame = 1
      local tFrame = (tEnd-tStart)/nFrame
      local nextIOt = tFrame
      local step = 1
      local tCurr = tStart
      local myDt = initDt
      local gcStatFile = io.open("gcMemHist.txt", "w")

      local tmSimStart = Time.clock()
      -- main simulation loop
      while true do
	 -- copy fields in case we need to take this step again
	 fieldDup:copy(field[1])

	 -- if needed adjust dt to hit tEnd exactly
	 if tCurr+myDt > tEnd then myDt = tEnd-tCurr end
	 log (string.format(" Taking step %5d at time %6g with dt %g\n", step, tCurr, myDt))
	 local status, dtSuggested = updateHyper(tCurr, tCurr+myDt) -- take a time-step

	 if not status then
	    -- updater failed, time-step too large
	    log (string.format(" ** Time step %g too large! Will retake with dt %g\n", myDt, dtSuggested))
	    myDt = dtSuggested
	    field[1]:copy(fieldDup)
	 else
	    calcDiagnostics(field[1], tCurr, tCurr+myDt)
	    -- write out data if needed
	    if tCurr+myDt > nextIOt or tCurr+myDt >= tEnd then
	       log (string.format(" Writing data at time %g (frame %d) ...\n\n", tCurr+myDt, frame))
	       fieldIo:write(field[1], string.format("field_%d.bp", frame), tCurr+myDt)
	       writeDiagnostics(frame, tCurr+myDt)
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

	 gcStatFile:write(string.format("%d\n", collectgarbage("count")))
      end -- end of time-step loop
      local tmSimEnd = Time.clock()

      local tmHyperSlvr = 0.0 -- total time in hyper solver
      for d = 1, ndim do
	 tmHyperSlvr = tmHyperSlvr+hyperSlvr[d].totalTime
      end
      local tmBC = 0.0 -- total time in BCs
      for _, bc in ipairs(boundaryConditions) do
	 tmBC = tmBC+bc.totalTime
      end
      log(string.format("Hyper solvers took %g sec\n", tmHyperSlvr))
      log(string.format("Diagnostics took %g sec\n", tmDiag))
      log(string.format("Boundary conditions took %g sec\n", tmBC))
      log(string.format("Main loop completed in %g sec\n\n", tmSimEnd-tmSimStart))
      log(date(false):fmt()); log("\n") -- time-stamp for sim end
   end
end

-- HyperEqn application object
local App = Proto()

-- methods
function App:init(tbl)
   self._runApplication = buildApplication(self, tbl)
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

-- add to table
M.App = App

return M
