-- Gkyl ------------------------------------------------------------------------
--
-- Euler solver on a Cartesian grid. Works in 1D, 2D and 3D.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local BoundaryCondition = require "Updater.BoundaryCondition"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local Eq = require "Eq"
local Logger = require "Lib.Logger"
local Time = require "Lib.Time"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"

-- For returning module table
local M = {}

-- Boundary condition objects
M.bcCopy = 1
M.bcWall = 2

-- top-level method to build simulation "run" method
local function buildSimulation(self, tbl)
   -- create logger
   local log = Logger {
      logToFile = tbl.logToFile and tbl.logToFile or false
   }

   -- function to warn user about default values
   local function warnDefault(varVal, varNm, default)
      if varVal then return varVal end
      log(string.format(" ** WARNING: %s not specified, assuming %s", varNm, tostring(default)))
      return default
   end

   log("Initializing EulerOnCartGrid simulation ...")
   local tmStart = Time.clock()

   -- basic parameters and checks
   local ndim = #tbl.lower -- simulation dimension
   assert(ndim == #tbl.upper, "upper should have exactly " .. ndim .. " entries")
   assert(ndim == #tbl.cells, "cells should have exactly " .. ndim .. " entries")

   -- Gas adiabatic constant
   local gasGamma = warnDefault(tbl.gasGamma, "gasGamma", 1.4)
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

   -- create decomposition and grid
   local decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      useShared = useShared
   }
   -- computational domain
   local grid = Grid.RectCart {
      lower = tbl.lower,
      upper = tbl.upper,
      cells = tbl.cells,
      periodicDirs = periodicDirs,
      decomposition = decomp,
   }

   -- allocate space for use in dimensional sweeps
   local field = {}
   for d = 1, 2 do -- we only two fields for dimensional split algorithm
      field[d] = DataStruct.Field {
	 onGrid = grid,
	 numComponents = 5,
	 ghost = {2, 2},
      }
   end
   -- in case we need to take time-step again
   local fieldDup = DataStruct.Field {
      onGrid = grid,
      numComponents = 5,
      ghost = {2, 2},
   }

   -- equation object: provides Reimann solver
   local eulerEqn = Eq.Euler {
      gasGamma = gasGamma
   }
   -- create solvers for updates in each direction
   local fluidSlvr = {}
   for d = 1, ndim do
      fluidSlvr[d] = Updater.WavePropagation {
	 onGrid = grid,
	 equation = eulerEqn,
	 limiter = limiter,
	 cfl = cfl,
	 updateDirections = {d}
      }
   end

   -- set flags to indicate which directions are periodic
   local isDirPeriodic = {false, false, false}
   for _, d in ipairs(periodicDirs) do isDirPeriodic[d] = true end

   -- objects to apply BCs
   local bcCopyAll = BoundaryCondition.Copy { components = {1, 2, 3, 4, 5} }
   local bcWallCopy = BoundaryCondition.Copy { components = {1, 5} }
   local bcWallZeroNormal = BoundaryCondition.ZeroNormal { components = {2, 3, 4} }

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
      if bcType == M.bcCopy then
	 table.insert(boundaryConditions, makeBcUpdater(dir, edge, {bcCopyAll}))
      elseif bcType == M.bcWall then
	 table.insert(boundaryConditions, makeBcUpdater(dir, edge, {bcWallCopy, bcWallZeroNormal}))
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
   end

   -- function to update fluids
   local function updateFluid(tCurr, t)
      local status, useLaxSolver = true, false
      local suggestedDt = GKYL_MAX_DOUBLE
      local fIdx = { {1,2}, {2,1}, {1,2} } -- for indexing inp/out fields
      
      -- update solution in each direction using dimensional splitting
      for d = 1, ndim do
	 local inpField, outField = field[fIdx[d][1]], field[fIdx[d][2]]
	 local myStatus, mySuggestedDt = fluidSlvr[d]:advance(tCurr, t-tCurr, {inpField}, {outField})
	 
	 status =  status and myStatus
	 suggestedDt = math.min(suggestedDt, mySuggestedDt)
	 if not status then
	    return status, suggestedDt, useLaxSolver
	 end
	 applyBc(outField, tCurr, t)
      end
       -- if solution is not already in field[1], copy for use in next time-step
      if fIdx[ndim][2] == 2 then
	 field[1]:copy(field[2])
      end
      return status, suggestedDt, useLaxSolver
   end

   -- create Adios object for field I/O
   local fieldIo = AdiosCartFieldIo { field[1]:elemType() }

   -- function to apply initial conditions
   local function init(fld)
      local xn = Lin.Vec(ndim)
      local indexer = fld:genIndexer()
      for idx in fld:localRangeIter() do
	 -- get cell-center coordinates
	 grid:setIndex(idx)
	 grid:cellCenter(xn)
	 -- evaluate supplied IC function
	 local rho, rhou, rhov, rhow, E = tbl.init(0.0, xn)

	 -- set values in cell
	 local fItr = fld:get(indexer(idx)) -- pointer to data in cell
	 fItr[1], fItr[2], fItr[3], fItr[4], fItr[5] = rho, rhou, rhov, rhow, E
      end
   end

   init(field[1]) -- initialize field
   applyBc(field[1], 0.0, 0.0) -- apply Bc to initial conditions

   -- write out ICs
   fieldIo:write(field[1], "fluid_0.bp", 0.0)

   local tmEnd = Time.clock()
   log(string.format("Initializing completed in %g sec\n", tmEnd-tmStart))

   -- return function that the main simulation loop
   return function (self)
      log("Starting main loop of EulerOnCartGrid simulation ...")
      local tmSimStart = Time.clock()


      local tStart, tEnd, nFrame = 0, tbl.tEnd, tbl.nFrame
      local initDt =  tbl.suggestedDt and tbl.suggestedDt or tEnd-tStart -- initial time-step
      local frame = 1
      local tFrame = (tEnd-tStart)/nFrame
      local nextIOt = tFrame
      local step = 1
      local tCurr = tStart
      local myDt = initDt

      -- main simulation loop
      while true do
	 -- copy fields in case we need to take this step again
	 fieldDup:copy(field[1])

	 -- if needed adjust dt to hit tEnd exactly
	 if tCurr+myDt > tEnd then myDt = tEnd-tCurr end

	 -- Take a time-step
	 log (string.format(" Taking step %5d at time %6g with dt %g", step, tCurr, myDt))
	 local status, dtSuggested, useLaxSolver = updateFluid(tCurr, tCurr+myDt)

	 if not status then
	    -- updater failed, time-step too large
	    log (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	    myDt = dtSuggested
	    field[1]:copy(fieldDup)
	 else
	    -- write out data if needed
	    if (tCurr+myDt > nextIOt or tCurr+myDt >= tEnd) then
	       log (string.format(" Writing data at time %g (frame %d) ...\n", tCurr+myDt, frame))
	       fieldIo:write(field[1], string.format("fluid_%d.bp", frame), tCurr+myDt)
	       frame = frame + 1
	       nextIOt = nextIOt + tFrame
	       step = 0
	    end
	    
	    tCurr = tCurr + myDt
	    myDt = dtSuggested
	    step = step + 1

	    -- check if done
	    if (tCurr >= tEnd) then
	       break
	    end
	 end 
      end -- end of time-step loop

      local tmSimEnd = Time.clock()
      local tmFluidSlvr = 0 -- total time in fluid solver
      for d = 1, ndim do
	 tmFluidSlvr = tmFluidSlvr+fluidSlvr[d].totalTime
      end
      log(string.format("Fluid solvers took %g sec", tmFluidSlvr))
      log(string.format("Main loop completed in %g sec\n", tmSimEnd-tmSimStart))
   end
end

-- Euler simulation object
local Sim = {}
-- constructor
function Sim:new(tbl)
   local self = setmetatable({}, Sim)
   self._runSimulation = buildSimulation(self, tbl)
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(Sim, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
Sim.__index = {
   run = function (self)
      return self:_runSimulation()
   end
}

-- add to table
M.Sim = Sim

return M
