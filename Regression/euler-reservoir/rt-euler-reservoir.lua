-- Gkyl ------------------------------------------------------------------------
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Updater = require "Updater"
local Euler = require "Eq.Euler"
local BoundaryCondition = require "Updater.BoundaryCondition"
local Time = require "Lib.Time"
local DecompRegionCalc = require "Lib.CartDecomp"
local Mpi = require "Comm.Mpi"

-- MPI rank we are on
local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)

local function log(msg)
   if rank == 0 then 
      io.write(msg); io.write("\n")
   end
end

gasGamma = 3.0

-- resolution and time-stepping
Lx = 1.0
NX = 128
cfl = 0.9
tStart = 0.0
tEnd = 50
nFrames = 9

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
decomp = DecompRegionCalc.CartProd { cuts = {1} }
-- computational domain
grid = Grid.RectCart {
   lower = {0.0},
   upper = {Lx},
   cells = {NX},
   periodicDirs = {},
   decomposition = decomp,
}

-- solution
fluid = DataStruct.Field {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}
-- final updated solution
fluidNew = DataStruct.Field {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
fluidDup = DataStruct.Field {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}

-----------------------
-- INITIAL CONDITION --
-----------------------
-- initial conditions
function init(field)
   local xlo, dx = grid:lower(1), grid:dx(1)
   
   local range, indexer = field:localRange(), field:genIndexer()
   for idx in range:colMajorIter() do
      local fItr = field:get(indexer(idx)) -- pointer to data in cell
      local x = xlo + idx[1]*dx
      if x < 0.5 then
	 fItr[1] = 1.0
	 fItr[5] = 1.0/(gasGamma-1)
      else
	 fItr[1] = 0.125
	 fItr[5] = 0.1/(gasGamma-1)
      end
      fItr[2] = 0.0
      fItr[3] = 0.0
      fItr[4] = 0.0
   end
end

------------------------
-- Boundary Condition --
------------------------

function applyPeriodic(fIn)
   local localRange = fIn:localRange()
   local indexer = fIn:genIndexer()

   for dir = 1,1 do
      local perpRange = localRange:shorten(dir) -- range orthogonal to 'dir'
      for idx in perpRange:colMajorIter() do
	 local idxGst, idxSkn = idx:copy(), idx:copy()

	 idxSkn[dir] = localRange:lower(dir)
	 idxGst[dir] = localRange:upper(dir)+1
	 local fGst, fSkn = fIn:get(indexer(idxGst)), fIn:get(indexer(idxSkn))
	 fGst[1] = fSkn[1]

	 idxSkn[dir] = localRange:upper(dir)
	 idxGst[dir] = localRange:lower(dir)-1
	 fGst, fSkn = fIn:get(indexer(idxGst)), fIn:get(indexer(idxSkn))
	 fGst[1] = fSkn[1]
      end
   end
end

-- boundary applicator objects for fluids and fields
bcConstL = BoundaryCondition.Const { components = {1,2,3,4,5},
				     values = {1.0, 0.0, 0.0, 0.0, 1.0/(gasGamma-1)} }
bcConstR = BoundaryCondition.Const { components = {1,2,3,4,5},
				     values = {0.125, 0.0, 0.0, 0.0, 0.1/(gasGamma-1)} }

bcLeft = Updater.Bc {
   onGrid = grid,
   boundaryConditions = {bcConstL},
   dir = 1,
   edge = "lower",
}
bcRight = Updater.Bc {
   onGrid = grid,
   boundaryConditions = {bcConstR},
   dir = 1,
   edge = "upper",
}

-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt)
   --applyPeriodic(fld)
   local bcList = {bcLeft, bcRight}
   for _, bc in ipairs(bcList) do
      bc:advance(tCurr, myDt, {}, {fld})
   end
end

----------------------
-- EQUATION SOLVERS --
----------------------
-- Burgers equations
eulerEqn = Euler { gasGamma = gasGamma }

-- ds solvers for regular Euler equations along X
fluidSlvrDir1 = Updater.WavePropagation {
   onGrid = grid,
   equation = eulerEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   updateDirections = {1}
}

-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t)
   local myStatus, useLaxSolver = true, false
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local xdirFluidSolver = {fluidSlvrDir1, {fluid}, {fluidNew}}

   -- X-direction updates
   for _, slvr in ipairs({xdirFluidSolver}) do
      local status, dtSuggested = slvr[1]:advance(tCurr, t-tCurr, slvr[2], slvr[3])
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end
   if (myStatus == false) then
      return myStatus, myDtSuggested, useLaxSolver
   end

   return myStatus, myDtSuggested, useLaxSolver
end

-- function to take one time-step with Euler solver
function solveTwoFluidSystem(tCurr, t)
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t)
   applyBc(fluidNew, tCurr, t-tCurr)
   return status, dtSuggested,useLaxSolver
end

----------------------------
-- TIME-STEPPING FUNCTION --
----------------------------
function runSimulation(tStart, tEnd, nFrames, initDt)
   local frame = 1
   local tFrame = (tEnd-tStart)/nFrames
   local nextIOt = tFrame
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested
   local useLaxSolver = false

   -- the grand loop 
   while true do
      -- copy fields in case we need to take this step again
      fluidDup:copy(fluid)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then myDt = tEnd-tCurr end

      -- advance fluids and fields
      if (useLaxSolver) then
        -- call Lax solver if positivity violated
        log (string.format(" Taking step %5d at time %6g with dt %g (using Lax solvers)", step, tCurr, myDt))
        status, dtSuggested = solveTwoFluidLaxSystem(tCurr, tCurr+myDt)
        useLaxSolver = false
      else
        log (string.format(" Taking step %5d at time %6g with dt %g", step, tCurr, myDt))
        status, dtSuggested, useLaxSolver = solveTwoFluidSystem(tCurr, tCurr+myDt)
      end

      if (status == false) then
        -- time-step too large
        log (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
        myDt = dtSuggested
        fluid:copy(fluidDup)
      elseif (useLaxSolver == true) then
        -- negative density/pressure occured
        log (string.format(" ** Negative pressure or density at %8g! Will retake step with Lax fluxes", tCurr+myDt))
        fluid:copy(fluidDup)
      else
        -- check if a nan occured
        if false then -- fluidNew:hasNan()
           log (string.format(" ** NaN occured at %g! Stopping simulation", tCurr))
           break
        end
        -- copy updated solution back
        fluid:copy(fluidNew)
     
        -- write out data
        if (tCurr+myDt > nextIOt or tCurr+myDt >= tEnd) then
           log (string.format(" Writing data at time %g (frame %d) ...\n", tCurr+myDt, frame))
	   fluid:write(string.format("fluid_%d.bp", frame), tCurr+myDt)
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
   return dtSuggested
end

----------------------------
-- RUNNING THE SIMULATION --
----------------------------
t1 = Time.clock()

-- setup initial condition
init(fluid)
applyBc(fluid, 0.0, 0.0)
fluid:write(string.format("fluid_%d.bp", 0), 0.0)
-- run simulation
initDt = 0.1
runSimulation(tStart, tEnd, nFrames, initDt)

-- Print timings
log(string.format("WavePropagation updaters took %g", fluidSlvrDir1.totalTime))
log(string.format("Simulation took %g secs\n", Time.clock()-t1))
