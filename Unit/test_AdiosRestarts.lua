-- Gkyl ------------------------------------------------------------------------
--
-- Test for ADIOS-based restarts.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Plasma      = require("App.PlasmaOnCartGrid").VlasovMaxwell
local lfs         = require "lfs"
local AdiosReader = require "Io.AdiosReader"
local Mpi         = require "Comm.Mpi"

local function log(msg)
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   if rank == 0 then
      print(msg)
   end
end

local verboseLog = function (msg) log(msg) end

local function twoStreamP2(endTime, frames, decomp, isRestart)
   if isRestart then
      GKYL_COMMANDS[1] = "restart"
   end
   -- This is the same two-stream p2 simulation as in Regression.
   knumber      = 0.5    -- Wave-number.
   elVTerm      = 0.2    -- Electron thermal velocity.
   vDrift       = 1.0    -- Drift velocity.
   perturbation = 1.0e-6 -- Distribution function perturbation.
   elVTerm      = 0.2    -- Electron thermal velocity.
   
   vlasovApp = Plasma.App {
      logToFile = false,
   
      tEnd        = endTime,            -- End time.
      nFrame      = frames,             -- Number of output frames.
      lower       = {-math.pi/knumber}, -- Configuration space lower left.
      upper       = { math.pi/knumber}, -- Configuration space upper right.
      cells       = {64},               -- Configuration space cells.
      basis       = "serendipity",      -- One of "serendipity" or "maximal-order".
      polyOrder   = 2,                  -- Polynomial order.
      timeStepper = "rk3",              -- One of "rk2" or "rk3".
   
      -- Decomposition for configuration space.
      decompCuts = decomp,   -- Cuts in each configuration direction.
      useShared  = false,    -- If to use shared memory.
   
      -- Boundary conditions for configuration space
      periodicDirs = {1}, -- Periodic directions.
   
      -- Electrons.
      elc = Plasma.Species {
         charge = -1.0, mass = 1.0,
         -- Velocity space grid.
         lower = {-6.0},
         upper = {6.0},
         cells = {32},
         -- Initial conditions.
         init = function (t, xn)
            local x, v  = xn[1], xn[2]
            local alpha = perturbation
            local k     = knumber
            local vt    = elVTerm
   
            local fv = 1/math.sqrt(8*math.pi*vt^2)*(math.exp(-(v-vDrift )^2/(2*vt^2))+math.exp(-(v+vDrift)^2/(2*vt^2)))
            return (1+alpha*math.cos(k*x))*fv
         end,
         evolve = true, -- Evolve species?
   
         diagnosticMoments = { "M0", "M1i", "M2" }
      },
   
      -- Field solver.
      field = Plasma.Field {
         epsilon0 = 1.0, mu0 = 1.0,
         init = function (t, xn)
            local alpha = perturbation
            local k     = knumber
            return -alpha*math.sin(k*xn[1])/k, 0.0, 0.0, 0.0, 0.0, 0.0
         end,
         evolve = true, -- Evolve field?
      },
   }
   -- Run application.
   vlasovApp:run()
end

-- (from Tool/runregression.lua) function to compare floats: the comparison is normalized to the
-- maximum value of the field being compared. Perhaps this is too
-- "coarse" but a direct comparison of floats is very tricky.
local function check_equal_numeric(expected, actual, maxVal)
   if maxVal < GKYL_MIN_DOUBLE then
      return math.max(expected-actual) > 10*GKYL_MIN_DOUBLE
   end
   if math.max(expected-actual)/maxVal > 1e-12 then
      return false
   end
   return true
end

-- (from Tool/runregression.lua) relative difference between two numbers (NOT SURE IF THIS IS BEST
-- WAY TO DO THINGS)
local function get_relative_numeric(expected, actual, maxVal)
   if maxVal < 1e-15 then
      return math.max(expected-actual)
   else
      return math.abs(expected-actual)/maxVal
   end
end

-- (from Tool/runregression.lua) Calculates maximum value in supplied field.
local function maxValueInField(fld)
   local maxVal = 0.0
   for i = 1, fld:size() do
      maxVal = math.max(maxVal, math.abs(fld[i]))
   end
   return maxVal
end

-- Function to compare files (from Tool/runregression.lua).
local function compareFiles(f1, f2)
   verboseLog(string.format("Comparing %s %s ...\n", f1, f2))
   if not lfs.attributes(f1) or not lfs.attributes(f2) then
      verboseLog(string.format(
                    " ... files %s and/or %s do not exist!\n", f1, f2))
      return false
   end

   local r1, r2 = AdiosReader.Reader(f1), AdiosReader.Reader(f2)

   local cmpPass = true
   local currMaxDiff = 0.0

   if r1:hasVar("CartGridField") and r2:hasVar("CartGridField") then
      -- compare CartField
      local d1, d2 = r1:getVar("CartGridField"):read(), r2:getVar("CartGridField"):read()

      if d1:size() ~= d2:size() then
         verboseLog(string.format(
                       " ... CartGridField in files %s and %s not the same size!\n", f1, f2))
         return false
      end

      local maxVal = maxValueInField(d1) -- maximum value (for numeric comparison)
      for i = 1, d1:size() do
         if check_equal_numeric(d1[i], d2[i], maxVal) == false then
            currMaxDiff = math.max(currMaxDiff, get_relative_numeric(d1[i], d2[i], maxVal))
            cmpPass = false
         end
      end
   elseif r1:hasVar("TimeMesh") and r2:hasVar("TimeMesh") then
      -- Compare DynVector
      local d1, d2 = r1:getVar("Data"):read(), r2:getVar("Data"):read()
      if d1:size() ~= d2:size() then
         verboseLog(string.format(
                       " ... DynVector in files %s and %s not the same size!\n", f1, f2))
         return false
      end

      local maxVal = maxValueInField(d1) -- maximum value (for numeric comparison)
      for i = 1, d1:size() do
         if check_equal_numeric(d1[i], d2[i], maxVal) == false then
            currMaxDiff = math.max(currMaxDiff, get_relative_numeric(d1[i], d2[i], maxVal))
            cmpPass = false
         end
      end
   end

   if cmpPass == false then
      verboseLog(string.format(" ... relative error in file %s is %g ...\n", f2, currMaxDiff))
   end

   r1:close(); r2:close()

   return cmpPass
end


-- Run one simulation until the end to produce reference data.
twoStreamP2(10.0, 2, {1}, false)
-- Copy over the moments and distribution function to reference files. 
fileName    = {GKYL_OUT_PREFIX .. "_elc_2.bp",
               GKYL_OUT_PREFIX .. "_elc_M0_2.bp",
               GKYL_OUT_PREFIX .. "_elc_M1i_2.bp",
               GKYL_OUT_PREFIX .. "_elc_M2_2.bp"}
fileNameRef = {GKYL_OUT_PREFIX .. "_elc_2ref.bp",
               GKYL_OUT_PREFIX .. "_elc_M0_2ref.bp",
               GKYL_OUT_PREFIX .. "_elc_M1i_2ref.bp",
               GKYL_OUT_PREFIX .. "_elc_M2_2ref.bp"}
for i,v in pairs(fileName) do
   os.execute("cp " .. fileName[i] .. " " .. fileNameRef[i])
end

-- Run first half of the simulation.
twoStreamP2(5.0,1,{1},false)
-- Restart and run the second half of the simulation.
twoStreamP2(10.0,2,{1},true)

print(' ')
-- Compare the outputs of this last run to reference data produced above.
compFlag      = true
areFilesEqual = true
for i,v in pairs(fileName) do
   compFlag = compareFiles(fileNameRef[i],fileName[i])
   areFilesEqual = areFilesEqual and compFlag
end

if not areFilesEqual then
   print(" --> ERROR <--- ")
   print(" Files are NOT equal => restarts are BROKEN ")
else
   print(" No errors in restarts detected.")
end
print(" ")

