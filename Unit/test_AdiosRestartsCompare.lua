-- Gkyl ------------------------------------------------------------------------
--
-- Compare files for test of ADIOS-based restarts.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

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

-- (from Tool/runregression.lua) Function to compare files (from Tool/runregression.lua).
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
   else
      verboseLog(string.format("     Identical CartGridField field data.\n"))
   end

   r1:close(); r2:close()

   return cmpPass
end

-- Names of latest and reference data files.
-- Note: simulationName is intended to be passed in script/command line.
fileName    = {simulationName .. "_elc_" .. tostring(compareFrame) .. ".bp",
               simulationName .. "_elc_M0_" .. tostring(compareFrame) .. ".bp",
               simulationName .. "_elc_M1i_" .. tostring(compareFrame) .. ".bp",
               simulationName .. "_elc_M2_" .. tostring(compareFrame) .. ".bp",
               simulationName .. "_field_" .. tostring(compareFrame) .. ".bp"}
fileNameRef = {simulationName .. "_elc_" .. tostring(compareFrame) .. "ref.bp",
               simulationName .. "_elc_M0_" .. tostring(compareFrame) .. "ref.bp",
               simulationName .. "_elc_M1i_" .. tostring(compareFrame) .. "ref.bp",
               simulationName .. "_elc_M2_" .. tostring(compareFrame) .. "ref.bp",
               simulationName .. "_field_" .. tostring(compareFrame) .. "ref.bp"}

print(' ')
-- Compare the outputs of this last run to reference data.
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

