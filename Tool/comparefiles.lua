-- Gkyl ------------------------------------------------------------------------
--
-- Compares two files and write out max/min differences
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local lfs = require "lfs"
local AdiosReader = require "Io.AdiosReader"
local Logger = require "Lib.Logger"
local argparse = require "Lib.argparse"

local log = Logger { logToFile = true }
local verboseLog = function (msg) end -- default no messages are written
local verboseLogger = function (msg) log(msg) end

-- function to compare floats: the comparison is normalized to the
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

-- relative difference between two numbers (NOT SURE IF THIS IS BEST
-- WAY TO DO THINGS)
local function get_relative_numeric(expected, actual, maxVal)
   if maxVal < 1e-15 then
      return math.max(expected-actual)
   else
      return math.abs(expected-actual)/maxVal
   end
end

-- calculates maximum value in supplied field
local function maxValueInField(fld)
   local maxVal = 0.0
   for i = 1, fld:size() do
      maxVal = math.max(maxVal, math.abs(fld[i]))
   end
   return maxVal
end

-- function to compare files
local function compareFiles(f1, f2)
   -- verboseLog(string.format("Comparing %s %s ...\n", f1, f2))   
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

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("help")
   :description [[Compare BP files 'a' and 'b']]
parser:option("-a --filea", "File 'a' to compare")
parser:option("-b --fileb", "File 'b' to compare")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS)

if (not args.filea) or (not args.fileb) then
   print("Must specify files to compare with -a and -b!")
else
   filea, fileb = args.filea, args.fileb
   local result = compareFiles(filea, fileb)
   if result then
      print(string.format("Files are same to numeric precision.", filea, fileb))
   else
      print(string.format("Files are different!", filea, fileb))
   end
end
