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

-- need to change this as it is keyed on input file name
GKYL_OUT_PREFIX = lfs.currentdir() .. "/" .. "comparefiles"

local log = Logger { logToFile = true }
local verboseLog = function (msg) end -- default no messages are written
local verboseLogger = function (msg) log(msg) end

-- function to compare floats: the comparison is normalized to the
-- maximum value of the field being compared. Perhaps this is too
-- "coarse" but a direct comparison of floats is very tricky.
local function check_equal_numeric(expected, actual, maxVal)
   if maxVal < GKYL_MIN_DOUBLE then
      return math.abs(expected-actual) > 10*GKYL_MIN_DOUBLE
   end
   if math.abs(expected-actual)/maxVal > 1e-12 then
      return false
   end
   return true
end

-- relative difference between two numbers (NOT SURE IF THIS IS BEST
-- WAY TO DO THINGS)
local function get_relative_numeric(expected, actual, maxVal)
   if maxVal < GKYL_MIN_DOUBLE then
      return math.abs(expected-actual)
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
   --verboseLog(string.format("Comparing %s %s ...\n", f1, f2))   
   if not lfs.attributes(f1) or not lfs.attributes(f2) then
      verboseLog(string.format(
		    " ... files %s and/or %s do not exist!\n", f1, f2))
      return false
   end
   
   local r1, r2 = AdiosReader.Reader(f1), AdiosReader.Reader(f2)

   -- check attribute
   local function checkVecAttr(attrNm)
      print(string.format("Checking attr %s in %s %s ...", attrNm, f1, f2))
      if not r1:hasAttr(attrNm) and not r2:hasAttr(attrNm) then
	 return true -- if both files have attribute missing, consider as pass
      end

      -- If both don't have it
      if not r1:hasAttr(attrNm) or not r2:hasAttr(attrNm) then
	 verboseLog(string.format(
		       " ... CartGridField attr %s not present in both files %s and %s!\n", attrNm, f1, f2))
	 return false
      end

      print(string.format("... comparing %s", attrNm))
      local r1_attrNm, r2_attrNm = r1:getAttr(attrNm):read(), r2:getAttr(attrNm):read()
      if #r1_attrNm ~= #r2_attrNm then
	 verboseLog(string.format(
		       " ... CartGridField attr %s in files %s and %s not the same size!\n", attrNm, f1, f2))
	 return false
      end
      for i = 1, #r1_attrNm do
	 if r1_attrNm[i] ~= r2_attrNm[i] then
	    verboseLog(string.format(
			  " ... CartGridField attr %s not the same files %s and %s not the same!\n", attrNm, f1, f2))
	    return false
	 end
      end
      return true
   end

   local cmpPass = true
   local currMaxDiff = 0.0
   
   if r1:hasVar("CartGridField") and r2:hasVar("CartGridField") then
      
      -- compare stable attributes (not all attributes can be compared)
      if not checkVecAttr("numCells") then return false end
      if not checkVecAttr("lowerBounds") then return false end
      if not checkVecAttr("upperBounds") then return false end
      if not checkVecAttr("basisType") then return false end
      if not checkVecAttr("polyOrder") then return false end
      
      -- compare CartField data
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
   -- Compare DynVector
   elseif (r1:hasVar("TimeMesh0") or r1:hasVar("TimeMesh")) and (r2:hasVar("TimeMesh0") or r2:hasVar("TimeMesh")) then
      local t1, t2
      local d1, d2
      
      -- Read from file 1
      -- If single dataset, read it into single array
      if r1:hasVar("TimeMesh") then 
         t1 = r1:getVar("TimeMesh"):read()
         d1 = r1:getVar("Data"):read()
      -- If multiple datasets, read and append into single array
      elseif r1:hasVar("TimeMesh0") then
         t1 = r1:getVar("TimeMesh0"):read()
         d1 = r1:getVar("Data0"):read()
         frNum = 1
         while r1:hasVar("TimeMesh"..frNum) do
            local t1N = r1:getVar("TimeMesh"..frNum):read()
            local d1N = r1:getVar("Data"..frNum):read()
            for i = 1, t1N:size() do
               t1:push(t1N[i])
            end
            for i = 1, d1N:size() do
               d1:push(d1N[i])
            end
            frNum = frNum + 1
         end
      end

      -- Read from file 2
      -- If single dataset, read it into single array
      if r2:hasVar("TimeMesh") then 
         t2 = r2:getVar("TimeMesh"):read()
         d2 = r2:getVar("Data"):read()
      -- If multiple datasets, read and append into single array
      elseif r2:hasVar("TimeMesh0") then
         t2 = r2:getVar("TimeMesh0"):read()
         d2 = r2:getVar("Data0"):read()
         frNum = 1
         while r2:hasVar("TimeMesh"..frNum) do
            local t2N = r2:getVar("TimeMesh"..frNum):read()
            local d2N = r2:getVar("Data"..frNum):read()
            for i = 1, t2N:size() do
               t2:push(t2N[i])
            end
            for i = 1, d2N:size() do
               d2:push(d2N[i])
            end
            frNum = frNum + 1
         end
      end

      -- Check equivalence
      if d1:size() ~= d2:size() then
         verboseLog(string.format(
          	  " ... DynVector in files %s and %s not the same size!\n", f1, f2))
         return false
      end
      
      local maxVal = math.max(maxValueInField(d1),maxValueInField(d2)) -- maximum value (for numeric comparison)
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
