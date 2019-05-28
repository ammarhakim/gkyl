-- Gkyl ------------------------------------------------------------------------
--
-- Runs regression tests and compares with accepted results.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosReader = require "Io.AdiosReader"
local Logger = require "Lib.Logger"
local Time = require "Lib.Time"
local argparse = require "Lib.argparse"
local date = require "xsys.date"
local lfs = require "lfs"
local lume = require "Lib.lume"
local sql = require "sqlite3"

local log = Logger { logToFile = true }
local verboseLog = function (msg) end -- default no messages are written
local verboseLogger = function (msg) log(msg) end

-- number of passed/failed and total tests
local numFailedTests = 0
local numPassedTests = 0
local numTotalTests = 0

-- flag to indicate if we are configuring system
local isConfiguring = false
-- global value to store config information
local configVals = nil
-- global list of tests to ignore
local ignoreTests = {}
-- list of MOAT regression
local moatTests = {}
-- Sqlite3 database for storing regression data
local sqlConn = nil
-- stored procedures to store data in tables
local insertRegressionDataProc = nil
local insertRegressionMetaProc = nil

-- UUID generator
math.randomseed( os.time() )
local function uuid()
   local random = math.random
   local template ='xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'
   return string.gsub(template, '[xy]', function (c)
			 local v = (c == 'x') and random(0, 0xf) or random(8, 0xb)
			 return string.format('%x', v)
   end)
end

-- generate ID for this run of regression system
local runID = uuid()
-- date when tests were run
local runDate = date(false):fmt("${iso}")

-- loads configuration file
local function loadConfigure(args)
   local f = loadfile("runregression.config.lua")
   if not f then
      log("Regression tests not configured! Run config command first\n")
      os.exit(1)
   end
   configVals = f()

   if args.verbose then
      verboseLog = verboseLogger -- set verbose logger
   end

   local g = loadfile("ignoretests.lua")
   if not args.all then
      ignoreTests = g()
   end

   local m = loadfile("moat.lua")
   moatTests = m()

   -- open connection to SQL DB
   sqlConn = sql.open(string.format("%s/regressiondb", configVals['results_dir']))
   -- SQL stored procedures to insert data into table
   insertRegressionDataProc = sqlConn:prepare [[
     insert into RegressionData values (
       ?, ?, ?, ?, ?
     )
   ]]
   insertRegressionMetaProc = sqlConn:prepare [[
     insert into RegressionMeta values (
       ?, ?, ?, ?, ?, ?, ?, ?
     )
   ]]
   
   return configVals
end

-- functions to insert data into the tables
local function insertRegressionData(id, nm, status, runtm, runlog)
   insertRegressionDataProc:reset():bind(
      id,
      nm,
      status,
      runtm,
      runlog):step()
end
local function insertRegressionMeta(id, tm, ntotal, npass, nfail)
   insertRegressionMetaProc:reset():bind(
      id,
      tm,
      GKYL_EXEC,
      GKYL_HG_CHANGESET,
      GKYL_BUILD_DATE,      
      ntotal,
      npass,
      nfail):step()
end

-- configure paths
local function configure(prefix, mpiExec)
   local mpiAttr = lfs.attributes(mpiExec)

   -- FOR NOW THIS IS DISABLED TILL I FIGURE OUT HOW TO DO PARALLEL
   -- TESTS (AH, 5/24/2019)
   if false then
      if mpiAttr and mpiAttr.mode == "file" and string.sub(mpiAttr.permissions, 3, 3) == "x" then
	 log(string.format("Setting MPIEXEC to %s ...\n", mpiExec))
      else
	 assert(false, string.format("MPI launch binary %s not found or is not an executable!", mpiExec))
      end
   end

   local prefixAttr = lfs.attributes(prefix)
   if prefixAttr and prefixAttr.mode == "directory" then
      log(string.format("Will write accepted results to %s/gkyl-results ...\n", prefix))
   else
      assert(false, string.format("Prefix %s is not a directory!", prefix))
   end

   -- create Sqlite3 database to store regression data
   log(string.format("Creating DB file %s/gkyl-results/regressiondb\n", prefix))
   local conn = sql.open(string.format("%s/gkyl-results/regressiondb", prefix))

   -- Following tables are created: RegressionMeta that stores date
   -- when test was run and how many passed or failed. RegressionData
   -- stores output from each test. GUIDs are shared between tables so
   -- one can get all tests and their data given the GUID
   conn:exec [[

    drop table if exists RegressionMeta;
    create table RegressionMeta (
      guid text,
      tstamp text,
      GKYL_EXEC text,
      GKYL_HG_CHANGESET text,
      GKYL_BUILD_DATE text,
      ntotal integer,
      npass integer,
      nfail integer
    );

    drop table if exists RegressionData;
    create table RegressionData (
      guid text,
      name text,
      status integer,
      runtime real,
      runlog text
    );
   ]]

   -- write information into config file
   local fn = io.open("runregression.config.lua", "w")
   fn:write("return {\n")
   fn:write(string.format("mpiExec = \"%s\",\n", mpiExec))
   fn:write(string.format("results_dir = \"%s/gkyl-results\",\n", prefix))
   fn:write("}")
   fn:close()
end

-- Walks down a directory tree recursively: returns a coroutine that
-- yields directory, file-name and attribute object
local function dirtree(dir)
   assert(dir and dir ~= "", "directory parameter is missing or empty")
   if string.sub(dir, -1) == "/" then
      dir = string.sub(dir, 1, -2)
   end

   local function yieldtree(dir)
      for fn in lfs.dir(dir) do
	 if fn ~= "." and fn ~= ".." then
	    local fullNm = dir .. "/" .. fn
	    local attr = lfs.attributes(fullNm)
	    coroutine.yield(dir, fn, attr)
	    if attr.mode == "directory" then
	       yieldtree(fullNm)
	    end
	 end
      end
   end
   return coroutine.wrap(function() yieldtree(dir) end)
end

-- creates directories if they don't exist
local function mkdir(path)
   local sep, pStr = package.config:sub(1,1), "/"
   for dir in path:gmatch("[^" .. sep .. "]+") do
      pStr = pStr .. dir .. sep
      lfs.mkdir(pStr)
   end
end

-- copies all files with specified extension to target directory
local function copyAllFiles(srcPath, ext, destPath)
   -- for now using os.execute to run the "cp" command. Perhaps this
   -- is not the best way to do this and one might want to use a more
   -- Lua-ish way
   local cp = "cp " .. srcPath .. "_*." .. ext .. " " .. destPath
   verboseLog(string.format("... executing %s ...\n", cp))
   os.execute(cp)
end

-- runs a single lua test
local function runLuaTest(test)
   log(string.format("\nRunning test %s ...\n", test))

   numTotalTests = numTotalTests+1

   local opts = {}
   -- open Lua file and check if first line is a command line
   local f = io.open(test, "r")
   local line1 = f:read()
   if string.find(line1, "--!") then
      local r = string.sub(line1, 4, -1)
      opts = loadstring("return " .. r)()
   end
   f:close()

   -- before running tests delete all old output
   local rm = "rm -f " .. string.sub(test, 1, -5) .. "_*.bp"
   os.execute(rm)

   -- construct command to run
   local runCmd = string.format("%s %s", GKYL_EXEC, test)
   if opts.numProc then
      runCmd = string.format("%s -n %d %s %s", configVals.mpiExec, opts.numProc, GKYL_EXEC, test)
   end

   -- run test
   local tmStart = Time.clock()

   local runlog = ""
   if opts.numProc then
      -- skip for now: NEED TO USE MPI_COMM_SPAWN FOR THIS, I THINK!
      log(string.format("**** NOT RUNNING PARALLEL TEST %s \n", runCmd))
   else
      local f = io.popen(runCmd, "r")
      for l in f:lines() do
	 runlog = runlog .. l .. "\n"
	 verboseLog(l.."\n")
      end
   end
   local runtm = Time.clock()-tmStart
   log(string.format("... completed in %g sec\n", runtm))
   return runtm, runlog
end

-- runs a single shell-script test
local function runShellTest(test)
   log(string.format("Running shell test %s\n", test))
end

-- runs a single lua unit test
local function runLuaUnitTest(test)
   log(string.format("\nRunning test %s ...\n", test))

   local runCmd = string.format("%s %s", GKYL_EXEC, test)
   local f = io.popen(runCmd, "r")
   local outPut = f:read("*a")
   if string.find(outPut, "FAILED") then
      log(string.format("... %s FAILED!\n", test))
      numFailedTests = numFailedTests+1
   else
      numPassedTests = numPassedTests+1
      log("... passed.\n")
   end
end

local function isLuaRegressionTest(fn)
   return string.find(fn, "/rt%-[^/]-%.lua$") ~= nil
end
local function isShellRegressionTest(fn)
   return string.find(fn, "/rt%-[^/]-%.sh$") ~= nil
end

local function isLuaUnitTest(fn)
   return string.find(fn, "/test_[^/]-%.lua$") ~= nil
end
local function isShellUnitTest(fn)
   return string.find(fn, "/test_[^/]-%.sh$") ~= nil
end

-- list of regression tests
local function list_tests(args)
   local luaRegTests, shellRegTests = {}, {}

   local function addTest(fn)
      if isLuaRegressionTest(fn) then
	 table.insert(luaRegTests, fn)
      elseif isShellRegressionTest(fn) then
	 table.insert(shellRegTests, fn)
      end
   end

   if args.moat then
      -- only MOAT regressions are to be run
      for _, t in ipairs(moatTests) do addTest(t) end
   elseif args.run_only then
      local a = lfs.attributes(args.run_only)
      if a.mode == "file" then
	 addTest(args.run_only)
      elseif a.mode == "directory" then
	 for dir, fn, _ in dirtree(args.run_only) do addTest(dir .. "/" .. fn) end
      end
   else
      for dir, fn, _ in dirtree(".") do addTest(dir .. "/" .. fn) end
   end

   -- this function is used to filter out tests that should not be
   -- run. Note the confusing name
   local function shouldRun(t)
      if lume.find(ignoreTests, t) then return false end
      return true
   end
   
   return lume.filter(luaRegTests, shouldRun), shellRegTests
end

-- list of unit tests
local function list_unit_tests(args)
   local luaUnitTests, shellUnitTests = {}, {}

   local function addTest(fn)
      if isLuaUnitTest(fn) then
	 table.insert(luaUnitTests, fn)
      elseif isShellUnitTest(fn) then
	 table.insert(shellUnitTests, fn)
      end
   end

   for dir, fn, _ in dirtree("../Unit") do addTest(dir .. "/" .. fn) end

   return luaUnitTests, shellUnitTests
end

-- function to handle "config" command
local function config_action(args, name)
   isConfiguring = true
   
   local prefix = args.config_prefix and args.config_prefix or
      os.getenv("HOME") .. "/gkylsoft"
   local mpiexec = args.config_mpiexec and  args.config_mpiexec or
      os.getenv("HOME") .. "/gkylsoft/openmpi/bin/mpiexec"

   configure(prefix, mpiexec)
end

-- function to handle "list" command
local function list_action(args, name)
   loadConfigure(args)
   
   local luaRegTests, shellRegTests = list_tests(args)
   lume.each(luaRegTests, print)
   lume.each(shellRegTests, print)
end

-- function to handle "create" sub-command of "run"
local function create_action(test)
   log(string.format("... creating accepted results for %s ...\n", test))
   local fullResultsDir = configVals.results_dir .. "/"
      .. string.sub(test, 1, -5) -- remove the initial ./ and last .lua

   verboseLog(string.format("... making directory %s ..\n", fullResultsDir))
   mkdir(fullResultsDir) -- recursively create all dirs as needed
   local srcPath = string.sub(test, 1, -5)
   copyAllFiles(srcPath, "bp", fullResultsDir)

   return -2
end

-- function to compare floats (NOT SURE IF THIS IS BEST WAY TO DO
-- THINGS)
local function check_equal_numeric(expected, actual)
   if math.abs(actual) < 1e-15 then
      if math.abs(expected-actual) > 1e-12 then
	 return false
      end
   else
      if math.abs(1-expected/actual) > 1e-12 then
	 return false
      end
   end
   return true
end

-- relative difference between two numbers (NOT SURE IF THIS IS BEST
-- WAY TO DO THINGS)
local function get_relative_numeric(expected, actual)
   if math.abs(actual) < 1e-15 then
      return math.abs(expected-actual)
   else
      return math.abs(1-expected/actual)
   end
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
      for i = 1, d1:size() do
	 if check_equal_numeric(d1[i], d2[i]) == false then
	    currMaxDiff = math.max(currMaxDiff, get_relative_numeric(d1[i], d2[i]))
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
      for i = 1, d1:size() do
	 if check_equal_numeric(d1[i], d2[i]) == false then
	    currMaxDiff = math.max(currMaxDiff, get_relative_numeric(d1[i], d2[i]))
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

-- function to handle "check" sub-command of "run"
local function check_action(test)
   local fullResultsDir = configVals.results_dir .. "/"
      .. string.sub(test, 1, -5) -- remove the initial ./ and last .lua

   local vloc = string.find(test, "/rt%-[^/]-%.lua$")
   local outDirName = string.sub(test, 1, vloc)
   -- the very strange looking pattern below puts an escape character
   -- (i.e. \) before characters that are treated as special by
   -- pattern matcher. (The final underscore is needed to avoid the
   -- regression system getting confused with other tests that have
   -- the same prefix)
   local testPrefix = string.gsub(
      string.sub(test, 3, -5), "(%W)", "%%%1") .. "_"

   local passed = true
   for fn in lfs.dir(outDirName) do
      local fullNm = outDirName .. fn
      local attr = lfs.attributes(fullNm)
      if attr.mode == "file" then
	 if string.find(fullNm, testPrefix) and string.sub(fullNm, -3, -1) == ".bp" then
	    local acceptedFileNm = fullResultsDir .. "/" .. string.sub(fullNm, vloc+1, -1)
	    local status = compareFiles(acceptedFileNm, fullNm)
	    passed = passed and status
	 end
      end
   end
   if passed then
      numPassedTests = numPassedTests+1
      log("... passed.\n")
   else
      numFailedTests = numFailedTests+1
      log(string.format("... %s FAILED!\n", test))
   end
   return passed and 1 or 0
end

-- function to handle "run" command
local function run_action(args, name)
   loadConfigure(args)
   local luaRegTests, shellRegTests = list_tests(args)

   -- function to run after simulation is finisihed
   local postRun = function(f) return -1 end
   if args.create then
      postRun = create_action
   elseif args.check then
      postRun = check_action
   end

   log("Running regression tests ...\n\n")
   local tmStart = Time.clock()

   local function trimname(nm)
      local s,e = string.find(nm, "./")
      if s and s==1 then return string.sub(nm, e+1) end
      return nm
   end

   -- run all tests
   for _, test in pairs(luaRegTests) do
      local runtm, runlog = runLuaTest(test)
      local status = postRun(test)
      -- insert data into SQL table
      insertRegressionData(runID, trimname(test), status, runtm, runlog)
   end

   log(string.format("All regression tests completed in %g secs\n", Time.clock()-tmStart))
end

-- function to handle "listunit" command
local function listunit_action(args, name)
   loadConfigure(args)
   
   local luaUnitTests, shellUnitTests = list_unit_tests(args)
   lume.each(luaUnitTests, print)
   lume.each(shellUnitTests, print)
end

-- function to handle "rununit" command
local function rununit_action(args, name)
   loadConfigure(args)
   local luaUnitTests, shellUnitTests = list_unit_tests(args)

   log("Running unit tests ...\n\n")
   local tmStart = Time.clock()
   lume.each(
      luaUnitTests, function (f) runLuaUnitTest(f) end)

   log(string.format("All unit tests completed in %g secs\n", Time.clock()-tmStart))
end

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("runregression")
   :require_command(true)
   :description [[
Run Gkyl regression tests. Regression tests are run and results are
compared to "accepted" results. For meaningful testing you first need
to configure the regression testing system and generate accepted
results. This can be done using the "config" and "run create"
commands. Obviously, this should be done only once.

Once accepted results are created or obtained from elsewhere,
regression testing can be  done using "run check" command.

When adding new tests it is that developer's responsibility to
ensure that the results produced are correct. Obviously, others
may not be able to determine just looking at output that the
results make sense.
]]

parser:flag("-v --verbose", "Print verbose messages as tests are run")
parser:flag("-a --all", "Run all tests, even those in ignoretests file")

-- "config" command
local c_conf = parser:command("config", "Configure regression tests")
   :action(config_action)

c_conf:option("-p --prefix", "Location to store accepted results")
   :target("config_prefix")
c_conf:option("-m --mpiexec", "Full path to MPI executable")
   :target("config_mpiexec")

-- "list" command
local c_list = parser:command("list", "List all regression tests")
   :action(list_action)
c_list:option("-r --run-only", "Only list this test or all tests in this directory")
c_list:flag("-m --moat", "Only list MOAT regression")

-- "run" tests
local c_run = parser:command("run")
   :description("Run regression tests.")
   :require_command(false)
   :action(run_action)
c_run:option("-r --run-only", "Only run this test or all tests in this directory")
c_run:flag("-m --moat", "Only run key MOAT regression")

-- check against accepted results
c_run:command("check", "Check results")

-- create accepted results
c_run:command("create", "Create accepted results")

-- "listunit" command
parser:command("listunit", "List all unit tests")
   :action(listunit_action)

-- "rununit" tests
local c_rununit = parser:command("rununit")
   :description("Run unit tests.")
   :require_command(false)
   :action(rununit_action)

-- parse command-line args (functionality encoded in command actions)
local _ = parser:parse(GKYL_COMMANDS)

if numPassedTests > 0 then
   log(string.format("\nTotal tests passed: %d\n", numPassedTests))
end
if numFailedTests > 0 then
   log(string.format("Total tests FAILED: %d\n", numFailedTests))
end

if numTotalTests>0 then
   insertRegressionMeta(runID, runDate, numTotalTests, numPassedTests, numFailedTests)
end
