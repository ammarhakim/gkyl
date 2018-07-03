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

local log = Logger { logToFile = true }

-- global value to store config information
local configVals = nil

-- loads configuration file
local function loadConfigure()
   local f = loadfile("runregression.config.lua")
   if not f then
      log("Regression tests not configured! Run config command first\n")
      os.exit(1)
   end
   configVals = f()
   return configVals
end

-- configure paths
local function configure(prefix, mpiExec)
   local mpiAttr = lfs.attributes(mpiExec)

   if mpiAttr and mpiAttr.mode == "file" and string.sub(mpiAttr.permissions, 3, 3) == "x" then
      log(string.format("Setting MPIEXEC to %s ...\n", mpiExec))
   else
      assert(false, string.format("MPI launch binary %s not found or is not an executable!", mpiExec))
   end

   local prefixAttr = lfs.attributes(prefix)
   if prefixAttr and prefixAttr.mode == "directory" then
      log(string.format("Will write accepted results to %s/gkyl-results ...\n", prefix))
   else
      assert(false, string.format("Prefix %s is not a directory!", prefix))
   end

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
   os.execute(cp)
end

-- runs a single lua test
local function runLuaTest(test)
   log(string.format("\nRunning test %s ...\n", test))

   local tmStart = Time.clock()
   local runCmd = string.format("%s %s", GKYL_EXEC, test)
   local f = io.popen(runCmd, "r")
   for _ in f:lines() do end -- silent output
   log(string.format("... completed in %g sec\n", Time.clock()-tmStart))
end

-- runs a single shell-script test
local function runShellTest(test)
   log(string.format("Running shell test %s\n", test))
end

local function isLuaRegressionTest(fn)
   return string.find(fn, "/rt%-[^/]-%.lua$") ~= nil 
end
local function isShellRegressionTest(fn)
   return string.find(fn, "/rt%-[^/]-%.sh$") ~= nil
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

   if args.run_only then
      addTest(args.run_only)
   else
      for dir, fn, _ in dirtree(".") do addTest(dir .. "/" .. fn) end
   end
   return luaRegTests, shellRegTests
end

-- function to handle "config" command
local function config_action(args, name)
   local prefix = args.config_prefix and args.config_prefix or
      os.getenv("HOME") .. "/gkylsoft"
   local mpiexec = args.config_mpiexec and  args.config_mpiexec or
      os.getenv("HOME") .. "/gkylsoft/openmpi/bin/mpiexec"

   configure(prefix, mpiexec)
end

-- function to handle "list" command
local function list_action(args, name)
   local luaRegTests, shellRegTests = list_tests(args)
   lume.each(luaRegTests, print)
   lume.each(shellRegTests, print)
end

-- function to handle "create" sub-command of "run"
local function create_action(test)
   log(string.format("... creating accepted results for %s ...\n", test))
   local fullResultsDir = configVals.results_dir .. "/"
      .. string.sub(test, 3, -5)
   mkdir(fullResultsDir) -- recursively create all dirs as needed
   local srcPath = string.sub(test, 1, -5)
   copyAllFiles(srcPath, "bp", fullResultsDir)
end

-- function to compare floats (NOT SURE IF THIS IS BEST WAY TO DO
-- THINGS)
local function check_equal_numeric(expected, actual)
   if math.abs(actual) < 1e-15 then
      if math.abs(expected-actual) > 5e-15 then
	 return false
      end
   else
      if math.abs(1-expected/actual) > 5e-14 then
	 return false
      end
   end
   return true
end

-- function to compare files
local function compareFiles(f1, f2)
   if not lfs.attributes(f1) or not lfs.attributes(f2) then
      return false
   end
   
   local r1, r2 = AdiosReader.Reader(f1), AdiosReader.Reader(f2)

   if r1:hasVar("CartGridField") and r2:hasVar("CartGridField") then
      -- compare CartField
      local d1, d2 = r1:getVar("CartGridField"):read(), r2:getVar("CartGridField"):read()

      if d1:size() ~= d2:size() then return false end
      for i = 1, d1:size() do
	 if check_equal_numeric(d1[i], d2[i]) == false then
	    return false
	 end
      end
   elseif r1:hasVar("TimeMesh") and r2:hasVar("TimeMesh") then
      -- Compare DynVector
      local d1, d2 = r1:getVar("Data"):read(), r2:getVar("Data"):read()
      if d1:size() ~= d2:size() then return false end
      for i = 1, d1:size() do
	 if check_equal_numeric(d1[i], d2[i]) == false then
	    return false
	 end
      end
   else
      -- skip file 
   end

   r1:close(); r2:close()
   
   return true
end

-- function to handle "check" sub-command of "run"
local function check_action(test)
   local fullResultsDir = configVals.results_dir .. "/"
      .. string.sub(test, 3, -5) -- remove the initial ./ and last .lua

   local vloc = string.find(test, "/rt%-[^/]-%.lua$")
   local outDirName = string.sub(test, 1, vloc)
   local testPrefix = string.gsub(
      string.sub(test, 3, -5), "(%W)", "%%%1")

   local passed = true
   for fn in lfs.dir(outDirName) do
      local fullNm = outDirName .. fn
      local attr = lfs.attributes(fullNm)
      if attr.mode == "file" then
	 if string.find(fullNm, testPrefix) and string.sub(fullNm, -3, -1) == ".bp" then
	    local acceptedFileNm = fullResultsDir .. "/" .. string.sub(fullNm, vloc+1, -1)
	    passed = passed and compareFiles(acceptedFileNm, fullNm)
	 end
      end
   end
   if passed then
      log("... passed.\n")
   else
      log("... FAILED!\n")
   end
end

-- function to handle "run" command
local function run_action(args, name)
   loadConfigure()
   local luaRegTests, shellRegTests = list_tests(args)

   -- function to run after simulation is finisihed
   local postRun = function(f) end
   if args.create then
      postRun = create_action
   elseif args.check then
      postRun = check_action
   end

   log("Running regression tests ...\n\n")
   local tmStart = Time.clock()
   
   lume.each(
      luaRegTests, function (f) runLuaTest(f); postRun(f) end)
   lume.each(
      shellRegTests, function (f) runShellTest(f); postRun(f) end)
   
   log(string.format("All tests completed in %g secs\n", Time.clock()-tmStart))
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

-- "run" tests
local c_run = parser:command("run")
   :description("Run regression tests.")
   :require_command(false)
   :action(run_action)

c_run:option("-r --run-only", "Only run this test")

-- check against accepted results
c_run:command("check", "Check results")

-- create accepted results
c_run:command("create", "Create accepted results")

-- parse command-line args (functionality encoded in command actions)
local args = parser:parse(GKYL_COMMANDS)
