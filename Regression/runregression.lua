-- Gkyl ------------------------------------------------------------------------
--
-- Runs regression tests and compares with accepted results.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

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

   -- write information in config file
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

-- runs a single lua test
local function runLuaTest(test)
   log(string.format("Running test %s ...\n", test))

   local tmStart = Time.clock()

   local runCmd = string.format("%s %s", GKYL_EXEC, test)
   local f = io.popen(runCmd, "r")
   for _ in f:lines() do
   end
   log(string.format("... completed in %g sec\n\n", Time.clock()-tmStart))
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

local function show(args)
   for v in pairs(args) do
      print(v, args[v])
   end
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
   -- check if top-level directory was created
   local attrs = lfs.attributes(configVals.results_dir)
   if not attrs then
      assert(lfs.mkdir(configVals.results_dir), string.format(
		"Unable to create results directory %s!", configVals.results_dir))
   end
end

-- function to handle "check" sub-command of "run"
local function check_action(test)
   print(string.format("Inside check for %s", test))
end

-- function to handle "run" command
local function run_action(args, name)
   loadConfigure()
   local luaRegTests, shellRegTests = list_tests(args)

   -- default post-run is to do nothing
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
