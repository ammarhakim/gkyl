-- Gkyl ------------------------------------------------------------------------
--
-- Walks down directory tree and runs all Lua and shell files with
-- prefix "rt-" through gkyl.
-- 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Logger = require "Lib.Logger"
local Time = require "Lib.Time"
local argparse = require "Lib.argparse"
local date = require "xsys.date"
local lfs = require "lfs"

local parser = argparse()
   :name("runregression")
   :require_command(false)
   :description [[
Run Gkyl regression tests. Regression tests are run and results
are compared to "accepted" results. For meaningful testing you
first need to configure the regression testing system and
generate accepted results. This can be done using the "config"
and "create" commands. Obviously, this should be done
only once.

When adding new tests it is that developer's responsibility to
ensure that the results produced are correct. Obviously, others
may not be able to determine just looking at output that the
results make sense.
]]

-- add options, flags and commands to CLI parser
parser:flag("-d --dont-compare", "Don't compare with accepted results")
parser:option("-r --run-only", "Only run specified test")

-- configure system
local conf = parser:command("config", "Configure regression tests")
conf:option("-p --prefix", "Location to store accepted results")
conf:option("-m --mpi-exec", "Full path to MPI executable")

-- create accepted results
parser:command("create", "Create accepted results")

local log = Logger { logToFile = true }

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
	    fullNm = dir .. "/" .. fn
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
   for l in f:lines() do
      -- silent output
   end
   log(string.format("... completed in %g sec\n\n", Time.clock()-tmStart))
   
   
end

-- runs a single shell-script test
local function runShellTest(test)
   log(string.format("Running shell test %s\n", test))
end

local function isLuaRegressionTest(fn)
   return string.find(fn, "^rt%-.-%.lua$") ~= nil 
end
local function isShellRegressionTest(fn)
   return string.find(fn, "^rt%-.-%.sh$") ~= nil
end

-- run all tests
local function runAll()
   log("Starting Gkyl regression tests ...\n")
   log(date(false):fmt() .. "\n\n")
   local tmStart = Time.clock()
   
   for dir, fn, attr in dirtree(".") do
      if isLuaRegressionTest(fn) then
	 runLuaTest(dir .. "/" .. fn)
      elseif isShellRegressionTest(fn) then
	 runShellTest(dir .. "/" .. fn)
      end
   end

   log(string.format("Completed regression tests in %g secs\n", Time.clock()-tmStart))
   log(date(false):fmt() .. "\n")
end

-- parse command-line args
local args = parser:parse(GKYL_COMMANDS)

runAll()
