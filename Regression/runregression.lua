-- Gkyl ------------------------------------------------------------------------
--
-- Walks down directory tree and runs all Lua and shell files with
-- prefix "rt-" through gkyl.
-- 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local lfs = require "lfs"
local Time = require "Lib.Time"
local date = require "xsys.date"
local Logger = require "Lib.Logger"

local log = Logger {
   logToFile = true
}

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

-- check if a file is a Lua-based regression test
local function isLuaRegressionTest(fn)
   if string.len(fn) > 3 then
      if string.sub(fn, 1, 3) == "rt-" then
	 if string.sub(fn, string.len(fn)-3, string.len(fn)) == ".lua" then
	    return true
	 end
      end
   end
   return false
end
-- check if a file is a shell-based regression test
local function isShellRegressionTest(fn)
   if string.len(fn) > 3 then
      if string.sub(fn, 1, 3) == "rt-" then
	 if string.sub(fn, string.len(fn)-2, string.len(fn)) == ".sh" then
	    return true
	 end
      end
   end
   return false
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

runAll()
