-- Gkyl ------------------------------------------------------------------------
--
-- Unit testing code
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local math = require "math"
local stats = { fail = 0, pass = 0 }

local function format_arg(arg)
   local argtype = type(arg)
   if argtype == "string" then
      return "'"..arg.."'"
   elseif argtype == "number" or argtype == "boolean" or argtype == "nil" then
      return tostring(arg)
   else
      return "["..tostring(arg).."]"
   end
end

local function assert_equal_numeric(expected, actual, msg)
   if math.abs(actual) < 1e-15 then
      if math.abs(expected-actual) > 5e-15 then
	 print( "** Assert_equal FAILED", msg, string.format("expected %s but was %s", format_arg(expected), format_arg(actual)) )
	 stats.fail = stats.fail+1
      else
	 stats.pass = stats.pass+1
      end
   else
      if math.abs(1-expected/actual) > 5e-14 then
	 print( "** Assert_equal FAILED", msg, string.format("expected %s but was %s", format_arg(expected), format_arg(actual)) )
	 stats.fail = stats.fail+1
      else
	 stats.pass = stats.pass+1
      end
   end
   return actual
end

local function assert_equal_string(expected, actual, msg)
   if expected ~= actual then
      print( "** Assert_equal FAILED", msg, string.format("expected %s but was %s", format_arg(expected), format_arg(actual)) )
      stats.fail = stats.fail+1
   else
      stats.pass = stats.pass+1
   end
   return actual
end

local function assert_equal_bool(expected, actual, msg)
   if expected ~= actual then
      print( "** Assert_equal FAILED", msg, string.format("expected %s but was %s", format_arg(expected), format_arg(actual)) )
      stats.fail = stats.fail+1
   else
      stats.pass = stats.pass+1
   end
   return actual
end

local function assert_equal(expected, actual, msg)
   if type(expected) == "number" then
      assert_equal_numeric(expected, actual, msg)
   elseif  type(expected) == "string" then
      assert_equal_string(expected, actual, msg)
   elseif  type(expected) == "boolean" then
      assert_equal_bool(expected, actual, msg)
   else
      io.write(string.format("Not performing test %s\n", msg))
   end
      
end

local function assert_close(expected, actual, tol, msg)
   if type(expected) ~= "number" then
      io.write(string.format("assert_close only works for numbers. Not performing test %s\n", msg))
   end

   if math.abs(actual) < 1e-15 then
      if math.abs(expected-actual) > tol then
	 print( "** Assert_equal FAILED", msg, string.format("expected %s but was %s; diff is greater than tol %s", format_arg(expected), format_arg(actual), format_arg(tol)) )
	 stats.fail = stats.fail+1
      else
	 stats.pass = stats.pass+1
      end
   else
      if math.abs(1-expected/actual) > tol/actual then
	 print( "** Assert_equal FAILED", msg, string.format("expected %s but was %s; diff is greater than tol %s", format_arg(expected), format_arg(actual), format_arg(tol)) )
	 stats.fail = stats.fail+1
      else
	 stats.pass = stats.pass+1
      end
   end
   return actual
   
end

return {
   assert_equal = assert_equal,
   assert_close = assert_close,
   stats = stats,
}
