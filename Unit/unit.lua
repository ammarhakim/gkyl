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

local function assert_equal(expected, actual, msg)
   if math.abs(actual) < 1e-15 then
      if math.abs(expected-actual) > 1e-15 then
	 print( "** Assert_equal FAILED", msg, string.format("expected %s but was %s", format_arg(expected), format_arg(actual)) )
	 stats.fail = stats.fail+1
      else
	 stats.pass = stats.pass+1
      end
   else
      if math.abs(1-expected/actual) > 1e-15 then
	 print( "** Assert_equal FAILED", msg, string.format("expected %s but was %s", format_arg(expected), format_arg(actual)) )
	 stats.fail = stats.fail+1
      else
	 stats.pass = stats.pass+1
      end
   end
   return actual
end

return {
   assert_equal = assert_equal,
   stats = stats,
}
