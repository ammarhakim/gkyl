-- Gkyl ------------------------------------------------------------------------
--
-- Test for range objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local xsys = require "xsys"

local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local assert_equal = Unit.assert_equal
local stats = Unit.stats

ffi.cdef [[
  double calcSum(int n, double *v);
]]

function test_1()
   local v = ffi.new(typeof("double[?]"), 10)
   for i = 1, 10 do
      v[i-1] = i
   end
   local sum = ffi.C.calcSum(10, v)
   assert_equal(55, sum, "Checking if external call to sum worked")
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
