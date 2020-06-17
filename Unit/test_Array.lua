-- Gkyl ------------------------------------------------------------------------
--
-- Test for range objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Array = require "DataStruct.Array"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local arr = Array.Array( {10, 20}, Array.double)

   assert_equal(2, arr.r, "Testing array rank")
   assert_equal(10, arr.d[0], "Testing array shape")
   assert_equal(20, arr.d[1], "Testing array shape")
   assert_equal(1, arr.c, "Use count")
   assert_equal(200, arr.n, "Total n")
   assert_equal(ffi.sizeof("double"), arr.s, "Size of element")
   local arrData = Array.dataPtr(arr, Array.double)
   assert_equal(true, arrData == arr.d+arr.r)

   local brr = arr:clone()
   assert_equal(2, brr.r, "Testing array rank")
   assert_equal(10, brr.d[0], "Testing array shape")
   assert_equal(20, brr.d[1], "Testing array shape")
   assert_equal(1, brr.c, "Use count")
   assert_equal(200, brr.n, "Total n")
   assert_equal(ffi.sizeof("double"), brr.s, "Size of element")
   local brrData = Array.dataPtr(brr, Array.double)
   assert_equal(true, brrData == brr.d+brr.r)
end

-- Run tests
test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
