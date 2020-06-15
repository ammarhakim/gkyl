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
   local a = Array.Array {10, 20}

   assert_equal(2, a.r, "Testing array rank")
   assert_equal(10, a.s[0], "Testing array shape")
   assert_equal(20, a.s[1], "Testing array shape")
   assert_equal(1, a.c, "Use count")
   assert_equal(200, a.N, "Total N")
   assert_equal(true, a.d == a:data(), "Checking if pointer addresses are same")
end

-- Run tests
test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
