-- Gkyl ------------------------------------------------------------------------
--
-- Test for range objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Range = require "Lib.Range"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local range = Range.Range({0, 0}, {1, 5})

   assert_equal(2, range:ndim(), "Checking dimension")

   assert_equal(0, range:lower(1), "Checking lower")
   assert_equal(0, range:lower(2), "Checking lower")

   assert_equal(1, range:upper(1), "Checking upper")
   assert_equal(5, range:upper(2), "Checking upper")

   assert_equal(2, range:shape(1), "Checking shape")
   assert_equal(6, range:shape(2), "Checking shape")
end

function test_2()
   local range = Range.Range({1, 1, 1}, {10, 20, 30})

   assert_equal(3, range:ndim(), "Checking dimension")

   assert_equal(1, range:lower(1), "Checking lower")
   assert_equal(1, range:lower(2), "Checking lower")
   assert_equal(1, range:lower(3), "Checking lower")   

   assert_equal(10, range:upper(1), "Checking upper")
   assert_equal(20, range:upper(2), "Checking upper")
   assert_equal(30, range:upper(3), "Checking upper")

   assert_equal(10, range:shape(1), "Checking shape")
   assert_equal(20, range:shape(2), "Checking shape")
   assert_equal(30, range:shape(3), "Checking shape")
end

-- Run tests
test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
