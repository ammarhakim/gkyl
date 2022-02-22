-- Gkyl ------------------------------------------------------------------------
--
-- Test for range objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local ZeroArray = require "DataStruct.ZeroArray"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local arr = ZeroArray.Array( ZeroArray.double, 3, 10)
   assert_equal(10, arr.size, "Testing array size")
   assert_equal(3, arr.ncomp, "Testing array ncomp")
   assert_equal(8, arr.elemsz, "Testing array elemsz")
end

-- Run tests
test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
