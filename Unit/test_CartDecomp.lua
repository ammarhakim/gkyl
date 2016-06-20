-- Gkyl ------------------------------------------------------------------------
--
-- Test for decomposition of Cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local DecompRegionCalc = require "Lib.CartDecomp"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local decomp = DecompRegionCalc.CartProd { cuts = {2, 3} }
   assert_equal(2, decomp:ndim(), "Checking ndim")
   
   assert_equal(2, decomp:cuts(1), "Checking cuts")
   assert_equal(3, decomp:cuts(2), "Checking cuts")

   decomp:decompose(Range.Range({1, 1}, {10, 10}))
   -- assert_equal(1, decomp:lower(1, 1), "Checking cuts")
   -- assert_equal(5, decomp:upper(1, 1), "Checking cuts")

   -- assert_equal(6, decomp:lower(1, 2), "Checking cuts")
   -- assert_equal(10, decomp:upper(1, 2), "Checking cuts")

   -- assert_equal(1, decomp:lower(2, 1), "Checking cuts")
   -- assert_equal(4, decomp:upper(2, 1), "Checking cuts")

   -- assert_equal(5, decomp:lower(2, 2), "Checking cuts")
   -- assert_equal(7, decomp:upper(2, 2), "Checking cuts")

   -- assert_equal(8, decomp:lower(2, 3), "Checking cuts")
   -- assert_equal(10, decomp:upper(2, 3), "Checking cuts")
end

-- Run tests
test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
