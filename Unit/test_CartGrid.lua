-- Gkyl ------------------------------------------------------------------------
--
-- Test for cartesian grid objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Grid = require "Grid"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   grid = Grid.CartGrid {
      cells = {10, 20}
   }

   assert_equal(2, grid:ndim(), "Checking NDIM")
   assert_equal(0.0, grid:lower(0), "Checking lower")
   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:upper(0), "Checking upper")
   assert_equal(1.0, grid:upper(1), "Checking upper")
   assert_equal(0.1, grid:dx(0), "Checking dx")
   assert_equal(0.05, grid:dx(1), "Checking dx")
end

function test_2()
   grid = Grid.CartGrid {
      lower = {0.0, 1.0},
      upper = {2.0, 5.0},
      cells = {10, 20}
   }

   assert_equal(2, grid:ndim(), "Checking NDIM")
   assert_equal(0.0, grid:lower(0), "Checking lower")
   assert_equal(1.0, grid:lower(1), "Checking lower")
   assert_equal(2.0, grid:upper(0), "Checking upper")
   assert_equal(5.0, grid:upper(1), "Checking upper")
   assert_equal(0.2, grid:dx(0), "Checking dx")
   assert_equal(0.2, grid:dx(1), "Checking dx")
end

test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\n**** FAILED %d tests", stats.fail))
else
   print("All tests passed!")
end
