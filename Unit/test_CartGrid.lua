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
   local grid = Grid.CartGrid {
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
   local grid = Grid.CartGrid {
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

   assert_equal(0.2*0.2, grid:cellVolume(), "Checking volume")
end

function test_3()
   local grid = Grid.CartGrid {
      lower = {0.0, 1.0, 2.0},
      upper = {2.0, 5.0, 10.0},
      cells = {10, 20, 40}
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(0.0, grid:lower(0), "Checking lower")
   assert_equal(1.0, grid:lower(1), "Checking lower")
   assert_equal(2.0, grid:lower(2), "Checking lower")   

   assert_equal(2.0, grid:upper(0), "Checking upper")
   assert_equal(5.0, grid:upper(1), "Checking upper")
   assert_equal(10.0, grid:upper(2), "Checking upper")   

   assert_equal(0.2, grid:dx(0), "Checking dx")
   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")

   assert_equal(0.2*0.2*0.2, grid:cellVolume(), "Checking volume")
end

function test_4()
   local grid = Grid.NonUniformCartGrid {
      cells = {10, 10}
   }

   assert_equal(2, grid:ndim(), "Checking NDIM")

   assert_equal(0.0, grid:lower(0), "Checking lower")
   assert_equal(0.0, grid:lower(1), "Checking lower")

   assert_equal(1.0, grid:upper(0), "Checking upper")
   assert_equal(1.0, grid:upper(1), "Checking upper")

   assert_equal(0.1, grid:dx(0), "Checking dx")
   assert_equal(0.1, grid:dx(1), "Checking dx")

   assert_equal(0.1*0.1, grid:cellVolume(), "Checking volume")   
end

function test_5()
   local grid = Grid.NonUniformCartGrid {   
      lower = {0.0, 1.0, 2.0},
      upper = {2.0, 5.0, 10.0},
      cells = {10, 20, 40}
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(0.0, grid:lower(0), "Checking lower")
   assert_equal(1.0, grid:lower(1), "Checking lower")
   assert_equal(2.0, grid:lower(2), "Checking lower")   

   assert_equal(2.0, grid:upper(0), "Checking upper")
   assert_equal(5.0, grid:upper(1), "Checking upper")
   assert_equal(10.0, grid:upper(2), "Checking upper")   

   assert_equal(0.2, grid:dx(0), "Checking dx")
   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")

   assert_equal(0.2*0.2*0.2, grid:cellVolume(), "Checking volume")
end

-- Run tests
test_1()
test_2()
test_3()
test_4()
test_5()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
