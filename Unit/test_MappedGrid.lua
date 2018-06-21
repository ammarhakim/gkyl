-- Gkyl ------------------------------------------------------------------------
--
-- Test for mapped-cartesian grid object
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local math = require("sci.math").generic

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local grid = Grid.MappedCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {10, 20},
      -- map computational to physical space
      mapc2p = function(xc)
	 return xc[1], xc[2]
      end
   }

   assert_equal(0.0, grid:lower(1), "Lower-left computational coordinate [1]")
   assert_equal(0.0, grid:lower(2), "Lower-left computational coordinate [2]")

   assert_equal(1.0, grid:upper(1), "Upper-right computational coordinate [1]")
   assert_equal(1.0, grid:upper(2), "Upper-right computational coordinate [2]")

   grid:setIndex({1,1})
   assert_equal(0.1, grid:dx(1), "Computational cell-space [1]")
   assert_equal(0.05, grid:dx(2), "Computational cell-space [2]")

   local xp = Lin.Vec(grid:ndim())
   -- check if mappings work
   xp[1], xp[2] = grid:mapc2p({0.5, 0.5})
   assert_equal(0.5, xp[1], "Checking mapping [1]")
   assert_equal(0.5, xp[2], "Checking mapping [2]")

   xp[1], xp[2] =  grid:mapc2p({0.75, 0.75})
   assert_equal(0.75, xp[1], "Checking mapping [1]")
   assert_equal(0.75, xp[2], "Checking mapping [2]")

   assert_equal(3, grid:numMetricElems(), "Checking size of metric")
   
   local g = Lin.Vec(grid:numMetricElems())
   -- compute metric quantities
   grid:calcMetric({0.5, 0.5}, g)
   
   assert_equal(1.0, g[1], "Checking metric [1,1]")
   assert_equal(0.0, g[2], "Checking metric [1,2]")
   assert_equal(1.0, g[3], "Checking metric [2,2]")
end

function test_2()
   -- lat-long coordinates on surface of unit sphere: first coordinate
   -- is lattitude and second longitude
   local grid = Grid.MappedCart {
      lower = {10*math.pi/180, 0.0},
      upper = {170*math.pi/180, 2*math.pi},
      cells = {64, 64},
      -- map computational to physical space
      mapc2p = function(xc)
	 local theta, phi = xc[1], xc[2] -- lattitude, longitude
	 return math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta)
      end
   }

   assert_equal(3, grid:rdim(), "Checking rdim")
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

