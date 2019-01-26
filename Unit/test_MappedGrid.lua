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

   local idx = Lin.IntVec(grid:ndim())

   grid:setIndex( idx:setValues {1,1} )
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

function test_salpha()

   local R0 = 0.5
   local r0 = 0.1*R0 -- minor radius of center of flux tube
   local R = R0 + r0 -- major radius of center of flux tube
   local dr = 0.01*R0
   local q0 = 2
   local shat = 1
   local grid = Grid.MappedCart {
      lower = {r0 - dr/2, -dr/2, -math.pi}, -- configuration space lower left
      upper = {r0 + dr/2,  dr/2, math.pi}, -- configuration space upper right
      cells = {8, 8, 8}, -- configuration space cells
         mapc2p = function(xc)
         local x, y, z = xc[1], xc[2], xc[3]
         local q = q0*(1+(x-r0)*shat/r0)
         local X = (R0 + x*math.cos(z))*math.cos(q*z - q0/r0*y)
         local Y = (R0 + x*math.cos(z))*math.sin(q*z - q0/r0*y)
         local Z = x*math.sin(z)
         return X, Y, Z
      end,
   }
   local g = Lin.Vec(grid:numMetricElems())
   local gContra = Lin.Vec(grid:numMetricElems())
   -- compute metric quantities

   xc_test = {1.0*r0, 0.0, math.pi/5}

   grid:calcMetric(xc_test, g)

   local g_xx = function(xc) 
      local x, y, z = xc[1], xc[2], xc[3]
      return 1 + z^2*(R0 + x*math.cos(z))^2*q0^2/r0^2*shat^2
   end

   assert_equal(g_xx(xc_test), g[1], "Checking metric g_xx")

   local g_xy = function(xc)
      local x, y, z = xc[1], xc[2], xc[3]
      return -q0^2/r0^2*z*(R0 + x*math.cos(z))^2*shat
   end

   assert_equal(g_xy(xc_test), g[2], "Checking metric g_xy")

   local g_yy = function(xc)
      local x, y, z = xc[1], xc[2], xc[3]
      return q0^2/r0^2*(R0 + x*math.cos(z))^2
   end

   assert_equal(g_yy(xc_test), g[4], "Checking metric g_yy")
      
   local g_xz = function(xc)
      local x, y, z = xc[1], xc[2], xc[3]
      return z*q0^2/r0*(R0 + x*math.cos(z))^2*shat*(1+(x-r0)*shat/r0)
   end

   assert_equal(g_xz(xc_test), g[3], "Checking metric g_xz")
      
   local g_yz = function(xc)
      local x, y, z = xc[1], xc[2], xc[3]
      return -q0^2/r0*(R0 + x*math.cos(z))^2*(1+(x-r0)*shat/r0)
   end

   assert_equal(g_yz(xc_test), g[5], "Checking metric g_yz")
      
   local g_zz = function(xc)
      local x, y, z = xc[1], xc[2], xc[3]
      return (R0 + x*math.cos(z))^2*q0^2*(1+(x-r0)*shat/r0)^2+x^2
   end

   assert_equal(g_zz(xc_test), g[6], "Checking metric g_zz")

end

-- Run tests
test_1()
test_2()
test_salpha()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end

