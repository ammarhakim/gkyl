-- Gkyl ------------------------------------------------------------------------
--
-- Test for cartesian grid objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local grid = Grid.RectCart {
      cells = {10, 20}
   }

   -- (just make sure setIndex() method works. For RectCart object
   -- setting index is not needed)
   idx = Lin.IntVec(grid:ndim())
   idx[1], idx[2] = 1, 1
   grid:setIndex(idx)
   
   assert_equal(2, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(0.0, grid:lower(2), "Checking lower")

   assert_equal(1.0, grid:upper(1), "Checking upper")
   assert_equal(1.0, grid:upper(2), "Checking upper")

   assert_equal(0.1, grid:dx(1), "Checking dx")
   assert_equal(0.05, grid:dx(2), "Checking dx")
end

function test_2()
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {2.0, 5.0},
      cells = {10, 20}
   }

   assert_equal(2, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")

   assert_equal(0.2*0.2, grid:cellVolume(), "Checking volume")
end

function test_3()
   local grid = Grid.RectCart {
      lower = {0.0, 1.0, 2.0},
      upper = {2.0, 5.0, 10.0},
      cells = {10, 20, 40}
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")
   assert_equal(2.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")
   assert_equal(10.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")
   assert_equal(0.2, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.2*0.2, grid:cellVolume(), "Checking volume")
end

function test_4()
   local grid = Grid.NonUniformRectCart {
      cells = {10, 10}
   }

   assert_equal(2, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(10, grid:numCells(2), "Checking numCells")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(0.0, grid:lower(2), "Checking lower")

   assert_equal(1.0, grid:upper(1), "Checking upper")
   assert_equal(1.0, grid:upper(2), "Checking upper")

   assert_equal(0.1, grid:dx(1), "Checking dx 1")
   assert_equal(0.1, grid:dx(2), "Checking dx 2")

   assert_equal(0.1*0.1, grid:cellVolume(), "Checking volume")   
end

function test_5()
   local grid = Grid.NonUniformRectCart {   
      lower = {0.0, 1.0, 2.0},
      upper = {2.0, 5.0, 10.0},
      cells = {10, 20, 40},
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")
   assert_equal(2.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")
   assert_equal(10.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")
   assert_equal(0.2, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.2*0.2, grid:cellVolume(), "Checking volume")
end

function test_5()
   local grid = Grid.NonUniformRectCart {   
      lower = {0.0, 1.0, 2.0},
      upper = {2.0, 5.0, 10.0},
      cells = {10, 20, 40},
      -- functions mapping computational space to physical space
      mappings = {
	 function (zeta)
	    return zeta
	 end,
	 function (zeta)
	    return zeta	    
	 end,
	 function (zeta)
	    return zeta
	 end,
      }
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")
   assert_equal(2.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")
   assert_equal(10.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")
   assert_equal(0.2, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.2*0.2, grid:cellVolume(), "Checking volume")
end

function test_6()
   local grid = Grid.NonUniformRectCart {   
      cells = {10, 20, 40},
      -- functions mapping computational space to physical space
      mappings = {
	 function (zeta)
	    return 2*zeta
	 end,
	 function (zeta)
	    return 2*zeta	    
	 end,
	 function (zeta)
	    return 2*zeta
	 end,
      }
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(0.0, grid:lower(2), "Checking lower")
   assert_equal(0.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(2.0, grid:upper(2), "Checking upper")
   assert_equal(2.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.1, grid:dx(2), "Checking dx")
   assert_equal(0.05, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.1*0.05, grid:cellVolume(), "Checking volume")
end

function test_7()
   local grid = Grid.NonUniformRectCart {   
      cells = {3},
      -- functions mapping computational space to physical space
      mappings = {
	 function (zeta)
	    return zeta*zeta
	 end,
      }
   }

   assert_equal(1, grid:ndim(), "Checking NDIM")

   assert_equal(3, grid:numCells(1), "Checking numCells")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:upper(1), "Checking upper")

   idx = Lin.IntVec(1) -- for indexing grid

   idx[1] = 1
   grid:setIndex(idx)
   assert_equal(1/9, grid:dx(1), "Checking dx")
   assert_equal(1/9, grid:cellVolume(), "Checking volume")

   idx[1] = 2
   grid:setIndex(idx)
   assert_equal(1/3, grid:dx(1), "Checking dx")
   assert_equal(1/3, grid:cellVolume(), "Checking volume")

   idx[1] = 3
   grid:setIndex(idx)
   assert_equal(5/9, grid:dx(1), "Checking dx")
   assert_equal(5/9, grid:cellVolume(), "Checking volume")
end

function test_8()
   local grid = Grid.NonUniformRectCart { cells = {3} }
   local xn = grid:nodeCoords(1)
   -- set nodes manually (this is the mapping zeta^2)
   xn[1] = 0.0
   xn[2] = 1/3*1/3
   xn[3] = 2/3*2/3
   xn[4] = 1*1
   -- done

   assert_equal(1, grid:ndim(), "Checking NDIM")

   assert_equal(3, grid:numCells(1), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:upper(1), "Checking upper")

   idx = Lin.IntVec(1) -- for indexing grid

   idx[1] = 1
   grid:setIndex(idx)
   assert_equal(1/9, grid:dx(1), "Checking dx")
   assert_equal(1/9, grid:cellVolume(), "Checking volume")

   idx[1] = 2
   grid:setIndex(idx)
   assert_equal(1/3, grid:dx(1), "Checking dx")
   assert_equal(1/3, grid:cellVolume(), "Checking volume")

   idx[1] = 3
   grid:setIndex(idx)
   assert_equal(5/9, grid:dx(1), "Checking dx")
   assert_equal(5/9, grid:cellVolume(), "Checking volume")   
end

-- Run tests
test_1()
test_2()
test_3()
test_4()
test_5()
test_6()
test_7()
test_8()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
