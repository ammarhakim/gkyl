-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to project on basis functions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"
local Lin        = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

function test_1d_1()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {4},
   }
   local basis = Basis.CartModalSerendipity { ndim = 1, polyOrder = 2 }
   local distf = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis(),
      ghost         = {0, 0},
   }
   local project = Updater.ProjectOnBasis {
      onGrid   = grid,
      basis    = basis,
      evaluate = function (t, xn)
	 return xn[1]
      end
   }

   -- Do projection.
   project:advance(0.0, {}, {distf})

   local idx     = Lin.IntVec(grid:ndim())
   local xc      = Lin.Vec(1)
   local indexer = distf:indexer()
   -- Check projection.
   for i = 1, grid:numCells(1) do
      grid:setIndex( idx:setValues {i} )
      grid:cellCenter(xc)
      local fItr = distf:get(indexer(i))
      assert_equal(xc[1], fItr[1]/math.sqrt(2), "Checking cell average")
      assert_equal(0.1020620726159657, fItr[2], "Checking slope")
      assert_equal(0.0, fItr[3], "Checking second-moment")
   end

end

function test_1d_2()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {4},
   }
   local basis = Basis.CartModalSerendipity { ndim = 1, polyOrder = 2 }
   local distf = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis(),
      ghost         = {0, 0},
   }
   local project = Updater.ProjectOnBasis {
      onGrid   = grid,
      basis    = basis,
      evaluate = function (t, xn)
	 return xn[1]^2
      end
   }

   -- Do projection.
   project:advance(0.0, {}, {distf})

   local idx     = Lin.IntVec(grid:ndim())
   local xc      = Lin.Vec(1)
   local indexer = distf:indexer()
   -- Check projection.
   for i = 1, grid:numCells(1) do
      grid:setIndex( idx:setValues {i} )
      grid:cellCenter(xc)
      local fItr = distf:get(indexer(i))
      assert_equal(1.414213562373095*xc[1]^2+0.007365695637359865, fItr[1], "Checking cell average")
      assert_equal(0.2041241452319315*xc[1], fItr[2], "Checking slope")
      --assert_equal(0.0065880784586841, fItr[3], "Checking second-moment " .. i)
   end

end

function test_2d()
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {4, 4},
   }
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 2 }
   local distf = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis(),
      ghost         = {0, 0},
   }
   local project = Updater.ProjectOnBasis {
      onGrid   = grid,
      basis    = basis,
      evaluate = function (t, xn)
	 return xn[1]
      end
   }

   -- Do projection.
   project:advance(0.0, {}, {distf})

   local idx     = Lin.IntVec(grid:ndim())
   local xc      = Lin.Vec(grid:ndim())
   local indexer = distf:indexer()
   -- Check projection.
   for i = 1, grid:numCells(1) do
      for j = 2, grid:numCells(2) do
	 grid:setIndex( idx:setValues {i, j} )
	 grid:cellCenter(xc)
	 local fItr = distf:get(indexer(i,j))
	 assert_equal(xc[1], fItr[1]/2, "Checking cell average")
	 assert_equal(math.sqrt(2)*0.1020620726159657, fItr[2], "Checking slope")
	 assert_equal(0.0, fItr[3], "Checking second-moment")
      end
   end
end

function test_2d_2()
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {4, 4},
   }
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 2 }
   local distf = DataStruct.Field {
      onGrid        = grid,
      numComponents = 2*basis:numBasis(),
      ghost         = {0, 0},
   }
   local project = Updater.ProjectOnBasis {
      onGrid   = grid,
      basis    = basis,
      evaluate = function (t, xn)
	 return xn[1], xn[1]
      end
   }

   -- Do projection.
   project:advance(0.0, {}, {distf})

   local idx     = Lin.IntVec(grid:ndim())
   local xc      = Lin.Vec(grid:ndim())
   local indexer = distf:indexer()
   -- Check projection.
   for i = 1, grid:numCells(1) do
      for j = 2, grid:numCells(2) do
	 grid:setIndex( idx:setValues {i, j} )
	 grid:cellCenter(xc)
	 local fItr = distf:get(indexer(i,j))

	 assert_equal(xc[1], fItr[1]/2, "Checking cell average")
	 assert_equal(math.sqrt(2)*0.1020620726159657, fItr[2], "Checking slope")
	 assert_equal(0.0, fItr[3], "Checking second-moment")

	 assert_equal(xc[1], fItr[9]/2, "Checking cell average")
	 assert_equal(math.sqrt(2)*0.1020620726159657, fItr[10], "Checking slope")
	 assert_equal(0.0, fItr[11], "Checking second-moment")

      end
   end
end

function test_3()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {4},
   }
   local basis = Basis.CartModalSerendipity { ndim = 1, polyOrder = 2 }
   local distf = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis(),
      ghost         = {0, 0},
   }
   local project = Updater.ProjectOnBasis {
      onGrid   = grid,
      basis    = basis,
      evaluate = function (t, xn)
	 return xn[1]
      end
   }

   -- Do projection.
   project:advance(0.0, {}, {distf})

   local idx     = Lin.IntVec(grid:ndim())
   local xc      = Lin.Vec(1)
   local indexer = distf:indexer()
   -- Check projection.
   for i = 1, grid:numCells(1) do
      grid:setIndex( idx:setValues {i} )
      grid:cellCenter(xc)
      local fItr = distf:get(indexer(i))
      assert_equal(xc[1], fItr[1]/math.sqrt(2), "Checking cell average")
      assert_equal(0.1020620726159657, fItr[2], "Checking slope")
      assert_equal(0.0, fItr[3], "Checking second-moment")
   end

   project:setFunc(function (t, xn) return 2*xn[1] end)

   -- Do projection.
   project:advance(0.0, {}, {distf})

   local idx     = Lin.IntVec(grid:ndim())
   local xc      = Lin.Vec(1)
   local indexer = distf:indexer()
   -- Check projection.
   for i = 1, grid:numCells(1) do
      grid:setIndex( idx:setValues {i} )
      grid:cellCenter(xc)
      local fItr = distf:get(indexer(i))
      assert_equal(2*xc[1], fItr[1]/math.sqrt(2), "Checking cell average")
      assert_equal(2*0.1020620726159657, fItr[2], "Checking slope")
      --assert_equal(2*0.0, fItr[3], "Checking second-moment")
   end   

end

-- Run tests.
test_1d_1()
test_1d_2()
test_2d()
test_2d_2()
test_3()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
