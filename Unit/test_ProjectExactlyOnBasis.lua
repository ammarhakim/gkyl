-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to project on basis functions
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Basis = require "Basis"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1d_1()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {4},
   }
   local basis = Basis.CartModalSerendipity { ndim = 1, polyOrder = 2 }
   local distf = DataStruct.Field {
      onGrid = grid,
      numComponents = basis:numBasis(),
      ghost = {0, 0},
   }
   local project = Updater.ProjectExactlyOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t, xn)
	 return xn[1]
      end
   }

   -- do projection
   project:advance(0.0, {}, {distf})

   local idx = Lin.IntVec(grid:ndim())
   local xc = Lin.Vec(1)
   local indexer = distf:indexer()
   -- check projection
   for i = 1, grid:numCells(1) do
      grid:setIndex( idx:setValues {i} )
      grid:cellCenter(xc)
      local fItr = distf:get(indexer(i))
      assert_equal(xc[1], fItr[1]/math.sqrt(2), "Checking cell average")
      assert_equal(0.1020620726159657, fItr[2], "Checking slope")
      assert_equal(0.0, fItr[3], "Checking second-moment")
   end

end
-- run tests
test_1d_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
