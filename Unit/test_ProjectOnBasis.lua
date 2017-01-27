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

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1d()
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
   local project = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t, xn)
	 return xn[1]
      end
   }
end

function test_2d()
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {4, 4},
   }
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 2 }
   local distf = DataStruct.Field {
      onGrid = grid,
      numComponents = basis:numBasis(),
      ghost = {0, 0},
   }
   local project = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function (t, xn)
	 return xn[1]
      end
   }
end

-- run tests
test_1d()
test_2d()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
