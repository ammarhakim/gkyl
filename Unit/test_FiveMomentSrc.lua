-- Gkyl ------------------------------------------------------------------------
--
-- Test for Euler equations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {4},
   }
   local elcFluid = DataStruct.Field {
      onGrid = grid,
      numComponents = 5,
      ghost = {1, 1},
   }
   local ionFluid = DataStruct.Field {
      onGrid = grid,
      numComponents = 5,
      ghost = {1, 1},
   }
   local em = DataStruct.Field {
      onGrid = grid,
      numComponents = 8,
      ghost = {1, 1}
   }
   local srcUpdater = Updater.FiveMomentSrc {
      onGrid = grid,
      numFluids = 2,
      charge = {-1.0, 1.0},
      mass = {1.0, 1/1836.2},
      epsilon0 = 1.0,
      scheme = "implicit", -- one of "implicit" or "explicit"
      hasStaticField = false, -- do we have static EM field?
      gravity = 0.0, -- gravitational force
      dir = 0.0, -- direction of force
   }
   srcUpdater:advance(0.0, 1.0, {}, {elcFluid, ionFluid, em})

end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
