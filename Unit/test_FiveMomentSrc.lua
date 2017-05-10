-- Gkyl ------------------------------------------------------------------------
--
-- Test for five-moment source updater
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

   local localRange = elcFluid:localRange()
   local elcIdxr, ionIdxr, emIdxr = elcFluid:indexer(), ionFluid:indexer(), em:indexer()   
   -- initialize fluids/fields

   for i = localRange:lower(1), localRange:upper(1) do
      local elcItr, ionItr, emItr = elcFluid:get(elcIdxr(i)), ionFluid:get(ionIdxr(i)), em:get(emIdxr(i))

      elcItr[1] = 2.0
      elcItr[2] = 2.0
      elcItr[3] = 0.0
      elcItr[3] = 0.0

      ionItr[1] = 1.0
      ionItr[2] = 0.0
      ionItr[3] = 0.0
      ionItr[3] = 0.0

      emItr[1] = 1.0
   end
   
   local srcUpdater = Updater.FiveMomentSrc {
      onGrid = grid,
      numFluids = 2,
      charge = {-1.0, 1.0},
      mass = {1.0, 1.0},
      epsilon0 = 1.0,
      scheme = "ssp-rk3", -- one of "ssp-rk3", "modified-boris"  or "time-centered"
      hasStaticField = false, -- do we have static EM field?
      gravity = 0.0, -- gravitational force
      dir = 0.0, -- direction of force
   }

   srcUpdater:advance(0.0, 1.0, {}, {elcFluid, ionFluid, em})

   for i = localRange:lower(1), localRange:upper(1) do
      local elcItr, ionItr, emItr = elcFluid:get(elcIdxr(i)), ionFluid:get(ionIdxr(i)), em:get(emIdxr(i))

      print("elc", elcItr[2], elcItr[3], elcItr[4])
      print("ion", ionItr[2], ionItr[3], ionItr[4])
   end
end

-- run tests
test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
