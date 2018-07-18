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

function test_1(scheme)
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

      elcItr[1] = 1
      elcItr[2] = 1
      elcItr[3] = 2
      elcItr[4] = 3
      elcItr[5] = 8.5

      ionItr[1] = 12
      ionItr[2] = 12
      ionItr[3] = 24
      ionItr[4] = 36
      ionItr[5] = 99

      emItr[1] = 4
      emItr[2] = 5
      emItr[3] = 6
      emItr[4] = 7
      emItr[5] = 8
      emItr[6] = 9

      emItr[7] = 3.33
      emItr[8] = 4.56
   end
   
   local srcUpdater = Updater.FiveMomentSrc {
      onGrid = grid,
      numFluids = 2,
      charge = {-1.0, 1.0},
      mass = {1.0, 4.0},
      epsilon0 = 1.0,
      elcErrorSpeedFactor = 1,
      mgnErrorSpeedFactor = 1,
      scheme = scheme,
   }
   srcUpdater:advance(0.0, 1.0, {}, {elcFluid, ionFluid, em})

   for i = localRange:lower(1), localRange:upper(1) do
      local elcItr, ionItr, emItr = elcFluid:get(elcIdxr(i)), ionFluid:get(ionIdxr(i)), em:get(emIdxr(i))
 
      print(string.format("  point %d", i))
      print(string.format("elc = {%g, %g, %g, %g, %g}", elcItr[1], elcItr[2], elcItr[3], elcItr[4], elcItr[5]))
      print(string.format("ion = {%g, %g, %g, %g, %g}", ionItr[1], ionItr[2], ionItr[3], ionItr[4], ionItr[5]))
      print(string.format("ele = {%g, %g, %g}", emItr[1], emItr[2], emItr[3]))
      print(string.format("mag = {%g, %g, %g}", emItr[4], emItr[5], emItr[6]))
      print(string.format("phiE, phiM = {%g, %g}", emItr[7], emItr[8]))
   end
   print("\n")
end

-- run tests
for _,scheme in ipairs({"time-centered", "ssp-rk3"}) do
  print(string.format("5m source update test with scheme: %s", scheme))
  test_1(scheme)
end

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
