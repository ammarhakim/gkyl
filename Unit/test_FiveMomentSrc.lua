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

function test_1(scheme, dt)
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
 
   local gasGamma = 5. / 3.
   local rhoe, vxe, vye, vze, pe = 1, 1, 2, 3, 1
   local rhoi, vxi, vyi, vzi, pi = 12, 1, 2, 3, 10
   local Ex, Ey, Ez = 4, 5, 6
   local Bx, By, Bz = 7, 8, 9

   local epsilon0 = 1
   local q = 1
   local me = 1
   local omega_pe = math.sqrt(rhoe / epsilon0 * q^2 / me^2)
   local B = math.sqrt(Bx^2 + By^2, Bz^2)
   local Omega_Ce = B * q / me

   print("1/omega_pe =", 1/omega_pe)
   print("1/Omega_Ce =", 1/Omega_Ce)

   for i = localRange:lower(1), localRange:upper(1) do
      local elcItr, ionItr, emItr = elcFluid:get(elcIdxr(i)), ionFluid:get(ionIdxr(i)), em:get(emIdxr(i))

      ee = 0.5 * (vxe^2 + vye^2 + vze^2) * rhoe + pe / (gasGamma - 1)
      elcItr[1] = rhoe
      elcItr[2] = rhoe * vxe
      elcItr[3] = rhoe * vye
      elcItr[4] = rhoe * vze
      elcItr[5] = ee

      ei = 0.5 * (vxi^2 + vyi^2 + vzi^2) * rhoi + pi / (gasGamma - 1)
      ionItr[1] = rhoi
      ionItr[2] = rhoi * vxi
      ionItr[3] = rhoi * vyi
      ionItr[4] = rhoi * vzi
      ionItr[5] = ei

      emItr[1] = Ex
      emItr[2] = Ey
      emItr[3] = Ez
      emItr[4] = Bx
      emItr[5] = By
      emItr[6] = Bz

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
   srcUpdater:advance(0.0, dt, {}, {elcFluid, ionFluid, em})

   for i = localRange:lower(1)+2, localRange:upper(1)-1 do
      local elcItr, ionItr, emItr = elcFluid:get(elcIdxr(i)), ionFluid:get(ionIdxr(i)), em:get(emIdxr(i))
 
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
  for _,dt in ipairs({1}) do
     print(string.format("5m source update test, scheme = %s, dt = %g", scheme, dt))
     test_1(scheme, dt)
  end
end

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
