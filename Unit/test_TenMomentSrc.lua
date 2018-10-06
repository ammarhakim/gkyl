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
      cells = {1},
   }
   local elcFluid = DataStruct.Field {
      onGrid = grid,
      numComponents = 10,
   }
   local ionFluid = DataStruct.Field {
      onGrid = grid,
      numComponents = 10,
   }
   local em = DataStruct.Field {
      onGrid = grid,
      numComponents = 8,
   }

   local localRange = elcFluid:localRange()
   local elcIdxr, ionIdxr, emIdxr = elcFluid:indexer(), ionFluid:indexer(), em:indexer()   
  -- initialize fluids/fields
 
   local gasGamma = 5. / 3.
   local epsilon0 = 1.
   local q = 1.
   local mi = 1.
   local me = 0.25
   local elcErrorSpeedFactor = 1
   local mgnErrorSpeedFactor = 1

   for i = localRange:lower(1), localRange:upper(1) do
      local elcItr, ionItr, emItr = elcFluid:get(
         elcIdxr(i)), ionFluid:get(ionIdxr(i)), em:get(emIdxr(i))

      local ne, vxe, vye, vze, pe = 1.1, 1.1, 2.1, 3.1, 1.1
      local ni, vxi, vyi, vzi, pi = 1.2, 1.2, 2.2, 3.2, 1.2
      local Ex, Ey, Ez = 1.3, 2.3, 3.3
      local Bx, By, Bz = 1.4, 2.4, 3.4
      local phiE, phiB = 0.5, 1.5

      local rhoe = ne * me
      local rhoi = ni * mi

      elcItr[1] = rhoe
      elcItr[2] = rhoe * vxe
      elcItr[3] = rhoe * vye
      elcItr[4] = rhoe * vze
      elcItr[5] = pe + vxe^2 * rhoe
      elcItr[6] = pe + vye^2 * rhoe
      elcItr[7] = pe + vze^2 * rhoe
      elcItr[8] = vxe * vye * rhoe
      elcItr[9] = vxe * vze * rhoe
      elcItr[10] = vye * vze * rhoe

      ionItr[1] = rhoi
      ionItr[2] = rhoi * vxi
      ionItr[3] = rhoi * vyi
      ionItr[4] = rhoi * vzi
      ionItr[5] = pi + vxi^2 * rhoi
      ionItr[6] = pi + vyi^2 * rhoi
      ionItr[7] = pi + vzi^2 * rhoi
      ionItr[8] = vxi * vyi * rhoi
      ionItr[9] = vxi * vzi * rhoi
      ionItr[10] = vyi * vzi * rhoi

      emItr[1] = Ex
      emItr[2] = Ey
      emItr[3] = Ez
      emItr[4] = Bx
      emItr[5] = By
      emItr[6] = Bz

      emItr[7] = phiE
      emItr[8] = phiB
   end
   
   local srcUpdater = Updater.TenMomentSrc {
      onGrid = grid,
      numFluids = 2,
      charge = {-q, q},
      mass = {me, mi},
      epsilon0 = epsilon0,
      elcErrorSpeedFactor = elcErrorSpeedFactor,
      mgnErrorSpeedFactor = mgnErrorSpeedFactor,
      scheme = scheme,
   }

   local printValues = function()
      for i = localRange:lower(1), localRange:upper(1) do
         local elcItr, ionItr, emItr = elcFluid:get(
           elcIdxr(i)), ionFluid:get(ionIdxr(i)), em:get(emIdxr(i))
    
         print(string.format("elc = {%g, %g, %g, %g,\n\t%g, %g, %g, %g, %g, %g}",
           elcItr[1], elcItr[2], elcItr[3], elcItr[4], elcItr[5],
           elcItr[6], elcItr[7], elcItr[8], elcItr[9], elcItr[10]))
         print(string.format("ion = {%g, %g, %g, %g,\n\t%g, %g, %g, %g, %g, %g}",
           ionItr[1], ionItr[2], ionItr[3], ionItr[4], ionItr[5],
           ionItr[6], ionItr[7], ionItr[8], ionItr[9], ionItr[10]))
         print(string.format("ele = {%g, %g, %g}",
           emItr[1], emItr[2], emItr[3]))
         print(string.format("mag = {%g, %g, %g}",
           emItr[4], emItr[5], emItr[6]))
         print(string.format("phiE, phiM = {%g, %g}", emItr[7], emItr[8]))
      end
      print("\n")
   end

   print("Before")
   printValues()
   srcUpdater:advance(0.0, dt, {}, {elcFluid, ionFluid, em})
   print("After")
   printValues()

end

-- run tests
for _,scheme in ipairs({"time-centered"}) do
  for _,dt in ipairs({0.1}) do
     print(string.format("1m source update test, scheme = %s, dt = %g", scheme, dt))
     test_1(scheme, dt)
  end
end

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
