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

   local localRange = elcFluid:localRange()
   local elcIdxr = elcFluid:indexer()   
 
   local gasGamma = 5. / 3.
   local epsilon0 = 1.
   local q = 1.
   local mi = 1.
   local me = 0.25

   for i = localRange:lower(1), localRange:upper(1) do
      local elcItr = elcFluid:get(elcIdxr(i))

      local ne, vxe, vye, vze, pe = 1.1, 1.1, 2.1, 3.1, 1.1
      local rhoe = ne * me

      elcItr[1] = rhoe
      elcItr[2] = rhoe * vxe
      elcItr[3] = rhoe * vye
      elcItr[4] = rhoe * vze
      elcItr[5] = pe + vxe^2 * rhoe
      elcItr[6] = vxe * vye * rhoe
      elcItr[7] = vxe * vze * rhoe
      elcItr[8] = pe + vye^2 * rhoe
      elcItr[9] = vye * vze * rhoe
      elcItr[10] = pe + vze^2 * rhoe
   end
   
   local srcUpdater = Updater.TenMomentRelax {
      onGrid = grid,
      charge = -1,
      mass = me,
      scheme = scheme,
      k = 1.,
   }

   local printValues = function()
      for i = localRange:lower(1), localRange:upper(1) do
         local elcItr = elcFluid:get(elcIdxr(i))
    
         print(string.format("elc = {%g, %g, %g, %g,\n\t%g, %g, %g, %g, %g, %g}",
           elcItr[1], elcItr[2], elcItr[3], elcItr[4], elcItr[5],
           elcItr[6], elcItr[7], elcItr[8], elcItr[9], elcItr[10]))
      end
      print("\n")
   end

   print("Before")
   printValues()
   srcUpdater:setupDtAndCflRate(dt, nil)
   srcUpdater:advance(0.0, {}, {elcFluid, ionFluid, em})
   print("After")
   printValues()

end

-- run tests
for _,scheme in ipairs({"explicit"}) do
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
