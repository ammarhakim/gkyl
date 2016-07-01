-- Gkyl
-- ------------------------------------------------------------------------
--
-- Test for boundary condition applicator objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Lin = require "Lib.Linalg"
local BoundaryCondition = require "Updater.BoundaryCondition"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local bcCopy = BoundaryCondition.Copy { components = {2, 3} }
   local qin, qout = Lin.Vec(5), Lin.Vec(5)

   for i = 1, 5 do qin[i] = i end
   for i = 1, 5 do qout[i] = 10 end
   bcCopy(1, 0.0, {}, qin, qout)

   assert_equal(10, qout[1],  "Checking Copy bc")
   assert_equal(qin[2], qout[2], "Checking Copy bc")
   assert_equal(qin[3], qout[3], "Checking Copy bc")
   assert_equal(10, qout[4],  "Checking Copy bc")
   assert_equal(10, qout[5],  "Checking Copy bc")
end

function test_2()
   local bcCopy = BoundaryCondition.Copy { components = {1, 2, 3, 4, 5} }
   local qin, qout = Lin.Vec(5), Lin.Vec(5)

   for i = 1, 5 do qin[i] = i end
   for i = 1, 5 do qout[i] = 10 end
   bcCopy(1, 0.0, {}, qin, qout)

   for i = 1, 5 do
      assert_equal(qin[i], qout[i], "Checking Copy bc")
   end
end

test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
