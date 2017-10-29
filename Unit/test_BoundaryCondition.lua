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

function test_3()
   local bcZeroNormal = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }
   local bcZeroTangent = BoundaryCondition.ZeroTangent { components = {4, 5, 6} }
   local qin, qout = Lin.Vec(6), Lin.Vec(6)

   for i = 1, 6 do qin[i] = i end
   for i = 1, 6 do qout[i] = 10 end
   bcZeroNormal(1, 0.0, {}, qin, qout)
   bcZeroTangent(1, 0.0, {}, qin, qout)      

   assert_equal(0.0, qout[1]+qin[1], "Testing zero-normal BC")
   assert_equal(qout[2], qin[2], "Testing zero-normal BC")
   assert_equal(qout[3], qin[3], "Testing zero-normal BC")

   assert_equal(qout[4], qin[4], "Testing zero-tangent BC")
   assert_equal(0, qout[5]+qin[5], "Testing zero-tangent BC")
   assert_equal(0, qout[6]+qin[6], "Testing zero-tangent BC")

   bcZeroNormal(2, 0.0, {}, qin, qout)
   bcZeroTangent(2, 0.0, {}, qin, qout)
   
   assert_equal(0.0, qout[2]+qin[2], "Testing zero-normal BC")
   assert_equal(qout[3], qin[3], "Testing zero-normal BC")
   assert_equal(qout[1], qin[1], "Testing zero-normal BC")

   assert_equal(qout[5], qin[5], "Testing zero-tangent BC")
   assert_equal(0, qout[4]+qin[4], "Testing zero-tangent BC")
   assert_equal(0, qout[6]+qin[6], "Testing zero-tangent BC")

   bcZeroNormal(3, 0.0, {}, qin, qout)
   bcZeroTangent(3, 0.0, {}, qin, qout)
   
   assert_equal(0.0, qout[3]+qin[3], "Testing zero-normal BC")
   assert_equal(qout[2], qin[2], "Testing zero-normal BC")
   assert_equal(qout[1], qin[1], "Testing zero-normal BC")

   assert_equal(qout[6], qin[6], "Testing zero-tangent BC")
   assert_equal(0, qout[4]+qin[4], "Testing zero-tangent BC")
   assert_equal(0, qout[5]+qin[5], "Testing zero-tangent BC")
end

function test_4()
   local bcConst = BoundaryCondition.Const { components = {2, 3}, values = {-1, -1} }
   local qin, qout = Lin.Vec(5), Lin.Vec(5)

   for i = 1, 5 do qin[i] = i end
   for i = 1, 5 do qout[i] = 10 end
   bcConst(1, 0.0, {}, qin, qout)

   assert_equal(10, qout[1],  "Checking Const bc")
   assert_equal(-1, qout[2], "Checking Const bc")
   assert_equal(-1, qout[3], "Checking Const bc")
   assert_equal(10, qout[4],  "Checking Const bc")
   assert_equal(10, qout[5],  "Checking Const bc")
end

test_1()
test_2()
test_3()
test_4()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
