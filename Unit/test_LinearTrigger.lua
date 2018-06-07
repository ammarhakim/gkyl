-- Gkyl ------------------------------------------------------------------------
--
-- Test for fields on cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local LinearTrigger = require "Lib.LinearTrigger"
local Unit = require "Unit"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local trig = LinearTrigger(0.0, 1.0, 10)

   assert_equal(true, trig(0.01), "Checking linear trigger 0.11")
   
   assert_equal(true, trig(0.11), "Checking linear trigger 0.11")
   assert_equal(false, trig(0.11), "Checking linear trigger 0.11")

   assert_equal(true, trig(0.21), "Checking linear trigger 0.21")
   assert_equal(false, trig(0.21), "Checking linear trigger 0.21")   

   assert_equal(true, trig(0.31), "Checking linear trigger 0.31")
   assert_equal(false, trig(0.31), "Checking linear trigger 0.31") 

   assert_equal(true, trig(0.41), "Checking linear trigger 0.41")
   assert_equal(false, trig(0.41), "Checking linear trigger 0.41")   

   assert_equal(true, trig(0.51), "Checking linear trigger 0.51")
   assert_equal(false, trig(0.51), "Checking linear trigger 0.51")   

   assert_equal(true, trig(0.61), "Checking linear trigger 0.61")
   assert_equal(false, trig(0.61), "Checking linear trigger 0.61")   

   assert_equal(true, trig(0.71), "Checking linear trigger 0.71")
   assert_equal(false, trig(0.71), "Checking linear trigger 0.71")   

   assert_equal(true, trig(0.81), "Checking linear trigger 0.81")
   assert_equal(false, trig(0.81), "Checking linear trigger 0.81")

   assert_equal(true, trig(0.91), "Checking linear trigger 0.91")
   assert_equal(false, trig(0.91), "Checking linear trigger 0.91")

   assert_equal(true, trig(1.01), "Checking linear trigger 1.01")
   assert_equal(false, trig(1.01), "Checking linear trigger 1.01")
end

function test_2()
   local trig = LinearTrigger(0.0, 1.5, 22)

   local nsteps = 25
   local dt = 1.5/nsteps

   local count = 0
   for s = 1, nsteps+1 do
      if trig(s*dt) then
	 count = count + 1
      end
   end
   assert_equal(23, count, "Checking if count is correct")
end

function test_3()
   local trig = LinearTrigger(0.0, 10.0, 10)
   assert_equal(true, trig(0.1), "Checking trigger")
end

test_1()
test_2()
test_3()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
