-- Gkyl ------------------------------------------------------------------------
--
-- Test for dynamic vectors
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local dynVec = DataStruct.DynVector { numComponents = 3 }
   assert_equal(dynVec:numComponents(), 3, "Testing number of components")

   dynVec:appendData(1.0, {1.5, 2.5, 3.5})

   assert_equal(dynVec:lastTime(), 1.0, "Checking last inserted time")
   local lv = dynVec:lastData()
   assert_equal(lv[1], 1.5, "Checking last inserted data")
   assert_equal(lv[2], 2.5, "Checking last inserted data")
   assert_equal(lv[3], 3.5, "Checking last inserted data")

   dynVec:appendData(2.0, {2.5, 3.5, 4.5})

   assert_equal(dynVec:lastTime(), 2.0, "Checking last inserted time")
   local lv = dynVec:lastData()
   assert_equal(lv[1], 2.5, "Checking last inserted data")
   assert_equal(lv[2], 3.5, "Checking last inserted data")
   assert_equal(lv[3], 4.5, "Checking last inserted data")
end

function test_2()
   local dynVec = DataStruct.DynVector { numComponents = 2 }
   assert_equal(dynVec:numComponents(), 2, "Testing number of components")

   for i = 1, 10 do
      dynVec:appendData(0.1*i, {2.5*i^2, 2.5*i^2+0.5})
   end
   dynVec:write("test_2.bp", 1.5)
end

test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
