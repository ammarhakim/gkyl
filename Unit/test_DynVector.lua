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
local Mpi = require "Comm.Mpi"

local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function log(msg)
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   if rank == 0 then
      print(msg)
   end
end

function allReduceOneInt(localv)
   local sendbuf, recvbuf = new("int[1]"), new("int[1]")
   sendbuf[0] = localv
   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[0]
end

function test_1()
   local dynVec = DataStruct.DynVector { numComponents = 3 }
   assert_equal(dynVec:numComponents(), 3, "Testing number of components")

   dynVec:appendData(1.0, {1.5, 2.5, 3.5})

   assert_equal(dynVec:lastTime(), 1.0, "Checking last inserted time")
   local tm, lv = dynVec:lastData()
   assert_equal(lv[1], 1.5, "Checking last inserted data")
   assert_equal(lv[2], 2.5, "Checking last inserted data")
   assert_equal(lv[3], 3.5, "Checking last inserted data")

   dynVec:appendData(2.0, {2.5, 3.5, 4.5})

   assert_equal(dynVec:lastTime(), 2.0, "Checking last inserted time")
   local tm, lv = dynVec:lastData()
   assert_equal(lv[1], 2.5, "Checking last inserted data")
   assert_equal(lv[2], 3.5, "Checking last inserted data")
   assert_equal(lv[3], 4.5, "Checking last inserted data")
end

function test_2()
   local dynVec = DataStruct.DynVector { numComponents = 2 }
   assert_equal(dynVec:numComponents(), 2, "Testing number of components")

   for i = 1, 5 do
      dynVec:appendData(0.1*i, {2.5*i^2, 2.5*i^2+0.5})
   end
   dynVec:write("test.bp", 1.5, 0)
   for i = 6, 10 do
      dynVec:appendData(0.1*i, {2.5*i^2, 2.5*i^2+0.5})
   end
   dynVec:write("test.bp", 1.5, 1)
end

function test_2r()
   local dynVec = DataStruct.DynVector { numComponents = 2 }
   assert_equal(dynVec:numComponents(), 2, "Testing number of components")

   local tmStamp, frNum = dynVec:read("test.bp")

   -- check contents
   assert_equal(10, dynVec:size(), "Checking size")

   local tmMesh = dynVec:timeMesh()
   for i = 1, dynVec:size() do
      assert_equal(0.1*i, tmMesh[i], "Checking time-mesh")
   end

   local dynData = dynVec:data()
   for i = 1, dynVec:size() do
      local v = dynData[i]
      assert_equal(2.5*i^2, v[1], "Checking contents[1]")
      assert_equal(2.5*i^2+0.5, v[2], "Checking contents[2]")
   end
end

function test_3()
   local dynVec = DataStruct.DynVector { numComponents = 2 }
   for i = 1, 10 do
      dynVec:appendData(0.1*i, {2.5*i^2, 2.5*i^2+0.5})
   end
   local tm, lv = dynVec:removeLast()
   assert_equal(tm, 0.1*10, "Testing removed time")
   assert_equal(lv[1], 2.5*10^2, "Testing removed value")
   assert_equal(lv[2], 2.5*10^2+0.5, "Testing removed value")

   local tm, lv = dynVec:lastData()
   assert_equal(tm, 0.1*9, "Testing last time")
   assert_equal(lv[1], 2.5*9^2, "Testing last value")
   assert_equal(lv[2], 2.5*9^2+0.5, "Testing last value")
end

function test_4()
   local dynVec = DataStruct.DynVector { numComponents = 2 }
   dynVec:appendData(0.1, {2.5, 3.5})

   local tm, lv = dynVec:lastData()
   assert_equal(0.1, tm, "Testing last time")
   assert_equal(2.5, lv[1], "Testing last value")
   assert_equal(3.5, lv[2], "Testing last value")

   dynVec:assignLastVal(2, {1.75, 1.25})
   tm, lv = dynVec:lastData()
   assert_equal(0.1, tm, "Testing last time")
   assert_equal(3.5, lv[1], "Testing last value after assignLastVal")
   assert_equal(2.5, lv[2], "Testing last value after assignLastVal")

   local dummyVec = DataStruct.DynVector { numComponents = 2 }
   dummyVec:appendData(0.1, {7.5, 4.5})

   for i = 2, 10 do
      dynVec:appendData(0.1*i, {2.5+i, 3.5+i})
   end
   dynVec:copyLast(dummyVec)
   tm, lv = dynVec:lastData()
   assert_equal(0.1, tm, "Testing last time")
   assert_equal(7.5, lv[1], "Testing last value after copyLast")
   assert_equal(4.5, lv[2], "Testing last value after copyLast")

   for i = 11, 20 do
      dynVec:appendData(0.1*i, {2.5+i, 3.5+i})
   end
   dynVec:accumulateLastOne(2.0,dummyVec)
   tm, lv = dynVec:lastData()
   assert_equal(2.0, tm, "Testing last time")
   assert_equal(22.5+7.5*2, lv[1], "Testing last value after accumulateLastOne")
   assert_equal(23.5+4.5*2, lv[2], "Testing last value after accumulateLastOne")

   dynVec:assignLast(0.2, dummyVec)
   tm, lv = dynVec:lastData()
   assert_equal(2.0, tm, "Testing last time")
   assert_equal(7.5*0.2, lv[1], "Testing last value after assignLast")
   assert_equal(4.5*0.2, lv[2], "Testing last value after assignLast")

   dynVec:combineLast(1.0, dummyVec, 2.0, dummyVec)
   tm, lv = dynVec:lastData()
   assert_equal(2.0, tm, "Testing last time")
   assert_equal(7.5*3.0, lv[1], "Testing last value after combineLast")
   assert_equal(4.5*3.0, lv[2], "Testing last value after combineLast")

   dynVec:accumulateLast(1.0, dummyVec, 2.0, dummyVec)
   tm, lv = dynVec:lastData()
   assert_equal(2.0, tm, "Testing last time")
   assert_equal(7.5*6.0, lv[1], "Testing last value after accumulateLast")
   assert_equal(4.5*6.0, lv[2], "Testing last value after accumulateLast")

   dynVec:assignLastTime(3.0)
   tm, lv = dynVec:lastData()
   assert_equal(3.0, tm, "Testing last time after assignLastTime")

   dynVec:appendLast(dummyVec)
   tm, lv = dynVec:lastData()
   assert_equal(0.1, tm, "Testing last time after appendLast")
   assert_equal(7.5, lv[1], "Testing last value after appendLast")
   assert_equal(4.5, lv[2], "Testing last value after appendLast")

end

test_1()
test_2()
test_2r()
test_3()
test_4()

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))   
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
