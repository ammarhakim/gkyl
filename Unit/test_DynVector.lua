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
   dynVec:write("test_1.bp", 1.5, 1)
   for i = 6, 10 do
      dynVec:appendData(0.1*i, {2.5*i^2, 2.5*i^2+0.5})
   end
   dynVec:write("test_2.bp", 1.5, 2)
end

function test_2r()
   local dynVec = DataStruct.DynVector { numComponents = 2 }
   assert_equal(dynVec:numComponents(), 2, "Testing number of components")

   local tmStamp, frNum = dynVec:read("test_1.bp")
   assert_equal(1.5, tmStamp, "Testing time-stamp")
   assert_equal(1, frNum, "Testing frame")

   -- check contents
   assert_equal(5, dynVec:size(), "Checking size")

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

test_1()
test_2()
test_2r()
test_3()

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))   
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
