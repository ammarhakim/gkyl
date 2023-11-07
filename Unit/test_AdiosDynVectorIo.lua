-- Gkyl ------------------------------------------------------------------------
--
-- Test for DynVector ADIOS IO.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit             = require "Unit"
local Mpi              = require "Comm.Mpi"
local Adios            = require "Io.Adios"
local AdiosDynVectorIo = require "Io.AdiosDynVectorIo"
local DataStruct       = require "DataStruct"
local Lin              = require "Lib.Linalg"

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function log(msg)
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   if rank == 0 then print(msg) end
end

function test_1w(comm)
   -- Test writing a single DynVector in one file.
   GKYL_ADIOS2_MPI = GKYL_ADIOS2_MPI or Adios.init_mpi(comm)

   local dynVec = DataStruct.DynVector {
      numComponents = 2,
      metaData      = {testNumber = 2}
   }
   assert_equal(dynVec:numComponents(), 2, "test_1w: Testing number of components")

   for i = 1, 5 do
      dynVec:appendData(0.1*i, {2.5*i^2, 2.5*i^2+0.5})
   end
   assert_equal(5, dynVec:size(), "test_1w: Checking size")
   dynVec:write("test_1w_0.bp", 1.5, "")  -- The "" here is so the variables in the file don't have 0 in the name.

   if GKYL_ADIOS2_MPI then Adios.finalize(GKYL_ADIOS2_MPI);  GKYL_ADIOS2_MPI = nil end
end

function test_1r(comm)
   -- Test reading a single DynVector in one file.
   GKYL_ADIOS2_MPI = GKYL_ADIOS2_MPI or Adios.init_mpi(comm)

   local dynVec = DataStruct.DynVector {
      numComponents = 2,
      metaData      = {testNumber = 2}
   }
   assert_equal(dynVec:numComponents(), 2, "test_1w: Testing number of components")

   local tmStamp, frNum = dynVec:read("test_1w_0.bp")

   -- check contents
   assert_equal(5, dynVec:size(), "test_1r: Checking size")

   local tmMesh = dynVec:timeMesh()
   for i = 1, dynVec:size() do
      assert_equal(0.1*i, tmMesh[i], "test_1r: Checking time-mesh")
   end

   local dynData = dynVec:data()
   for i = 1, dynVec:size() do
      local v = dynData[i]
      assert_equal(2.5*i^2, v[1], "test_1r: Checking contents[1]")
      assert_equal(2.5*i^2+0.5, v[2], "test_1r: Checking contents[2]")
   end

   if GKYL_ADIOS2_MPI then Adios.finalize(GKYL_ADIOS2_MPI);  GKYL_ADIOS2_MPI = nil end
end

function test_2w(comm)
   -- Test writing two separate DynVectors in two separate files.
   GKYL_ADIOS2_MPI = GKYL_ADIOS2_MPI or Adios.init_mpi(comm)

   local dynVecA = DataStruct.DynVector {
      numComponents = 2,
      metaData      = {testNumber = 2}
   }
   assert_equal(dynVecA:numComponents(), 2, "test_2w: Testing number of components")

   for i = 1, 5 do
      dynVecA:appendData(0.1*i, {2.5*i^2, 2.5*i^2+0.5})
   end
   assert_equal(5, dynVecA:size(), "test_2w: Checking size")
   dynVecA:write("test_2w_0.bp", 1.5, 0)

   local dynVecB = DataStruct.DynVector {
      numComponents = 3,
      metaData      = {testNumber = 2}
   }
   for i = 6, 11 do
      dynVecB:appendData(0.1*i, {i, 2.5*i^2, 2.5*i^2+0.5})
   end
   assert_equal(6, dynVecB:size(), "test_2w: Checking size")
   dynVecB:write("test_2w_1.bp", 1.5, 1)

   if GKYL_ADIOS2_MPI then Adios.finalize(GKYL_ADIOS2_MPI);  GKYL_ADIOS2_MPI = nil end
end

function test_2r(comm)
   -- Test reading two separate DynVectors from two separate files.
   GKYL_ADIOS2_MPI = GKYL_ADIOS2_MPI or Adios.init_mpi(comm)

   local dynVecA = DataStruct.DynVector {
      numComponents = 2,
      metaData      = {testNumber = 2}
   }

   local tmStamp, frNum = dynVecA:read("test_2w_0.bp")

   -- check contents
   assert_equal(5, dynVecA:size(), "test_2r: Checking A size")

   local tmMesh = dynVecA:timeMesh()
   for i = 1, dynVecA:size() do
      assert_equal(0.1*i, tmMesh[i], "test_2r: Checking A time-mesh")
   end

   local dynData = dynVecA:data()
   for i = 1, dynVecA:size() do
      local v = dynData[i]
      assert_equal(2.5*i^2, v[1], "test_2r: Checking A contents[1]")
      assert_equal(2.5*i^2+0.5, v[2], "test_2r: Checking A contents[2]")
   end

   -- Now read and check the other DynVector.
   local dynVecB = DataStruct.DynVector {
      numComponents = 3,
      metaData      = {testNumber = 2}
   }

   local tmStamp, frNum = dynVecB:read("test_2w_1.bp")

   -- check contents
   assert_equal(6, dynVecB:size(), "test_2r: Checking B size")

   local tmMesh = dynVecB:timeMesh()
   for i = 1, dynVecB:size() do
      assert_equal(0.1*(i+5), tmMesh[i], "test_2r: Checking B time-mesh")
   end

   local dynData = dynVecB:data()
   for i = 1, dynVecB:size() do
      local v = dynData[i]
      assert_equal((i+5), v[1], "test_2r: Checking B contents[1]")
      assert_equal(2.5*(i+5)^2, v[2], "test_2r: Checking B contents[2]")
      assert_equal(2.5*(i+5)^2+0.5, v[3], "test_2r: Checking B contents[3]")
   end

   if GKYL_ADIOS2_MPI then Adios.finalize(GKYL_ADIOS2_MPI);  GKYL_ADIOS2_MPI = nil end
end

function test_3w(comm)
   -- Test writing a DynVector to the same file multiple times.
   GKYL_ADIOS2_MPI = GKYL_ADIOS2_MPI or Adios.init_mpi(comm)

   local dynVec = DataStruct.DynVector {
      numComponents = 2,
      metaData      = {testNumber = 2}
   }
   assert_equal(dynVec:numComponents(), 2, "test_3w: Testing number of components")

   for i = 1, 5 do
      dynVec:appendData(0.1*i, {2.5*i^2, 2.5*i^2+0.5})
   end
   assert_equal(5, dynVec:size(), "test_3w: Checking size")
   dynVec:write("test_3w_0.bp", 1.5, 0)  -- The "" here is so the variables in the file don't have 0 in the name.

   for i = 6, 9 do
      dynVec:appendData(0.1*i, {1.5*i^2, 1.5*i^2+0.5})
   end
   assert_equal(4, dynVec:size(), "test_3w: Checking size")
   dynVec:write("test_3w_0.bp", 3., 1)  -- The "" here is so the variables in the file don't have 0 in the name.

   for i = 10, 12 do
      dynVec:appendData(0.1*i, {0.5*i^2, 0.5*i^2+0.5})
   end
   assert_equal(3, dynVec:size(), "test_3w: Checking size")
   dynVec:write("test_3w_0.bp", 4.5, 2)  -- The "" here is so the variables in the file don't have 0 in the name.

   if GKYL_ADIOS2_MPI then Adios.finalize(GKYL_ADIOS2_MPI);  GKYL_ADIOS2_MPI = nil end
end

function test_3r(comm)
   -- Test reading a single DynVector written to one file multiple times.
   GKYL_ADIOS2_MPI = GKYL_ADIOS2_MPI or Adios.init_mpi(comm)

   local dynVec = DataStruct.DynVector {
      numComponents = 2,
      metaData      = {testNumber = 2}
   }

   local tmStamp, frNum = dynVec:read("test_3w_0.bp")

   -- check contents
   assert_equal(5+4+3, dynVec:size(), "test_3r: Checking size")

   local tmMesh = dynVec:timeMesh()
   for i = 1, dynVec:size() do
      assert_equal(0.1*i, tmMesh[i], "test_3r: Checking time-mesh")
   end

   local dynData = dynVec:data()
   for i = 1, dynVec:size() do
      local v = dynData[i]
      if i < 6 then
         assert_equal(2.5*i^2, v[1], "test_3r: Checking contents[1]")
         assert_equal(2.5*i^2+0.5, v[2], "test_3r: Checking contents[2]")
      elseif i < 10 then
         assert_equal(1.5*i^2, v[1], "test_3r: Checking contents[1]")
         assert_equal(1.5*i^2+0.5, v[2], "test_3r: Checking contents[2]")
      else
         assert_equal(0.5*i^2, v[1], "test_3r: Checking contents[1]")
         assert_equal(0.5*i^2+0.5, v[2], "test_3r: Checking contents[2]")
      end
   end

   if GKYL_ADIOS2_MPI then Adios.finalize(GKYL_ADIOS2_MPI);  GKYL_ADIOS2_MPI = nil end
end

function allReduceOneInt(localv)
   local sendbuf, recvbuf = Lin.IntVec(1), Lin.IntVec(1)
   sendbuf[1] = localv
   Mpi.Allreduce(sendbuf:data(), recvbuf:data(), 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[1]
end

-- Run tests
test_1w(Mpi.COMM_WORLD)
test_1r(Mpi.COMM_WORLD)
test_2w(Mpi.COMM_WORLD)
test_2r(Mpi.COMM_WORLD)
test_3w(Mpi.COMM_WORLD)
test_3r(Mpi.COMM_WORLD)

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
