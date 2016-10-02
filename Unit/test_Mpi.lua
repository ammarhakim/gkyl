-- Gkyl ------------------------------------------------------------------------
--
-- Test for MPI wrappers
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local DecompRegionCalc = require "Lib.CartDecomp"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"

local ffi  = require "ffi"
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

function test_1(comm)
   local sz = Mpi.Comm_size(comm)
   local rank = Mpi.Comm_rank(comm)
   local sendbuf = new("int[1]")
   sendbuf[0] = rank

   local recvbuf = new("int[1]")
   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.SUM, comm)
   assert_equal(recvbuf[0], sz*(sz-1)/2, "Checking allReduce sum")

   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.MAX, comm)
   assert_equal(recvbuf[0], sz-1, "Checking allReduce max")

   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.MIN, comm)
   assert_equal(recvbuf[0], 0, "Checking allReduce min")

   Mpi.Barrier(comm)
end

function test_2(comm)
   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
      
   local nz = 100
   local vIn, vOut = new("double[?]", nz), new("double[?]", nz)

   if rank == 0 then
      -- send data from rank 0 to all other ranks
      for i = 0, nz-1 do
	 vIn[i] = i
      end
      for dest = 1, sz-1 do
	 Mpi.Send(vIn, nz, Mpi.DOUBLE, dest, 22, comm)
      end      
   else
      -- recv stuff from rank 0
      Mpi.Recv(vOut, nz, Mpi.DOUBLE, 0, 22, comm, nil)
      for i = 0, nz-1 do
	 assert_equal(vOut[i], i, "Checking recv data")
      end
   end
   Mpi.Barrier(comm)
end

function test_3(comm)
   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
   local tag = 42
   local nz = 100
   local vIn, vOut = new("double[?]", nz), new("double[?]", nz)
   local status = Mpi.Status()

   if rank == 0 then
      for i = 0, nz-1 do vIn[i] = i end
      -- send data from rank 0 to all other ranks      
      for dest = 1, sz-1 do
	 Mpi.Send(vIn, nz, Mpi.DOUBLE, dest, tag, comm)
      end
   else
      -- recv stuff from rank 0
      Mpi.Recv(vOut, 100, Mpi.DOUBLE, 0, tag, comm, status)
      local count = Mpi.Get_count(status, Mpi.DOUBLE)
      assert_equal(count, 100, "Checking if correct number of elements were recv-ed")
      for i = 0, nz-1 do
	 assert_equal(vOut[i], i, "Checking recv data")
      end
   end
   Mpi.Barrier(comm)
end

-- Run tests
test_1(Mpi.COMM_WORLD)
test_2(Mpi.COMM_WORLD)
test_3(Mpi.COMM_WORLD)

function allReduceOneInt(localv)
   local sendbuf, recvbuf = new("int[1]"), new("int[1]")
   sendbuf[0] = localv
   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[0]
end

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))   
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end

