-- Gkyl ------------------------------------------------------------------------
--
-- Test for MPI wrappers
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local Alloc = require "Lib.Alloc"

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

function test_0(comm)
   local grp = Mpi.Comm_group(comm)
   
   Mpi.Barrier(comm)   
end

-- Comm and Groups
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

   -- test communicators
   local myComm = Mpi.Comm_dup(comm)
   assert_equal(true, ffi.istype(typeof("MPI_Comm[1]"), myComm), "Checking new comm type")
   local mySz = Mpi.Comm_size(myComm)
   assert_equal(sz, mySz, "Checking if duplicated communicator has same number of elements")

   -- test groups
   local grp = Mpi.Comm_group(myComm)
   local gSz = Mpi.Group_size(grp)
   assert_equal(sz, gSz, "Checking size of group")
   local gRank = Mpi.Group_rank(grp)
   assert_equal(rank, gRank, "Checking rank of group")

   local gComm = Mpi.Comm_create(myComm, grp)
   assert_equal(sz, Mpi.Comm_size(gComm), "Checking if new communicator has correct size")

   Mpi.Barrier(comm)
end

-- Send/Recv
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

-- Send/Recv
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

-- Split_comm
function test_4(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 4 then
      log("Test for Split_comm not run as number of procs not exactly 4")
      return
   end

   local ranks = Lin.IntVec(2)
   ranks[1] = 0; ranks[2] = 1;
   local newComm = Mpi.Split_comm(comm, ranks)

   if rnk < 2 then
      local newSz = Mpi.Comm_size(newComm)
      assert_equal(2, newSz, "Testing new comminucator size")
      assert_equal(true, Mpi.Is_comm_valid(newComm), "Checking if comm is valid")
   else
      assert_equal(false, Mpi.Is_comm_valid(newComm), "Checking if comm is valid")
   end
   Mpi.Barrier(comm)
end

-- MPI_Bcast tests
function test_5(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz < 2 then
      log("Test for MPI_Bcast not run as number of procs not 2 or more")
      return
   end
   
   local vals = new("double [?]", 10)
   if rnk == 1 then
      for i = 0, 9 do
	 vals[i] = i+0.5
      end
   end
   -- send from rank 1 to all ranks
   Mpi.Bcast(vals, 10, Mpi.DOUBLE, 1, comm)
   -- check 
   for i = 0, 9 do
      assert_equal(i+0.5, vals[i], "Testing Bcast data")
   end
   Mpi.Barrier(comm)
end

-- Non-blocking recv
function test_6(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for non-blocking calls not run as number of procs not exactly 2")
      return
   end
   
   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
      
   local nz = 1000
   local vIn, vOut = Alloc.Double(nz), Alloc.Double(nz)
   for i = 1, nz do
      vIn[i] = (rank+1)*(i+0.5)
   end

   -- Send message from rank 0 -> rank 1
   if rank == 0 then
      Mpi.Send(vIn:data(), nz, Mpi.DOUBLE, 1, 42, comm)
   end   
   if rank == 1 then
      local request = Mpi.Irecv(vOut:data(), nz, Mpi.DOUBLE, 0, 42, comm)
      local status = Mpi.Status()
      Mpi.Wait(request, status)
      local count = Mpi.Get_count(status, Mpi.DOUBLE)

      assert_equal(42, status.TAG, "Checking tag")
      assert_equal(nz, count, "Checking if correct number of elements were recv-ed")
      for i = 1, nz do
	 assert_equal(i+0.5, vOut[i], "Checking recv-ed data on rank 1")
      end
   end

   -- Send message from rank 1 -> rank 0
   if rank == 1 then
      Mpi.Send(vIn:data(), nz, Mpi.DOUBLE, 0, 42, comm)
   end      
   if rank == 0 then
      local request = Mpi.Irecv(vOut:data(), nz, Mpi.DOUBLE, 1, 42, comm)
      local status = Mpi.Status()
      Mpi.Wait(request, status)
      local count = Mpi.Get_count(status, Mpi.DOUBLE)

      assert_equal(42, status.TAG, "Checking tag")
      assert_equal(nz, count, "Checking if correct number of elements were recv-ed")
      for i = 1, nz do
	 assert_equal(2*(i+0.5), vOut[i], "Checking recv-ed data on rank 0")
      end      
   end
   
   Mpi.Barrier(comm)
end

-- SHM
function test_7(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 4 then
      log("Test for SHM calls not run as number of procs not exactly 4")
      return
   end

   local shmComm = Mpi.Comm_split_type(comm, Mpi.COMM_TYPE_SHARED, 0, Mpi.INFO_NULL)
   assert_equal(4, Mpi.Comm_size(shmComm), "Checking number of ranks")
end   

-- Run tests
test_0(Mpi.COMM_WORLD)
test_1(Mpi.COMM_WORLD)
test_2(Mpi.COMM_WORLD)
test_3(Mpi.COMM_WORLD)
test_4(Mpi.COMM_WORLD)
test_5(Mpi.COMM_WORLD)
test_6(Mpi.COMM_WORLD)
test_7(Mpi.COMM_WORLD)

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

