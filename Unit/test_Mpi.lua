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
local Range = require "Lib.Range"

local cuda
local cuAlloc
if GKYL_HAVE_CUDA then
   cuda = require "Cuda.RunTime"
   cuAlloc = require "Cuda.Alloc"
end

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
   assert_equal(true, Mpi.Is_comm_valid(comm))
   
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

   local ranks = Lin.IntVec(2) --creating an int vector with 2 elements
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

   -- Test broadcasting strings.
   -- MF 2021/05/05: as currently implemented it should probably only be used for strings
   --                with length>0. It occasionally seg faults with empty strings.
   -- MF 2021/05/25: below we use len+1 because without the +1 I was seeing random failures.
   -- MF 2021/06/17: we have removed the +1 and +2 altogether. Lua strings are not NULL terminated.
   -- MF 2021/06/22: even without the +1,+2, a seg fault was encountered in Stellar but not Frontera. Conclusion: we don't yet know how to broadcast strings reliably.
   local myStr = "myRank".. rnk
   local Cstr = new("char [?]", string.len(myStr))
   ffi.copy(Cstr, myStr)
   Mpi.Bcast(Cstr, string.len(myStr), Mpi.CHAR, 1, comm)
   myStr = ffi.string(Cstr)
   assert_equal("myRank1", myStr, "Testing Bcast string")
end

-- Non-blocking recv
function test_6a(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for non-blocking recv not run as number of procs not exactly 2")
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
      local request = Mpi.Request()
      local err = Mpi.Irecv(vOut:data(), nz, Mpi.DOUBLE, 0, 42, comm, request)
      local status = Mpi.Status()
      Mpi.Wait(request, status)
      local count = Mpi.Get_count(status, Mpi.DOUBLE)

      assert_equal(42, Mpi.getStatusTAG(status), "Checking tag")
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
      local request = Mpi.Request()
      local err = Mpi.Irecv(vOut:data(), nz, Mpi.DOUBLE, 1, 42, comm, request)
      local status = Mpi.Status()
      Mpi.Wait(request, status)
      local count = Mpi.Get_count(status, Mpi.DOUBLE)

      assert_equal(42, Mpi.getStatusTAG(status), "Checking tag")
      assert_equal(nz, count, "Checking if correct number of elements were recv-ed")
      for i = 1, nz do
	 assert_equal(2*(i+0.5), vOut[i], "Checking recv-ed data on rank 0")
      end      
   end
   
   Mpi.Barrier(comm)
end

-- Non-blocking send and recv
function test_6b(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for non-blocking send/recv not run as number of procs not exactly 2")
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
      local request = Mpi.Request()
      local err = Mpi.Isend(vIn:data(), nz, Mpi.DOUBLE, 1, 42, comm, request)
      local status = Mpi.Status()
      Mpi.Wait(request, status)
      local count = Mpi.Get_count(status, Mpi.DOUBLE)
      assert_equal(42, Mpi.getStatusTAG(status), "rank 0: Checking tag")
      assert_equal(nz, count, "rank 0: Checking if correct number of elements were sent")
   end   
   if rank == 1 then
      local request = Mpi.Request()
      local err = Mpi.Irecv(vOut:data(), nz, Mpi.DOUBLE, 0, 42, comm, request)
      local status = Mpi.Status()
      Mpi.Wait(request, status)
      local count = Mpi.Get_count(status, Mpi.DOUBLE)

      assert_equal(42, Mpi.getStatusTAG(status), "rank 1: Checking tag")
      assert_equal(nz, count, "rank 1: Checking if correct number of elements were recv-ed")
      for i = 1, nz do
	 assert_equal(i+0.5, vOut[i], "rank 1: Checking recv-ed data")
      end
   end

   -- Send message from rank 1 -> rank 0
   if rank == 1 then
      local request = Mpi.Request()
      local err = Mpi.Isend(vIn:data(), nz, Mpi.DOUBLE, 0, 42, comm, request)
      local status = Mpi.Status()
      Mpi.Wait(request, status)
      local count = Mpi.Get_count(status, Mpi.DOUBLE)

      assert_equal(42, Mpi.getStatusTAG(status), "rank 1: Checking tag")
      assert_equal(nz, count, "rank 1: Checking if correct number of elements were sent")
   end      
   if rank == 0 then
      local request = Mpi.Request()
      local err = Mpi.Irecv(vOut:data(), nz, Mpi.DOUBLE, 1, 42, comm, request)
      local status = Mpi.Status()
      Mpi.Wait(request, status)
      local count = Mpi.Get_count(status, Mpi.DOUBLE)

      assert_equal(42, Mpi.getStatusTAG(status), "rank 0: Checking tag")
      assert_equal(nz, count, "rank 0: Checking if correct number of elements were recv-ed")
      for i = 1, nz do
	 assert_equal(2*(i+0.5), vOut[i], "rank 0: Checking recv-ed data")
      end      
   end
   
   Mpi.Barrier(comm)
end

-- Waitall with non-blocking send and recv
function test_6c(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for non-blocking send/recv not run as number of procs not exactly 2")
      return
   end
   
   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
      
   local nz0, nz1 = 1000, 600
   local vIn0, vOut0 = Alloc.Double(nz0), Alloc.Double(nz0)
   local vIn1, vOut1 = Alloc.Double(nz1), Alloc.Double(nz1)
   for i = 1, nz0 do
      vIn0[i] = (rank+1)*(i+0.5)
   end
   for i = 1, nz1 do
      vIn1[i] = (rank+2)*(i+0.1)
   end

   -- Send message from rank 0 -> rank 1
   local requests = Mpi.Request(2)
   if rank == 0 then
      local err
      err = Mpi.Isend(vIn0:data(), nz0, Mpi.DOUBLE, 1, 32, comm, requests.req+0)
      err = Mpi.Isend(vIn1:data(), nz1, Mpi.DOUBLE, 1, 42, comm, requests.req+1)
   end   
   if rank == 1 then
      local err
      err = Mpi.Irecv(vOut0:data(), nz0, Mpi.DOUBLE, 0, 32, comm, requests.req+0)
      err = Mpi.Irecv(vOut1:data(), nz1, Mpi.DOUBLE, 0, 42, comm, requests.req+1)
   end
   local status = Mpi.Status(2)
   local mpierr = Mpi.Waitall(2, requests, status)

   local count = {Mpi.Get_count(status, Mpi.DOUBLE, 0),
                  Mpi.Get_count(status, Mpi.DOUBLE, 1)}
   assert_equal(32, Mpi.getStatusTAG(status,0), "Checking tag")
   assert_equal(42, Mpi.getStatusTAG(status,1), "Checking tag")
   assert_equal(nz0, count[1], "Checking if correct number of elements were sent/received")
   assert_equal(nz1, count[2], "Checking if correct number of elements were sent/received")
   if rank == 1 then
      for i = 1, nz0 do
	 assert_equal((0+1)*(i+0.5), vOut0[i], "rank 1: Checking recv-ed data0")
      end
      for i = 1, nz1 do
	 assert_equal((0+2)*(i+0.1), vOut1[i], "rank 1: Checking recv-ed data1")
      end
   end

   -- Send message from rank 1 -> rank 0
   if rank == 1 then
      local err
      err = Mpi.Isend(vIn0:data(), nz0, Mpi.DOUBLE, 0, 32, comm, requests.req+0)
      err = Mpi.Isend(vIn1:data(), nz1, Mpi.DOUBLE, 0, 42, comm, requests.req+1)
   end      
   if rank == 0 then
      local err
      err = Mpi.Irecv(vOut0:data(), nz0, Mpi.DOUBLE, 1, 32, comm, requests.req+0)
      err = Mpi.Irecv(vOut1:data(), nz1, Mpi.DOUBLE, 1, 42, comm, requests.req+1)
   end
   local mpierr = Mpi.Waitall(2, requests, status)

   local count = {Mpi.Get_count(status, Mpi.DOUBLE, 0),
                  Mpi.Get_count(status, Mpi.DOUBLE, 1)}
   assert_equal(32, Mpi.getStatusTAG(status,0), "Checking tag")
   assert_equal(42, Mpi.getStatusTAG(status,1), "Checking tag")
   assert_equal(nz0, count[1], "Checking if correct number of elements were sent/received")
   assert_equal(nz1, count[2], "Checking if correct number of elements were sent/received")
   if rank == 0 then
      for i = 1, nz0 do
	 assert_equal((1+1)*(i+0.5), vOut0[i], "rank 0: Checking recv-ed data0")
      end
      for i = 1, nz1 do
	 assert_equal((1+2)*(i+0.1), vOut1[i], "rank 1: Checking recv-ed data1")
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

   Mpi.Barrier(comm)
end

-- Group_translate_ranks
function test_8(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 4 then
      log("Test for Group_translate_ranks not run as number of procs not exactly 4")
      return
   end

   local shmComm = Mpi.Comm_split_type(comm, Mpi.COMM_TYPE_SHARED, 0, Mpi.INFO_NULL)

   local worldGrp = Mpi.Comm_group(comm)
   local shmGrp = Mpi.Comm_group(shmComm)

   local worldRanks = Lin.IntVec(sz)
   for d = 1, sz do worldRanks[d] = d-1 end -- ranks are indexed from 0
   local shmRanks = Mpi.Group_translate_ranks(worldGrp, worldRanks, shmGrp)

   for d = 1, #shmRanks do
      assert_equal(worldRanks[d], shmRanks[d], "Checking rank mapping")
   end

   Mpi.Barrier(comm)
end

-- Datatypes
function test_9(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)

   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for MPI_Datatype not run as number of procs not exactly 2")
      return
   end   

   local nz, nsend = 100, 20
   local buff = Alloc.Double(nz)

   -- '"datatype" for send/recv. This is not really a data-type but a
   -- shape that tells what portion of buffers to copy. The
   -- "contiguous" type is simplest, allowing one to copy a series of 
   -- memory locations.
   local dType = Mpi.Type_contiguous(nsend, Mpi.DOUBLE)
   Mpi.Type_commit(dType)

   if rnk == 0 then
      for i = 1, nz do buff[i] = i+0.5 end
      Mpi.Send(buff:data(), 1, dType, 1, 42, comm)
   end
   if rnk == 1 then
      Mpi.Recv(buff:data(), 1, dType, 0, 42, comm, nil)
      for i = 1, nsend do
	 assert_equal(buff[i], i+0.5, "Checking send/recv via MPI_Datatype")
      end
      for i = nsend+1, nz do
	 assert_equal(0.0, buff[i], "Checking send/recv via MPI_Datatype")
      end
   end

   if rnk == 0 then
      for i = 1, nz do buff[i] = i+10.5 end
      Mpi.Send(buff:data()+5, 1, dType, 1, 42, comm)
   end
   if rnk == 1 then
      Mpi.Recv(buff:data(), 1, dType, 0, 42, comm, nil)
      for i = 1, nsend do
   	 assert_equal(i+10.5+5, buff[i], "Checking send/recv via MPI_Datatype")
      end
      for i = nsend+1, nz do
   	 assert_equal(0.0, buff[i], "Checking send/recv via MPI_Datatype")
      end      
   end   
   Mpi.Type_free(dType)
end

-- Datatypes (2D test)
function test_10(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)

   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for MPI_Datatype (test_10) not run as number of procs not exactly 2")
      return
   end

   local range = Range.Range({1, 1}, {10, 20})
   local indexer = Range.makeRowMajorIndexer(range)
   
   local nz = range:volume()
   local buff = Alloc.Double(nz)

   -- Construct '"datatype" for send/recv. This is not really a
   -- data-type but a shape that tells what portion of buffers to
   -- copy.
   local dTypeX, dTypeY = nil, nil

   local distance = indexer(range:upper(1), range:lower(2))
      - indexer(range:lower(1), range:lower(2))
   if distance == range:upper(1)-range:lower(1) then
      dTypeX = Mpi.Type_vector(range:shape(2), 1, range:shape(1), Mpi.DOUBLE)
      dTypeY = Mpi.Type_contiguous(range:shape(1), Mpi.DOUBLE)
   else
      dTypeX = Mpi.Type_contiguous(range:shape(2), Mpi.DOUBLE) 
      dTypeY = Mpi.Type_vector(range:shape(1), 1, range:shape(2), Mpi.DOUBLE)
   end
   Mpi.Type_commit(dTypeX)
   Mpi.Type_commit(dTypeY)

   if rnk == 0 then
      for i = range:lower(1), range:upper(1) do
	 for j = range:lower(2), range:upper(2) do
	    buff[indexer(i,j)] = i+20*j+0.5
	 end
      end
   end

   -- X direction
   if rnk == 0 then
      -- left edge
      Mpi.Send(buff:data(), 1, dTypeX, 1, 42, comm)
      -- right edge
      local ix, iy = range:upper(1), range:lower(2)
      Mpi.Send(buff:data()+indexer(ix,iy)-1, 1, dTypeX, 1, 43, comm)
   else
      -- left edge
      Mpi.Recv(buff:data(), 1, dTypeX, 0, 42, comm, nil)
      -- right edge
      local ix, iy = range:upper(1), range:lower(2)
      Mpi.Recv(buff:data()+indexer(ix,iy)-1, 1, dTypeX, 0, 43, comm, nil)

      local i = range:lower(1)
      for j = range:lower(2), range:upper(2) do
	 assert_equal(i+20*j+0.5, buff[indexer(i,j)], "Checking X-direction send/recv")
      end

      i = range:upper(1)
      for j = range:lower(2), range:upper(2) do
	 assert_equal(i+20*j+0.5, buff[indexer(i,j)], "Checking X-direction send/recv")
      end      
   end

   -- Y direction
   if rnk == 0 then
      -- bottom edge
      Mpi.Send(buff:data(), 1, dTypeY, 1, 44, comm)
      -- top edge
      local ix, iy = range:lower(1), range:upper(2)
      Mpi.Send(buff:data()+indexer(ix,iy)-1, 1, dTypeY, 1, 45, comm)
   else
      -- bottom edge
      Mpi.Recv(buff:data(), 1, dTypeY, 0, 44, comm, nil)
      -- top edge
      local ix, iy = range:lower(1), range:upper(2)
      Mpi.Recv(buff:data()+indexer(ix,iy)-1, 1, dTypeY, 0, 45, comm, nil)

      local j = range:lower(2)
      for i = range:lower(1), range:upper(1) do
	 assert_equal(i+20*j+0.5, buff[indexer(i,j)], "Checking Y-direction send/recv")
      end

      j = range:upper(2)
      for i = range:lower(1), range:upper(1) do
	 assert_equal(i+20*j+0.5, buff[indexer(i,j)], "Checking Y-direction send/recv")
      end      
   end   
end

-- Datatypes (2D test with Type_indexed)
function test_11(comm)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)

   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for MPI_Datatype (test_11) not run as number of procs not exactly 2")
      return
   end

   local range = Range.Range({1, 1}, {10, 20})
   local indexer = Range.makeRowMajorIndexer(range)
   
   local nz = range:volume()
   local buff = Alloc.Double(nz)

   -- Construct '"datatype" for send/recv. This is not really a
   -- data-type but a shape that tells what portion of buffers to
   -- copy.
   local dTypeX, dTypeY = nil, nil

   local distance = indexer(range:upper(1), range:lower(2))
      - indexer(range:lower(1), range:lower(2))
   if distance == range:upper(1)-range:lower(1) then
      -- i-indices are ctt
      local xbl, xdp = Lib.IntVec(range:shape(2)), Lib.IntVec(range:shape(2))
      for i = 1, #xbl do xbl[i] = 1 end
      xdp[1] = 0
      for i = 2, #xdp do xdp[i] = xdp[i-1]+range:shape(1) end
      dTypeX = Mpi.Type_indexed(xbl, xdp, Mpi.DOUBLE)

      local ybl, ydp = Lin.IntVec(1), Lin.IntVec(1)
      ybl[1] = range:shape(1)
      ydp[1] = 0
      dTypeY = Mpi.Type_indexed(ybl, ydp, Mpi.DOUBLE)
   else
      -- j-indices are ctt
      local xbl, xdp = Lin.IntVec(1), Lin.IntVec(1)
      xbl[1] = range:shape(2)
      xdp[1] = 0
      dTypeX = Mpi.Type_indexed(xbl, xdp, Mpi.DOUBLE)

      local ybl, ydp = Lin.IntVec(range:shape(1)), Lin.IntVec(range:shape(1))
      for i = 1, #ybl do ybl[i] = 1 end
      ydp[1] = 0
      for i = 2, #ydp do ydp[i] = ydp[i-1]+range:shape(2) end
      dTypeY = Mpi.Type_indexed(ybl, ydp, Mpi.DOUBLE)
   end
   Mpi.Type_commit(dTypeX)
   Mpi.Type_commit(dTypeY)

   if rnk == 0 then
      for i = range:lower(1), range:upper(1) do
	 for j = range:lower(2), range:upper(2) do
	    buff[indexer(i,j)] = i+20*j+0.5
	 end
      end
   end

   -- X direction
   if rnk == 0 then
      -- left edge
      Mpi.Send(buff:data(), 1, dTypeX, 1, 42, comm)
      -- right edge
      local ix, iy = range:upper(1), range:lower(2)
      Mpi.Send(buff:data()+indexer(ix,iy)-1, 1, dTypeX, 1, 43, comm)
   else
      -- left edge
      Mpi.Recv(buff:data(), 1, dTypeX, 0, 42, comm, nil)
      -- right edge
      local ix, iy = range:upper(1), range:lower(2)
      Mpi.Recv(buff:data()+indexer(ix,iy)-1, 1, dTypeX, 0, 43, comm, nil)

      local i = range:lower(1)
      for j = range:lower(2), range:upper(2) do
	 assert_equal(i+20*j+0.5, buff[indexer(i,j)], "Checking X-direction send/recv")
      end

      i = range:upper(1)
      for j = range:lower(2), range:upper(2) do
	 assert_equal(i+20*j+0.5, buff[indexer(i,j)], "Checking X-direction send/recv")
      end      
   end

   -- Y direction
   if rnk == 0 then
      -- bottom edge
      Mpi.Send(buff:data(), 1, dTypeY, 1, 44, comm)
      -- top edge
      local ix, iy = range:lower(1), range:upper(2)
      Mpi.Send(buff:data()+indexer(ix,iy)-1, 1, dTypeY, 1, 45, comm)
   else
      -- bottom edge
      Mpi.Recv(buff:data(), 1, dTypeY, 0, 44, comm, nil)
      -- top edge
      local ix, iy = range:lower(1), range:upper(2)
      Mpi.Recv(buff:data()+indexer(ix,iy)-1, 1, dTypeY, 0, 45, comm, nil)

      local j = range:lower(2)
      for i = range:lower(1), range:upper(1) do
	 assert_equal(i+20*j+0.5, buff[indexer(i,j)], "Checking Y-direction send/recv")
      end

      j = range:upper(2)
      for i = range:lower(1), range:upper(1) do
	 assert_equal(i+20*j+0.5, buff[indexer(i,j)], "Checking Y-direction send/recv")
      end      
   end   
end

function test_12(comm, nlayer, ordering)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)

   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for MPI_Datatype (test_12) not run as number of procs not exactly 2")
      return
   end

   local range = Range.Range({1, 1}, {10, 20})
   local dTypeX = Mpi.createDataTypeFromRange(1, range, nlayer, 1, ordering, Mpi.DOUBLE)
   local dTypeY = Mpi.createDataTypeFromRange(2, range, nlayer, 1, ordering, Mpi.DOUBLE)

   local indexer = range:indexer(ordering)
   local nz = range:volume()
   local buff = Alloc.Double(nz)   
   if rnk == 0 then
      for i = range:lower(1), range:upper(1) do
	 for j = range:lower(2), range:upper(2) do
	    buff[indexer(i,j)] = i+20*j+0.5
	 end
      end
   end

   -- X direction
   if rnk == 0 then
      -- left edge
      Mpi.Send(buff:data(), 1, dTypeX, 1, 42, comm)
      -- right edge
      local ix, iy = range:upper(1)-nlayer+1, range:lower(2)
      Mpi.Send(buff:data()+indexer(ix,iy)-1, 1, dTypeX, 1, 43, comm)
   else
      -- left edge
      Mpi.Recv(buff:data(), 1, dTypeX, 0, 42, comm, nil)
      -- right edge
      local ix, iy = range:upper(1)-nlayer+1, range:lower(2)
      Mpi.Recv(buff:data()+indexer(ix,iy)-1, 1, dTypeX, 0, 43, comm, nil)

      for n = 1, nlayer do
	 local i = range:lower(1)+n-1
	 for j = range:lower(2), range:upper(2) do
	    assert_equal(i+20*j+0.5, buff[indexer(i,j)],
			 string.format("Checking lower X-direction send/recv (%d,%d)", i,j))
	 end
      end

      for n = 1, nlayer do
      	 local i = range:upper(1)-n+1
      	 for j = range:lower(2), range:upper(2) do
      	    assert_equal(i+20*j+0.5, buff[indexer(i,j)],
      			 string.format("Checking upper X-direction send/recv (%d,%d)", i,j))
      	 end
      end
   end

   -- Y direction
   if rnk == 0 then
      -- bottom edge
      Mpi.Send(buff:data(), 1, dTypeY, 1, 44, comm)
      -- top edge
      local ix, iy = range:lower(1), range:upper(2)-nlayer+1
      Mpi.Send(buff:data()+indexer(ix,iy)-1, 1, dTypeY, 1, 45, comm)
   else
      -- bottom edge
      Mpi.Recv(buff:data(), 1, dTypeY, 0, 44, comm, nil)
      -- top edge
      local ix, iy = range:lower(1), range:upper(2)-nlayer+1
      Mpi.Recv(buff:data()+indexer(ix,iy)-1, 1, dTypeY, 0, 45, comm, nil)

      for n = 1, nlayer do
	 local j = range:lower(2)+n-1
	 for i = range:lower(1), range:upper(1) do
	    assert_equal(i+20*j+0.5, buff[indexer(i,j)], "Checking Y-direction send/recv")
	 end
      end

      for n = 1, nlayer do
	 local j = range:upper(2)-n+1
	 for i = range:lower(1), range:upper(1) do
	    assert_equal(i+20*j+0.5, buff[indexer(i,j)], "Checking Y-direction send/recv")
	 end
      end
   end
   
   Mpi.Barrier(comm)
end

function test_13(comm, nlayer, numComponents, ordering)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)

   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for MPI_Datatype (test_13) not run as number of procs not exactly 2")
      return
   end

   local range = Range.Range({1, 1}, {10, 20})
   local dTypeX = Mpi.createDataTypeFromRange(1, range, nlayer, numComponents, ordering, Mpi.DOUBLE)
   local dTypeY = Mpi.createDataTypeFromRange(2, range, nlayer, numComponents, ordering, Mpi.DOUBLE)

   local nz = range:volume()*numComponents
   local buff = Alloc.Double(nz)

   local indexer = range:indexer(ordering)   
   local function cidx(i,j,c)
      return (indexer(i,j)-1)*numComponents+c
   end
   
   if rnk == 0 then
      for i = range:lower(1), range:upper(1) do
   	 for j = range:lower(2), range:upper(2) do
   	    for k = 1, numComponents do
   	       buff[cidx(i,j,k)] = i+20*j+0.5+1000*k
   	    end
   	 end
      end
   end

   -- X direction
   if rnk == 0 then
      -- left edge
      Mpi.Send(buff:data(), 1, dTypeX, 1, 42, comm)
      -- right edge
      local ix, iy = range:upper(1)-nlayer+1, range:lower(2)
      Mpi.Send(buff:data()+cidx(ix,iy,1)-1, 1, dTypeX, 1, 43, comm)
   else
      -- left edge
      Mpi.Recv(buff:data(), 1, dTypeX, 0, 42, comm, nil)
      -- right edge
      local ix, iy = range:upper(1)-nlayer+1, range:lower(2)
      Mpi.Recv(buff:data()+cidx(ix,iy,1)-1, 1, dTypeX, 0, 43, comm, nil)

      for n = 1, nlayer do
   	 local i = range:lower(1)+n-1
   	 for j = range:lower(2), range:upper(2) do
   	    for k = 1, numComponents do
   	       assert_equal(i+20*j+0.5+1000*k, buff[cidx(i,j,k)],
   			    string.format("Checking lower X-direction send/recv (%d,%d;%d [%d])",
					  i,j,k,cidx(i,j,k)))
   	    end
   	 end
      end

      for n = 1, nlayer do
      	 local i = range:upper(1)-n+1
      	 for j = range:lower(2), range:upper(2) do
   	    for k = 1, numComponents do	    
   	       assert_equal(i+20*j+0.5+1000*k, buff[cidx(i,j,k)],
   			    string.format("Checking upper X-direction send/recv (%d,%d;%d [%d])",
					  i,j,k,cidx(i,j,k)))
   	    end
      	 end
      end
   end

   -- Y direction
   if rnk == 0 then
      -- bottom edge
      Mpi.Send(buff:data(), 1, dTypeY, 1, 44, comm)
      -- top edge
      local ix, iy = range:lower(1), range:upper(2)-nlayer+1
      Mpi.Send(buff:data()+cidx(ix,iy,1)-1, 1, dTypeY, 1, 45, comm)
   else
      -- bottom edge
      Mpi.Recv(buff:data(), 1, dTypeY, 0, 44, comm, nil)
      -- top edge
      local ix, iy = range:lower(1), range:upper(2)-nlayer+1
      Mpi.Recv(buff:data()+cidx(ix,iy,1)-1, 1, dTypeY, 0, 45, comm, nil)

      for n = 1, nlayer do
   	 local j = range:lower(2)+n-1
   	 for i = range:lower(1), range:upper(1) do
	    for k = 1, numComponents do
	       assert_equal(i+20*j+0.5+1000*k, buff[cidx(i,j,k)], "Checking Y-direction send/recv")
	    end
   	 end
      end

      for n = 1, nlayer do
   	 local j = range:upper(2)-n+1
   	 for i = range:lower(1), range:upper(1) do
	    for k = 1, numComponents do
	       assert_equal(i+20*j+0.5+1000*k, buff[cidx(i,j,k)], "Checking Y-direction send/recv")
	    end
   	 end
      end
   end
end

function test_14(comm, nlayer, numComponents, ordering)
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)

   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for MPI_Datatype (test_14) not run as number of procs not exactly 2")
      return
   end

   local range = Range.Range({1, 1}, {10, 20})
   local rangeX = range:shorten(1, nlayer)
   local rangeY = range:shorten(2, nlayer)
   
   local dTypeX = Mpi.createDataTypeFromRangeAndSubRange(rangeX, range, numComponents, ordering, Mpi.DOUBLE)
   local dTypeY = Mpi.createDataTypeFromRangeAndSubRange(rangeY, range, numComponents, ordering, Mpi.DOUBLE)

   local nz = range:volume()*numComponents
   local buff = Alloc.Double(nz)

   local indexer = range:indexer(ordering)   
   local function cidx(i,j,c)
      return (indexer(i,j)-1)*numComponents+c
   end
   
   if rnk == 0 then
      for i = range:lower(1), range:upper(1) do
   	 for j = range:lower(2), range:upper(2) do
   	    for k = 1, numComponents do
   	       buff[cidx(i,j,k)] = i+20*j+0.5+1000*k
   	    end
   	 end
      end
   end

   -- X direction
   if rnk == 0 then
      -- left edge
      Mpi.Send(buff:data(), 1, dTypeX, 1, 42, comm)
      -- right edge
      local ix, iy = range:upper(1)-nlayer+1, range:lower(2)
      Mpi.Send(buff:data()+cidx(ix,iy,1)-1, 1, dTypeX, 1, 43, comm)
   else
      -- left edge
      Mpi.Recv(buff:data(), 1, dTypeX, 0, 42, comm, nil)
      -- right edge
      local ix, iy = range:upper(1)-nlayer+1, range:lower(2)
      Mpi.Recv(buff:data()+cidx(ix,iy,1)-1, 1, dTypeX, 0, 43, comm, nil)

      for n = 1, nlayer do
   	 local i = range:lower(1)+n-1
   	 for j = range:lower(2), range:upper(2) do
   	    for k = 1, numComponents do
   	       assert_equal(i+20*j+0.5+1000*k, buff[cidx(i,j,k)],
   			    string.format("Checking lower X-direction send/recv (%d,%d;%d [%d])",
					  i,j,k,cidx(i,j,k)))
   	    end
   	 end
      end

      for n = 1, nlayer do
      	 local i = range:upper(1)-n+1
      	 for j = range:lower(2), range:upper(2) do
   	    for k = 1, numComponents do	    
   	       assert_equal(i+20*j+0.5+1000*k, buff[cidx(i,j,k)],
   			    string.format("Checking upper X-direction send/recv (%d,%d;%d [%d])",
					  i,j,k,cidx(i,j,k)))
   	    end
      	 end
      end
   end

   -- Y direction
   if rnk == 0 then
      -- bottom edge
      Mpi.Send(buff:data(), 1, dTypeY, 1, 44, comm)
      -- top edge
      local ix, iy = range:lower(1), range:upper(2)-nlayer+1
      Mpi.Send(buff:data()+cidx(ix,iy,1)-1, 1, dTypeY, 1, 45, comm)
   else
      -- bottom edge
      Mpi.Recv(buff:data(), 1, dTypeY, 0, 44, comm, nil)
      -- top edge
      local ix, iy = range:lower(1), range:upper(2)-nlayer+1
      Mpi.Recv(buff:data()+cidx(ix,iy,1)-1, 1, dTypeY, 0, 45, comm, nil)

      for n = 1, nlayer do
   	 local j = range:lower(2)+n-1
   	 for i = range:lower(1), range:upper(1) do
	    for k = 1, numComponents do
	       assert_equal(i+20*j+0.5+1000*k, buff[cidx(i,j,k)], "Checking Y-direction send/recv")
	    end
   	 end
      end

      for n = 1, nlayer do
   	 local j = range:upper(2)-n+1
   	 for i = range:lower(1), range:upper(1) do
	    for k = 1, numComponents do
	       assert_equal(i+20*j+0.5+1000*k, buff[cidx(i,j,k)], "Checking Y-direction send/recv")
	    end
   	 end
      end
   end   
end

function test_15(comm, nlayer, numComponents, ordering)
   -- tests MPI_Datatype API with a sub-region of a given region
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)

   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for MPI_Datatype (test_15) not run as number of procs not exactly 2")
      return
   end

   local range = Range.Range({1, 1}, {10, 20})
   local rangeX = Range.Range(
      {range:lower(1), range:lower(2)+nlayer},
      {range:lower(1), range:upper(2)-nlayer})
   local rangeY = Range.Range(
      {range:lower(1)+nlayer, range:lower(2)},
      {range:upper(1)-nlayer, range:lower(2)})
   
   local dTypeX = Mpi.createDataTypeFromRangeAndSubRange(rangeX, range, numComponents, ordering, Mpi.DOUBLE)
   --local dTypeY = Mpi.createDataTypeFromRangeAndSubRange(rangeY, range, numComponents, ordering, Mpi.DOUBLE)

   local nz = range:volume()*numComponents
   local buff = Alloc.Double(nz)

   local indexer = range:indexer(ordering)   
   local function cidx(i,j,c)
      return (indexer(i,j)-1)*numComponents+c
   end
   
   if rnk == 0 then
      for i = range:lower(1), range:upper(1) do
   	 for j = range:lower(2), range:upper(2) do
   	    for k = 1, numComponents do
   	       buff[cidx(i,j,k)] = i+20*j+0.5+1000*k
   	    end
   	 end
      end
   end

   -- X direction
   if rnk == 0 then
      -- send box
      local loc = cidx(rangeX:lower(1), rangeX:lower(2), 0)
      Mpi.Send(buff:data()+loc, 1, dTypeX, 1, 42, comm)
   else
      -- recv box
      local loc = cidx(rangeX:lower(1), rangeX:lower(2), 0)      
      Mpi.Recv(buff:data()+loc, 1, dTypeX, 0, 42, comm, nil)

      for i = rangeX:lower(1), rangeX:upper(1) do
	 for j = rangeX:lower(2), rangeX:upper(2) do
	    for k = 1, numComponents do
	       assert_equal(i+20*j+0.5+1000*k, buff[cidx(i,j,k)],
			    string.format("Checking lower X-direction send/recv (%d,%d;%d [%d])",
					  i,j,k,cidx(i,j,k)))
	    end
	 end
      end
   end
end

function test_16(comm, nlayer, numComponents, ordering)
   -- tests MPI_Datatype API with a sub-region of a given region
   local sz = Mpi.Comm_size(comm)
   local rnk = Mpi.Comm_rank(comm)

   local rnk = Mpi.Comm_rank(comm)
   if sz ~= 2 then
      log("Test for MPI_Datatype (test_16) not run as number of procs not exactly 2")
      return
   end

   local range = Range.Range({1, 1, 1}, {10, 20, 30})
   local rangeX = Range.Range(
      {range:lower(1), range:lower(2)+nlayer, range:lower(3)+nlayer},
      {range:lower(1), range:upper(2)-nlayer, range:upper(3)-nlayer})
   
   local dTypeX = Mpi.createDataTypeFromRangeAndSubRange(rangeX, range, numComponents, ordering, Mpi.DOUBLE)

   local nz = range:volume()*numComponents
   local buff = Alloc.Double(nz)

   local indexer = range:indexer(ordering)   
   local function cidx(i,j,k,c)
      return (indexer(i,j,k)-1)*numComponents+c
   end
   
   if rnk == 0 then
      for i = range:lower(1), range:upper(1) do
   	 for j = range:lower(2), range:upper(2) do
	    for k = range:lower(3), range:upper(3) do
	       for c = 1, numComponents do
		  buff[cidx(i,j,k,c)] = i+20*j+10*k+0.5+1000*c
	       end
	    end
	 end
      end
   end

   -- X direction
   if rnk == 0 then
      -- send box
      local loc = cidx(rangeX:lower(1), rangeX:lower(2), rangeX:lower(3), 0)
      Mpi.Send(buff:data()+loc, 1, dTypeX, 1, 42, comm)
   else
      -- recv box
      local loc = cidx(rangeX:lower(1), rangeX:lower(2), rangeX:lower(3), 0)
      Mpi.Recv(buff:data()+loc, 1, dTypeX, 0, 42, comm, nil)

      for i = rangeX:lower(1), rangeX:upper(1) do
	 for j = rangeX:lower(2), rangeX:upper(2) do
	    for k = rangeX:lower(3), rangeX:upper(3) do
	       for c = 1, numComponents do
		  assert_equal(i+20*j+10*k+0.5+1000*c, buff[cidx(i,j,k,c)],
			       string.format("Checking lower X-direction send/recv (%d,%d,%d;%d [%d])",
					     i,j,k,c,cidx(i,j,k,c)))
	       end
	    end
	 end
      end
   end
end

-- This test is a repeat of test 2, but communicating device data using CUDA-aware MPI.
-- We copy a vIn array into rank 0 (one GPU) and the vOut array onto every other rank (every other GPU).
-- We then check that this communication has occurred correctly and the bandwidth of this communication.
function test_17(comm)
   if not GKYL_HAVE_CUDA then return end

   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
      
   local nz = 100
   -- allocate memory on host
   local vIn, vOut = Alloc.Double(nz), Alloc.Double(nz)
   -- Allocate memory on device.
   local d_vIn, d_vOut = cuAlloc.Double(nz), cuAlloc.Double(nz)

   if rank == 0 then
      -- send data from rank 0 to all other ranks
      for i = 1, vIn:size() do
	 vIn[i] = i
      end
      -- Copy vIn from host memory (now contains non-zero values).
      local err = d_vIn:copyHostToDevice(vIn)
      assert_equal(cuda.Success, err, "Checking if Memcpy worked")
      for dest = 1, sz-1 do
	 Mpi.Send(d_vIn:data(), nz, Mpi.DOUBLE, dest, 22, comm)
      end      
   else
      -- copy vOut from host memory.
      local err = d_vOut:copyHostToDevice(vOut)
      assert_equal(cuda.Success, err, "Checking if Memcpy worked")
      -- recv stuff from rank 0
      Mpi.Recv(d_vOut:data(), nz, Mpi.DOUBLE, 0, 22, comm, nil)

      -- Copy from device back to host memory to check communication.
      local err = d_vOut:copyDeviceToHost(vOut)
      assert_equal(cuda.Success, err, "Checking if Memcpy worked")
      for i = 1, vOut:size() do
	 assert_equal(vOut[i], i, "Checking recv data")
      end
   end
   Mpi.Barrier(comm)
end

-- Test MPI_Cart routines.
function test_18(comm)
   assert_equal(true, Mpi.Is_comm_valid(comm))

   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
   if sz ~= 4 then
      log("Test of MPI_Cart not run as number of procs not exactly 4")
      return
   end

   local cuts = {2,2}
   local isDirPeriodic = {false, true}

   -- Creation of Cartesian communicator.
   local commNumDims = #cuts
   local commDims = Lin.IntVec(commNumDims)
   local intIsDirPeriodic = Lin.IntVec(commNumDims)
   for d = 1,commNumDims do
      commDims[d] = cuts[d]
      intIsDirPeriodic[d] = isDirPeriodic[d] and 1 or 0
   end
   local reorder = 0
   local commCart = Mpi.Cart_create(comm, commNumDims, commDims, intIsDirPeriodic, reorder)
   assert_equal(sz, Mpi.Comm_size(commCart), "Checking if Cart communicator has correct size")

   -- Obtain Cartesian coordinates of this rank within Cartesian comm.
   -- Test assuming a row-major order distribution of ranks.
   local commCartCoords = Mpi.Cart_coords(commCart, rank, commNumDims)
   assert_equal(math.floor(rank/cuts[1]), commCartCoords[1], "Checking Cart 0th-coordinate")
   assert_equal(rank % cuts[2], commCartCoords[2], "Checking Cart 1st-coordinate")

   -- Retreive information about Cartesian communicator.
   local commDimsNew, intIsDirPeriodicNew, commCartCoordsNew = Mpi.Cart_get(commCart, commNumDims)
   for d = 1,commNumDims do
      assert_equal(commDims[d], commDimsNew[d], "Checking Cart_get dimensions")
      assert_equal(intIsDirPeriodic[d], intIsDirPeriodicNew[d], "Checking Cart_get isDirPeriodic")
      assert_equal(commCartCoords[d], commCartCoordsNew[d], "Checking Cart_get commCartCoords")
   end

   -- Given the Cartesian coordinates, obtain the linear rank.
   local rankNew = Mpi.Cart_rank(commCart, commCartCoords)
   assert_equal(rank, rankNew, "Checking Cart_rank")

   -- Obtain the source/destination for send-receive along each positive direction.
   local srcRankP, destRankP, disp = {}, {}, 1
   srcRankP[1], destRankP[1] = Mpi.Cart_shift(commCart, 0, disp)
   srcRankP[2], destRankP[2] = Mpi.Cart_shift(commCart, 1, disp)
   if commCartCoords[1] == 0 then
      assert_equal(Mpi.PROC_NULL, srcRankP[1], "Checking Cart_shift in 0th positive direction (src a)")
      assert_equal(rank+cuts[1], destRankP[1], "Checking Cart_shift in 0th positive direction (dest a)")
   else
      assert_equal(rank-cuts[1], srcRankP[1], "Checking Cart_shift in 0th positive direction (src b)")
      assert_equal(Mpi.PROC_NULL, destRankP[1], "Checking Cart_shift in 0th positive direction (dest b)")
   end
   if commCartCoords[2] == 0 then
      assert_equal(rank+1, srcRankP[2], "Checking Cart_shift in 1st positive direction (src a)")
      assert_equal(rank+1, destRankP[2], "Checking Cart_shift in 1st positive direction (dest a)")
   else
      assert_equal(rank-1, srcRankP[2], "Checking Cart_shift in 1st positive direction (src b)")
      assert_equal(rank-1, destRankP[2], "Checking Cart_shift in 1st positive direction (dest b)")
   end

   -- Obtain the source/destination for send-receive along each negative direction.
   local srcRankM, destRankM, disp = {}, {}, -1
   srcRankM[1], destRankM[1] = Mpi.Cart_shift(commCart, 0, disp)
   srcRankM[2], destRankM[2] = Mpi.Cart_shift(commCart, 1, disp)
   if commCartCoords[1] == 1 then
      assert_equal(Mpi.PROC_NULL, srcRankM[1], "Checking Cart_shift in 0th negative direction (src a)")
      assert_equal(rank-cuts[1], destRankM[1], "Checking Cart_shift in 0th negative direction (dest a)")
   else
      assert_equal(rank+cuts[1], srcRankM[1], "Checking Cart_shift in 0th negative direction (src b)")
      assert_equal(Mpi.PROC_NULL, destRankM[1], "Checking Cart_shift in 0th negative direction (dest b)")
   end
   if commCartCoords[2] == 1 then
      assert_equal(rank-1, srcRankM[2], "Checking Cart_shift in 1st negative direction (src a)")
      assert_equal(rank-1, destRankM[2], "Checking Cart_shift in 1st negative direction (dest a)")
   else
      assert_equal(rank+1, srcRankM[2], "Checking Cart_shift in 1st negative direction (src b)")
      assert_equal(rank+1, destRankM[2], "Checking Cart_shift in 1st negative direction (dest b)")
   end

   -- Create sub-communicators along each direction.
   local commCart1D = {}
   local keepDir = Lin.IntVec(commNumDims)
   keepDir[1], keepDir[2] = 1, 0
   commCart1D[1] = Mpi.Cart_sub(commCart, keepDir)
   keepDir[1], keepDir[2] = 0, 1
   commCart1D[2] = Mpi.Cart_sub(commCart, keepDir)

   local cc1d = Mpi.Cartdim_get(commCart1D[1])
   local cc2d = Mpi.Cartdim_get(commCart1D[2])

   assert_equal(1, cc1d, "Checking if dimension of sub-comm is correct")
   assert_equal(1, cc2d, "Checking if dimension of sub-comm is correct")

   -- Obtain the dimensionality of the Cartesian communicator. 
   local commCartNumDims = Mpi.Cartdim_get(commCart)
   assert_equal(commNumDims, commCartNumDims, "Checking Cartdim_get")

   Mpi.Barrier(comm)
end

function test_19(comm)
   -- Test graph creation.
   -- The graph is connected as follows:
   -- 0 -> 1
   --   -> 3
   -- 1 -> 2
   -- 2 -> 3
   -- 3 -> 2
   -- 4 -> 1
   --   -> 2
   --   -> 3

   assert_equal(true, Mpi.Is_comm_valid(comm))

   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
   if sz ~= 5 then
      log("Test of MPI_Dist_graph_create_adjacent (test_19) not run as number of procs not exactly 5")
      return
   end

   local indeg, outdeg
   local src, dest
   if rank==0 then
      indeg  = 0
      src    = Lin.IntVec(0)
      outdeg = 2
      dest   = Lin.IntVec(2)
      dest[1], dest[2] = 1, 3
   elseif rank==1 then
      indeg  = 2
      src    = Lin.IntVec(2)
      src[1], src[2] = 0, 4
      outdeg = 1
      dest   = Lin.IntVec(1)
      dest[1] = 2
   elseif rank==2 then
      indeg  = 3
      src    = Lin.IntVec(3)
      src[1], src[2], src[3] = 1, 3, 4
      outdeg = 1
      dest   = Lin.IntVec(1)
      dest[1] = 3
   elseif rank==3 then
      indeg  = 3
      src    = Lin.IntVec(3)
      src[1], src[2], src[3] = 0, 2, 4
      outdeg = 1
      dest   = Lin.IntVec(1)
      dest[1] = 2
   elseif rank==4 then
      indeg  = 0
      src    = Lin.IntVec(0)
      outdeg = 3
      dest   = Lin.IntVec(3)
      dest[1], dest[2], dest[3] = 1, 2, 3
   end
   local srcw, destw = Lin.IntVec(indeg), Lin.IntVec(outdeg)
   for d = 1, indeg do srcw[d] = 1 end
   for d = 1, outdeg do destw[d] = 1 end

   local reorder = 0
   local graphComm = Mpi.Dist_graph_create_adjacent(comm, indeg, src, srcw:data(),
                                                    outdeg, dest, destw:data(), Mpi.INFO_NULL, reorder)
--   -- Could also use UNWEIGHTED since the weights are not being used.
--   local graphComm = Mpi.Dist_graph_create_adjacent(comm, indeg, src, Mpi.UNWEIGHTED,
--                                                    outdeg, dest, Mpi.UNWEIGHTED, Mpi.INFO_NULL, reorder)
   Mpi.Barrier(comm)

   -- Check graph:
   local indegree, outdegree, weighted = Mpi.Dist_graph_neighbors_count(graphComm)

   assert_equal(indeg, indegree, "test_19: indeg != indegree.")
   assert_equal(outdeg, outdegree, "test_19: indeg != indegree.")

   local sources, sourceweights, destinations, destinationweights = Mpi.Dist_graph_neighbors(graphComm, 4, 4)
   for d = 1, #src do
      assert_equal(sources[d], src[d], "test_19: src != sources")
      assert_equal(sourceweights[d], srcw[d], "test_19: srcw != sourceweights")
   end
   for d = 1, #dest do
      assert_equal(destinations[d], dest[d], "test_19: dest != destinations")
      assert_equal(destinationweights[d], destw[d], "test_19: destw != destinationweights")
   end

   Mpi.Barrier(comm)
end

function test_20(comm)
   -- Test Neighbor_allgather.

   assert_equal(true, Mpi.Is_comm_valid(comm))

   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
   if sz ~= 4 then
      log("Test of MPI_Neighbor_allgather not run as number of procs not exactly 4")
      return
   end

   -- Use the graph
   -- 2 -> 1
   -- 3 -> 1

   local data = Lin.Vec(2)
   data[1], data[2] = -math.pi*rank, math.pi*rank+1.

   local indeg, outdeg
   local src, dest
   if rank==0 then
      indeg, outdeg = 0, 0
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
   elseif rank==1 then
      indeg, outdeg = 2, 0
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      src[1], src[2] = 2, 3
   elseif rank==2 then
      indeg, outdeg = 0, 1
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      dest[1] = 1
   elseif rank==3 then
      indeg, outdeg = 0, 1
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      dest[1] = 1
   end

   local srcw, destw = Lin.IntVec(indeg), Lin.IntVec(outdeg)
   for d = 1, indeg do srcw[d] = 1 end
   for d = 1, outdeg do destw[d] = 1 end

   local reorder = 0
   local graphComm = Mpi.Dist_graph_create_adjacent(comm, indeg, src, srcw:data(),
                                                    outdeg, dest, destw:data(), Mpi.INFO_NULL, reorder)
--   -- Could also use UNWEIGHTED since the weights are not being used.
--   local graphComm = Mpi.Dist_graph_create_adjacent(comm, indeg, src, Mpi.UNWEIGHTED,
--                                                    outdeg, dest, Mpi.UNWEIGHTED, Mpi.INFO_NULL, reorder)

   Mpi.Barrier(comm)

   local dataGlobal = Lin.Vec(#data*2)

   Mpi.Neighbor_allgather(data:data(), #data, Mpi.DOUBLE, dataGlobal:data(), #data, Mpi.DOUBLE, graphComm)

   Mpi.Barrier(comm)

   if rank==1 then
      assert_equal(dataGlobal[1], -math.pi*2,    "test_20: dataGlobal[1] in rank=1 is erroneous.")
      assert_equal(dataGlobal[2],  math.pi*2+1., "test_20: dataGlobal[2] in rank=1 is erroneous.")
      assert_equal(dataGlobal[3], -math.pi*3,    "test_20: dataGlobal[3] in rank=1 is erroneous.")
      assert_equal(dataGlobal[4],  math.pi*3+1., "test_20: dataGlobal[4] in rank=1 is erroneous.")
   else
      assert_equal(dataGlobal[1], 0., "test_20: dataGlobal[1] is erroneous.")
      assert_equal(dataGlobal[2], 0., "test_20: dataGlobal[2] is erroneous.")
      assert_equal(dataGlobal[3], 0., "test_20: dataGlobal[3] is erroneous.")
      assert_equal(dataGlobal[4], 0., "test_20: dataGlobal[4] is erroneous.")
   end
end

function test_21(comm)
   -- Test Neighbor_allgather with an MPI derived datatype.

   assert_equal(true, Mpi.Is_comm_valid(comm))

   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
   if sz ~= 4 then
      log("Test of MPI_Neighbor_allgather (test_21) not run as number of procs not exactly 4")
      return
   end

   -- Use the graph
   -- 2 -> 0
   -- 3 -> 0
   -- 2 -> 1
   -- 3 -> 1

   local indeg, outdeg
   local src, dest
   if rank==0 then
      indeg, outdeg = 2, 0
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      src[1], src[2] = 2, 3
   elseif rank==1 then
      indeg, outdeg = 2, 0
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      src[1], src[2] = 2, 3
   elseif rank==2 then
      indeg, outdeg = 0, 2
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      dest[1], dest[2] = 0, 1
   elseif rank==3 then
      indeg, outdeg = 0, 2
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      dest[1], dest[2] = 0, 1
   end

   local srcw, destw = Lin.IntVec(indeg), Lin.IntVec(outdeg)
   for d = 1, indeg do srcw[d] = 1 end
   for d = 1, outdeg do destw[d] = 1 end

   local reorder = 0
   local graphComm = Mpi.Dist_graph_create_adjacent(comm, indeg, src, srcw:data(),
                                                    outdeg, dest, destw:data(), Mpi.INFO_NULL, reorder)
--   -- Could also use UNWEIGHTED since the weights are not being used.
--   local graphComm = Mpi.Dist_graph_create_adjacent(comm, indeg, src, Mpi.UNWEIGHTED,
--                                                    outdeg, dest, Mpi.UNWEIGHTED, Mpi.INFO_NULL, reorder)

   Mpi.Barrier(comm)

   local data = Lin.Vec(7)
   for i = 1, #data do data[i] = math.pi*rank+i end

   -- Have ranks 0 and 1 gather elements 2-3 and 5-6 from ranks 2 and 3:
   -- Construct a vector MPI datatype for this.
   local count, blocklength, stride = 2, 2, 3
   local sendType = Mpi.Type_vector(count, blocklength, stride, Mpi.DOUBLE)
   Mpi.Type_commit(sendType)

   -- The receiver will put them all in a compressed buffer.
   local count, blocklength, stride = 2, 2, 2
   local recvType = Mpi.Type_vector(count, blocklength, stride, Mpi.DOUBLE)
   Mpi.Type_commit(recvType)

   local dataGlobal = Lin.Vec(count*blocklength*2)

   Mpi.Neighbor_allgather(data:data()+1, 1, sendType, dataGlobal:data(), 1, recvType, graphComm)

   Mpi.Barrier(comm)

   if rank==0 or rank==1 then
      assert_equal(math.pi*2+2., dataGlobal[1], "test_21: dataGlobal[1] in rank=0,1 is erroneous.")
      assert_equal(math.pi*2+3., dataGlobal[2], "test_21: dataGlobal[2] in rank=0,1 is erroneous.")
      assert_equal(math.pi*2+5., dataGlobal[3], "test_21: dataGlobal[3] in rank=0,1 is erroneous.")
      assert_equal(math.pi*2+6., dataGlobal[4], "test_21: dataGlobal[4] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+2., dataGlobal[5], "test_21: dataGlobal[5] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+3., dataGlobal[6], "test_21: dataGlobal[6] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+5., dataGlobal[7], "test_21: dataGlobal[7] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+6., dataGlobal[8], "test_21: dataGlobal[8] in rank=0,1 is erroneous.")
   else
      for i = 1, #dataGlobal do
        assert_equal(dataGlobal[i], 0., "test_21: dataGlobal is erroneous in rank=2,3.")
      end
   end

   Mpi.Type_free(sendType)
   Mpi.Type_free(recvType)
end

function test_22(comm)
   -- Test Neighbor_alltoall.

   assert_equal(true, Mpi.Is_comm_valid(comm))

   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
   if sz ~= 4 then
      log("Test of MPI_Neighbor_alltoall (test_22) not run as number of procs not exactly 4")
      return
   end

   -- Use the hourglass graph
   -- 2----3   0 -> 1
   --  \  /    1 -> 2
   --   \/     2 -> 3
   --   /\     3 -> 0
   --  /  \
   -- 0----1

   local indeg, outdeg = 1, 1
   local src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
   if rank==0 then
      src[1], dest[1] = 3, 1
   elseif rank==1 then
      src[1], dest[1] = 0, 2
   elseif rank==2 then
      src[1], dest[1] = 1, 3
   elseif rank==3 then
      src[1], dest[1] = 2, 0
   end
   local reorder = 0
   local graphComm = Mpi.Dist_graph_create_adjacent(comm, indeg, src, Mpi.UNWEIGHTED,
                                                    outdeg, dest, Mpi.UNWEIGHTED, Mpi.INFO_NULL, reorder)
   Mpi.Barrier(comm)

   local numE = 6
   local data = Lin.Vec(numE)  -- Assume first and last elements are "ghost" elements.
   for i = 2, #data-1 do data[i] = math.pi*rank+i end

   -- Send right "skin" element and place it in left "ghost".
   Mpi.Neighbor_alltoall(data:data()+numE-2, 1, Mpi.DOUBLE, data:data(), 1, Mpi.DOUBLE, graphComm)

   Mpi.Barrier(comm)

   assert_equal(math.pi*rank+2., data[2], "test_22: data[2] is erroneous.")
   assert_equal(math.pi*rank+3., data[3], "test_22: data[3] is erroneous.")
   assert_equal(math.pi*rank+4., data[4], "test_22: data[4] is erroneous.")
   assert_equal(math.pi*rank+5., data[5], "test_22: data[5] is erroneous.")
   if rank==0 then
      assert_equal(math.pi*3+5., data[1], "test_22: data[1] is erroneous.")
   else
      assert_equal(math.pi*(rank-1)+5., data[1], "test_22: data[1] is erroneous.")
   end
end

function test_23(comm)
   -- Test Neighbor_alltoallv.

   assert_equal(true, Mpi.Is_comm_valid(comm))

   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
   if sz ~= 4 then
      log("Test of MPI_Neighbor_alltoallv (test_23) not run as number of procs not exactly 4")
      return
   end

   -- Use the graph
   -- 2 -> 0
   -- 3 -> 0
   -- 2 -> 1
   -- 3 -> 1
   -- And perform the an operation similar to that in test_21.
   local indeg, outdeg
   local src, dest
   if rank==0 then
      indeg, outdeg = 2, 0
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      src[1], src[2] = 2, 3
   elseif rank==1 then
      indeg, outdeg = 2, 0
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      src[1], src[2] = 2, 3
   elseif rank==2 then
      indeg, outdeg = 0, 2
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      dest[1], dest[2] = 0, 1
   elseif rank==3 then
      indeg, outdeg = 0, 2
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      dest[1], dest[2] = 0, 1
   end
   local reorder = 0
   local graphComm = Mpi.Dist_graph_create_adjacent(comm, indeg, src, Mpi.UNWEIGHTED,
                                                    outdeg, dest, Mpi.UNWEIGHTED, Mpi.INFO_NULL, reorder)


   local data = Lin.Vec(7)
   for i = 1, #data do data[i] = math.pi*rank+i end
   data[1] = -1*(rank+1)
   data[4] = -4*(rank+1)

   -- Have ranks 0 and 1 gather elements 1-3 and 4-6 (in 2 blocks) from ranks 2 and 3:
   -- Construct a vector MPI datatype for this.
   local count, blocklength, stride = 2, 3, 3
   local sendType = Mpi.Type_vector(count, blocklength, stride, Mpi.DOUBLE)
   Mpi.Type_commit(sendType)

   -- The receiver will put them all in a compressed buffer.
   local count, blocklength, stride = 2, 3, 3
   local recvType = Mpi.Type_vector(count, blocklength, stride, Mpi.DOUBLE)
   Mpi.Type_commit(recvType)

   local dataGlobal = Lin.Vec(count*blocklength*2)

   local sendCounts, recvCounts = Lin.IntVec(2), Lin.IntVec(2)
   local sendDispls, recvDispls = Lin.IntVec(2), Lin.IntVec(2)

   sendCounts[1], sendCounts[2] = 1, 1
   recvCounts[1], recvCounts[2] = 1, 1
   sendDispls[1], sendDispls[2] = 0, 0
   -- MF 03/10/2022: Documentation/Examples for combining alltoallv with MPIDDs are terrible.
   -- They led me to believe that recvDispls[2] below should be count*blocklenth,
   -- but that leads to a seg fault. It works with 1, which I interpret to mean that
   -- the displacement should be in terms of some size of the MPIDD. Not totally sure.
   recvDispls[1], recvDispls[2] = 0, 1
   Mpi.Neighbor_alltoallv(data:data(), sendCounts:data(), sendDispls:data(), sendType,
                          dataGlobal:data(), recvCounts:data(), recvDispls:data(), recvType, graphComm)

--   -- Could also do it without MPI data types as follows:
--   sendCounts[1], sendCounts[2] = 6, 6
--   recvCounts[1], recvCounts[2] = 6, 6
--   sendDispls[1], sendDispls[2] = 0, 0
--   recvDispls[1], recvDispls[2] = 0, count*blocklength
--   Mpi.Neighbor_alltoallv(data:data(), sendCounts:data(), sendDispls:data(), Mpi.DOUBLE,
--                          dataGlobal:data(), recvCounts:data(), recvDispls:data(), Mpi.DOUBLE, graphComm)

   Mpi.Barrier(comm)

   if rank==0 or rank==1 then
      assert_equal(-3          , dataGlobal[1 ], "test_23: dataGlobal[1 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*2+2., dataGlobal[2 ], "test_23: dataGlobal[2 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*2+3., dataGlobal[3 ], "test_23: dataGlobal[3 ] in rank=0,1 is erroneous.")
      assert_equal(-12         , dataGlobal[4 ], "test_23: dataGlobal[4 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*2+5., dataGlobal[5 ], "test_23: dataGlobal[5 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*2+6., dataGlobal[6 ], "test_23: dataGlobal[6 ] in rank=0,1 is erroneous.")
      assert_equal(-4          , dataGlobal[7 ], "test_23: dataGlobal[7 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+2., dataGlobal[8 ], "test_23: dataGlobal[8 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+3., dataGlobal[9 ], "test_23: dataGlobal[9 ] in rank=0,1 is erroneous.")
      assert_equal(-16         , dataGlobal[10], "test_23: dataGlobal[10] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+5., dataGlobal[11], "test_23: dataGlobal[11] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+6., dataGlobal[12], "test_23: dataGlobal[12] in rank=0,1 is erroneous.")
   else
      for i = 1, #dataGlobal do
        assert_equal(dataGlobal[i], 0., "test_23: dataGlobal is erroneous in rank=2,3.")
      end
   end

end

function test_24(comm)
   -- Test Neighbor_alltoallw.

   assert_equal(true, Mpi.Is_comm_valid(comm))

   local rank = Mpi.Comm_rank(comm)
   local sz = Mpi.Comm_size(comm)
   if sz ~= 4 then
      log("Test of MPI_Neighbor_alltoallw (test_24) not run as number of procs not exactly 4")
      return
   end

   -- Use the graph
   -- 2 -> 0
   -- 3 -> 0
   -- 2 -> 1
   -- 3 -> 1
   -- And perform the an operation similar to that in test_21.
   local indeg, outdeg
   local src, dest
   if rank==0 then
      indeg, outdeg = 2, 0
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      src[1], src[2] = 2, 3
   elseif rank==1 then
      indeg, outdeg = 2, 0
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      src[1], src[2] = 2, 3
   elseif rank==2 then
      indeg, outdeg = 0, 2
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      dest[1], dest[2] = 0, 1
   elseif rank==3 then
      indeg, outdeg = 0, 2
      src, dest     = Lin.IntVec(indeg), Lin.IntVec(outdeg)
      dest[1], dest[2] = 0, 1
   end
   local srcw, destw = Lin.IntVec(indeg), Lin.IntVec(outdeg)
   for d = 1, indeg do srcw[d] = 1 end
   for d = 1, outdeg do destw[d] = 1 end
   local reorder = 0
   local graphComm = Mpi.Dist_graph_create_adjacent(comm, indeg, src, srcw:data(),
                                                    outdeg, dest, destw:data(), Mpi.INFO_NULL, reorder)
--   -- Could also use UNWEIGHTED since the weights are not being used.
--   local graphComm = Mpi.Dist_graph_create_adjacent(comm, indeg, src, Mpi.UNWEIGHTED,
--                                                    outdeg, dest, Mpi.UNWEIGHTED, Mpi.INFO_NULL, reorder)


   local data = Lin.Vec(7)
   for i = 1, #data do data[i] = math.pi*rank+i end
   data[1] = -1*(rank+1)
   data[4] = -4*(rank+1)

   -- We wish ranks 0 and 1 to end with an array that has
   -- { 0, data[1] from rank 2, data[2] from rank 2, data[4] from rank 2, data[5] from rank 2,
   -- { 0, data[1] from rank 3, data[2] from rank 3, data[4] from rank 3, data[5] from rank 3}
   -- We will use MPI_Type_vectors for this.
   local dataGlobal = Lin.Vec(10)
   for i = 1, 10 do dataGlobal[i] = 0. end

   -- Construct the send MPI datatype.
   local count, blocklength, stride = 2, 2, 3
   local sendType = Mpi.MPI_Datatype_vec(2)
   sendType[0] = Mpi.Type_vector(count, blocklength, stride, Mpi.DOUBLE)[0]
   sendType[1] = Mpi.Type_vector(count, blocklength, stride, Mpi.DOUBLE)[0]
   Mpi.Type_commit(sendType)
   Mpi.Type_commit(sendType+1)

   -- The receiver will put them all in a compressed buffer, with one empty
   -- element in between what's received from rank 2 and rank 3.
   local count, blocklength, stride = 2, 2, 2
   local recvType = Mpi.MPI_Datatype_vec(2)
   recvType[0] = Mpi.Type_vector(count, blocklength, stride, Mpi.DOUBLE)[0]
   recvType[1] = Mpi.Type_vector(count, blocklength, stride, Mpi.DOUBLE)[0]
   Mpi.Type_commit(recvType)
   Mpi.Type_commit(recvType+1)

   -- Counts =1 since we are using MPI derived data types.
   local sendCounts, recvCounts = Lin.IntVec(2), Lin.IntVec(2)
   sendCounts[1], sendCounts[2] = 1, 1
   recvCounts[1], recvCounts[2] = 1, 1

   local sendDispls, recvDispls = Mpi.MPI_Aint_vec(2), Mpi.MPI_Aint_vec(2)
   -- Send displacements are 0 since we'll pass globalData[1] as the send
   -- buffer. Receive displacements are 0 for the first 4 elements received
   -- but 5*sizeof(double) for the next 4 elements.
   sendDispls[0], sendDispls[0] = 0, 0
   recvDispls[0], recvDispls[1] = 0, 5*sizeof(Mpi.DOUBLE)
   Mpi.Neighbor_alltoallw(data:data()+1, sendCounts:data(), sendDispls, sendType,
                          dataGlobal:data()+1, recvCounts:data(), recvDispls, recvType, graphComm)
---- ......................................
--   -- There's another way to set the displacements. One uses MPI_Bottom as the
--   -- buffer, and places the desired addresses in the displacements.
--   sendDispls[0] = Mpi.Get_address(data:data()+1)[0]
--   sendDispls[1] = Mpi.Get_address(data:data()+1)[0]
--   recvDispls[0] = Mpi.Get_address(dataGlobal:data()+1)[0]
--   recvDispls[1] = Mpi.Get_address(dataGlobal:data()+6)[0]
--   Mpi.Neighbor_alltoallw(Mpi.BOTTOM, sendCounts:data(), sendDispls, sendType,
--                          Mpi.BOTTOM, recvCounts:data(), recvDispls, recvType, graphComm)
---- ......................................

   Mpi.Barrier(comm)

   if rank==0 or rank==1 then
      assert_equal(0.          , dataGlobal[1 ], "test_24: dataGlobal[1 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*2+2., dataGlobal[2 ], "test_24: dataGlobal[2 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*2+3., dataGlobal[3 ], "test_24: dataGlobal[3 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*2+5., dataGlobal[4 ], "test_24: dataGlobal[4 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*2+6., dataGlobal[5 ], "test_24: dataGlobal[5 ] in rank=0,1 is erroneous.")
      assert_equal(0.          , dataGlobal[6 ], "test_24: dataGlobal[6 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+2., dataGlobal[7 ], "test_24: dataGlobal[7 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+3., dataGlobal[8 ], "test_24: dataGlobal[8 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+5., dataGlobal[9 ], "test_24: dataGlobal[9 ] in rank=0,1 is erroneous.")
      assert_equal(math.pi*3+6., dataGlobal[10], "test_24: dataGlobal[10] in rank=0,1 is erroneous.")
   else
      for i = 1, #dataGlobal do
        assert_equal(dataGlobal[i], 0., "test_24: dataGlobal is erroneous in rank=2,3.")
      end
   end

end

function test_25(comm)
   -- Test MPI_Comm_create_group
   assert_equal(true, Mpi.Is_comm_valid(comm))

   local sz, rank = Mpi.Comm_size(comm), Mpi.Comm_rank(comm)

   if sz ~= 5 then
      log("Test of MPI_Comm_create_group (test_25) not run as number of procs not exactly 5")
      return
   end

   -- Create a communicator with even ranks, and one with odd ranks.
   local commGroup = Mpi.Comm_group(comm)
   local subComm, subGroup
   local subGroupRanks, tag
   if rank % 2 == 0 then
      subGroupRanks = Lin.IntVec(3)
      subGroupRanks[1], subGroupRanks[2], subGroupRanks[3] = 0, 2, 4
      tag = 0
   else
      subGroupRanks = Lin.IntVec(2)
      subGroupRanks[1], subGroupRanks[2] = 1, 3
      tag = 1
   end
   subGroup = Mpi.Group_incl(commGroup, #subGroupRanks, subGroupRanks:data());
   subComm = Mpi.Comm_create_group(comm, subGroup, tag);
--   -- In this example we could've instead called these outside of the if-statement:
--   subGroup = Mpi.Group_incl(commGroup, #subGroupRanks, subGroupRanks:data());
--   subComm = Mpi.Comm_create_group(comm, subGroup, 1);

   local data, red = Lin.Vec(1), Lin.Vec(1)
   data[1] = 0.1+rank

   Mpi.Allreduce(data:data(), red:data(), 1, Mpi.DOUBLE, Mpi.SUM, subComm)

   if rank % 2 == 0 then
      assert_equal(3*0.1+6, red[1], "test_25: Checking even allReduce sum")
   else
      assert_equal(2*0.1+4, red[1], "test_25: Checking odd allReduce sum")
   end

   Mpi.Barrier(comm)
end

function test_26(comm)
   local sz = Mpi.Comm_size(comm)
   if sz < 2 then
      log("Test 19 for Mpi.Allgather not run as the number of processes is less than 2")
      return
   end

   -- Get the current process rank
   local rank = Mpi.Comm_rank(comm)

   -- Prepare test data (unique value per rank)
   local nz = 100
   local sendBuf = Alloc.Double(nz)
   local recvBuf = Alloc.Double(nz * sz)
   for i = 1, nz do
      sendBuf[i] = rank + (i * 0.1)  -- Unique value for each rank
   end

   -- Perform the Mpi.Allgather operation
   Mpi.Allgather(sendBuf:data(), nz, Mpi.DOUBLE, recvBuf:data(), nz, Mpi.DOUBLE, comm)

   -- Verify the correctness of the gathered data
   for r = 0, sz - 1 do
      for i = 1, nz do
         local expectedValue = r + (i * 0.1)  -- Expected value for each rank
         local index = r * nz + i
         assert_equal(expectedValue, recvBuf[index], "Checking gathered data")
      end
   end

   -- Barrier synchronization to ensure all processes complete the test
   Mpi.Barrier(comm)
end

-- Run tests
test_0(Mpi.COMM_WORLD)
test_1(Mpi.COMM_WORLD)
test_2(Mpi.COMM_WORLD)
test_3(Mpi.COMM_WORLD)
test_4(Mpi.COMM_WORLD)
test_5(Mpi.COMM_WORLD)
test_6a(Mpi.COMM_WORLD)
test_6b(Mpi.COMM_WORLD)
test_6c(Mpi.COMM_WORLD)
test_7(Mpi.COMM_WORLD)
test_8(Mpi.COMM_WORLD)
test_9(Mpi.COMM_WORLD)
test_10(Mpi.COMM_WORLD)
test_11(Mpi.COMM_WORLD)

test_12(Mpi.COMM_WORLD, 1, Range.rowMajor)
test_12(Mpi.COMM_WORLD, 1, Range.colMajor)
test_12(Mpi.COMM_WORLD, 2, Range.rowMajor)
test_12(Mpi.COMM_WORLD, 2, Range.colMajor)
test_12(Mpi.COMM_WORLD, 3, Range.rowMajor)
test_12(Mpi.COMM_WORLD, 3, Range.colMajor)

test_13(Mpi.COMM_WORLD, 1, 2, Range.rowMajor)

test_14(Mpi.COMM_WORLD, 1, 2, Range.colMajor)
test_14(Mpi.COMM_WORLD, 1, 2, Range.rowMajor)

test_15(Mpi.COMM_WORLD, 0, 1, Range.rowMajor)
test_15(Mpi.COMM_WORLD, 0, 2, Range.rowMajor)
test_15(Mpi.COMM_WORLD, 0, 1, Range.colMajor)
test_15(Mpi.COMM_WORLD, 0, 2, Range.colMajor)

test_15(Mpi.COMM_WORLD, 1, 1, Range.rowMajor)
test_15(Mpi.COMM_WORLD, 1, 2, Range.rowMajor)
test_15(Mpi.COMM_WORLD, 1, 1, Range.colMajor)
test_15(Mpi.COMM_WORLD, 1, 2, Range.colMajor)

test_15(Mpi.COMM_WORLD, 2, 1, Range.rowMajor)
test_15(Mpi.COMM_WORLD, 2, 2, Range.rowMajor)
test_15(Mpi.COMM_WORLD, 2, 1, Range.colMajor)
test_15(Mpi.COMM_WORLD, 2, 2, Range.colMajor)

test_16(Mpi.COMM_WORLD, 0, 1, Range.rowMajor)
test_16(Mpi.COMM_WORLD, 0, 2, Range.rowMajor)
test_16(Mpi.COMM_WORLD, 0, 1, Range.colMajor)
test_16(Mpi.COMM_WORLD, 0, 2, Range.colMajor)

test_16(Mpi.COMM_WORLD, 1, 1, Range.rowMajor)
test_16(Mpi.COMM_WORLD, 1, 2, Range.rowMajor)
test_16(Mpi.COMM_WORLD, 1, 1, Range.colMajor)
test_16(Mpi.COMM_WORLD, 1, 2, Range.colMajor)

test_16(Mpi.COMM_WORLD, 2, 1, Range.rowMajor)
test_16(Mpi.COMM_WORLD, 2, 2, Range.rowMajor)
test_16(Mpi.COMM_WORLD, 2, 1, Range.colMajor)
test_16(Mpi.COMM_WORLD, 2, 2, Range.colMajor)

test_17(Mpi.COMM_WORLD)
test_18(Mpi.COMM_WORLD)

test_19(Mpi.COMM_WORLD)
test_20(Mpi.COMM_WORLD)
test_21(Mpi.COMM_WORLD)
test_22(Mpi.COMM_WORLD)
test_23(Mpi.COMM_WORLD)
test_24(Mpi.COMM_WORLD)
test_25(Mpi.COMM_WORLD)

test_26(Mpi.COMM_WORLD)

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

