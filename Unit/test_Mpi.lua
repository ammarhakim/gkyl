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
      	 i = range:upper(1)-n+1
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
      	 i = range:upper(1)-n+1
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
      log("Test for MPI_Datatype (test_13) not run as number of procs not exactly 2")
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
      	 i = range:upper(1)-n+1
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

-- Run tests
test_0(Mpi.COMM_WORLD)
test_1(Mpi.COMM_WORLD)
test_2(Mpi.COMM_WORLD)
test_3(Mpi.COMM_WORLD)
test_4(Mpi.COMM_WORLD)
test_5(Mpi.COMM_WORLD)
test_6(Mpi.COMM_WORLD)
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

