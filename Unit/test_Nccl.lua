-- Gkyl ------------------------------------------------------------------------
--
-- Test for MPI wrappers
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Don't do anything if we were not built with CUDA.
if GKYL_HAVE_CUDA == false then
   print("**** Can't run CUDA tests without CUDA enabled GPUs!")
   return 0
end

local Unit  = require "Unit"
local Mpi   = require "Comm.Mpi"
local Alloc = require "Lib.Alloc"
local Nccl  = require "Comm.Nccl"

local cuda
local cuAlloc
if GKYL_HAVE_CUDA then
   cuda    = require "Cuda.RunTime"
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
   if rank == 0 then print(msg) end
end

-- ncclComm and ncclAllReduce.
function test_1(comm)
   assert_equal(true, Mpi.Is_comm_valid(comm))

   local comm_sz = Mpi.Comm_size(comm)
   local rank = Mpi.Comm_rank(comm)

   local num_devices, _ = cuda.GetDeviceCount() 
   local local_rank = rank % num_devices

   -- Picking a GPU based on localRank, allocate device buffers.
   local err = cuda.SetDevice(local_rank)

   -- Initialize some data to be reduced.
   local data_len = 32*1024*1024;
   local type_sz = ffi.sizeof("double")
   local h_sendbuff, h_recvbuff = Alloc.Double(data_len), Alloc.Double(data_len)
   local d_sendbuff, d_recvbuff = cuda.Malloc(type_sz*data_len), cuda.Malloc(type_sz*data_len)
   for i = 1, data_len do h_sendbuff[i] = (rank+1)*i end
   local err = cuda.Memcpy(d_sendbuff, h_sendbuff:data(), h_sendbuff:size()*type_sz, cuda.MemcpyHostToDevice)
   assert_equal(cuda.Success, err, "Checking if Memcpy worked")

   -- Get NCCL unique ID at rank 0 and broadcast it to all others.
   local ncclId = Nccl.UniqueId()
   if rank == 0 then
      local _ = Nccl.GetUniqueId(ncclId)
   end
   Mpi.Bcast(ncclId, sizeof(ncclId), Mpi.BYTE, 0, comm)

   -- Create a new CUDA stream (needed by NCCL).
   local custream = cuda.StreamCreate()

   -- Initializing NCCL.
   local devComm = Nccl.Comm()
   local _ = Nccl.CommInitRank(devComm, comm_sz, ncclId, rank)

   -- Reduce using NCCL.
   local _ Nccl.AllReduce(d_sendbuff, d_recvbuff, data_len, Nccl.Double, Nccl.Sum,
        devComm, custream)

   -- Completing NCCL operation by synchronizing on the CUDA stream.
   local _ = cuda.StreamSynchronize(custream)

   -- Copy data to host and check results.
   local err = cuda.Memcpy(h_recvbuff:data(), d_recvbuff, h_recvbuff:size()*type_sz, cuda.MemcpyDeviceToHost)
   for i = 1, data_len do
      local res = 0
      for r = 0, comm_sz-1 do res = res + (r+1)*i end
      assert_equal(res, h_recvbuff[i], "test_1: Checking Allreduce result.")
   end

   Mpi.Barrier(comm)

   cuda.Free(d_sendbuff); cuda.Free(d_recvbuff)
end

-- ncclSend and ncclRecv.
function test_2(comm)
   assert_equal(true, Mpi.Is_comm_valid(comm))

   local comm_sz = Mpi.Comm_size(comm)
   local rank = Mpi.Comm_rank(comm)
   if comm_sz ~= 2 then
      log("test_2 not run as number of procs not 2 or more")
      return
   end

   local num_devices, _ = cuda.GetDeviceCount() 
   local local_rank = rank % num_devices

   -- Picking a GPU based on localRank, allocate device buffers.
   local err = cuda.SetDevice(local_rank)

   -- Initialize some data to be reduced.
   local data_len = 32*1024;
   local type_sz = ffi.sizeof("double")
   local h_sendbuff, h_recvbuff = Alloc.Double(data_len), Alloc.Double(data_len)
   local d_sendbuff, d_recvbuff = cuda.Malloc(type_sz*data_len), cuda.Malloc(type_sz*data_len)

   -- Get NCCL unique ID at rank 0 and broadcast it to all others.
   local ncclId = Nccl.UniqueId()
   if rank == 0 then
      local _ = Nccl.GetUniqueId(ncclId)
   end
   Mpi.Bcast(ncclId, sizeof(ncclId), Mpi.BYTE, 0, comm)

   -- Create a new CUDA stream (needed by NCCL).
   local custream = cuda.StreamCreate()

   -- Initializing NCCL.
   local devComm = Nccl.Comm()
   local _ = Nccl.CommInitRank(devComm, comm_sz, ncclId, rank)

   -- Communicate data from rank 0 to rank 1.
   if rank == 0 then
      for i = 1, data_len do h_sendbuff[i] = (rank+1)*i end
      local err = cuda.Memcpy(d_sendbuff, h_sendbuff:data(), h_sendbuff:size()*type_sz, cuda.MemcpyHostToDevice)
      assert_equal(cuda.Success, err, "Checking if Memcpy worked")
   end
   if rank == 1 then
      local err = Nccl.Recv(d_recvbuff, data_len, Nccl.Double, 0, devComm, custream)
   end
   if rank == 0 then
      local err = Nccl.Send(d_sendbuff, data_len, Nccl.Double, 1, devComm, custream)
   end

   -- Completing NCCL operation by synchronizing on the CUDA stream.
   local _ = cuda.StreamSynchronize(custream)

   -- Copy data to host and check results.
   if rank == 1 then
     local err = cuda.Memcpy(h_recvbuff:data(), d_recvbuff, h_recvbuff:size()*type_sz, cuda.MemcpyDeviceToHost)
      for i = 1, data_len do
         assert_equal((0+1)*i, h_recvbuff[i], "test_2: Checking 0 -> 1 send/recv result.")
      end
   end

   -- Communicate data from rank 1 to rank 0.
   if rank == 1 then
      for i = 1, data_len do h_sendbuff[i] = (rank+1)*i end
      local err = cuda.Memcpy(d_sendbuff, h_sendbuff:data(), h_sendbuff:size()*type_sz, cuda.MemcpyHostToDevice)
      assert_equal(cuda.Success, err, "Checking if Memcpy worked")
   end
   if rank == 0 then
      local err = Nccl.Recv(d_recvbuff, data_len, Nccl.Double, 1, devComm, custream)
   end
   if rank == 1 then
      local err = Nccl.Send(d_sendbuff, data_len, Nccl.Double, 0, devComm, custream)
   end

   -- Completing NCCL operation by synchronizing on the CUDA stream.
   local _ = cuda.StreamSynchronize(custream)

   -- Copy data to host and check results.
   if rank == 0 then
     local err = cuda.Memcpy(h_recvbuff:data(), d_recvbuff, h_recvbuff:size()*type_sz, cuda.MemcpyDeviceToHost)
      for i = 1, data_len do
         assert_equal((1+1)*i, h_recvbuff[i], "test_2: Checking 1 -> 0 send/recv result.")
      end
   end

   Mpi.Barrier(comm)

   cuda.Free(d_sendbuff); cuda.Free(d_recvbuff)
end

-- ncclSend and ncclRecv with a nonblocking comm.
function test_3(comm)
   assert_equal(true, Mpi.Is_comm_valid(comm))

   local comm_sz = Mpi.Comm_size(comm)
   local rank = Mpi.Comm_rank(comm)
   if comm_sz ~= 2 then
      log("test_3 not run as number of procs not 2 or more")
      return
   end

   local num_devices, _ = cuda.GetDeviceCount() 
   local local_rank = rank % num_devices

   -- Picking a GPU based on localRank, allocate device buffers.
   local err = cuda.SetDevice(local_rank)

   -- Initialize some data to be reduced.
   local data_len = 32*1024;
   local type_sz = ffi.sizeof("double")
   local h_sendbuff, h_recvbuff = Alloc.Double(data_len), Alloc.Double(data_len)
   local d_sendbuff, d_recvbuff = cuda.Malloc(type_sz*data_len), cuda.Malloc(type_sz*data_len)

   -- Get NCCL unique ID at rank 0 and broadcast it to all others.
   local ncclId = Nccl.UniqueId()
   if rank == 0 then
      local _ = Nccl.GetUniqueId(ncclId)
   end
   Mpi.Bcast(ncclId, sizeof(ncclId), Mpi.BYTE, 0, comm)

   -- Create a new CUDA stream (needed by NCCL).
   local custream = cuda.StreamCreate()

   -- Initializing NCCL.
   local commConfig = Nccl.Config()
   local devComm = Nccl.Comm()
   local _ = Nccl.CommInitRankConfig(devComm, comm_sz, ncclId, rank, commConfig)

   -- Communicate data from rank 0 to rank 1.
   if rank == 0 then
      for i = 1, data_len do h_sendbuff[i] = (rank+1)*i end
      local err = cuda.Memcpy(d_sendbuff, h_sendbuff:data(), h_sendbuff:size()*type_sz, cuda.MemcpyHostToDevice)
      assert_equal(cuda.Success, err, "Checking if Memcpy worked")
   end
   if rank == 1 then
      local err = Nccl.Recv(d_recvbuff, data_len, Nccl.Double, 0, devComm, custream)
   end
   if rank == 0 then
      local err = Nccl.Send(d_sendbuff, data_len, Nccl.Double, 1, devComm, custream)
   end

   local devResult = Nccl.Result()
   local _ = Nccl.CommGetAsyncError(devComm, devResult)
   while (devResult[0] == Nccl.InProgress) do
      local _ = Nccl.CommGetAsyncError(devComm, devResult)
   end
   -- Completing NCCL operation by synchronizing on the CUDA stream.
   local _ = cuda.StreamSynchronize(custream)

   -- Copy data to host and check results.
   if rank == 1 then
     local err = cuda.Memcpy(h_recvbuff:data(), d_recvbuff, h_recvbuff:size()*type_sz, cuda.MemcpyDeviceToHost)
      for i = 1, data_len do
         assert_equal((0+1)*i, h_recvbuff[i], "test_3: Checking 0 -> 1 send/recv result.")
      end
   end

   -- Communicate data from rank 1 to rank 0.
   if rank == 1 then
      for i = 1, data_len do h_sendbuff[i] = (rank+1)*i end
      local err = cuda.Memcpy(d_sendbuff, h_sendbuff:data(), h_sendbuff:size()*type_sz, cuda.MemcpyHostToDevice)
      assert_equal(cuda.Success, err, "Checking if Memcpy worked")
   end
   if rank == 0 then
      local err = Nccl.Recv(d_recvbuff, data_len, Nccl.Double, 1, devComm, custream)
   end
   if rank == 1 then
      local err = Nccl.Send(d_sendbuff, data_len, Nccl.Double, 0, devComm, custream)
   end

   local _ = Nccl.CommGetAsyncError(devComm, devResult)
   while (devResult[0] == Nccl.InProgress) do
      local _ = Nccl.CommGetAsyncError(devComm, devResult)
   end
   -- Completing NCCL operation by synchronizing on the CUDA stream.
   local _ = cuda.StreamSynchronize(custream)

   -- Copy data to host and check results.
   if rank == 0 then
     local err = cuda.Memcpy(h_recvbuff:data(), d_recvbuff, h_recvbuff:size()*type_sz, cuda.MemcpyDeviceToHost)
      for i = 1, data_len do
         assert_equal((1+1)*i, h_recvbuff[i], "test_3: Checking 1 -> 0 send/recv result.")
      end
   end

   Mpi.Barrier(comm)

   cuda.Free(d_sendbuff); cuda.Free(d_recvbuff)
end

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
