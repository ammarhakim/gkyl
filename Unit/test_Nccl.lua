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
local Grid             = require "Grid"
local DecompRegionCalc = require "Lib.CartDecomp"
local DataStruct       = require "DataStruct"

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

-- ncclSend and ncclRecv a whole CartField.
function test_4(comm)

   local defaultCommSize = 2
   local comm_sz = Mpi.Comm_size(comm)
   if comm_sz ~= defaultCommSize then
      log(string.format("Not running test_4 as numProcs not exactly %d",defaultCommSize))
      return
   end

   local worldRank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   local numCutsConf = 1

   local confColor = math.floor(worldRank/numCutsConf)
   local confComm  = Mpi.Comm_split(Mpi.COMM_WORLD, confColor, worldRank)
   local speciesColor = worldRank % numCutsConf
   local speciesComm  = Mpi.Comm_split(Mpi.COMM_WORLD, speciesColor, worldRank)
   local confRank = Mpi.Comm_rank(confComm)
   local speciesRank = Mpi.Comm_rank(speciesComm)

   local rank = numCutsConf==1 and speciesRank or confRank
   local comm = numCutsConf==1 and speciesComm or confComm

   local decomp = DecompRegionCalc.CartProd { cuts = {numCutsConf}, comm = confComm }
   local grid = Grid.RectCart {
      lower = {0.0},  cells = {12},
      upper = {1.0},  decomposition = decomp,
   }
   local field0 = DataStruct.Field {
      onGrid        = grid,
      numComponents = 3,
      ghost         = {1, 1},
   }
   local field1 = DataStruct.Field {
      onGrid        = grid,
      numComponents = 3,
      ghost         = {1, 1},
   }
   if speciesRank==0 then field0:clear(0.3) end
   if speciesRank==1 then field1:clear(1.5) end

   local speciesComm_sz = Mpi.Comm_size(speciesComm)
   local num_devices, _ = cuda.GetDeviceCount() 
   local local_rank = speciesRank % num_devices

   -- Picking a GPU based on localRank, allocate device buffers.
   local err = cuda.SetDevice(local_rank)

   -- Get NCCL unique ID at rank 0 and broadcast it to all others.
   local ncclId = Nccl.UniqueId()
   if speciesRank == 0 then
      local _ = Nccl.GetUniqueId(ncclId)
   end
   Mpi.Bcast(ncclId, sizeof(ncclId), Mpi.BYTE, 0, speciesComm)

   -- Create a new CUDA stream (needed by NCCL).
   local custream = cuda.StreamCreate()

   -- Initializing NCCL.
   local devComm = Nccl.Comm()
   local _ = Nccl.CommInitRank(devComm, speciesComm_sz, ncclId, speciesRank)

   -- Communicate data from rank 0 to rank 1.
   if speciesRank == 0 then
      local err = Nccl.Recv(field1:deviceDataPointer(), field1:size(), Nccl.Double, 1, devComm, custream)
      local err = Nccl.Send(field0:deviceDataPointer(), field0:size(), Nccl.Double, 1, devComm, custream)
   end
   if speciesRank == 1 then
      local err = Nccl.Send(field1:deviceDataPointer(), field1:size(), Nccl.Double, 0, devComm, custream)
      local err = Nccl.Recv(field0:deviceDataPointer(), field0:size(), Nccl.Double, 0, devComm, custream)
   end

   -- Completing NCCL operation by synchronizing on the CUDA stream.
   local _ = cuda.StreamSynchronize(custream)

   -- Copy data to host and check results.
   field1:copyDeviceToHost()
   field0:copyDeviceToHost()

   local field = speciesRank==0 and field1 or field0
   local value = speciesRank==0 and 1.5 or 0.3
   local localRange = field:localRange()
   local indexer = field:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      assert_equal(value, fitr[1], "test_4: Checking field value")
      assert_equal(value, fitr[2], "test_4: Checking field value")
      assert_equal(value, fitr[3], "test_4: Checking field value")
   end

   Mpi.Barrier(comm)
end

-- ncclAllReduce a whole CartField.
function test_5(comm)

   local defaultCommSize = 2
   local comm_sz = Mpi.Comm_size(comm)
   if comm_sz ~= defaultCommSize then
      log(string.format("Not running test_4 as numProcs not exactly %d",defaultCommSize))
      return
   end

   local worldRank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   local numCutsConf = 1

   local confColor = math.floor(worldRank/numCutsConf)
   local confComm  = Mpi.Comm_split(Mpi.COMM_WORLD, confColor, worldRank)
   local speciesColor = worldRank % numCutsConf
   local speciesComm  = Mpi.Comm_split(Mpi.COMM_WORLD, speciesColor, worldRank)
   local confRank = Mpi.Comm_rank(confComm)
   local speciesRank = Mpi.Comm_rank(speciesComm)

   local rank = numCutsConf==1 and speciesRank or confRank
   local comm = numCutsConf==1 and speciesComm or confComm

   local decomp = DecompRegionCalc.CartProd { cuts = {numCutsConf}, comm = confComm }
   local grid = Grid.RectCart {
      lower = {0.0},  cells = {12},
      upper = {1.0},  decomposition = decomp,
   }
   local field0 = DataStruct.Field {
      onGrid        = grid,
      numComponents = 3,
      ghost         = {1, 1},
   }
   local field1 = DataStruct.Field {
      onGrid        = grid,
      numComponents = 3,
      ghost         = {1, 1},
   }
   local localRange = field0:localRange()
   local indexer = field0:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field0:get(indexer(i))
      for k = 1, field0:numComponents() do
         fitr[k] = speciesRank*0.3 + k*i
      end
   end
   field0:copyHostToDevice()

   local speciesComm_sz = Mpi.Comm_size(speciesComm)
   local num_devices, _ = cuda.GetDeviceCount() 
   local local_rank = speciesRank % num_devices

   -- Picking a GPU based on localRank, allocate device buffers.
   local err = cuda.SetDevice(local_rank)

   -- Get NCCL unique ID at rank 0 and broadcast it to all others.
   local ncclId = Nccl.UniqueId()
   if speciesRank == 0 then
      local _ = Nccl.GetUniqueId(ncclId)
   end
   Mpi.Bcast(ncclId, sizeof(ncclId), Mpi.BYTE, 0, speciesComm)

   -- Create a new CUDA stream (needed by NCCL).
   local custream = cuda.StreamCreate()

   -- Initializing NCCL.
   local devComm = Nccl.Comm()
   local _ = Nccl.CommInitRank(devComm, speciesComm_sz, ncclId, speciesRank)

   -- Reduce using NCCL.
   local _ Nccl.AllReduce(field0:deviceDataPointer(), field1:deviceDataPointer(),
        field0:size(), Nccl.Double, Nccl.Sum, devComm, custream)

   -- Completing NCCL operation by synchronizing on the CUDA stream.
   local _ = cuda.StreamSynchronize(custream)

   -- Copy data to host and check results.
   field1:copyDeviceToHost()

   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field1:get(indexer(i))
      for k = 1, field0:numComponents() do
         assert_equal(0*0.3 + k*i + 1*0.3 + k*i, fitr[k], "test_5: Checking field value")
      end
   end

   Mpi.Barrier(comm)
end

test_1(Mpi.COMM_WORLD)
test_2(Mpi.COMM_WORLD)
test_3(Mpi.COMM_WORLD)
test_4(Mpi.COMM_WORLD)
test_5(Mpi.COMM_WORLD)


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
