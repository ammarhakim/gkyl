-- Gkyl ------------------------------------------------------------------------
--
-- Commmunication manager, for handling MPI and NCCL communications
-- at the same time.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local xsys             = require "xsys"
local DecompRegionCalc = require "Lib.CartDecomp"
local ffi              = require "ffi"
local lume             = require "Lib.lume"
local Mpi              = require "Comm.Mpi"
local Proto            = require "Lib.Proto"
local sizeof           = xsys.from(ffi, "sizeof")
local xsys             = require "xsys"
local ZeroArray        = require "DataStruct.ZeroArray"
local Lin              = require "Lib.Linalg"
local cuda, Nccl
if GKYL_HAVE_CUDA then
   cuda = require "Cuda.RunTime"
   Nccl = require "Comm.Nccl"
end

local Messenger = Proto()  -- Shell class.

function Messenger:init(tbl)
   -- Create configuration space and species decomposition and communicators.
   local cells         = tbl.cells
   local numSpecies    = tbl.numSpecies
   self.decompCutsConf = tbl.decompCutsConf
   local parallelizeSpecies = xsys.pickBool(tbl.parallelizeSpecies, false)

   local cdim = #cells

   -- Check the configuration space parallel decomposition and
   -- establish the decomposition of the species.
   self.worldSize, self.worldRank = Mpi.Comm_size(Mpi.COMM_WORLD), Mpi.Comm_rank(Mpi.COMM_WORLD)
   local numCutsConf
   if self.decompCutsConf then
      assert(cdim == #self.decompCutsConf, "Messenger: decompCuts should have exactly " .. cdim .. " entries")
      numCutsConf = lume.reduce(self.decompCutsConf, function(a,b) return a*b end)
      -- Calculate the species decomposition.
      self.decompCutsSpecies = numSpecies==1 and 1 or (parallelizeSpecies and self.worldSize/numCutsConf or 1)
      assert(self.worldSize == self.decompCutsSpecies*numCutsConf, "Messenger: decompCuts*decompCutSpecies must equal the number of MPI processes")
   else
      -- If not specified, construct a decompCuts automatically for configuration space.
      if parallelizeSpecies then
         local dofs = {}
         for d = 1, #cells do table.insert(dofs, cells[d]) end
         table.insert(dofs, numSpecies)
         decompCutsWspec = DecompRegionCalc.makeCuts(#dofs, self.worldSize, dofs)
         self.decompCutsConf = {}
         for d = 1, #cells do self.decompCutsConf[d] = decompCutsWspec[d] end
         self.decompCutsSpecies = decompCutsWspec[cdim+1]
      else
         self.decompCutsConf = DecompRegionCalc.makeCuts(cdim, self.worldSize, cells)
         self.decompCutsSpecies = 1
      end
      numCutsConf = lume.reduce(self.decompCutsConf, function(a,b) return a*b end)
   end

   -- Create species and domain communicators.
   local confColor = math.floor(self.worldRank/numCutsConf)
   self.confComm   = Mpi.Comm_split(Mpi.COMM_WORLD, confColor, self.worldRank)
   local speciesColor = self.worldRank % numCutsConf
   self.speciesComm   = Mpi.Comm_split(Mpi.COMM_WORLD, speciesColor, self.worldRank)
   self.confRank, self.speciesRank = Mpi.Comm_rank(self.confComm), Mpi.Comm_rank(self.speciesComm)

   -- Buffers for CartField sync, one for each type of
   -- CartField communicated. They get assigned later.
   self.syncBufs = {}

   -- Initiage GPU comms if needed, and select Allreduce, Send, Isend and Irecv functions.
   if GKYL_USE_GPU then
      self:initCommsDevice()

      self.defaultComms = {world=Mpi.COMM_WORLD, conf=self.confComm_dev, species=self.speciesComm_dev}
      self.reduceOps    = {max = Nccl.Max, min = Nccl.Min, sum = Nccl.Sum}
      self.commTypes    = {double=Nccl.Double, float=Nccl.Float, int=Nccl.Int}

      -- Keep a list of available streams.
      self.twoPointCommStreams = {
        {inUseBy = -1, stream = self.ncclStream}
      }

      self.newRequestStatusFunc = function()
         local reqstat = {srcdest = 0,  tag = -1,  streamIdx = 0,}
         return reqstat
      end

      -- A NCCL group is needed to send multiple messages concurrently, e.g. electrons
      -- sending their moments to ions, and ions sending their moments to electrons.
      self.startCommGroupFunc = function() local ncclErr = Nccl.GroupStart() end

      self.endCommGroupFunc = function() local ncclErr = Nccl.GroupEnd() end

      self.AllreduceFunc = function(dataIn, dataOut, numValuesOut, ncclDataType, ncclOp, comm)
         Nccl.AllReduce(dataIn, dataOut, numValuesOut, ncclDataType, ncclOp, comm, self.ncclStream)
         local _ = Nccl.CommGetAsyncError(comm, self.ncclResult)
         while (self.ncclResult[0] == Nccl.InProgress) do
            local _ = Nccl.CommGetAsyncError(comm, self.ncclResult)
         end
         -- Completing NCCL operation by synchronizing on the CUDA stream.
         local _ = cuda.StreamSynchronize(self.ncclStream)
      end
      
      self.AllreduceByCellFunc = function(fldIn, fldOut, ncclOp, comm)
         Nccl.AllReduce(fldIn:deviceDataPointer(), fldOut:deviceDataPointer(), fldIn:size(),
            self:getCommDataType(fldIn:elemType()), ncclOp, comm, self.ncclStream)
         local _ = Nccl.CommGetAsyncError(comm, self.ncclResult)
         while (self.ncclResult[0] == Nccl.InProgress) do
            local _ = Nccl.CommGetAsyncError(comm, self.ncclResult)
         end
         -- Completing NCCL operation by synchronizing on the CUDA stream.
         local _ = cuda.StreamSynchronize(self.ncclStream)
      end
      
      self.AllgatherFunc = function(fldIn, fldOut, comm)
         Nccl.AllGather(fldIn:deviceDataPointer(), fldOut:deviceDataPointer(), fldIn:size(),
            self:getCommDataType(fldIn:elemType()), comm, self.ncclStream)
         local _ = Nccl.CommGetAsyncError(comm, self.ncclResult)
         while (self.ncclResult[0] == Nccl.InProgress) do
            local _ = Nccl.CommGetAsyncError(comm, self.ncclResult)
         end
         -- Completing NCCL operation by synchronizing on the CUDA stream.
         local _ = cuda.StreamSynchronize(self.ncclStream)
      end

      self.SendCartFieldFunc = function(fld, dest, tag, comm)
         Nccl.Send(fld:deviceDataPointer(), fld:size(), self:getCommDataType(fld:elemType()), dest, comm, self.ncclStream)
      end

      self.IrecvCartFieldFunc = function(fld, src, tag, comm, req)
         local idx, stream
         if req then
            -- Find an available stream.
            local streamTbl
            streamTbl, idx = lume.match(self.twoPointCommStreams, function(e) return (e.inUseBy==-1) or (e.inUseBy==tag) end)
            if idx == nil then
                -- No available streams. Create a new one.
               table.insert(self.twoPointCommStreams, {inUseBy=tag, stream=cuda.StreamCreate()})
               idx, stream = #self.twoPointCommStreams, self.twoPointCommStreams[#self.twoPointCommStreams].stream
            else
               stream = streamTbl.stream
            end
            req.srcdest, req.tag, req.streamIdx = src, tag, idx
         else
            idx, stream = 1, self.ncclStream
         end

         Nccl.Recv(fld:deviceDataPointer(), fld:size(), self:getCommDataType(fld:elemType()), src, comm, self.ncclStream)

         self.twoPointCommStreams[idx].inUseBy = tag
      end

      self.IsendCartFieldFunc = function(fld, dest, tag, comm, req)
         local idx, stream
         if req then
            -- Find an available stream.
            local streamTbl
            streamTbl, idx = lume.match(self.twoPointCommStreams, function(e) return (e.inUseBy==-1) or (e.inUseBy==tag) end)
            if idx == nil then
                -- No available streams. Create a new one.
               table.insert(self.twoPointCommStreams, {inUseBy=tag, stream=cuda.StreamCreate()})
               idx, stream = #self.twoPointCommStreams, self.twoPointCommStreams[#self.twoPointCommStreams].stream
            else
               stream = streamTbl.stream
            end
            req.srcdest, req.tag, req.streamIdx = src, tag, idx
         else
            idx, stream = 1, self.ncclStream
         end

         Nccl.Send(fld:deviceDataPointer(), fld:size(), self:getCommDataType(fld:elemType()), dest, comm, self.ncclStream)

         self.twoPointCommStreams[idx].inUseBy = tag
      end

      self.WaitFunc = function(req, stat, comm)
         local _ = Nccl.CommGetAsyncError(comm, self.ncclResult)
         while (self.ncclResult[0] == Nccl.InProgress) do
            local _ = Nccl.CommGetAsyncError(comm, self.ncclResult)
         end
         -- Completing NCCL operation by synchronizing on the CUDA stream.
         local stream
         if req then
            stream = self.twoPointCommStreams[req.streamIdx].stream
            self.twoPointCommStreams[req.streamIdx].inUseBy = -1  -- Free stream for use by others.
         else
            stream = self.ncclStream
         end
         local _ = cuda.StreamSynchronize(self.ncclStream)
      end

      self.syncCartFieldFunc = function(fld, comm)
         Messenger["syncCartFieldNCCL"](self, fld, comm)
      end

      self.syncPeriodicCartFieldFunc = function(fld, comm, dirs)
         Messenger["syncPeriodicCartFieldNCCL"](self, fld, comm, dirs)
      end
   else
      self.defaultComms = {world=Mpi.COMM_WORLD, conf=self.confComm, species=self.speciesComm}
      self.reduceOps    = {max = Mpi.MAX, min = Mpi.MIN, sum = Mpi.SUM}
      self.commTypes    = {double=Mpi.DOUBLE, float=Mpi.FLOAT, int=Mpi.INT}

      self.newRequestStatusFunc = function()
         local reqstat = Mpi.RequestStatus()
         return reqstat
      end

      self.startCommGroupFunc = function() end

      self.endCommGroupFunc = function() end

      self.AllgatherFunc = function(fldIn, fldOut, comm)
         Mpi.Allgather(fldIn:dataPointer(), fldIn:size(), fldIn:elemCommType(),
                       fldOut:dataPointer(), fldIn:size(), fldOut:elemCommType(), comm)
      end
      self.AllreduceFunc = function(dataIn, dataOut, numValuesOut, mpiDataType, mpiOp, comm)
         Mpi.Allreduce(dataIn, dataOut, numValuesOut, mpiDataType, mpiOp, comm)
      end
      self.AllreduceByCellFunc = function(fldIn, fldOut, mpiOp, comm)
         Mpi.Allreduce(fldIn:dataPointer(), fldOut:dataPointer(),
            fldOut:size(), fldOut:elemCommType(), mpiOp, comm)
      end
      self.SendCartFieldFunc = function(fld, dest, tag, comm)
         Mpi.Send(fld:dataPointer(), fld:size(), fld:elemCommType(), dest, tag, comm)
      end

      self.IrecvCartFieldFunc = function(fld, src, tag, comm, req)
         Mpi.Irecv(fld:dataPointer(), fld:size(), fld:elemCommType(), src, tag, comm, req)
      end
      self.IsendCartFieldFunc = function(fld, dest, tag, comm, req)
         Mpi.Isend(fld:dataPointer(), fld:size(), fld:elemCommType(), dest, tag, comm, req)
      end

      self.WaitFunc = function(req, stat, comm)
         local _ = Mpi.Wait(req, stat)
      end

      self.syncReqStat = Mpi.RequestStatus()

      self.syncCartFieldFunc = function(fld, comm)
         Messenger["syncCartFieldMPI"](self, fld, comm)
      end

      self.syncPeriodicCartFieldFunc = function(fld, comm, dirs)
         Messenger["syncPeriodicCartFieldMPI"](self, fld, comm, dirs)
      end
   end

   -- Configuration space decomp object.
   self.decompConf = DecompRegionCalc.CartProd {
      cuts  = self.decompCutsConf,
      comms = {host=self.confComm, device=self.confComm_dev},
   }

end

function Messenger:newCommDevice(commIn)
   local rank, commSz = Mpi.Comm_rank(commIn), Mpi.Comm_size(commIn)

   -- Get NCCL unique ID at rank 0 and broadcast it to all others.
   local ncclId = Nccl.UniqueId()
   if rank == 0 then local _ = Nccl.GetUniqueId(ncclId) end
   Mpi.Bcast(ncclId, sizeof(ncclId), Mpi.BYTE, 0, commIn)

   -- NCCL Result object needed to query status of a communicator.
   local ncclResult, commConfig, commIn_dev = Nccl.Result(), Nccl.Config(), Nccl.Comm()
   commConfig[0].blocking = 0;   -- Nonblocking comm.
   local _ = Nccl.CommInitRankConfig(commIn_dev, commSz, ncclId, rank, commConfig)
   local _ = Nccl.CommGetAsyncError(commIn_dev, ncclResult)
   while (ncclResult[0] == Nccl.InProgress) do
      local _ = Nccl.CommGetAsyncError(commIn_dev, ncclResult)
   end

   return commIn_dev, ncclId
end

function Messenger:initCommsDevice()
   -- If Cuda/NCCL are present create device communicators for conf/species comms.
   local numDevices, _ = cuda.GetDeviceCount()
   -- Assume number of devices must be the same as number of MPI processes.
   local deviceRank = self.worldRank % numDevices
--   -- MF 2023/02/04: This assert causes trouble. Logic may be wrong for multi-node runs.
--   assert(self.worldSize == numDevices, string.format("Messenger: number of GPUs (%d) must match number of MPI processes (%d).", numDevices, self.worldSize))
   local _ = cuda.SetDevice(deviceRank)   -- Picking a GPU based on localRank

   -- Get NCCL unique ID at rank 0 and broadcast it to all others.
   self.ncclIdConf = Nccl.UniqueId()
   if self.confRank == 0 then local _ = Nccl.GetUniqueId(self.ncclIdConf) end
   Mpi.Bcast(self.ncclIdConf, sizeof(self.ncclIdConf), Mpi.BYTE, 0, self.confComm)

   self.ncclIdSpecies = Nccl.UniqueId()
   if self.speciesRank == 0 then local _ = Nccl.GetUniqueId(self.ncclIdSpecies) end
   Mpi.Bcast(self.ncclIdSpecies, sizeof(self.ncclIdSpecies), Mpi.BYTE, 0, self.speciesComm)
   
   -- NCCL Result object needed to query status of a communicator.
   self.ncclResult = Nccl.Result()

   local confComm_sz    = Mpi.Comm_size(self.confComm)
   local speciesComm_sz = Mpi.Comm_size(self.speciesComm)
   
   -- Initializing NCCL comms.
   local commConfig = Nccl.Config()
   commConfig[0].blocking = 0;   -- Nonblocking comm.

   self.confComm_dev = Nccl.Comm()
   local _ = Nccl.CommInitRankConfig(self.confComm_dev, confComm_sz, self.ncclIdConf, self.confRank, commConfig)
   local _ = Nccl.CommGetAsyncError(self.confComm_dev, self.ncclResult)
   while (self.ncclResult[0] == Nccl.InProgress) do
      local _ = Nccl.CommGetAsyncError(self.confComm_dev, self.ncclResult)
   end

   self.speciesComm_dev = Nccl.Comm()
   local _ = Nccl.CommInitRankConfig(self.speciesComm_dev, speciesComm_sz, self.ncclIdSpecies, self.speciesRank, commConfig)
   local _ = Nccl.CommGetAsyncError(self.speciesComm_dev, self.ncclResult)
   while (self.ncclResult[0] == Nccl.InProgress) do
      local _ = Nccl.CommGetAsyncError(self.speciesComm_dev, self.ncclResult)
   end

   -- Create a new CUDA stream (needed by NCCL).
   self.ncclStream = cuda.StreamCreate()
end

function Messenger:createSubComms(grid)
   -- Create subcommunicators along some dimensions.
   local decompRange  = grid:decomposedRange()
   local cutsIdxr     = decompRange:cutsIndexer()
   local cutsRange    = decompRange:cutsRange()
   local subdomainIdx = {} -- Grid ID on this processor.
   decompRange:cutsInvIndexer()(grid:subGridId(), subdomainIdx)
   local confComm = self:getConfComm_host()
   local group = Mpi.Comm_group(confComm)

   local ndim = grid:ndim()
   -- Create comm along the last dimension.
   local zNumRanks = decompRange:cuts(ndim)
   local zCommGroupRanks = Lin.IntVec(zNumRanks)
   local j = 0
   for idx in cutsRange:colMajorIter() do  -- MPI processes are column-major ordered.
      local sameXY = true
      for d = ndim-1,1,-1 do sameXY = sameXY and idx[d] == subdomainIdx[d] end
      if sameXY then
         j = j+1;  zCommGroupRanks[j] = cutsIdxr(idx)-1
      end
   end
   local zGroup   = Mpi.Group_incl(group, zNumRanks, zCommGroupRanks:data());
   local tag = ndim > 1 and subdomainIdx[1]-1 or 0
   for d = ndim-2,1,-1 do 
      tag = tag + (subdomainIdx[d+1]-1)*decompRange:cuts(d)
   end
   self.confCommZ = Mpi.Comm_create_group(confComm, zGroup, tag);
   Mpi.Group_free(zGroup)

   -- Create a x-y communicator.
   local xyNumRanks = 1
   for d = ndim-1,1,-1 do xyNumRanks = xyNumRanks*decompRange:cuts(d) end
   local xyCommGroupRanks = Lin.IntVec(xyNumRanks)
   local j = 0
   for idx in cutsRange:colMajorIter() do  -- MPI processes are column-major ordered.
      if idx[ndim] == subdomainIdx[ndim] then
         j = j+1;  xyCommGroupRanks[j] = cutsIdxr(idx)-1
      end
   end
   local xyGroup   = Mpi.Group_incl(group, xyNumRanks, xyCommGroupRanks:data());
   local tag       = subdomainIdx[ndim]-1
   self.confCommXY = Mpi.Comm_create_group(confComm, xyGroup, tag);
   Mpi.Group_free(xyGroup)

   self.subComms = {z = self.confCommZ, xy = self.confCommXY}

   if GKYL_USE_GPU then
      -- Create NCCL comms.
      self.confCommZ_dev, self.ncclIdConfZ   = self:newNCCLcomm(self.confCommZ)
      self.confCommXY_dev, self.ncclIdConfXY = self:newNCCLcomm(self.confCommXY)
      self.defaultComms["z"]  = self.confCommZ_dev
      self.defaultComms["xy"] = self.confCommXY_dev
      self.subComms_dev       = {z = self.confCommZ_dev, xy = self.confCommXY_dev}
      self.defaultSubComms    = self.subComms_dev
   else
      self.defaultComms["z"]  = self.confCommZ
      self.defaultComms["xy"] = self.confCommXY
      self.defaultSubComms    = self.subComms
   end

   Mpi.Group_free(group)
end

function Messenger:newNCCLcomm(mpiComm)
   local commConfig = Nccl.Config()
   commConfig[0].blocking = 0;   -- Nonblocking comm.

   -- Get NCCL unique ID at rank 0 and broadcast it to all others.
   local ncclId = Nccl.UniqueId()
   local mpiRank, mpiCommSize = Mpi.Comm_rank(mpiComm), Mpi.Comm_size(mpiComm)
   if mpiRank == 0 then local _ = Nccl.GetUniqueId(ncclId) end
   Mpi.Bcast(ncclId, sizeof(ncclId), Mpi.BYTE, 0, mpiComm)

   local ncclComm = Nccl.Comm()
   local _ = Nccl.CommInitRankConfig(ncclComm, mpiCommSize, ncclId, mpiRank, commConfig)
   local _ = Nccl.CommGetAsyncError(ncclComm, self.ncclResult)
   while (self.ncclResult[0] == Nccl.InProgress) do
      local _ = Nccl.CommGetAsyncError(ncclComm, self.ncclResult)
   end

   return ncclComm, ncclId
end

function Messenger:chooseSyncBuf(fld)
   local nc = fld:numComponents()
   local bufVol = {recv=fld._syncRecvBufVol, send=fld._syncSendBufVol}
   local bufIdx = nil
   -- Look through existing buffers and see if one
   -- has the right shape/size. Otherwise create one.
   for i = 1, #self.syncBufs do
      local recvnc, sendnc = self.syncBufs[i].recv:get_ncomp(),  self.syncBufs[i].send:get_ncomp()
      local recvsz, sendsz = self.syncBufs[i].recv:get_size(),  self.syncBufs[i].send:get_size()
      if recvnc>=nc and sendnc>=nc and recvsz>=bufVol.recv and sendsz>=bufVol.send then
         bufIdx = i
         break
      end
   end
   if bufIdx == nil then
      bufIdx = #self.syncBufs+1
      self.syncBufs[bufIdx] = {recv=ZeroArray.Array(ZeroArray.double, nc, bufVol.recv, fld.useDevice),
                               send=ZeroArray.Array(ZeroArray.double, nc, bufVol.send, fld.useDevice),} 
   end
   return bufIdx
end

function Messenger:syncCartFieldMPI(fld, comm)
   -- No need to do anything if communicator is not valid.
   if not Mpi.Is_comm_valid(comm) then return end
   -- Immediately return if nothing to sync.
   if fld._lowerGhost == 0 and fld._upperGhost == 0 then return end
   -- Immediately return if there is no decomp.
   if fld:grid():decomposedRange():numSubDomains() == 1 then return end

   local nc = fld:numComponents()
   local bufIdx = self:chooseSyncBuf(fld)
   local recvBuf, sendBuf = self.syncBufs[bufIdx].recv, self.syncBufs[bufIdx].send
   local recvPtr, sendPtr = recvBuf:data(), sendBuf:data()
   local dataType = self:getCommDataType(fld:elemType())
   local myId     = fld:grid():subGridId() -- Grid ID on this processor.
   local neigIds  = fld._decompNeigh:neighborData(myId) -- List of neighbors.
   local tag      = 42 -- Communicator tag for ghost (non-periodic) messages.
   for _, rank in ipairs(neigIds) do
      local recvRgn = fld._syncRecvRng[rank]
      local recvSz  = nc*recvRgn:volume()
      -- Post recv for expected data into the recv buffer.
      local _ = Mpi.Irecv(recvPtr, recvSz, dataType, rank-1, tag, comm, self.syncReqStat)
      -- Copy data to send buffer, and send it.
      local sendRgn = fld._syncSendRng[rank]
      local sendSz  = nc*sendRgn:volume()
      fld:copyRangeToBuffer(sendRgn, sendBuf:data()) -- copy from skin cells.
      Mpi.Send(sendPtr, sendSz, dataType, rank-1, tag, comm)
      -- Complete the recv and copy data from buffer to field.
      local _ = Mpi.Wait(self.syncReqStat, self.syncReqStat)
      fld:copyRangeFromBuffer(recvRgn, recvBuf:data()) -- copy data into ghost cells.
   end
end

function Messenger:syncCartFieldNCCL(fld, comm)
   -- Immediately return if nothing to sync.
   if fld._lowerGhost == 0 and fld._upperGhost == 0 then return end
   -- Immediately return if there is no decomp.
   if fld:grid():decomposedRange():numSubDomains() == 1 then return end

   local nc = fld:numComponents()
   local bufIdx = self:chooseSyncBuf(fld)
   local recvBuf, sendBuf = self.syncBufs[bufIdx].recv, self.syncBufs[bufIdx].send
   local recvPtr, sendPtr = recvBuf:data(), sendBuf:data()
   local dataType = self:getCommDataType(fld:elemType())
   local myId     = fld:grid():subGridId() -- Grid ID on this processor.
   local neigIds  = fld._decompNeigh:neighborData(myId) -- List of neighbors.
   for _, rank in ipairs(neigIds) do
      local isLo = myId < rank
      for sI = 0, 1 do
         if isLo then
            -- Post recv for expected data into the recv buffer.
            local recvRgn = fld._syncRecvRng[rank]
            local recvSz  = nc*recvRgn:volume()
            Nccl.Recv(recvPtr, recvSz, dataType, rank-1, comm, self.ncclStream)
            -- Complete the recv and copy data from buffer to field.
            self:Wait(req, stat, comm)
            fld:copyRangeFromBuffer(recvRgn, recvBuf:data()) -- copy data into ghost cells.
         else
            -- Copy data to send buffer, and send it.
            local sendRgn = fld._syncSendRng[rank]
            local sendSz  = nc*sendRgn:volume()
            fld:copyRangeToBuffer(sendRgn, sendBuf:data()) -- copy from skin cells.
            Nccl.Send(sendPtr, sendSz, dataType, rank-1, comm, self.ncclStream)
            -- Complete the send.
            self:Wait(req, stat, comm)
         end
	 isLo = not isLo
      end
   end
end

function Messenger:syncPeriodicCartFieldMPI(fld, comm, dirs)
   -- No need to do anything if communicator is not valid.
   if not Mpi.Is_comm_valid(comm) then return end
   -- Immediately return if nothing to sync.
   if fld._lowerGhost == 0 and fld._upperGhost == 0 then return end
   -- If there is no decomp do copies, not comm.
   if fld:grid():decomposedRange():numSubDomains() == 1 and not fld._syncCorners then
      return fld:periodicCopy()
   end

   local nc = fld:numComponents()
   local bufIdx = self:chooseSyncBuf(fld)
   local recvBuf, sendBuf = self.syncBufs[bufIdx].recv, self.syncBufs[bufIdx].send
   local recvPtr, sendPtr = recvBuf:data(), sendBuf:data()
   local dataType = self:getCommDataType(fld:elemType())
   local tag      = 32 -- Communicator tag for periodic messages.
   local syncDirs = dirs or fld:grid():getPeriodicDirs()
   local onBoundary = fld._onBoundary
   for _, dir in ipairs(syncDirs) do
      if onBoundary[dir] then
         if fld:grid():cuts(dir) == 1 then
            fld:periodicCopyInDir(dir)
         else
            local rank    = fld._syncPerNeigh[dir]
            local recvRgn = fld._syncPerRecvRng[dir]
            local recvSz  = nc*recvRgn:volume()
            -- Post recv for expected data into the recv buffer.
            local _ = Mpi.Irecv(recvPtr, recvSz, dataType, rank-1, tag, comm, self.syncReqStat)
            -- Copy data to send buffer, and send it.
            local sendRgn = fld._syncPerSendRng[dir]
            local sendSz  = nc*sendRgn:volume()
            fld:copyRangeToBuffer(sendRgn, sendBuf:data()) -- copy from skin cells.
            Mpi.Send(sendPtr, sendSz, dataType, rank-1, tag, comm)
            -- Complete the recv and copy data from buffer to field.
            local _ = Mpi.Wait(self.syncReqStat, self.syncReqStat)
            fld:copyRangeFromBuffer(recvRgn, recvBuf:data()) -- copy data into ghost cells.
         end
      end
   end
end

function Messenger:syncPeriodicCartFieldNCCL(fld, comm, dirs)
   -- Immediately return if nothing to sync.
   if fld._lowerGhost == 0 and fld._upperGhost == 0 then return end
   -- If there is no decomp do copies, not comm.
   if fld:grid():decomposedRange():numSubDomains() == 1 and not fld._syncCorners then
      return fld:periodicCopy()
   end

   local nc = fld:numComponents()
   local bufIdx = self:chooseSyncBuf(fld)
   local recvBuf, sendBuf = self.syncBufs[bufIdx].recv, self.syncBufs[bufIdx].send
   local recvPtr, sendPtr = recvBuf:data(), sendBuf:data()
   local dataType = self:getCommDataType(fld:elemType())
   local syncDirs = dirs or fld:grid():getPeriodicDirs()
   local myId     = fld:grid():subGridId() -- Grid ID on this processor.
   local onBoundary = fld._onBoundary
   for _, dir in ipairs(syncDirs) do
      if onBoundary[dir] then
         if fld:grid():cuts(dir) == 1 then
            fld:periodicCopyInDir(dir)
         else
            local rank = fld._syncPerNeigh[dir]
            local isLo = myId < rank

            for sI = 0, 1 do
               if isLo then
                  local recvRgn = fld._syncPerRecvRng[dir]
                  local recvSz  = nc*recvRgn:volume()
                  -- Post recv for expected data into the recv buffer.
                  Nccl.Recv(recvPtr, recvSz, dataType, rank-1, comm, self.ncclStream)
                  -- Complete the recv and copy data from buffer to field.
                  self:Wait(req, stat, comm)
                  fld:copyRangeFromBuffer(recvRgn, recvBuf:data()) -- copy data into ghost cells.
               else
                  -- Copy data to send buffer, and send it.
                  local sendRgn = fld._syncPerSendRng[dir]
                  local sendSz  = nc*sendRgn:volume()
                  fld:copyRangeToBuffer(sendRgn, sendBuf:data()) -- copy from skin cells.
                  Nccl.Send(sendPtr, sendSz, dataType, rank-1, comm, self.ncclStream)
                  -- Complete the send.
                  self:Wait(req, stat, comm)
               end
               isLo = not isLo
            end
         end
      end
   end
end

function Messenger:getComms()       return self.defaultComms end
function Messenger:getWorldComm()   return self.defaultComms["world"] end
function Messenger:getConfComm()    return self.defaultComms["conf"] end
function Messenger:getSpeciesComm() return self.defaultComms["species"] end
function Messenger:getSubComms()    return self.defaultSubComms end

function Messenger:getConfComm_host()    return self.confComm end
function Messenger:getSpeciesComm_host() return self.speciesComm end
function Messenger:getSubComms_host()       return self.subComms end
function Messenger:getConfComm_device()    return self.confComm_dev end
function Messenger:getSpeciesComm_device() return self.speciesComm_dev end
function Messenger:getSubComms_device()       return self.subComms_dev end

function Messenger:getRanks() return {world=self.worldRank, conf=self.confRank, species=self.speciesRank} end
function Messenger:getWorldRank() return self.worldRank end
function Messenger:getConfRank() return self.confRank end
function Messenger:getSpeciesRank() return self.speciesRank end

function Messenger:getConfDecompCuts() return self.decompCutsConf end
function Messenger:getSpeciesDecompCuts() return self.decompCutsSpecies end

function Messenger:getConfDecomp() return self.decompConf end

function Messenger:getCommDataType(dataType)
   local dataCommType = type(dataType)=="string" and self.commTypes[dataType]
                       or (ffi.istype("double",dataType) and self.commTypes["double"] or
                           (ffi.istype("int",dataType) and self.commTypes["int"] or
                            (ffi.istype("float",dataType) and self.commTypes["float"] or
		             assert(false,"Messenger: type not supported"))))
   return dataCommType
end

function Messenger:newRequestStatus()
   return self.newRequestStatusFunc()
end

function Messenger:startCommGroup()
   return self.startCommGroupFunc()
end

function Messenger:endCommGroup()
   return self.endCommGroupFunc()
end

function Messenger:Allreduce(dataIn, dataOut, numValuesOut, dataType, op, comm)
   self.AllreduceFunc(dataIn, dataOut, numValuesOut, self.commTypes[dataType], self.reduceOps[op], comm)
end

function Messenger:AllreduceByCell(fieldIn, fieldOut, op, comm)
   self.AllreduceByCellFunc(fieldIn, fieldOut, self.reduceOps[op], comm)
end

function Messenger:Allgather(localField, globalField, comm)
   self.AllgatherFunc(localField, globalField, comm)
end

function Messenger:SendCartField(fld, dest, tag, comm)
   self.SendCartFieldFunc(fld, dest, tag, comm)
end

function Messenger:IsendCartField(fld, dest, tag, comm, req)
   self.IsendCartFieldFunc(fld, dest, tag, comm, req)
end

function Messenger:IrecvCartField(fld, src, tag, comm, req)
   self.IrecvCartFieldFunc(fld, src, tag, comm, req)
end

function Messenger:Wait(req, stat, comm) self.WaitFunc(req, stat, comm) end

function Messenger:syncCartField(fld, comm)
   self.syncCartFieldFunc(fld, comm)
end

function Messenger:syncPeriodicCartField(fld, comm, dirs)
   self.syncPeriodicCartFieldFunc(fld, comm, dirs)
end

function Messenger:func() end

return Messenger
