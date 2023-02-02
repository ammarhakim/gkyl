-- Gkyl ------------------------------------------------------------------------
--
-- Multi-component fields on cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

-- Gkyl libraries.
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Alloc            = require "Lib.Alloc"
local CartDecompNeigh  = require "Lib.CartDecompNeigh"
local Grid             = require "Grid.RectCart"
local Lin              = require "Lib.Linalg"
local LinearDecomp     = require "Lib.LinearDecomp"
local Mpi              = require "Comm.Mpi"
local Range            = require "Lib.Range"
local ZeroArray        = require "DataStruct.ZeroArray"
local lume             = require "Lib.lume"

-- Load CUDA allocators (or dummy when CUDA is not found).
local cuda = nil
local cuAlloc = require "Cuda.AllocDummy"
if GKYL_HAVE_CUDA then
   cuAlloc = require "Cuda.Alloc"
   cuda    = require "Cuda.RunTime"
end

local RangeVec = Lin.new_vec_ct(ffi.typeof("struct gkyl_range"))

-- C interfaces
ffi.cdef [[
    // s: start index. nv: number of values.
    void gkylCartFieldAbs(unsigned s, unsigned nv, double *out);

/**
 * Create ghost and skin sub-ranges given parent range. The skin and
 * ghost ranges are sub-ranges of the parent range and DO NOT include
 * corners.
 *
 * @param skin On output, skin range
 * @param ghost On outout, ghost range
 * @param dir Direction in which skin/ghost are computed
 * @param edge Edge on which skin/ghost are computed
 * @param parent Range for which skin/ghost are computed
 * @param nghost Number of ghost cells in 'dir' are nghost[dir]
 */
void gkyl_skin_ghost_ranges(struct gkyl_range *skin, struct gkyl_range *ghost,
  int dir, int edge, const struct gkyl_range *parent, const int *nghost);
]]

-- Local definitions.
local rowMajLayout, colMajLayout = Range.rowMajor, Range.colMajor -- data layout
local indexerMakerFuncs = {} -- list of functions that make indexers
indexerMakerFuncs[rowMajLayout] = Range.makeRowMajorIndexer
indexerMakerFuncs[colMajLayout] = Range.makeColMajorIndexer
-- Default layout.
local defaultLayout = rowMajLayout

local genIndexerMakerFuncs = {} -- List of functions that make generic indexers.
genIndexerMakerFuncs[rowMajLayout] = Range.makeRowMajorGenIndexer
genIndexerMakerFuncs[colMajLayout] = Range.makeColMajorGenIndexer

-- Helper function to check if two fields are compatible.
local function field_compatible(y, x)
   return y:localRange() == x:localRange() and y:numComponents() == x:numComponents()
end
local function field_compatible_ext(y, x)
   return y:localExtRange() == x:localExtRange() and y:numComponents() == x:numComponents()
end
-- Helper function to check if two fields have the same range.
-- Useful when manipulating two fields with different number of components, but same range.
local function field_check_range(y, x)
   return y:localRange() == x:localRange()
end

-- Turn numbers in a table into strings and concatenate them.
local tblToStr = function(tblIn)
   local strOut = ""
   for _, v in ipairs(tblIn) do 
      strOut = v<0 and (strOut .. math.abs(v) .. 1) or (strOut .. math.abs(v) .. 2) 
   end
   return strOut
end

-- Field accessor object: allows access to field values in cell.
local function new_field_comp_ct(elct)
   local field_comp_mf = {
      data = function(self)
	 return self._cdata
      end
   }
   local field_comp_mt = {
      __index = function(self, k)
	 if type(k) == "number" then
	    return self._cdata[k-1]
	 else
	    return field_comp_mf[k]
	 end
      end,
      __newindex = function(self, k, v)
	 self._cdata[k-1] = v
      end,
   }
   return metatype(typeof("struct { int numComponents; $* _cdata; }", elct), field_comp_mt)
end


-- A function to create constructors for Field objects
local function Field_meta_ctor(elct)
   -- ctor for component data
   local fcompct = new_field_comp_ct(elct)
   -- ctor for creating vector of element types
   local ElemVec = Lin.new_vec_ct(elct)   

   local elctSize = sizeof(elct)
   local elctMinValue, elctMaxValue = 0, 0
   
   -- Meta-data for type
   local isNumberType = false   
   local elctCommType = nil
   if ffi.istype(new(elct), new("double")) then
      elctCommType = Mpi.DOUBLE
      isNumberType = true
      elctMinValue, elctMaxValue = GKYL_MIN_DOUBLE, GKYL_MAX_DOUBLE
   elseif ffi.istype(new(elct), new("float")) then
      elctCommType = Mpi.FLOAT
      isNumberType = true
      elctMinValue, elctMaxValue = GKYL_MIN_FLOAT, GKYL_MAX_FLOAT
   elseif ffi.istype(new(elct), new("int")) then
      elctCommType = Mpi.INT
      isNumberType = true
      elctMinValue, elctMaxValue = GKYL_MIN_INT, GKYL_MAX_INT
   elseif ffi.istype(new(elct), new("long")) then
      elctCommType = Mpi.LONG
      isNumberType = true
      elctMinValue, elctMaxValue = GKYL_MIN_LONG, GKYL_MAX_LONG
   else
      elctCommType = Mpi.BYTE -- by default, send stuff as byte array
   end

   -- functions for regular and device memory allocations
   local allocFunc = Alloc.Alloc_meta_ctor(elct)
   local allocCudaFunc = cuAlloc.Alloc_meta_ctor(elct, false) -- don't used managed memory
   
   -- allocator for use in applications
   local function allocatorFunc(comm, numElem)
      return allocFunc(numElem)
   end
   -- allocator for use in memory on device
   local function deviceAllocatorFunc(comm, numElem)
      return allocCudaFunc(numElem)
   end

   -- Binary operation functions and reduce MPI types (used in reduce method).
   local binOpFuncs = {
      max = function(a,b) return math.max(a,b) end,
      min = function(a,b) return math.min(a,b) end,
      sum = function(a,b) return a+b end
   }
   local binOpFlags = {min = 1, max = 2, sum = 3}
   local reduceOpsMPI = {max = Mpi.MAX, min = Mpi.MIN, sum = Mpi.SUM}
   local reduceInitialVal = {max = elctMinValue, min = elctMaxValue , sum = 0.0}
   
   -- Make constructor for Field.
   local Field = {}
   function Field:new(tbl)
      local self = setmetatable({}, Field)

      -- Read data from input table.
      local grid  = tbl.onGrid
      local nc    = tbl.numComponents and tbl.numComponents or 1 -- Default numComponents=1.
      local ghost = tbl.ghost and tbl.ghost or {0, 0} -- No ghost cells by default.

      self._syncCorners      = xsys.pickBool(tbl.syncCorners, false) -- Don't sync corners by default.
      self._syncPeriodicDirs = xsys.pickBool(tbl.syncPeriodicDirs, true) -- Sync periodic BCs by default.

      -- Local and global ranges.
      local globalRange = grid:globalRange()
      local localRange  = grid:localRange()

      local nodeComm = grid:commSet().nodeComm

      -- Allocator function.
      local allocator = allocatorFunc
      
      local sz = localRange:extend(ghost[1], ghost[2]):volume()*nc -- Amount of data in field.
      -- Setup object.
      self._grid = grid
      self._ndim = grid:ndim()
      self._lowerGhost, self._upperGhost = ghost[1], ghost[2]
      self._numComponents = nc
      self._size = sz

      -- Underlying ZeroArray data structure, which handles memory allocation
      self._zero = ZeroArray.Array(ZeroArray.double, self._numComponents, self._size/self._numComponents, 0)
      -- get pointer to data in ZeroArray
      self._data = self._zero:data()

      -- Create a device copy if needed.
      if xsys.pickBool(tbl.useDevice, GKYL_USE_GPU) then
         self.useDevice = 1
         self._zeroDevice = ZeroArray.Array(ZeroArray.double, self._numComponents, self._size/self._numComponents, 1)
         self._devAllocData = self._zeroDevice:data()

         self._zeroForOps = self._zeroDevice
      else
         self._zeroDevice = nil
         self._devAllocData = nil

         self._zeroForOps = self._zero
      end

      -- Extended ranges include ghost cells.
      self._globalExtRange = globalRange:extend(self._lowerGhost, self._upperGhost)
      self._localExtRange  = localRange:extend(self._lowerGhost, self._upperGhost)

      -- re-initialize localRange and globalRange as sub-ranges of localExtRange and globalExtRange
      self._localRange  = self._localExtRange:subRange(localRange:lowerAsVec(), localRange:upperAsVec())
      self._globalRange = self._globalExtRange:subRange(globalRange:lowerAsVec(), globalRange:upperAsVec())

      self._skinRgnUpper = RangeVec(self._ndim)
      self._skinRgnLower = RangeVec(self._ndim)
      self._ghostRgnUpper = RangeVec(self._ndim)
      self._ghostRgnLower = RangeVec(self._ndim)
      self._lowerGhostVec = Lin.IntVec(self._ndim)
      self._upperGhostVec = Lin.IntVec(self._ndim)
      local lowerEdge, upperEdge = 0, 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
      for d=1, self._ndim do
         self._lowerGhostVec[d] = self._lowerGhost
         self._upperGhostVec[d] = self._upperGhost
      end
      for d=1, self._ndim do
         ffiC.gkyl_skin_ghost_ranges(self._skinRgnLower[d], self._ghostRgnLower[d], d-1, lowerEdge, self._localExtRange, self._lowerGhostVec:data())
         ffiC.gkyl_skin_ghost_ranges(self._skinRgnUpper[d], self._ghostRgnUpper[d], d-1, upperEdge, self._localExtRange, self._upperGhostVec:data())
      end

      -- All real-cell edges.
      self._localEdgeRange = self._localRange:extend(1, 0) -- Or (1, 0)?

      -- All cell-cell edges, including those of a ghost cell.
      self._localExtEdgeRange = self._localRange:extend(
	 self._lowerGhost-1, self._upperGhost)

      -- Local and (MPI) global values of a reduction (reduce method).
      if self.useDevice then
	 local on_gpu = 1
         self.localReductionVal = ZeroArray.Array(ZeroArray.double, self._numComponents, 1, on_gpu)
         self.globalReductionVal = ZeroArray.Array(ZeroArray.double, self._numComponents, 1)
         self.localReductionVal_h = ZeroArray.Array(ZeroArray.double, self._numComponents, 1)
      else
         self.localReductionVal = ZeroArray.Array(ZeroArray.double, self._numComponents, 1)
         self.globalReductionVal = ZeroArray.Array(ZeroArray.double, self._numComponents, 1)
	 self.localReductionVal_h = self.localReductionVal
      end
      
      self._layout = defaultLayout -- Default layout is column-major.
      if tbl.layout then
         self._layout = tbl.layout=="row-major" and rowMajLayout or colMajLayout
      end

      -- Store start index and size handled by local rank for local and extended range.
      self._localStartIdx, self._localNumBump       = self._localRange:lowerAsVec(), self._localRange:volume()
      self._localExtStartIdx, self._localExtNumBump = self._localExtRange:lowerAsVec(), self._localExtRange:volume()

      -- Compute communication neighbors.
      self._decompNeigh = CartDecompNeigh(grid:decomposedRange())
      if self._syncCorners then
         self._decompNeigh:calcAllCommNeigh(ghost[1], ghost[2])
      else
         self._decompNeigh:calcFaceCommNeigh(ghost[1], ghost[2])
      end

      -- Pre-create MPI DataTypes for send/recv calls when doing ghost-cell
      -- sync(). Using MPI DataTypes, we do not require temporary buffers
      -- for send/recv.
      -- Also pre-create the location in memory required so that we know 
      -- what parts of the data structure are being sent and received to.
      self._sendMPIDataType, self._recvMPIDataType = {}, {}
      self._sendMPILoc, self._recvMPILoc = {}, {}
      self._recvMPIReq = {}
      local localExtRange   = self._localExtRange
      local indexer         = self:genIndexer()
      local decomposedRange = self._grid:decomposedRange()
      local myId            = self._grid:subGridId() -- Grid ID on this processor.
      local neigIds         = self._decompNeigh:neighborData(myId) -- List of neighbors.

      for _, sendId in ipairs(neigIds) do
         local neighRgn = decomposedRange:subDomain(sendId)
         local sendRgn  = self._localRange:intersect(
            neighRgn:extend(self._lowerGhost, self._upperGhost))
         local idx      = sendRgn:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self._sendMPILoc[sendId]      = (indexer(idx)-1)*self._numComponents
         self._sendMPIDataType[sendId] = Mpi.createDataTypeFromRangeAndSubRange(
            sendRgn, localExtRange, self._numComponents, self._layout, elctCommType)
      end

      for _, recvId in ipairs(neigIds) do
         local neighRgn = decomposedRange:subDomain(recvId)
         local recvRgn  = self._localExtRange:intersect(neighRgn)
         local idx      = recvRgn:lowerAsVec()
         -- set idx to starting point of region you want to recv
         self._recvMPILoc[recvId]      = (indexer(idx)-1)*self._numComponents
         self._recvMPIDataType[recvId] = Mpi.createDataTypeFromRangeAndSubRange(
            recvRgn, localExtRange, self._numComponents, self._layout, elctCommType)
         self._recvMPIReq[recvId]      = Mpi.Request()
      end

      -- Create MPI DataTypes for periodic directions.
      -- Also store location in memory required for sending/receiving periodic data.
      self._sendLowerPerMPIDataType, self._recvLowerPerMPIDataType = {}, {}
      self._sendUpperPerMPIDataType, self._recvUpperPerMPIDataType = {}, {}
      self._sendLowerPerMPILoc, self._recvLowerPerMPILoc = {}, {}
      self._sendUpperPerMPILoc, self._recvUpperPerMPILoc = {}, {}
      self._recvLowerPerMPIReq = {}
      self._recvUpperPerMPIReq = {}

      -- Create buffers for periodic copy if Mpi.Comm_size(nodeComm) = 1.
      for dir = 1, self._ndim do
         if self._grid:cuts(dir) == 1 then
            self._lowerPeriodicBuff, self._upperPeriodicBuff = {}, {}
            break
         end
      end

      -- Following loop creates Datatypes for periodic directions.
      -- This is complicated as one needs to treat lower -> upper
      -- transfers differently than upper -> lower as the number of
      -- ghost cells may be different on each lower/upper side. (AHH)
      for dir = 1, self._ndim do
         -- set up periodic-sync Datatypes for all dirs, in case we want to change periodicDirs later
         if self._lowerGhost > 0 and self._upperGhost > 0 then
            local skelIds = decomposedRange:boundarySubDomainIds(dir)
            for i = 1, #skelIds do
               local loId, upId = skelIds[i].lower, skelIds[i].upper

               -- Only create if we are on proper ranks.
               -- Note that if the node communicator has rank size of 1, then we can access all the 
               -- memory needed for periodic boundary conditions and we do not need MPI Datatypes.
               if myId == loId then
                  local rgnSend = decomposedRange:subDomain(loId):lowerSkin(dir, self._upperGhost)
                  if self._grid:cuts(dir) == 1 then
                     local szSend = rgnSend:volume()*self._numComponents
                     self._lowerPeriodicBuff[dir] = ZeroArray.Array(ZeroArray.double, self._numComponents, rgnSend:volume(), self.useDevice)
                  end
                  local idx = rgnSend:lowerAsVec()
                  -- Set idx to starting point of region you want to recv.
                  self._sendLowerPerMPILoc[dir]      = (indexer(idx)-1)*self._numComponents
                  self._sendLowerPerMPIDataType[dir] = Mpi.createDataTypeFromRangeAndSubRange(
                     rgnSend, localExtRange, self._numComponents, self._layout, elctCommType)
                  
                  local rgnRecv = decomposedRange:subDomain(loId):lowerGhost(dir, self._lowerGhost)
                  local idx     = rgnRecv:lowerAsVec()
                  -- Set idx to starting point of region you want to recv.
                  self._recvLowerPerMPILoc[dir]      = (indexer(idx)-1)*self._numComponents
                  self._recvLowerPerMPIDataType[dir] = Mpi.createDataTypeFromRangeAndSubRange(
                     rgnRecv, localExtRange, self._numComponents, self._layout, elctCommType)
                  self._recvLowerPerMPIReq[dir]      = Mpi.Request()
               end
               if myId == upId then
                  local rgnSend = decomposedRange:subDomain(upId):upperSkin(dir, self._lowerGhost)
                  if self._grid:cuts(dir) == 1 then
                     local szSend = rgnSend:volume()*self._numComponents
                     self._upperPeriodicBuff[dir] = ZeroArray.Array(ZeroArray.double, self._numComponents, rgnSend:volume(), self.useDevice)
                  end
                  local idx = rgnSend:lowerAsVec()
                  -- Set idx to starting point of region you want to recv.
                  self._sendUpperPerMPILoc[dir]      = (indexer(idx)-1)*self._numComponents
                  self._sendUpperPerMPIDataType[dir] = Mpi.createDataTypeFromRangeAndSubRange(
                     rgnSend, localExtRange, self._numComponents, self._layout, elctCommType)
                  
                  local rgnRecv = decomposedRange:subDomain(upId):upperGhost(dir, self._upperGhost)
                  local idx     = rgnRecv:lowerAsVec()
                  -- Set idx to starting point of region you want to recv.
                  self._recvUpperPerMPILoc[dir]      = (indexer(idx)-1)*self._numComponents
                  self._recvUpperPerMPIDataType[dir] = Mpi.createDataTypeFromRangeAndSubRange(
                     rgnRecv, localExtRange, self._numComponents, self._layout, elctCommType)
                  self._recvUpperPerMPIReq[dir]      = Mpi.Request()
               end	       
            end
         end
      end
      

      -- Get info for syncs without MPI datatypes. We'll loop over blocks and post a send/recv pair ourselves.
      self._syncSendLoc, self._syncRecvLoc = {}, {}
      self._syncSendBlockSz, self._syncRecvBlockSz = {}, {}
      self._syncSendBlockOff, self._syncRecvBlockOff = {}, {}
      self._syncSendRng, self._syncRecvRng = {}, {}
      self._syncSendBufVol, self._syncRecvBufVol = 0, 0
      for _, sendId in ipairs(neigIds) do
         local neighRgn = decomposedRange:subDomain(sendId)
         local sendRgn  = self._localRange:intersect(
            neighRgn:extend(self._lowerGhost, self._upperGhost))
         local idx      = sendRgn:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self._syncSendLoc[sendId] = (indexer(idx)-1)*self._numComponents
         self._syncSendBlockSz[sendId], self._syncSendBlockOff[sendId] = Mpi.createBlockInfoFromRangeAndSubRange(
            sendRgn, localExtRange, self._numComponents, self._layout)
         self._syncSendRng[sendId] = self._localExtRange:subRange(sendRgn:lowerAsVec(), sendRgn:upperAsVec())
         self._syncSendBufVol = math.max(self._syncSendBufVol, sendRgn:volume())
      end
      for _, recvId in ipairs(neigIds) do
         local neighRgn = decomposedRange:subDomain(recvId)
         local recvRgn  = self._localExtRange:intersect(neighRgn)
         local idx      = recvRgn:lowerAsVec()
         -- set idx to starting point of region you want to recv
         self._syncRecvLoc[recvId] = (indexer(idx)-1)*self._numComponents
         self._syncRecvBlockSz[recvId], self._syncRecvBlockOff[recvId] = Mpi.createBlockInfoFromRangeAndSubRange(
            recvRgn, localExtRange, self._numComponents, self._layout)
         self._syncRecvRng[recvId] = self._localExtRange:subRange(recvRgn:lowerAsVec(), recvRgn:upperAsVec())
         self._syncRecvBufVol = math.max(self._syncRecvBufVol, recvRgn:volume())
      end

      -- Get info for periodic syncs without MPI datatypes.
      self._syncPerNeigh = {}
      self._syncPerSendLoc,      self._syncPerRecvLoc = {}, {}
      self._syncPerSendBlockSz,  self._syncPerRecvBlockSz = {}, {}
      self._syncPerSendBlockOff, self._syncPerRecvBlockOff = {}, {}
      self._syncPerRecvRng,      self._syncPerRecvRng = {}, {}
      self._syncPerRecvBufVol,   self._syncPerRecvBufVol = 0, 0
      self._syncPerSendRng,      self._syncPerRecvRng = {}, {}
      self._syncPerSendBufVol,   self._syncPerRecvBufVol = 0, 0
      for dir = 1, self._ndim do
         -- set up periodic-sync Datatypes for all dirs, in case we want to change periodicDirs later
         if self._lowerGhost > 0 and self._upperGhost > 0 and decomposedRange:numSubDomains() > 1 then
            local skelIds = decomposedRange:boundarySubDomainIds(dir)
            for i = 1, #skelIds do
               local loId, upId = skelIds[i].lower, skelIds[i].upper
               -- Only create if we are on proper ranks.
               -- Note that if the node communicator has rank size of 1, then we can access all the
               -- memory needed for periodic boundary conditions and no communication is needed.
               local rgnSend, rgnRecv, oppId
               if myId == loId and self._grid:cuts(dir) > 1 then
                  rgnSend = decomposedRange:subDomain(loId):lowerSkin(dir, self._upperGhost)
                  rgnRecv = decomposedRange:subDomain(loId):lowerGhost(dir, self._lowerGhost)
                  oppId = upId
               end
               if myId == upId and self._grid:cuts(dir) > 1 then
                  rgnSend = decomposedRange:subDomain(upId):upperSkin(dir, self._lowerGhost)
                  rgnRecv = decomposedRange:subDomain(upId):upperGhost(dir, self._upperGhost)
                  oppId = loId
               end
               if oppId ~= nil then
                  self._syncPerNeigh[dir] = oppId
                  self._syncPerSendRng[dir] = self._localExtRange:subRange(rgnSend:lowerAsVec(), rgnSend:upperAsVec())
                  self._syncPerSendBufVol = math.max(self._syncPerSendBufVol, rgnSend:volume())
                  local idx = rgnSend:lowerAsVec()
                  -- Set idx to starting point of region you want to send.
                  self._syncPerSendLoc[dir] = (indexer(idx)-1)*self._numComponents
                  self._syncPerSendBlockSz[dir], self._syncPerSendBlockOff[dir] = Mpi.createBlockInfoFromRangeAndSubRange(
                     rgnSend, localExtRange, self._numComponents, self._layout)

                  self._syncPerRecvRng[dir] = self._localExtRange:subRange(rgnRecv:lowerAsVec(), rgnRecv:upperAsVec())
                  self._syncPerRecvBufVol = math.max(self._syncPerRecvBufVol, rgnRecv:volume())
                  local idx = rgnRecv:lowerAsVec()
                  -- Set idx to starting point of region you want to recv.
                  self._syncPerRecvLoc[dir] = (indexer(idx)-1)*self._numComponents
                  self._syncPerRecvBlockSz[dir], self._syncPerRecvBlockOff[dir] = Mpi.createBlockInfoFromRangeAndSubRange(
                     rgnRecv, localExtRange, self._numComponents, self._layout)
               end
            end
         end
      end

      local structAny = function(key,structIn, idIn)
         for i = 1, #structIn do if structIn[i][key] == idIn then return true end end
         return false
      end
      local tblAbs = function(tblIn)
         local tblOut = lume.clone(tblIn)
         for i, v in ipairs(tblOut) do tblOut[i] = math.abs(v) end
         return tblOut
      end
      -- For corner sync, decomposedRange provides all the corner neighbors assuming all
      -- directions are periodic. Because some directions may not be periodic, not all of
      -- those corner neighbors need to communicate. Remove such pairs.
      self._cornersToSync = {}
      for dir = 1, self._ndim do
         self._cornersToSync[dir] = {}
         if grid:isDirPeriodic(dir) then
            local corIds = decomposedRange:boundarySubDomainCornerIds(dir)
            for i, bD in ipairs(corIds) do   -- Loop over lower boundary subdomains.
               self._cornersToSync[dir][i] = {}
               for j, dC in ipairs(bD) do   -- Loop over corners.
                  local loId, upId, corDirs = dC.lower, dC.upper, dC.dirs

                  -- Check if this subdomain abutts a boundary in another periodic dimension.
                  local syncThisCorner = true
                  for dI = 2,#corDirs do
                     local oDir    = math.abs(corDirs[dI])   -- Signs used to differentiate corners. See CartDecomp.
                     local skelIds = decomposedRange:boundarySubDomainIds(oDir)
                     if structAny("lower",skelIds,loId) then
                        -- This subdomain is on lower boundary of other dimension. If
                        -- other direction is not periodic then don't sync this corner.
                        if corDirs[dI]<0 and not grid:isDirPeriodic(oDir) then syncThisCorner=false end
                     end
                     if structAny("upper",skelIds,loId) then
                        -- This subdomain is on upper boundary of other dimension. If
                        -- other direction is not periodic then don't sync this corner.
                        if corDirs[dI]>0 and not grid:isDirPeriodic(oDir) then syncThisCorner=false end
                     end
                     -- Subdomains not on the boundary in the other direction do sync this corner. 
                  end

                  if syncThisCorner then
                     table.insert(self._cornersToSync[dir][i],{lower=loId, upper=upId, dirs=corDirs})
                  end

               end
            end
         end
      end

      -- Create MPI DataTypes for syncing corners in periodic directions.
      -- Also store location in memory required for sending/receiving periodic data.
      -- Note: please understand the code for the face-sync first. It'll help in understanding corner sync.
      self._sendLowerCornerPerMPIDataType, self._recvLowerCornerPerMPIDataType = {}, {}
      self._sendUpperCornerPerMPIDataType, self._recvUpperCornerPerMPIDataType = {}, {}
      self._sendLowerCornerPerMPILoc, self._recvLowerCornerPerMPILoc = {}, {}
      self._sendUpperCornerPerMPILoc, self._recvUpperCornerPerMPILoc = {}, {}
      self._recvLowerCornerPerMPIReq = {}
      self._recvUpperCornerPerMPIReq = {}
      for dir = 1, self._ndim do
         self._sendLowerCornerPerMPILoc[dir]     , self._sendUpperCornerPerMPILoc[dir]      = {}, {}
         self._sendLowerCornerPerMPIDataType[dir], self._sendUpperCornerPerMPIDataType[dir] = {}, {}
         self._recvLowerCornerPerMPILoc[dir]     , self._recvUpperCornerPerMPILoc[dir]      = {}, {}
         self._recvLowerCornerPerMPIDataType[dir], self._recvUpperCornerPerMPIDataType[dir] = {}, {}
         self._recvLowerCornerPerMPIReq[dir]     , self._recvUpperCornerPerMPIReq[dir]      = {}, {}
         if grid:isDirPeriodic(dir) then
            local cTs = self._cornersToSync[dir]
            for dI, bD in ipairs(cTs) do   -- Loop over lower boundary subdomains.
               for cI, dC in ipairs(bD) do   -- Loop over corners.
                  local loId, upId, corDirs = dC.lower, dC.upper, dC.dirs
                  if myId == loId then
                     local rgnSend = decomposedRange:subDomain(loId):lowerSkin(dir, self._upperGhost)
                     for dI = 2,#corDirs do
                        local oDir = math.abs(corDirs[dI])   -- Signs used to differentiate corners. See CartDecomp.
                        if corDirs[dI]<0 then
                           rgnSend = rgnSend:shorten(oDir, self._upperGhost)
                        else
                           rgnSend = rgnSend:shortenFromBelow(oDir, self._upperGhost)
                        end
                     end
                     local idx = rgnSend:lowerAsVec()
                     -- Set idx to starting point of region you want to recv.
                     table.insert(self._sendLowerCornerPerMPILoc[dir], (indexer(idx)-1)*self._numComponents)
                     table.insert(self._sendLowerCornerPerMPIDataType[dir], Mpi.createDataTypeFromRangeAndSubRange(
                        rgnSend, localExtRange, self._numComponents, self._layout, elctCommType))

                     local rgnRecv = decomposedRange:subDomain(loId):extendDirs(tblAbs(corDirs),self._lowerGhost,self._upperGhost)
                     local rgnRecv = rgnRecv:shorten(dir, self._lowerGhost)
                     for dI = 2,#corDirs do
                        local oDir = math.abs(corDirs[dI])   -- Signs used to differentiate corners. See CartDecomp.
                        if corDirs[dI]<0 then
                           rgnRecv = rgnRecv:shorten(oDir, self._upperGhost)
                        else
                           rgnRecv = rgnRecv:shortenFromBelow(oDir, self._upperGhost)
                        end
                     end
                     local idx = rgnRecv:lowerAsVec()
                     -- Set idx to starting point of region you want to recv.
                     table.insert(self._recvLowerCornerPerMPILoc[dir], (indexer(idx)-1)*self._numComponents)
                     table.insert(self._recvLowerCornerPerMPIDataType[dir], Mpi.createDataTypeFromRangeAndSubRange(
                        rgnRecv, localExtRange, self._numComponents, self._layout, elctCommType))
                     table.insert(self._recvLowerCornerPerMPIReq[dir], Mpi.Request())
                  end
                  if myId == upId then
                     local rgnSend = decomposedRange:subDomain(upId):upperSkin(dir, self._lowerGhost)
                     for dI = 2,#corDirs do
                        local oDir = math.abs(corDirs[dI])   -- Signs used to differentiate corners. See CartDecomp.
                        if corDirs[dI]<0 then
                           rgnSend = rgnSend:shortenFromBelow(oDir, self._upperGhost)
                        else
                           rgnSend = rgnSend:shorten(oDir, self._upperGhost)
                        end
                     end
                     local idx = rgnSend:lowerAsVec()
                     -- Set idx to starting point of region you want to recv.
                     table.insert(self._sendUpperCornerPerMPILoc[dir], (indexer(idx)-1)*self._numComponents)
                     table.insert(self._sendUpperCornerPerMPIDataType[dir], Mpi.createDataTypeFromRangeAndSubRange(
                        rgnSend, localExtRange, self._numComponents, self._layout, elctCommType))

                     local rgnRecv = decomposedRange:subDomain(upId):extendDirs(tblAbs(corDirs),self._lowerGhost,self._upperGhost)
                     local rgnRecv = rgnRecv:shortenFromBelow(dir, self._upperGhost)
                     for dI = 2,#corDirs do
                        local oDir = math.abs(corDirs[dI])   -- Signs used to differentiate corners. See CartDecomp.
                        if corDirs[dI]<0 then
                           rgnRecv = rgnRecv:shortenFromBelow(oDir, self._upperGhost)
                        else
                           rgnRecv = rgnRecv:shorten(oDir, self._upperGhost)
                        end
                     end
                     local idx = rgnRecv:lowerAsVec()
                     -- Set idx to starting point of region you want to recv.
                     table.insert(self._recvUpperCornerPerMPILoc[dir], (indexer(idx)-1)*self._numComponents)
                     table.insert(self._recvUpperCornerPerMPIDataType[dir], Mpi.createDataTypeFromRangeAndSubRange(
                        rgnRecv, localExtRange, self._numComponents, self._layout, elctCommType))
                     table.insert(self._recvUpperCornerPerMPIReq[dir], Mpi.Request())
                  end
               end
            end
         end
      end

      -- Create IO object.
      self._adiosIo = AdiosCartFieldIo {
         elemType = elct,
         metaData = tbl.metaData,
      }
      -- tag to identify basis used to set this field
      self._metaData = tbl.metaData
      self._basisId  = "none"

      return self
   end
   setmetatable(Field, { __call = function (self, o) return self.new(self, o) end })

   -- Set callable methods.
   Field.__index = {
      elemType = function (self)
	 return elct
      end,
      elemSize = function (self)
	 return elctSize
      end,      
      elemCommType = function (self)
	 return elctCommType
      end,
      ndim = function (self)
	 return self._ndim
      end,
      grid = function (self)
	 return self._grid
      end,
      numComponents = function (self)
	 return self._numComponents
      end,
      hasCuDev = function (self)
         if self._zeroDevice then return self._zeroDevice:is_cu_dev() else return false end
      end,
      setHostOps = function (self)
         self._zeroForOps = self._zero
      end,
      setDeviceOps = function (self)
         self._zeroForOps = self._zeroDevice
      end,
      copy = function (self, fIn)
	 self._zeroForOps:copy(fIn._zeroForOps)
      end,
      deviceCopy = function (self, fIn)
	 self._zeroDevice:copy(fIn._zeroDevice)
      end,
      copyHostToDevice = function (self)
         if self._zeroDevice then self._zeroDevice:copy(self._zero) end
      end,
      copyDeviceToHost = function (self)
         if self._zeroDevice then self._zero:copy(self._zeroDevice) end
      end,
      copyDeviceToHostAsync = function (self)
	 self._zero:copyAsync(self._zeroDevice)
      end,
      copyHostToDeviceAsync = function (self)
	 self._zero:copyAsync(self._zeroDevice)
      end,
      deviceDataPointer = function (self)
	 return self._devAllocData
      end,
      dataPointer = function (self)
	 return self._data
      end,
      clear = function (self, val)
         self._zeroForOps:clear(val)
      end,
      fill = function (self, k, fc)
	 local loc = (k-1)*self._numComponents -- (k-1) as k is 1-based index	 
	 fc._cdata = self._data+loc
      end,
      _localLower = function (self)
	 return 0
      end,
      _localShape = function (self)
	 return self._localExtRange:volume()*self:numComponents()
      end,
      _assign = function(self, fact, fld)
         self._zeroForOps:set(fact, fld._zeroForOps)
      end,
      _assignRange = function(self, fact, fld, rng)
         self._zeroForOps:setRange(fact, fld._zeroForOps, rng)
      end,
      -- assignOffsetOneFld assumes that one of the input or output fields have fewer components than the other.
      --   a) nCompOut > nCompIn: assigns all of the input field w/ part of the output field, the (0-based)
      --                          offset indicating which is the 1st component in the output field to assign.
      --   b) nCompOut < nCompIn: assigns nCompOut components of the input field onto the output field, the
      --                          (0-based) offset indicates which is the first component in the input field to read.
      -- It assumes that the components being summed are continuous.
      _assignOffsetOneFld = function(self, fact, fld, compStart)
	 assert(field_check_range(self, fld),
		"CartField:assignOffsetOneFld: Can only assign fields with the same range")
	 assert(type(fact) == "number",
		"CartField:assignOffsetOneFld: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:assignOffsetOneFld: Fields should have same layout for sums to make sense")
         self._zeroForOps:setOffset(fact, fld._zeroForOps, compStart)
      end,
      _assignOffsetOneFldRange = function(self, fact, fld, compStart, rng)
	 assert(field_check_range(self, fld),
		"CartField:assignOffsetOneFld: Can only assign fields with the same range")
	 assert(type(fact) == "number",
		"CartField:assignOffsetOneFld: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:assignOffsetOneFld: Fields should have same layout for sums to make sense")
         self._zeroForOps:setOffsetRange(fact, fld._zeroForOps, compStart, rng)
      end,
      _accumulateOneFld = function(self, fact, fld)
         self._zeroForOps:accumulate(fact, fld._zeroForOps)
      end,
      _accumulateOneFldRange = function(self, fact, fld, rng)
         self._zeroForOps:accumulateRange(fact, fld._zeroForOps, rng)
      end,
      -- accumulateOffsetOneFld assumes that one of the input or output fields have fewer components than the other.
      --   a) nCompOut > nCompIn: accumulates all of the input field w/ part of the output field, the (0-based)
      --                          offset indicating which is the 1st component in the output field to accumulate to.
      --   b) nCompOut < nCompIn: accumulates nCompIn components of the input field onto the output field, the
      --                          (0-based) offset indicates which is the first component in the input field to accumulate. 
      -- It assumes that the components being summed are continuous.
      _accumulateOffsetOneFld = function(self, fact, fld, compStart)
	 assert(field_check_range(self, fld),
		"CartField:accumulateOffsetOneFld: Can only accumulate fields with the same range")
	 assert(type(fact) == "number",
		"CartField:accumulateOffsetOneFld: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:accumulateOffsetOneFld: Fields should have same layout for sums to make sense")
         self._zeroForOps:accumulateOffset(fact, fld._zeroForOps, compStart)
      end,
      _accumulateOffsetOneFldRange = function(self, fact, fld, compStart, rng)
	 assert(field_check_range(self, fld),
		"CartField:accumulateOffsetOneFld: Can only accumulate fields with the same range")
	 assert(type(fact) == "number",
		"CartField:accumulateOffsetOneFld: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:accumulateOffsetOneFld: Fields should have same layout for sums to make sense")
         self._zeroForOps:accumulateOffsetRange(fact, fld._zeroForOps, compStart, rng)
      end,
      accumulate = isNumberType and
	 function (self, c1, fld1, ...)
	    local args = {...} -- Package up rest of args as table.
	    local nFlds = #args/2
	    self:_accumulateOneFld(c1, fld1) -- Accumulate first field.
	    for i = 1, nFlds do -- Accumulate rest of the fields.
	       self:_accumulateOneFld(args[2*i-1], args[2*i])
	    end
	 end or
	 function (self, c1, fld1, ...)
	    assert(false, "CartField:accumulate: Accumulate only works on numeric fields")
	 end,
      accumulateOffset = isNumberType and
	 function (self, c1, fld1, compStart1, ...)
	    local args = {...} -- Package up rest of args as table.
	    local nFlds = #args/3
	    self:_accumulateOffsetOneFld(c1, fld1, compStart1) -- Accumulate first field.
	    for i = 1, nFlds do -- Accumulate rest of the fields
	       self:_accumulateOffsetOneFld(args[3*i-2], args[3*i-1], args[3*i])
	    end
	 end or
	 function (self, c1, fld1, compStart1, ...)
	    assert(false, "CartField:accumulateOffset: Accumulate only works on numeric fields")
	 end,
      combine = isNumberType and
         function (self, c1, fld1, ...)
            local args = {...} -- Package up rest of args as table.
            local nFlds = #args/2
            self:_assign(c1, fld1) -- Assign first field.
            for i = 1, nFlds do -- Accumulate rest of the fields.
               self:_accumulateOneFld(args[2*i-1], args[2*i])
            end
         end or
         function (self, c1, fld1, ...)
            assert(false, "CartField:combine: Combine only works on numeric fields")
         end,
      combineOffset = isNumberType and
	 function (self, c1, fld1, compStart1, ...)
            local args = {...} -- Package up rest of args as table.
            local nFlds = #args/3
            local notAssigned = {}
            for i = 1, self:numComponents() do table.insert(notAssigned,true) end   -- Boolean indicates if already assigned.
            self:_assignOffsetOneFld(c1, fld1, compStart1) -- Assign first field.
            notAssigned[compStart1+1] = false
            for i = 1, nFlds do -- Accumulate rest of the fields.
               local cOff = args[3*i]
               if notAssigned[cOff+1] then
                  self:_assignOffsetOneFld(args[3*i-2], args[3*i-1], cOff)
                  notAssigned[cOff+1] = false
               else
                  self:_accumulateOffsetOneFld(args[3*i-2], args[3*i-1], cOff)
               end
            end
         end or
         function (self, c1, fld1, ...)
            assert(false, "CartField:combineOffset: Combine only works on numeric fields")
         end,
      accumulateRange = isNumberType and
	 function (self, c1, fld1, ...)
	    local args = {...} -- Package up rest of args as table.
	    local nFlds = (#args-1)/2
	    local rng = args[#args]
	    self:_accumulateOneFldRange(c1, fld1, rng) -- Accumulate first field.
	    for i = 1, nFlds do -- Accumulate rest of the fields.
	       self:_accumulateOneFldRange(args[2*i-1], args[2*i], rng)
	    end
	 end or
	 function (self, c1, fld1, ...)
	    assert(false, "CartField:accumulate: Accumulate only works on numeric fields")
	 end,
      accumulateOffsetRange = isNumberType and
	 function (self, c1, fld1, compStart1, ...)
	    local args = {...} -- Package up rest of args as table.
	    local nFlds = #args/3
	    local rng = args[#args]
	    self:_accumulateOffsetOneFldRange(c1, fld1, compStart1, rng) -- Accumulate first field.
	    for i = 1, nFlds do -- Accumulate rest of the fields
	       self:_accumulateOffsetOneFldRange(args[3*i-2], args[3*i-1], args[3*i], rng)
	    end
	 end or
	 function (self, c1, fld1, compStart1, ...)
	    assert(false, "CartField:accumulateOffset: Accumulate only works on numeric fields")
	 end,
      combineRange = isNumberType and
         function (self, c1, fld1, ...)
            local args = {...} -- Package up rest of args as table.
            local nFlds = #args/2
	    local rng = args[#args]
            self:_assignRange(c1, fld1, rng) -- Assign first field.
            for i = 1, nFlds do -- Accumulate rest of the fields.
               self:_accumulateOneFldRange(args[2*i-1], args[2*i], rng)
            end
         end or
         function (self, c1, fld1, ...)
            assert(false, "CartField:combine: Combine only works on numeric fields")
         end,
      combineOffsetRange = isNumberType and
	 function (self, c1, fld1, compStart1, ...)
            local args = {...} -- Package up rest of args as table.
            local nFlds = #args/3
	    local rng = args[#args]
            local notAssigned = {}
            for i = 1, self:numComponents() do table.insert(notAssigned,true) end   -- Boolean indicates if already assigned.
            self:_assignOffsetOneFldRange(c1, fld1, compStart1, rng) -- Assign first field.
            notAssigned[compStart1+1] = false
            for i = 1, nFlds do -- Accumulate rest of the fields.
               local cOff = args[3*i]
               if notAssigned[cOff+1] then
                  self:_assignOffsetOneFldRange(args[3*i-2], args[3*i-1], cOff, rng)
                  notAssigned[cOff+1] = false
               else
                  self:_accumulateOffsetOneFldRange(args[3*i-2], args[3*i-1], cOff, rng)
               end
            end
         end or
         function (self, c1, fld1, ...)
            assert(false, "CartField:combineOffset: Combine only works on numeric fields")
         end,
      scale = isNumberType and
	 function (self, fact)
            self._zeroForOps:scale(fact)
	 end or
	 function (self, fact)
	    assert(false, "CartField:scale: Scale only works on numeric fields")
	 end,
      scaleByCell = function (self, factByCell)
         self._zeroForOps:scale_by_cell(factByCell._zeroForOps)
      end,
      abs = isNumberType and
         function (self)
            ffiC.gkylCartFieldAbs(self:_localLower(), self:_localShape(), self._data)
	 end or
	 function (self, fact)
	    assert(false, "CartField:abs: Abs only works on numeric fields")
	 end,
      defaultLayout = function (self)
	 if defaultLayout == rowMajLayout then
	    return "row-major"
	 end
	 return "col-major"
      end,
      layout = function (self)
	 if self._layout == rowMajLayout then
	    return "row-major"
	 end
	 return "col-major"
      end,
      lowerGhost = function (self)
	 return self._lowerGhost
      end,
      upperGhost = function (self)
	 return self._upperGhost
      end,
      localRange = function (self)
	 return self._localRange
      end,
      localExtRange = function (self) -- includes ghost cells
	 return self._localExtRange
      end,      
      localEdgeRange = function (self)
	 return self._localEdgeRange
      end,      
      localExtEdgeRange = function (self)
	 return self._localExtEdgeRange
      end,      
      globalRange = function (self)
	 return self._globalRange
      end,
      globalExtRange = function (self) -- includes ghost cells
	 return self._globalExtRange
      end,
      localRangeIter = function (self)
	 if self._layout == rowMajLayout then
	    return self._localRange:rowMajorIter(self._localStartIdx, self._localNumBump)
	 end
	 return self._localRange:colMajorIter(self._localStartIdx, self._localNumBump)
      end,
      localExtRangeIter = function (self) -- includes ghost cells
	 local lext = self:localRange():extend(self:lowerGhost(), self:upperGhost())
	 if self._layout == rowMajLayout then
	    return lext:rowMajorIter(self._localExtStartIdx, self._localExtNumBump)
	 end
	 return lext:colMajorIter(self._localExtStartIdx, self._localExtNumBump)
      end,      
      size = function (self)
	 return self._size
      end,
      indexer = function (self) -- linear indexer taking (i,j,...)
	 return indexerMakerFuncs[self._layout](self:localExtRange())
      end,
      genIndexer = function (self) -- linear indexer taking indices as a vector
	 return genIndexerMakerFuncs[self._layout](self:localExtRange())
      end,
      get = function (self, k) -- k is an integer returned by a linear indexer
	 local loc = (k-1)*self._numComponents -- (k-1) as k is 1-based index
	 return fcompct(self._numComponents, self._data+loc)
      end,
      getDataPtrAt = function (self, k) -- k is an integer returned by a linear indexer
	 local loc = (k-1)*self._numComponents -- (k-1) as k is 1-based index
	 return self._data+loc
      end,
      write = function (self, fName, tmStamp, frNum, writeGhost)
	 self._adiosIo:write(self, fName, tmStamp, frNum, writeGhost)
      end,
      read = function (self, fName) --> time-stamp, frame-number
	 return self._adiosIo:read(self, fName)
      end,
      sync = function (self, syncPeriodicDirs_)
         local syncPeriodicDirs = xsys.pickBool(syncPeriodicDirs_, true)
--         self._field_sync(self, self._zeroForOps:data())
         local mess = self._grid:getMessenger()
         mess:syncCartField(self, mess:getConfComm())
	 if self._syncPeriodicDirs and syncPeriodicDirs then
	    --self._field_periodic_sync(self, self._zeroForOps:data())
            mess:syncPeriodicCartField(self, mess:getConfComm())
	 end
      end,
      -- This method is an alternative function for applying periodic boundary
      -- conditions when Mpi.Comm_size(nodeComm) = 1 and we do not need to call
      -- Send/Recv to copy the skin cell data into ghost cells.
      periodicCopy = function (self, syncPeriodicDirs_)
         local syncPeriodicDirs = xsys.pickBool(syncPeriodicDirs_, true)
	 if self._syncPeriodicDirs and syncPeriodicDirs then
            self._field_periodic_copy(self)
	 end
      end,
      periodicCopyInDir = function (self, dir)
         self._field_periodic_copy_indir(self, dir)
      end,
      setBasisId = function(self, basisId)
         self._basisId = basisId
      end,
      getBasisId = function(self)
         return self._basisId
      end,
      getMetaData = function(self)
         return self._metaData 
      end,
      compatible = function(self, fld)
         return field_compatible(self, fld)
      end,
      checkRange = function(self, fld)
         return field_check_range(self, fld)
      end,
      reduce = isNumberType and
	 function(self, opIn)
            self._zeroForOps:reduceRange(self.localReductionVal:data(), opIn, self._localRange)
	    if self.useDevice then self.localReductionVal_h:copy(self.localReductionVal) end
	      
	    -- Input 'opIn' must be one of the binary operations in binOpFuncs.
	    local grid = self._grid
	    local localVal = {}
            Mpi.Allreduce(self.localReductionVal_h:data(), self.globalReductionVal:data(),
               self._numComponents, elctCommType, reduceOpsMPI[opIn], grid:commSet().comm)

            --self.localReductionVal:copy(self.globalReductionVal)
            for k = 1, self._numComponents do localVal[k] = self.globalReductionVal:data()[k-1] end
            return localVal
            --return self.localReductionVal:data()
	 end or
	 function (self, opIn)
	    assert(false, "CartField:reduce: Reduce only works on numeric fields")
	 end,
      reduceByCell = function(self, opIn, comm, localFldIn)
         -- Reduce the CartField 'localFldIn' across communicator 'comm',
         -- and put the result in this CartField.
         assert(field_compatible_ext(self, localFldIn))
         Mpi.Allreduce(localFldIn:dataPointer(), self._data,
            self:size(), elctCommType, reduceOpsMPI[opIn], comm)
      end,
      _copy_from_field_region = function (self, rgn, data)
	 self._zeroForOps:copy_to_buffer(data:data(), rgn)
      end,
      _copy_to_field_region = function (self, rgn, data)
	 self._zeroForOps:copy_from_buffer(data:data(), rgn)
      end,
      _field_sync = function (self, dataPtr)
         local comm = self._grid:commSet().nodeComm -- communicator to use
         if not Mpi.Is_comm_valid(comm) then
            return -- no need to do anything if communicator is not valid
         end
         -- immediately return if nothing to sync
         if self._lowerGhost == 0 and self._upperGhost == 0 then return end

         -- immediately return if there is no decomp
         if self._grid:decomposedRange():numSubDomains() == 1 then return end
        
         -- Steps: (1) Post non-blocking recv requests. (2) Do
         -- blocking sends, (3) Complete recv and copy data into ghost
         -- cells.
        
         local myId    = self._grid:subGridId() -- grid ID on this processor
         local neigIds = self._decompNeigh:neighborData(myId) -- list of neighbors
         local tag     = 42 -- Communicator tag for regular (non-periodic) messages
         -- post a non-blocking recv request
         for _, recvId in ipairs(neigIds) do
            local dataType = self._recvMPIDataType[recvId]
            local loc      = self._recvMPILoc[recvId]
            local recvReq  = self._recvMPIReq[recvId]
            -- recv data: (its from recvId-1 as MPI ranks are zero indexed)
            local _ = Mpi.Irecv(dataPtr+loc, 1, dataType, recvId-1, tag, comm, recvReq)
         end
         
         -- Do a blocking send (does not really block as recv requests
         -- are already posted).
         for _, sendId in ipairs(neigIds) do
            local dataType = self._sendMPIDataType[sendId]
            local loc      = self._sendMPILoc[sendId]
            -- Send data: (its to sendId-1 as MPI ranks are zero indexed).
            Mpi.Send(dataPtr+loc, 1, dataType, sendId-1, tag, comm)
         end
        
         -- Complete recv.
         -- Since MPI DataTypes eliminate the need for buffers, 
         -- all we have to do is wait for non-blocking receives to finish.
         for _, recvId in ipairs(neigIds) do local _ = Mpi.Wait(self._recvMPIReq[recvId], nil) end
      end,
      _field_periodic_sync = function (self, dataPtr)
         local comm = self._grid:commSet().nodeComm -- Communicator to use.
         if not Mpi.Is_comm_valid(comm) then
            return -- No need to do anything if communicator is not valid
         end
         
         -- Immediately return if nothing to sync.
         if self._lowerGhost == 0 and self._upperGhost == 0 then return end
         
         local grid = self._grid
         
         -- Steps: (1) Post non-blocking recv requests. (2) Do
         -- blocking sends, (3) Complete recv and copy data into ghost
         -- cells.
         
         local decomposedRange = self._grid:decomposedRange()
         if decomposedRange:numSubDomains() == 1 and not self._syncCorners then
            return self:periodicCopy()
         end
         local myId       = self._grid:subGridId() -- Grid ID on this processor.
         local basePerTag = 53 -- Tag for periodic BCs.
         local cornerBasePerTag = 70 -- Tag for periodic corner sync.
         
         -- Note on tags: Each MPI message (send/recv pair) must have
         -- a unique tag. With periodic BCs it is possible that a rank
         -- may send/recv more than one message and hence some way to
         -- distinguishing various messages is needed. The
         -- non-periodic messages are all tagged 42. The periodic
         -- messages have a base tag of 53, with up->lo tags being
         -- 53+dir+10, while lo->up tags being 53+dir. As at most we
         -- will have 6 directions, this generates enough unique tags
         -- to do periodic communication safely.
         -- However we must also account for corner syncs. In 6D there
         -- are 716 such "corners" that may need to be sync-ed. We will
         -- thus say that up->lo tags have 70+tn+800, while lo->up tags 
         -- have 70+tn, where tn is the corner tag number. One way to make
         -- these tag numbers identical is to have them be composed of the
         -- ID of the participating ranks and the corner directions.
         
         -- Post non-blocking recv requests for periodic directions.
         for dir = 1, self._ndim do
            if grid:isDirPeriodic(dir) then
               local skelIds = decomposedRange:boundarySubDomainIds(dir)
               for i = 1, #skelIds do
                  local loId, upId = skelIds[i].lower, skelIds[i].upper
         
                  if myId == loId then
                     local loTag    = basePerTag+dir+10
                     local dataType = self._recvLowerPerMPIDataType[dir]
                     local loc      = self._recvLowerPerMPILoc[dir]
                     local req      = self._recvLowerPerMPIReq[dir]
                     local _        = Mpi.Irecv(dataPtr+loc, 1, dataType,
                                                upId-1, loTag, comm, req)
                  end
                  if myId == upId then
                     local upTag    = basePerTag+dir
                     local dataType = self._recvUpperPerMPIDataType[dir]
                     local loc      = self._recvUpperPerMPILoc[dir]
                     local req      = self._recvUpperPerMPIReq[dir]
                     local _        = Mpi.Irecv(dataPtr+loc, 1, dataType,
                                                loId-1, upTag, comm, req)
                  end
               end

               if self._syncCorners then
                  local cTs = self._cornersToSync[dir]
                  local ccLo, ccUp = 0, 0
                  for bI, bD in ipairs(cTs) do   -- Loop over lower boundary subdomains.
                     for _, dC in ipairs(bD) do   -- Loop over corners.
                        local loId, upId, corDirs = dC.lower, dC.upper, dC.dirs
                        if myId == loId then
                           ccLo = ccLo+1
                           local loTag    = cornerBasePerTag+tonumber(loId..upId..tblToStr(corDirs)..1)
                           local dataType = self._recvLowerCornerPerMPIDataType[dir][ccLo]
                           local loc      = self._recvLowerCornerPerMPILoc[dir][ccLo]
                           local req      = self._recvLowerCornerPerMPIReq[dir][ccLo]
                           local _        = Mpi.Irecv(dataPtr+loc, 1, dataType,
                                                      upId-1, loTag, comm, req)
                        end
                        if myId == upId then
                           ccUp = ccUp+1
                           local upTag    = cornerBasePerTag+tonumber(loId..upId..tblToStr(corDirs)..2)
                           local dataType = self._recvUpperCornerPerMPIDataType[dir][ccUp]
                           local loc      = self._recvUpperCornerPerMPILoc[dir][ccUp]
                           local req      = self._recvUpperCornerPerMPIReq[dir][ccUp]
                           local _        = Mpi.Irecv(dataPtr+loc, 1, dataType,
                                                      loId-1, upTag, comm, req)
                        end
                     end
                  end
               end
            end
         end
         
         -- Do a blocking send for periodic directions (does not
         -- really block as recv requests are already posted).
         for dir = 1, self._ndim do
            if grid:isDirPeriodic(dir) then
               local skelIds = decomposedRange:boundarySubDomainIds(dir)
               for i = 1, #skelIds do
                  local loId, upId = skelIds[i].lower, skelIds[i].upper
         
                  if myId == loId then
                     local loTag    = basePerTag+dir -- This must match recv tag posted above.
                     local dataType = self._sendLowerPerMPIDataType[dir]
                     local loc      = self._sendLowerPerMPILoc[dir]
                     Mpi.Send(dataPtr+loc, 1, dataType, upId-1, loTag, comm)
                  end
                  if myId == upId then
                     local upTag    = basePerTag+dir+10 -- This must match recv tag posted above.
                     local dataType = self._sendUpperPerMPIDataType[dir]
                     local loc      = self._sendUpperPerMPILoc[dir]
                     Mpi.Send(dataPtr+loc, 1, dataType, loId-1, upTag, comm)
                  end
               end

               if self._syncCorners then
                  local cTs = self._cornersToSync[dir]
                  local ccLo, ccUp = 0, 0
                  for _, bD in ipairs(cTs) do   -- Loop over lower boundary subdomains.
                     for _, dC in ipairs(bD) do   -- Loop over corners.
                        local loId, upId, corDirs = dC.lower, dC.upper, dC.dirs
                        if myId == loId then
                           ccUp = ccUp+1
                           local loTag    = cornerBasePerTag+tonumber(loId..upId..tblToStr(corDirs)..2)
                           local dataType = self._sendLowerCornerPerMPIDataType[dir][ccUp]
                           local loc      = self._sendLowerCornerPerMPILoc[dir][ccUp]
                           Mpi.Send(dataPtr+loc, 1, dataType, upId-1, loTag, comm)
                        end
                        if myId == upId then
                           ccLo = ccLo+1
                           local upTag    = cornerBasePerTag+tonumber(loId..upId..tblToStr(corDirs)..1)
                           local dataType = self._sendUpperCornerPerMPIDataType[dir][ccLo]
                           local loc      = self._sendUpperCornerPerMPILoc[dir][ccLo]
                           Mpi.Send(dataPtr+loc, 1, dataType, loId-1, upTag, comm)
                        end
                     end
                  end
               end
            end
         end
         
         -- Complete recv for periodic directions.
         -- Since MPI DataTypes eliminate the need for buffers,
         -- all we have to do is wait for non-blocking receives to finish.
         for dir = 1, self._ndim do
            if grid:isDirPeriodic(dir) then
               local skelIds = decomposedRange:boundarySubDomainIds(dir)
               for i = 1, #skelIds do
                  local loId, upId = skelIds[i].lower, skelIds[i].upper
                  if myId == loId then local _ = Mpi.Wait(self._recvLowerPerMPIReq[dir], nil) end
                  if myId == upId then local _ = Mpi.Wait(self._recvUpperPerMPIReq[dir], nil) end
               end

               if self._syncCorners then
                  local cTs = self._cornersToSync[dir]
                  local ccLo, ccUp = 0, 0
                  for dI, bD in ipairs(cTs) do   -- Loop over lower boundary subdomains.
                     for _, dC in ipairs(bD) do   -- Loop over corners.
                        local loId, upId = dC.lower, dC.upper
                        if myId == loId then
                           ccLo = ccLo+1
                           local _ = Mpi.Wait(self._recvLowerCornerPerMPIReq[dir][ccLo], nil)
                        end
                        if myId == upId then
                           ccUp = ccUp+1
                           local _ = Mpi.Wait(self._recvUpperCornerPerMPIReq[dir][ccUp], nil)
                        end
                     end
                  end
               end
            end
         end
      end,
      _field_periodic_copy_indir = function (self, dir)
         -- First get region for skin cells for upper ghost region (the lower skin cells).
         local skinRgnLower = self._skinRgnLower[dir]
         local skinRgnUpper = self._skinRgnUpper[dir]
         local ghostRgnLower = self._ghostRgnLower[dir]
         local ghostRgnUpper = self._ghostRgnUpper[dir]

         local periodicBuffUpper = self._upperPeriodicBuff[dir]
         -- Copy skin cells into temporary buffer.
         self:_copy_from_field_region(skinRgnUpper, periodicBuffUpper)
         -- Get region for looping over upper ghost cells and copy lower skin cells into upper ghost cells.
         self:_copy_to_field_region(ghostRgnLower, periodicBuffUpper)

	 -- Now do the same, but for the skin cells for the lower ghost region (the upper skin cells).
         local periodicBuffLower = self._lowerPeriodicBuff[dir]
         self:_copy_from_field_region(skinRgnLower, periodicBuffLower)
         self:_copy_to_field_region(ghostRgnUpper, periodicBuffLower)
      end,
      _field_periodic_copy = function (self)
         local grid = self._grid
         for dir = 1, self._ndim do
            if grid:isDirPeriodic(dir) then self:_field_periodic_copy_indir(dir) end
         end
      end,
   }
   
   return Field
end

return {
   new_field_ct = Field_meta_ctor,
   Field = Field_meta_ctor(typeof("double")),
}
