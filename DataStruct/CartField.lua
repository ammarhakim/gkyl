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
local AllocShared      = require "Lib.AllocShared"
local CartDecompNeigh  = require "Lib.CartDecompNeigh"
local Grid             = require "Grid.RectCart"
local Lin              = require "Lib.Linalg"
local LinearDecomp     = require "Lib.LinearDecomp"
local Mpi              = require "Comm.Mpi"
local Range            = require "Lib.Range"
local lume             = require "Lib.lume"

-- Load CUDA allocators (or dummy when CUDA is not found).
local cuda = nil
local cuAlloc = require "Cuda.AllocDummy"
if GKYL_HAVE_CUDA then
   cuAlloc = require "Cuda.Alloc"
   cuda    = require "Cuda.RunTime"
end

-- C interfaces
ffi.cdef [[
    typedef struct {
        int ndim;
        int elemSize;
        int numComponents;
        GkylRange_t *localRange, *localExtRange;
        GkylRange_t *localEdgeRange, *localExtEdgeRange;
        GkylRange_t *globalRange, *globalExtRange;
        GkylRectCart_t *grid;
        double *_data;
    } GkylCartField_t;

    // s: start index. nv: number of values.
    void gkylCartFieldAccumulate(unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldAssign(unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldScale(unsigned s, unsigned nv, double fact, double *out);
    void gkylCartFieldScaleByCell(unsigned s, unsigned nv, unsigned ncomp, double *fact, double *out);
    void gkylCartFieldAbs(unsigned s, unsigned nv, double *out);
    void gkylCopyFromField(double *data, double *f, unsigned numComponents, unsigned c);
    void gkylCopyToField(double *f, double *data, unsigned numComponents, unsigned c);
    void gkylCartFieldAssignAll(unsigned s, unsigned nv, double val, double *out);

    // sInp/sOut: start index for input/output fields. nCells: number of cells being looped over. 
    // compStart: starting component for offset. nCompInp/nCompOut: input/output field's number of components.
    void gkylCartFieldAssignOffset(unsigned sInp, unsigned sOut, unsigned nCells, unsigned compStart, unsigned nCompInp, unsigned nCompOut, double fact, const double *inp, double *out);
    void gkylCartFieldAccumulateOffset(unsigned sInp, unsigned sOut, unsigned nCells, unsigned compStart, unsigned nCompInp, unsigned nCompOut, double fact, const double *inp, double *out);
    void gkylCartFieldAssignOffset(unsigned sInp, unsigned sOut, unsigned nCells, unsigned compStart, unsigned nCompInp, unsigned nCompOut, double fact, const double *inp, double *out);

    void gkylCartFieldDeviceAccumulate(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldDeviceAccumulateOffset(int numBlocks, int numThreads, unsigned sInp, unsigned sOut, unsigned nCells, unsigned compStart, unsigned nCompInp, unsigned nCompOut, double fact, const double *inp, double *out);
    void gkylCartFieldDeviceAssign(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldDeviceScale(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, double *out);
    void gkylCartFieldDeviceAbs(int numBlocks, int numThreads, unsigned s, unsigned nv, double *out);

    // copy periodic boundary conditions when using only one GPU
    void gkylDevicePeriodicCopy(int numBlocks, int numThreads, GkylRange_t *rangeSkin, GkylRange_t *rangeGhost, GkylCartField_t *f, unsigned numComponents);

    // Assign all elements to specified value.
    void gkylCartFieldDeviceAssignAll(int numBlocks, int numThreads, unsigned s, unsigned nv, double val, double *out);
]]

if GKYL_HAVE_CUDA then
   ffi.cdef [[
    // Reduction down to a single value (e.g. min, max, sum).
    void reductionBlocksAndThreads(GkDeviceProp *prop, int numElements, int maxBlocks,
                                   int maxThreads, int &blocks, int &threads);
    typedef double (*redBinOpFunc_t)(double a, double b);
    typedef struct {
      double initValue;
      redBinOpFunc_t reduceFunc;
    } baseReduceOp_t; 
    redBinOpFunc_t getRedMinFuncFromDevice();
    redBinOpFunc_t getRedMaxFuncFromDevice();
    redBinOpFunc_t getRedSumFuncFromDevice();
    void gkylCartFieldDeviceReduce(baseReduceOp_t *redOp, int numCellsTot, int numComponents, int numBlocks, int numThreads, int maxBlocks, int maxThreads, GkDeviceProp *prop, GkylCartField_t *fIn, double *blockOut, double *intermediate, double *out);

   ]]
end

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
-- Helper function to check if two fields have the same range.
-- Useful when manipulating two fields with different number of components, but same range.
local function field_check_range(y, x)
   return y:localRange() == x:localRange()
end

-- Return local start and num times to bump.
local function getStartAndBump(self, decomp)
   if self._layout == colMajLayout then
      return decomp:colStartIndex(self._shmIndex), decomp:shape(self._shmIndex)
   end
   return decomp:rowStartIndex(self._shmIndex), decomp:shape(self._shmIndex)
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

   -- functions for regular, shared and device memory allocations
   local allocFunc = Alloc.Alloc_meta_ctor(elct)
   local allocSharedFunc = AllocShared.AllocShared_meta_ctor(elct)
   local allocCudaFunc = cuAlloc.Alloc_meta_ctor(elct, false) -- don't used managed memory
   
   -- allocator for use in non-shared applications
   local function allocatorFunc(comm, numElem)
      return allocFunc(numElem)
   end
   -- allocator for use in shared applications
   local function sharedAllocatorFunc(comm, numElem)
      return allocSharedFunc(comm, numElem)
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

      -- Various communicators for use in shared allocator.
      local nodeComm = grid:commSet().nodeComm
      local shmComm  = grid:commSet().sharedComm

      -- Allocator function.
      local allocator = grid:isShared() and sharedAllocatorFunc or allocatorFunc
      
      -- Allocate memory: this is NOT managed by the LuaJIT GC, allowing fields to be arbitrarly large.
      local sz = localRange:extend(ghost[1], ghost[2]):volume()*nc -- Amount of data in field.
      self._allocData = allocator(shmComm, sz) -- Store this so it does not vanish under us.
      self._data      = self._allocData:data() -- Pointer to data.
      
      -- Setup object.
      self._grid = grid
      self._ndim = grid:ndim()
      self._lowerGhost, self._upperGhost = ghost[1], ghost[2]
      self._numComponents = nc
      self._size = sz

      self._globalRange = globalRange
      self._globalExtRange = self._globalRange:extend(
	 self._lowerGhost, self._upperGhost)
      
      self._localRange = localRange
      self._localExtRange = self._localRange:extend(
	 self._lowerGhost, self._upperGhost)

      -- All real-cell edges.
      self._localEdgeRange = self._localRange:extend(1, 0) -- Or (1, 0)?

      -- All cell-cell edges, including those of a ghost cell.
      self._localExtEdgeRange = self._localRange:extend(
	 self._lowerGhost-1, self._upperGhost)

      -- Local and (MPI) global values of a reduction (reduce method).
      self.localReductionVal  = ElemVec(self._numComponents)
      self.globalReductionVal = ElemVec(self._numComponents)

      -- Create a device copy is needed.
      local createDeviceCopy = xsys.pickBool(tbl.createDeviceCopy, GKYL_USE_DEVICE)
      if createDeviceCopy then
         -- Allocate device memory.
	 self._devAllocData = deviceAllocatorFunc(shmComm, sz)
         -- Package data and info into struct on device.
         local f = ffi.new("GkylCartField_t")
         local sz = sizeof("GkylCartField_t")
         f.ndim = self._ndim
         f.elemSize = elctSize
         f.numComponents = self._numComponents
         f._data = self._devAllocData:data()
         f.localRange = Range.copyHostToDevice(self._localRange)
         f.localExtRange = Range.copyHostToDevice(self._localExtRange)
         f.localEdgeRange = Range.copyHostToDevice(self._localEdgeRange)
         f.localExtEdgeRange = Range.copyHostToDevice(self._localExtEdgeRange)
         f.globalRange = Range.copyHostToDevice(self._globalRange)
         f.globalExtRange = Range.copyHostToDevice(self._globalExtRange)
         f.grid = self._grid._onDevice
         self._onDevice, err = cuda.Malloc(sz)
         cuda.Memcpy(self._onDevice, f, sz, cuda.MemcpyHostToDevice)

         local devNum, _     = cuda.GetDevice()
         self.deviceProps, _ = cuda.GetDeviceProperties(devNum)
         -- Establish number of blocks and threads/block for deviceReduce, and allocate memory.
         self.reduceBlocksMAX  = 64
         self.reduceThreadsMAX = GKYL_DEFAULT_NUM_THREADS
         local numBlocksC, numThreadsC = Alloc.Int(1), Alloc.Int(1)
         ffiC.reductionBlocksAndThreads(self.deviceProps,self._localRange:volume(),self.reduceBlocksMAX,
                                        self.reduceThreadsMAX,numBlocksC:data(),numThreadsC:data());
         self.reduceBlocks  = numBlocksC[1]
         self.reduceThreads = numThreadsC[1]
        
         numBlocksC:delete()
         numThreadsC:delete()
         self.d_blockRed, self.d_intermediateRed = deviceAllocatorFunc(shmComm, self.reduceBlocks), deviceAllocatorFunc(shmComm, self.reduceBlocks)
         -- Create reduction operator on host, and copy to device.
         local redOp           = {}
         for k, v in pairs(reduceInitialVal) do
           redOp[k]            = ffi.new("baseReduceOp_t")
           redOp[k].initValue  = v
         end
         redOp["min"].reduceFunc = ffi.C.getRedMinFuncFromDevice()
         redOp["max"].reduceFunc = ffi.C.getRedMaxFuncFromDevice()
         redOp["sum"].reduceFunc = ffi.C.getRedSumFuncFromDevice()
         local sz = ffi.sizeof("baseReduceOp_t")
         self.d_redOp = {min = cuda.Malloc(sz), max = cuda.Malloc(sz), sum = cuda.Malloc(sz)}
         for k, _ in pairs(reduceInitialVal) do
            err = cuda.Memcpy(self.d_redOp[k], redOp[k], sz, cuda.MemcpyHostToDevice)
         end
      end
      if not GKYL_HAVE_CUDA then self._devAllocData = nil end
      
      self._layout = defaultLayout -- Default layout is column-major.
      if tbl.layout then
         self._layout = tbl.layout=="row-major" and rowMajLayout or colMajLayout
      end

      self._shmIndex = Mpi.Comm_rank(shmComm)+1 -- Our local index on SHM comm (one more than rank).

      -- Construct linear decomposition of various ranges.
      self._localRangeDecomp = LinearDecomp.LinearDecompRange {
         range    = localRange,
         numSplit = Mpi.Comm_size(shmComm)
      }
      self._localExtRangeDecomp = LinearDecomp.LinearDecompRange {
         range    = localRange:extend(self._lowerGhost, self._upperGhost),
         numSplit = Mpi.Comm_size(shmComm)
      }

      -- Store start index and size handled by local SHM-rank for local and extended range.
      self._localStartIdx, self._localNumBump       = getStartAndBump(self, self._localRangeDecomp)
      self._localExtStartIdx, self._localExtNumBump = getStartAndBump(self, self._localExtRangeDecomp)

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
      local localExtRange   = self._localExtRange
      local indexer         = self:genIndexer()
      local decomposedRange = self._grid:decomposedRange()
      local myId            = self._grid:subGridId() -- Grid ID on this processor.
      local neigIds         = self._decompNeigh:neighborData(myId) -- List of neighbors.

      for _, sendId in ipairs(neigIds) do
         local neighRgn = decomposedRange:subDomain(sendId)
         local sendRgn  = localRange:intersect(
            neighRgn:extend(self._lowerGhost, self._upperGhost))
         local idx      = sendRgn:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self._sendMPILoc[sendId]      = (indexer(idx)-1)*self._numComponents
         self._sendMPIDataType[sendId] = Mpi.createDataTypeFromRangeAndSubRange(
            sendRgn, localExtRange, self._numComponents, self._layout, elctCommType)
      end

      for _, recvId in ipairs(neigIds) do
         local neighRgn = decomposedRange:subDomain(recvId)
         local recvRgn  = localExtRange:intersect(neighRgn)
         local idx      = recvRgn:lowerAsVec()
         -- set idx to starting point of region you want to recv
         self._recvMPILoc[recvId]      = (indexer(idx)-1)*self._numComponents
         self._recvMPIDataType[recvId] = Mpi.createDataTypeFromRangeAndSubRange(
            recvRgn, localExtRange, self._numComponents, self._layout, elctCommType)
      end

      -- Create MPI DataTypes for periodic directions.
      -- Also store location in memory required for sending/receiving periodic data.
      self._sendLowerPerMPIDataType, self._recvLowerPerMPIDataType = {}, {}
      self._sendUpperPerMPIDataType, self._recvUpperPerMPIDataType = {}, {}
      self._sendLowerPerMPILoc, self._recvLowerPerMPILoc = {}, {}
      self._sendUpperPerMPILoc, self._recvUpperPerMPILoc = {}, {}

      -- Create buffers for periodic copy if Mpi.Comm_size(nodeComm) = 1.
      -- Note that since nodeComm is only valid on shmComm = 0, need to check whether this is a valid comm.
      if Mpi.Comm_size(shmComm) == 1 and Mpi.Comm_size(nodeComm) == 1 then
         self._lowerPeriodicBuff, self._upperPeriodicBuff = {}, {}
      end

      -- Following loop creates Datatypes for periodic directions.
      -- This is complicated as one needs to treat lower -> upper
      -- transfers differently than upper -> lower as the number of
      -- ghost cells may be different on each lower/upper side. (AHH)
      for dir = 1, self._ndim do
         if grid:isDirPeriodic(dir) then
            local skelIds = decomposedRange:boundarySubDomainIds(dir)
            for i = 1, #skelIds do
               local loId, upId = skelIds[i].lower, skelIds[i].upper

               -- Only create if we are on proper ranks.
               -- Note that if the node communicator has rank size of 1, then we can access all the 
               -- memory needed for periodic boundary conditions and we do not need MPI Datatypes.
               if myId == loId then
                  local rgnSend = decomposedRange:subDomain(loId):lowerSkin(dir, self._upperGhost)
                  if Mpi.Comm_size(shmComm) == 1 and Mpi.Comm_size(nodeComm) == 1 then
                     local szSend = rgnSend:volume()*self._numComponents
                     self._lowerPeriodicBuff[dir] = allocator(shmComm, szSend)
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
               end
               if myId == upId then
                  local rgnSend = decomposedRange:subDomain(upId):upperSkin(dir, self._lowerGhost)
                  if Mpi.Comm_size(shmComm) == 1 and Mpi.Comm_size(nodeComm) == 1 then
                     local szSend = rgnSend:volume()*self._numComponents
                     self._upperPeriodicBuff[dir] = allocator(shmComm, szSend)
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
      for dir = 1, self._ndim do
         self._sendLowerCornerPerMPILoc[dir]     , self._sendUpperCornerPerMPILoc[dir]      = {}, {}
         self._sendLowerCornerPerMPIDataType[dir], self._sendUpperCornerPerMPIDataType[dir] = {}, {}
         self._recvLowerCornerPerMPILoc[dir]     , self._recvUpperCornerPerMPILoc[dir]      = {}, {}
         self._recvLowerCornerPerMPIDataType[dir], self._recvUpperCornerPerMPIDataType[dir] = {}, {}
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
      ndim = function (self)
	 return self._ndim
      end,
      grid = function (self)
	 return self._grid
      end,
      numComponents = function (self)
	 return self._numComponents
      end,
      copy = function (self, fIn)
	 self:_assign(1.0, fIn)
      end,
      copyOffset = function (self, fIn, compStart)
        self:_assignOffset(1.0, fIn, compStart)
      end,
      deviceCopy = function (self, fIn)
         if self._devAllocData then
	    self:_deviceAssign(1.0, fIn)
	 end
      end,
      copyHostToDevice = function (self)
	 if self._devAllocData then
	    return self._devAllocData:copyHostToDevice(self._allocData)
	 end
	 return 0
      end,
      copyDeviceToHost = function (self)
	 if self._devAllocData then
	    return self._devAllocData:copyDeviceToHost(self._allocData)
	 end
	 return 0
      end,
      deviceData = function (self)
	 return self._devAllocData
      end,
      deviceDataPointer = function (self)
	 return self._devAllocData:data()
      end,
      dataPointer = function (self)
	 return self._allocData:data()
      end,
      clear = function (self, val)
	 ffiC.gkylCartFieldAssignAll(self:_localLower(), self:_localShape(), val, self._data)
      end,
      deviceClear = function (self, val)
         if self._devAllocData then
	    local numThreads = GKYL_DEFAULT_NUM_THREADS
	    local shape = self:_localShape()
	    local numBlocks = math.floor(shape/numThreads)+1
	    ffiC.gkylCartFieldDeviceAssignAll(numBlocks, numThreads, self:_localLower(), self:_localShape(), val, self:deviceDataPointer())
	 end
      end,
      fill = function (self, k, fc)
	 local loc = (k-1)*self._numComponents -- (k-1) as k is 1-based index	 
	 fc._cdata = self._data+loc
      end,
      _localLower = function (self)
	 return (self._localExtRangeDecomp:lower(self._shmIndex)-1)*self:numComponents()
      end,
      _localShape = function (self)
	 return self._localExtRangeDecomp:shape(self._shmIndex)*self:numComponents()
      end,
      _assign = function(self, fact, fld)
	 assert(field_compatible(self, fld), "CartField:assign: Can only assign compatible fields")
	 assert(type(fact) == "number", "CartField:assign: Factor not a number")

	 ffiC.gkylCartFieldAssign(self:_localLower(), self:_localShape(), fact, fld._data, self._data)
      end,
      _assignOffset = function(self, fact, fld, compStart)
   assert(field_check_range(self, fld),
    "CartField:assignOffset: Can only assign fields with the same range")
   assert(type(fact) == "number",
    "CartField:assignOffset: Factor not a number")
         assert(self:layout() == fld:layout(),
    "CartField:assignOffset: Fields should have same layout for assignment to make sense")

         -- Get number of cells for outer loop.
         -- We do not need to use an indexer since we are simply accumulating cell-wise a subset of the components.
         local numCells = self:_localShape()/self:numComponents()
	 ffiC.gkylCartFieldAssignOffset(fld:_localLower(), self:_localLower(), numCells, compStart, fld:numComponents(), self:numComponents(), fact, fld._data, self._data)
      end,
      _deviceAssign = function(self, fact, fld)
	 assert(field_compatible(self, fld), "CartField:deviceAssign: Can only accumulate compatible fields")
	 assert(type(fact) == "number", "CartField:deviceAssign: Factor not a number")

	 local numThreads = GKYL_DEFAULT_NUM_THREADS
	 local shape = self:_localShape()
	 local numBlocks = math.floor(shape/numThreads)+1
	 ffiC.gkylCartFieldDeviceAssign(numBlocks, numThreads, self:_localLower(), self:_localShape(), fact, fld:deviceDataPointer(), self:deviceDataPointer())
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

         -- Get number of cells for outer loop.
         -- We do not need to use an indexer since we are simply accumulating cell-wise a subset of the components.
         local numCells = self:_localShape()/self:numComponents()
	 ffiC.gkylCartFieldAssignOffset(fld:_localLower(), self:_localLower(), numCells, compStart, fld:numComponents(), self:numComponents(), fact, fld._data, self._data)
      end,
      _accumulateOneFld = function(self, fact, fld)
	 assert(field_compatible(self, fld),
		"CartField:accumulateOneFld: Can only accumulate compatible fields")
	 assert(type(fact) == "number",
		"CartField:accumulateOneFld: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:accumulateOneFld: Fields should have same layout for sums to make sense")

	 ffiC.gkylCartFieldAccumulate(self:_localLower(), self:_localShape(), fact, fld._data, self._data)
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

         -- Get number of cells for outer loop.
         -- We do not need to use an indexer since we are simply accumulating cell-wise a subset of the components.
         local numCells = self:_localShape()/self:numComponents()
	 ffiC.gkylCartFieldAccumulateOffset(fld:_localLower(), self:_localLower(), numCells, compStart, fld:numComponents(), self:numComponents(), fact, fld._data, self._data)
      end,
      _deviceAccumulateOneFld = function(self, fact, fld)
	 assert(field_compatible(self, fld),
		"CartField:deviceAccumulateOneFld: Can only accumulate compatible fields")
	 assert(type(fact) == "number",
		"CartField:deviceAccumulateOneFld: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:deviceAccumulateOneFld: Fields should have same layout for sums to make sense")

	 local numThreads = GKYL_DEFAULT_NUM_THREADS
	 local shape = self:_localShape()
	 local numBlocks = math.floor(shape/numThreads)+1
	 ffiC.gkylCartFieldDeviceAccumulate(numBlocks, numThreads, self:_localLower(), self:_localShape(), fact, fld:deviceDataPointer(), self:deviceDataPointer())
      end,
      _deviceAccumulateOffsetOneFld = function(self, fact, fld, compStart)
	 assert(field_check_range(self, fld),
		"CartField:accumulateOffset: Can only accumulate fields with the same range")
	 assert(type(fact) == "number",
		"CartField:accumulateOffset: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:accumulateOffset: Fields should have same layout for sums to make sense")

         -- Get number of cells for outer loop.
         -- We do not need to use an indexer since we are simply accumulating cell-wise a subset of the components.
         local numCells = self:_localShape()/self:numComponents()

         local numThreads = GKYL_DEFAULT_NUM_THREADS
         local numBlocks = math.floor(numCells/numThreads)+1
	 ffiC.gkylCartFieldDeviceAccumulateOffset(numBlocks, numThreads, fld:_localLower(), self:_localLower(), numCells, compStart, fld:numComponents(), self:numComponents(), fact, fld:deviceDataPointer(), self:deviceDataPointer())
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
      deviceAccumulate = isNumberType and
	 function (self, c1, fld1, ...)
	    if self._devAllocData then
	       local args = {...} -- Package up rest of args as table.
	       local nFlds = #args/2
	       self:_deviceAccumulateOneFld(c1, fld1) -- Accumulate first field.
	       for i = 1, nFlds do -- Accumulate rest of the fields.
	          self:_deviceAccumulateOneFld(args[2*i-1], args[2*i])
	       end
	    end
	 end or
	 function (self, c1, fld1, ...)
	    assert(false, "CartField:deviceAccumulate: Accumulate only works on numeric fields")
	 end,
      deviceAccumulateOffset = isNumberType and
	 function (self, c1, fld1, compStart1, ...)
            if self._devAllocData then
	       local args = {...} -- Package up rest of args as table.
	       local nFlds = #args/3
	       self:_deviceAccumulateOffsetOneFld(c1, fld1, compStart1) -- Accumulate first field.
	       for i = 1, nFlds do -- Accumulate rest of the fields.
	          self:_deviceAccumulateOffsetOneFld(args[3*i-2], args[3*i-1], args[3*i])
	       end
            end
	 end or
	 function (self, c1, fld1, compStart1, ...)
	    assert(false, "CartField:deviceAccumulateOffset: Accumulate only works on numeric fields")
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
      deviceCombine = isNumberType and
	 function (self, c1, fld1, ...)
	    if self._devAllocData then
	       local args = {...} -- Package up rest of args as table.
	       local nFlds = #args/2
	       self:_deviceAssign(c1, fld1) -- Assign first field.
	       for i = 1, nFlds do -- Accumulate rest of the fields.
	          self:_deviceAccumulateOneFld(args[2*i-1], args[2*i])
	       end
	    end
	 end or
	 function (self, c1, fld1, ...)
	    assert(false, "CartField:deviceCombine: Combine only works on numeric fields")
	 end,
      scale = isNumberType and
	 function (self, fact)
	    ffiC.gkylCartFieldScale(self:_localLower(), self:_localShape(), fact, self._data)
	 end or
	 function (self, fact)
	    assert(false, "CartField:scale: Scale only works on numeric fields")
	 end,
      deviceScale = isNumberType and
	 function (self, fact)
	    if self._devAllocData then
	       local numThreads = GKYL_DEFAULT_NUM_THREADS
	       local shape = self:_localShape()
	       local numBlocks = math.floor(shape/numThreads)+1
	       ffiC.gkylCartFieldDeviceScale(numBlocks, numThreads,  self:_localLower(), self:_localShape(), fact, self:deviceDataPointer())
	    end
	 end or
	 function (self, fact)
	    assert(false, "CartField:deviceScale: Scale only works on numeric fields")
	 end,
      scaleByCell = function (self, factByCell)
         assert(factByCell:numComponents() == 1, "CartField:scaleByCell: scalar must be a 1-component field")
         assert(factByCell:localRange() == self:localRange() and factByCell:layout() == self:layout(), "CartField:scaleByCell: scalar and field must be compatible")
            
	 ffiC.gkylCartFieldScaleByCell(self:_localLower(), self:_localShape(), self:numComponents(), factByCell._data, self._data)
      end,
      abs = isNumberType and
         function (self)
            ffiC.gkylCartFieldAbs(self:_localLower(), self:_localShape(), self._data)
	 end or
	 function (self, fact)
	    assert(false, "CartField:abs: Abs only works on numeric fields")
	 end,
      deviceAbs = isNumberType and
	 function (self)
	    if self._devAllocData then
	       local numThreads = GKYL_DEFAULT_NUM_THREADS
	       local shape = self:_localShape() 
	       local numBlocks = math.floor(shape/numThreads)+1
	       ffiC.gkylCartFieldDeviceAbs(numBlocks, numThreads,  self:_localLower(), self:_localShape(), self:deviceDataPointer())
	    end
	 end or
	 function (self, fact)
	    assert(false, "CartField:deviceAbs: Abs only works on numeric fields")
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
	 -- this barrier is needed as when using MPI-SHM some
	 -- processors will get to the sync method before others
         -- this is especially troublesome in the RK combine step
	 Mpi.Barrier(self._grid:commSet().sharedComm)
	 self._field_sync(self, self:dataPointer())
	 if self._syncPeriodicDirs and syncPeriodicDirs then
	    self._field_periodic_sync(self, self:dataPointer())
	 end
	 -- this barrier is needed as when using MPI-SHM some
	 -- processors will not participate in sync()
	 Mpi.Barrier(self._grid:commSet().sharedComm)
      end,
      deviceSync = function (self, syncPeriodicDirs_)
         local syncPeriodicDirs = xsys.pickBool(syncPeriodicDirs_, true)
	 -- this barrier is needed as when using MPI-SHM some
	 -- processors will get to the sync method before others
         -- this is especially troublesome in the RK combine step
	 Mpi.Barrier(self._grid:commSet().sharedComm)
	 self._field_sync(self, self:deviceDataPointer())
	 if self._syncPeriodicDirs and syncPeriodicDirs then
	    self._field_periodic_sync(self, self:deviceDataPointer())
	 end
	 -- this barrier is needed as when using MPI-SHM some
	 -- processors will not participate in sync()
	 Mpi.Barrier(self._grid:commSet().sharedComm)
      end,
      -- This method is an alternative function for applying periodic boundary
      -- conditions when Mpi.Comm_size(nodeComm) = 1 and we do not need to call
      -- Send/Recv to copy the skin cell data into ghost cells.
      periodicCopy = function (self, syncPeriodicDirs_)
         local syncPeriodicDirs = xsys.pickBool(syncPeriodicDirs_, true)
	 -- this barrier is needed as when using MPI-SHM some
	 -- processors could apply periodic boundary conditions before others
         -- this is especially troublesome in the RK combine step
	 Mpi.Barrier(self._grid:commSet().sharedComm)
	 if self._syncPeriodicDirs and syncPeriodicDirs then
            self._field_periodic_copy(self)
	 end
	 -- This barrier is needed as when using MPI-SHM some
	 -- processors will not participate in applying periodic
         -- boundary conditions
	 Mpi.Barrier(self._grid:commSet().sharedComm)
      end,
      devicePeriodicCopy = function (self, syncPeriodicDirs_)
         local syncPeriodicDirs = xsys.pickBool(syncPeriodicDirs_, true)
         if self._syncPeriodicDirs and syncPeriodicDirs then
            local numThreads = GKYL_DEFAULT_NUM_THREADS

            local grid = self._grid
            local decomposedRange = grid:decomposedRange()
            for dir = 1, self._ndim do
               if grid:isDirPeriodic(dir) then
                  -- First get region for skin cells for upper ghost region (the lower skin cells).
                  local skinRgnUpper = decomposedRange:subDomain(1):lowerSkin(dir, self._upperGhost)
                  local ghostRgnUpper = decomposedRange:subDomain(1):upperGhost(dir, self._upperGhost)

                  -- Both skinRgnUpper and ghostRgnUpper have the same shape, so pick one to set shape.
                  -- and thus the number of blocks and threads
                  local shape = skinRgnUpper:shape(self._shmIndex)
                  local numBlocks = math.floor(shape/numThreads)+1

                  -- Copy lower skin cells into upper ghost cells.
                  ffiC.gkylDevicePeriodicCopy(numBlocks, numThreads, skinRgnUpper, ghostRgnUpper, self._onDevice, self._numComponents)

	          -- Now do the same, but for the skin cells for the lower ghost region (the upper skin cells).
                  local skinRgnLower = decomposedRange:subDomain(1):upperSkin(dir, self._lowerGhost)
                  local ghostRgnLower = decomposedRange:subDomain(1):lowerGhost(dir, self._lowerGhost)

                  -- skinRgnLower and ghostRgnLower have the same shape, and potentially a different shape than skin/ghostRgnUpper.
                  -- Pick one to set shape and thus the number of blocks and threads.
                  shape = skinRgnLower:shape(self._shmIndex)
                  numBlocks = math.floor(shape/numThreads)+1

                  -- Copy lower skin cells into upper ghost cells.
                  ffiC.gkylDevicePeriodicCopy(numBlocks, numThreads, skinRgnLower, ghostRgnLower, self._onDevice, self._numComponents)
               end
            end
         end
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
	    -- Input 'opIn' must be one of the binary operations in binOpFuncs.
	    local grid = self._grid
	    local tId = grid:subGridSharedId() -- Local thread ID.
	    local localRangeDecomp = LinearDecomp.LinearDecompRange {
	       range = self._localRange, numSplit = grid:numSharedProcs() }
	    local indexer = self:genIndexer()
	    local itr = self:get(1)
	    
	    local localVal = {}
	    for k = 1, self._numComponents do localVal[k] = reduceInitialVal[opIn] end
	    for idx in localRangeDecomp:rowMajorIter(tId) do
	       self:fill(indexer(idx), itr)
	       for k = 1, self._numComponents do
		  localVal[k] = binOpFuncs[opIn](localVal[k], itr:data()[k-1])
	       end
	    end

	    for k = 1, self._numComponents do self.localReductionVal[k] = localVal[k] end
	    Mpi.Allreduce(self.localReductionVal:data(), self.globalReductionVal:data(),
			  self._numComponents, elctCommType, reduceOpsMPI[opIn], grid:commSet().comm)

	    for k = 1, self._numComponents do localVal[k] = self.globalReductionVal[k] end
            return localVal
	 end or
	 function (self, opIn)
	    assert(false, "CartField:reduce: Reduce only works on numeric fields")
	 end,
      deviceReduce = isNumberType and
	 function(self, opIn, d_reduction)
	    -- Input 'opIn' must be one of the binary operations in binOpFuncs.
            ffi.C.gkylCartFieldDeviceReduce(self.d_redOp[opIn],self._localRange:volume(),self._numComponents,self.reduceBlocks,self.reduceThreads,self.reduceBlocksMAX,self.reduceThreadsMAX,
               self.deviceProps,self._onDevice,self.d_blockRed:data(),self.d_intermediateRed:data(),d_reduction:data())
	 end or
	 function (self, opIn, d_reduction)
	    assert(false, "CartField:deviceReduce: Reduce only works on numeric fields.")
	 end,
      _copy_from_field_region = function (self, rgn, data)
	 local indexer = self:genIndexer()
	 local c = 0
         local fitr = self:get(1)
	 for idx in rgn:rowMajorIter() do
	    self:fill(indexer(idx), fitr)
            ffiC.gkylCopyFromField(data:data(), fitr:data(), self._numComponents, c)
            c = c + self._numComponents
	 end
      end,
      _copy_to_field_region = function (self, rgn, data)
	 local indexer = self:genIndexer()
	 local c = 0
         local fitr = self:get(1)
	 for idx in rgn:rowMajorIter() do
	    self:fill(indexer(idx), fitr)
            ffiC.gkylCopyToField(fitr:data(), data:data(), self._numComponents, c)
            c = c + self._numComponents
	 end
      end,
      _field_sync = function (self, dataPtr)
         local comm = self._grid:commSet().nodeComm -- communicator to use
         if not Mpi.Is_comm_valid(comm) then
            return -- no need to do anything if communicator is not valid
         end
         -- immediately return if nothing to sync
         if self._lowerGhost == 0 and self._upperGhost == 0 then return end
        
         -- Steps: (1) Post non-blocking recv requests. (2) Do
         -- blocking sends, (3) Complete recv and copy data into ghost
         -- cells.
        
         local myId    = self._grid:subGridId() -- grid ID on this processor
         local neigIds = self._decompNeigh:neighborData(myId) -- list of neighbors
         local tag     = 42 -- Communicator tag for regular (non-periodic) messages
         local recvReq = {} -- list of recv requests
         -- post a non-blocking recv request
         for _, recvId in ipairs(neigIds) do
            local dataType = self._recvMPIDataType[recvId]
            local loc      = self._recvMPILoc[recvId]
            -- recv data: (its from recvId-1 as MPI ranks are zero indexed)
            recvReq[recvId] = Mpi.Irecv(dataPtr+loc, 1, dataType, recvId-1, tag, comm)
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
         for _, recvId in ipairs(neigIds) do Mpi.Wait(recvReq[recvId], nil) end
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
         
         local recvUpperReq, recvLowerReq = {}, {}
         local recvUpperCornerReq, recvLowerCornerReq = {}, {}
         -- Post non-blocking recv requests for periodic directions.
         for dir = 1, self._ndim do
            recvUpperCornerReq[dir], recvLowerCornerReq[dir] = {}, {}
            if grid:isDirPeriodic(dir) then
               local skelIds = decomposedRange:boundarySubDomainIds(dir)
               for i = 1, #skelIds do
                  local loId, upId = skelIds[i].lower, skelIds[i].upper
         
                  if myId == loId then
                     local loTag       = basePerTag+dir+10
                     local dataType    = self._recvLowerPerMPIDataType[dir]
                     local loc         = self._recvLowerPerMPILoc[dir]
                     recvLowerReq[dir] = Mpi.Irecv(dataPtr+loc, 1, dataType,
                                                   upId-1, loTag, comm)
                  end
                  if myId == upId then
                     local upTag       = basePerTag+dir
                     local dataType    = self._recvUpperPerMPIDataType[dir]
                     local loc         = self._recvUpperPerMPILoc[dir]
                     recvUpperReq[dir] = Mpi.Irecv(dataPtr+loc, 1, dataType,
                                                   loId-1, upTag, comm)
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
                           recvLowerCornerReq[dir][ccLo] = Mpi.Irecv(dataPtr+loc, 1, dataType,
                                                                     upId-1, loTag, comm)
                        end
                        if myId == upId then
                           ccUp = ccUp+1
                           local upTag    = cornerBasePerTag+tonumber(loId..upId..tblToStr(corDirs)..2)
                           local dataType = self._recvUpperCornerPerMPIDataType[dir][ccUp]
                           local loc      = self._recvUpperCornerPerMPILoc[dir][ccUp]
                           recvUpperCornerReq[dir][ccUp] = Mpi.Irecv(dataPtr+loc, 1, dataType,
                                                                     loId-1, upTag, comm)
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
                  if myId == loId then Mpi.Wait(recvLowerReq[dir], nil) end
                  if myId == upId then Mpi.Wait(recvUpperReq[dir], nil) end
               end

               if self._syncCorners then
                  local cTs = self._cornersToSync[dir]
                  local ccLo, ccUp = 0, 0
                  for dI, bD in ipairs(cTs) do   -- Loop over lower boundary subdomains.
                     for _, dC in ipairs(bD) do   -- Loop over corners.
                        local loId, upId = dC.lower, dC.upper
                        if myId == loId then
                           ccLo = ccLo+1
                           Mpi.Wait(recvLowerCornerReq[dir][ccLo], nil)
                        end
                        if myId == upId then
                           ccUp = ccUp+1
                           Mpi.Wait(recvUpperCornerReq[dir][ccUp], nil)
                        end
                     end
                  end
               end
            end
         end
      end,
      _field_periodic_copy = function (self)
         local grid = self._grid
         local decomposedRange = grid:decomposedRange()
         for dir = 1, self._ndim do
            if grid:isDirPeriodic(dir) then
               -- First get region for skin cells for upper ghost region (the lower skin cells).
               local skinRgnUpper = decomposedRange:subDomain(1):lowerSkin(dir, self._upperGhost)
               local periodicBuffUpper = self._upperPeriodicBuff[dir]
               -- Copy skin cells into temporary buffer.
               self:_copy_from_field_region(skinRgnUpper, periodicBuffUpper)
               -- Get region for looping over upper ghost cells and copy lower skin cells into upper ghost cells.
               local ghostRgnUpper = decomposedRange:subDomain(1):upperGhost(dir, self._upperGhost)
               self:_copy_to_field_region(ghostRgnUpper, periodicBuffUpper)

	       -- Now do the same, but for the skin cells for the lower ghost region (the upper skin cells).
               local skinRgnLower = decomposedRange:subDomain(1):upperSkin(dir, self._lowerGhost)
               local periodicBuffLower = self._lowerPeriodicBuff[dir]
               self:_copy_from_field_region(skinRgnLower, periodicBuffLower)
               local ghostRgnLower = decomposedRange:subDomain(1):lowerGhost(dir, self._lowerGhost)
               self:_copy_to_field_region(ghostRgnLower, periodicBuffLower)
            end
         end
      end,
   }
   
   return Field
end

return {
   new_field_ct = Field_meta_ctor,
   Field = Field_meta_ctor(typeof("double")),
}
