-- Gkyl ------------------------------------------------------------------------
--
-- Multi-component fields on cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

-- Gkyl libraries
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Alloc = require "Lib.Alloc"
local AllocShared = require "Lib.AllocShared"
local CartDecompNeigh = require "Lib.CartDecompNeigh"
local Grid = require "Grid.RectCart"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Mpi = require "Comm.Mpi"
local Range = require "Lib.Range"

-- load CUDA allocators (or dummy when CUDA is not found)
local cuda = nil
local cuAlloc = require "Cuda.AllocDummy"
if GKYL_HAVE_CUDA then
   cuAlloc = require "Cuda.Alloc"
   cuda = require "Cuda.RunTime"
end

-- C interfaces
ffi.cdef [[
    // s: start index. nv: number of values
    void gkylCartFieldAccumulate(unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldAssign(unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldScale(unsigned s, unsigned nv, double fact, double *out);
    void gkylCartFieldScaleByCell(unsigned s, unsigned nv, unsigned ncomp, double *fact, double *out);
    void gkylCartFieldAbs(unsigned s, unsigned nv, double *out);
    void gkylCopyFromField(double *data, double *f, unsigned numComponents, unsigned c);
    void gkylCopyToField(double *f, double *data, unsigned numComponents, unsigned c);
    void gkylCartFieldAssignAll(unsigned s, unsigned nv, double val, double *out);

    void gkylCartFieldDeviceAccumulate(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldDeviceAssign(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldDeviceScale(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, double *out);
    void gkylCartFieldDeviceAbs(int numBlocks, int numThreads, unsigned s, unsigned nv, double *out);

    // copy component data from/to field
    void gkylCopyFromFieldDevice(int numBlocks, int numThreads, double *data, double *f, unsigned numComponents, unsigned c);
    void gkylCopyToFieldDevice(int numBlocks, int numThreads, double *f, double *data, unsigned numComponents, unsigned c);

    // Assign all elements to specified value.
    void gkylCartFieldDeviceAssignAll(int numBlocks, int numThreads, unsigned s, unsigned nv, double val, double *out);

    typedef struct {
        int ndim;
        int elemSize;
        int numComponents;
        GkylRange_t *localRange, *localExtRange;
        GkylRange_t *globalRange, *globalExtRange;
        GkylRectCart_t *grid;
        double *_data;
    } GkylCartField_t;

]]

if GKYL_HAVE_CUDA then
   ffi.cdef [[
    // Reduction down to a single value (e.g. min, max, sum).
    void reductionBlocksAndThreads(GkDeviceProp *prop, int numElements, int maxBlocks,
                                   int maxThreads, int &blocks, int &threads);

    void gkylCartFieldDeviceReduce(const int reduceOp, int numCellsTot, int numBlocks, int numThreads, int maxBlocks, int maxThreads,
       GkDeviceProp *prop, GkylCartField_t *fIn, double *blockOut, double *intermediate, double *out);

   ]]
end

-- Local definitions
local rowMajLayout, colMajLayout = Range.rowMajor, Range.colMajor -- data layout
local indexerMakerFuncs = {} -- list of functions that make indexers
indexerMakerFuncs[rowMajLayout] = Range.makeRowMajorIndexer
indexerMakerFuncs[colMajLayout] = Range.makeColMajorIndexer
-- Default layout
local defaultLayout = rowMajLayout

local genIndexerMakerFuncs = {} -- list of functions that make generic indexers
genIndexerMakerFuncs[rowMajLayout] = Range.makeRowMajorGenIndexer
genIndexerMakerFuncs[colMajLayout] = Range.makeColMajorGenIndexer

-- helper to check if two field are compatible
local function field_compatible(y, x)
   return y:localRange() == x:localRange() and y:numComponents() == x:numComponents()
end

-- return local start and num times to bump
local function getStartAndBump(self, decomp)
   if self._layout == colMajLayout then
      return decomp:colStartIndex(self._shmIndex), decomp:shape(self._shmIndex)
   end
   return decomp:rowStartIndex(self._shmIndex), decomp:shape(self._shmIndex)
end

-- Field accessor object: allows access to field values in cell
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
   local reduceInitialVal = {max = elctMinValue, min = elctMaxValue , sum = 0}
   
   -- make constructor for Field
   local Field = {}
   function Field:new(tbl)
      local self = setmetatable({}, Field)

      -- read data from input table
      local grid = tbl.onGrid
      local nc = tbl.numComponents and tbl.numComponents or 1 -- default numComponents=1
      local ghost = tbl.ghost and tbl.ghost or {0, 0} -- No ghost cells by default

      local syncCorners = xsys.pickBool(tbl.syncCorners, false) -- don't sync corners by default
      self._syncPeriodicDirs = xsys.pickBool(tbl.syncPeriodicDirs, true) -- sync periodic BCs by default

      -- local and global ranges
      local globalRange = grid:globalRange()
      local localRange = grid:localRange()

      -- various communicators for use in shared allocator
      local shmComm = grid:commSet().sharedComm

      -- allocator function
      local allocator = grid:isShared() and sharedAllocatorFunc or allocatorFunc
      
      -- allocate memory: this is NOT managed by the LuaJIT GC, allowing fields to be arbitrarly large
      local sz = localRange:extend(ghost[1], ghost[2]):volume()*nc -- amount of data in field
      self._allocData = allocator(shmComm, sz) -- store this so it does not vanish under us
      self._data = self._allocData:data() -- pointer to data

      -- for number types fill it with zeros (for others, the
      -- assumption is that users will initialize themselves)
      if isNumberType then self._allocData:fill(0) end
      
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

      -- Local and (MPI) global values of a reduction (reduce method).
      self.localReductionVal  = ElemVec(1)
      self.globalReductionVal = ElemVec(1)

      -- Create device memory if needed.
      local createDeviceCopy = xsys.pickBool(tbl.createDeviceCopy, false) -- by default, no device mem allocated
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
         self.d_blockRed, self.d_intermediateRed = cuAlloc.Double(numBlocks), cuAlloc.Double(numBlocks)
      end
      if not GKYL_HAVE_CUDA then self._devAllocData = nil end
      
      self._layout = defaultLayout -- default layout is column-major
      if tbl.layout then
	 if tbl.layout == "row-major" then
	    self._layout = rowMajLayout
	 else
	    self._layout = colMajLayout
	 end
      end

      self._shmIndex = Mpi.Comm_rank(shmComm)+1 -- our local index on SHM comm (one more than rank)

      -- construct linear decomposition of various ranges
      self._localRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = localRange,
	 numSplit = Mpi.Comm_size(shmComm)
      }
      self._localExtRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = localRange:extend(self._lowerGhost, self._upperGhost),
	 numSplit = Mpi.Comm_size(shmComm)
      }

      -- store start index and size handled by local SHM-rank for local and extended range
      self._localStartIdx, self._localNumBump = getStartAndBump(self, self._localRangeDecomp)
      self._localExtStartIdx, self._localExtNumBump = getStartAndBump(self, self._localExtRangeDecomp)

      -- compute communication neighbors
      self._decompNeigh = CartDecompNeigh(grid:decomposedRange())
      if syncCorners then
	 self._decompNeigh:calcAllCommNeigh(ghost[1], ghost[2])
      else
	 self._decompNeigh:calcFaceCommNeigh(ghost[1], ghost[2])
      end

      -- pre-create MPI DataTypes for send/recv calls when doing ghost-cell
      -- sync(). Using MPI DataTypes, we do not require temporary buffers
      -- for send/recv.
      -- also pre-create the location in memory required so that we know 
      -- what parts of the data structure are being sent and received to
      self._sendMPIDataType, self._recvMPIDataType = {}, {}
      self._sendMPILoc, self._recvMPILoc = {}, {}
      local localExtRange = self._localExtRange
      local indexer = self:genIndexer()
      local decomposedRange = self._grid:decomposedRange()
      local myId = self._grid:subGridId() -- grid ID on this processor
      local neigIds = self._decompNeigh:neighborData(myId) -- list of neighbors

      for _, sendId in ipairs(neigIds) do
	 local neighRgn = decomposedRange:subDomain(sendId)
	 local sendRgn = localRange:intersect(
	    neighRgn:extend(self._lowerGhost, self._upperGhost))
         local idx = sendRgn:lowerAsVec()
         -- set idx to starting point of region you want to recv
         self._sendMPILoc[sendId] = (indexer(idx)-1)*self._numComponents
         self._sendMPIDataType[sendId] = Mpi.createDataTypeFromRangeAndSubRange(
	    sendRgn, localExtRange, self._numComponents, self._layout, elctCommType)
      end

      for _, recvId in ipairs(neigIds) do
	 local neighRgn = decomposedRange:subDomain(recvId)
	 local recvRgn = localExtRange:intersect(neighRgn)
         local idx = recvRgn:lowerAsVec()
         -- set idx to starting point of region you want to recv
         self._recvMPILoc[recvId] = (indexer(idx)-1)*self._numComponents
         self._recvMPIDataType[recvId] = Mpi.createDataTypeFromRangeAndSubRange(
	    recvRgn, localExtRange, self._numComponents, self._layout, elctCommType)
      end

      -- create MPI DataTypes for periodic directions
      -- also store location in memory required for sending/receiving periodic data
      self._sendLowerPerMPIDataType, self._recvLowerPerMPIDataType = {}, {}
      self._sendUpperPerMPIDataType, self._recvUpperPerMPIDataType = {}, {}
      self._sendLowerPerMPILoc, self._recvLowerPerMPILoc = {}, {}
      self._sendUpperPerMPILoc, self._recvUpperPerMPILoc = {}, {}

      -- Following loop creates Datatypes for periodic
      -- directions. This is complicated as one needs to treat lower
      -- -> upper transfers differently than upper -> lower as the
      -- number of ghost cells may be different on each lower/upper
      -- side. (AHH)
      for dir = 1, self._ndim do
	 if grid:isDirPeriodic(dir) then
	    local skelIds = decomposedRange:boundarySubDomainIds(dir)
	    for i = 1, #skelIds do
	       local loId, upId = skelIds[i].lower, skelIds[i].upper

	       -- only create if we are on proper ranks
	       if myId == loId then
		  local rgnSend = decomposedRange:subDomain(loId):lowerSkin(dir, self._upperGhost)
		  local idx = rgnSend:lowerAsVec()
		  -- set idx to starting point of region you want to recv
		  self._sendLowerPerMPILoc[dir] = (indexer(idx)-1)*self._numComponents
                  self._sendLowerPerMPIDataType[dir] = Mpi.createDataTypeFromRangeAndSubRange(
		     rgnSend, localExtRange, self._numComponents, self._layout, elctCommType)
		  
		  local rgnRecv = decomposedRange:subDomain(loId):lowerGhost(dir, self._lowerGhost)
		  local idx = rgnRecv:lowerAsVec()
		  -- set idx to starting point of region you want to recv
		  self._recvLowerPerMPILoc[dir] = (indexer(idx)-1)*self._numComponents
                  self._recvLowerPerMPIDataType[dir] = Mpi.createDataTypeFromRangeAndSubRange(
		     rgnRecv, localExtRange, self._numComponents, self._layout, elctCommType)
	       end
	       if myId == upId then
		  local rgnSend = decomposedRange:subDomain(upId):upperSkin(dir, self._lowerGhost)
		  local idx = rgnSend:lowerAsVec()
		  -- set idx to starting point of region you want to recv
		  self._sendUpperPerMPILoc[dir] = (indexer(idx)-1)*self._numComponents
                  self._sendUpperPerMPIDataType[dir] = Mpi.createDataTypeFromRangeAndSubRange(
		     rgnSend, localExtRange, self._numComponents, self._layout, elctCommType)

		  local rgnRecv = decomposedRange:subDomain(upId):upperGhost(dir, self._upperGhost)
		  local idx = rgnRecv:lowerAsVec()
		  -- set idx to starting point of region you want to recv
		  self._recvUpperPerMPILoc[dir] = (indexer(idx)-1)*self._numComponents
                  self._recvUpperPerMPIDataType[dir] = Mpi.createDataTypeFromRangeAndSubRange(
		     rgnRecv, localExtRange, self._numComponents, self._layout, elctCommType)
	       end	       
	    end
	 end
      end
      -- create IO object
      self._adiosIo = AdiosCartFieldIo {
	 elemType = elct,
	 metaData = tbl.metaData,
      }
      -- tag to identify basis used to set this field
      self._metaData = tbl.metaData
      self._basisId = "none"

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
	    local shape = self._localExtRangeDecomp:shape(self._shmIndex)
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
	 assert(field_compatible(self, fld), "CartField:combine: Can only accumulate compatible fields")
	 assert(type(fact) == "number", "CartField:combine: Factor not a number")

	 ffiC.gkylCartFieldAssign(self:_localLower(), self:_localShape(), fact, fld._data, self._data)
      end,
      _deviceAssign = function(self, fact, fld)
	 assert(field_compatible(self, fld), "CartField:combine: Can only accumulate compatible fields")
	 assert(type(fact) == "number", "CartField:combine: Factor not a number")

	 local numThreads = GKYL_DEFAULT_NUM_THREADS
	 local shape = self._localExtRangeDecomp:shape(self._shmIndex)
	 local numBlocks = math.floor(shape/numThreads)+1
	 ffiC.gkylCartFieldDeviceAssign(numBlocks, numThreads, self:_localLower(), self:_localShape(), fact, fld:deviceDataPointer(), self:deviceDataPointer())
      end,
      _accumulateOneFld = function(self, fact, fld)
	 assert(field_compatible(self, fld),
		"CartField:accumulate/combine: Can only accumulate/combine compatible fields")
	 assert(type(fact) == "number",
		"CartField:accumulate/combine: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:accumulate/combine: Fields should have same layout for sums to make sense")

	 ffiC.gkylCartFieldAccumulate(self:_localLower(), self:_localShape(), fact, fld._data, self._data)
      end,
      _deviceAccumulateOneFld = function(self, fact, fld)
	 assert(field_compatible(self, fld),
		"CartField:accumulate/combine: Can only accumulate/combine compatible fields")
	 assert(type(fact) == "number",
		"CartField:accumulate/combine: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:accumulate/combine: Fields should have same layout for sums to make sense")

	 local numThreads = GKYL_DEFAULT_NUM_THREADS
	 local shape = self._localExtRangeDecomp:shape(self._shmIndex)
	 local numBlocks = math.floor(shape/numThreads)+1
	 ffiC.gkylCartFieldDeviceAccumulate(numBlocks, numThreads, self:_localLower(), self:_localShape(), fact, fld:deviceDataPointer(), self:deviceDataPointer())
      end,
      accumulate = isNumberType and
	 function (self, c1, fld1, ...)
	    local args = {...} -- package up rest of args as table
	    local nFlds = #args/2
	    self:_accumulateOneFld(c1, fld1) -- accumulate first field
	    for i = 1, nFlds do -- accumulate rest of the fields
	       self:_accumulateOneFld(args[2*i-1], args[2*i])
	    end
	 end or
	 function (self, c1, fld1, ...)
	    assert(false, "CartField:accumulate: Accumulate only works on numeric fields")
	 end,
      deviceAccumulate = isNumberType and
	 function (self, c1, fld1, ...)
	    if self._devAllocData then
	       local args = {...} -- package up rest of args as table
	       local nFlds = #args/2
	       self:_deviceAccumulateOneFld(c1, fld1) -- accumulate first field
	       for i = 1, nFlds do -- accumulate rest of the fields
	          self:_deviceAccumulateOneFld(args[2*i-1], args[2*i])
	       end
	    end
	 end or
	 function (self, c1, fld1, ...)
	    assert(false, "CartField:accumulate: Accumulate only works on numeric fields")
	 end,
      combine = isNumberType and
         function (self, c1, fld1, ...)
            local args = {...} -- package up rest of args as table
            local nFlds = #args/2
            self:_assign(c1, fld1) -- assign first field
            for i = 1, nFlds do -- accumulate rest of the fields
               self:_accumulateOneFld(args[2*i-1], args[2*i])
            end
         end or
         function (self, c1, fld1, ...)
            assert(false, "CartField:combine: Combine only works on numeric fields")
         end,
      deviceCombine = isNumberType and
	 function (self, c1, fld1, ...)
	    if self._devAllocData then
	       local args = {...} -- package up rest of args as table
	       local nFlds = #args/2
	       self:_deviceAssign(c1, fld1) -- assign first field
	       for i = 1, nFlds do -- accumulate rest of the fields
	          self:_deviceAccumulateOneFld(args[2*i-1], args[2*i])
	       end
	    end
	 end or
	 function (self, c1, fld1, ...)
	    assert(false, "CartField:combine: Combine only works on numeric fields")
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
	       local shape = self._localExtRangeDecomp:shape(self._shmIndex)
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
	       local shape = self._localExtRangeDecomp:shape(self._shmIndex)
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
      write = function (self, fName, tmStamp, frNum, writeSkin)
	 self._adiosIo:write(self, fName, tmStamp, frNum, writeSkin)
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
      reduce = isNumberType and
	 function(self, opIn)
	    -- Input 'opIn' must be one of the binary operations in binOpFuncs.
	    local grid = self._grid
	    local tId = grid:subGridSharedId() -- Local thread ID.
	    local localRangeDecomp = LinearDecomp.LinearDecompRange {
	       range = self._localRange, numSplit = grid:numSharedProcs() }
	    local indexer = self:genIndexer()
	    local itr = self:get(1)
	    
	    local localVal = reduceInitialVal[opIn]
	    for idx in localRangeDecomp:rowMajorIter(tId) do
	       self:fill(indexer(idx), itr)
	       for k = 0, self._numComponents-1 do
		  localVal = binOpFuncs[opIn](localVal, itr:data()[k])
	       end
	    end

	    self.localReductionVal[1] = localVal
	    Mpi.Allreduce(self.localReductionVal:data(), self.globalReductionVal:data(),
			  1, elctCommType, reduceOpsMPI[opIn], grid:commSet().comm)
	    return self.globalReductionVal[1]
	 end or
	 function (self, opIn)
	    assert(false, "CartField:reduce: Reduce only works on numeric fields")
	 end,
      deviceReduce = isNumberType and
	 function(self, opIn, d_reduction)
	    -- Input 'opIn' must be one of the binary operations in binOpFuncs.
            ffi.C.gkylCartFieldDeviceReduce(binOpFlags[opIn],self._localRange:volume(),self.reduceBlocks,self.reduceThreads,self.reduceBlocksMAX,self.reduceThreadsMAX,
               self.deviceProps,self._onDevice,self.d_blockRed:data(),self.d_intermediateRed:data(),d_reduction:data())
	 end or
	 function (self, opIn)
	    assert(false, "CartField:deviceReduce: Reduce only works on numeric fields")
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
	 -- cells

	 local myId = self._grid:subGridId() -- grid ID on this processor
	 local neigIds = self._decompNeigh:neighborData(myId) -- list of neighbors
	 local tag = 42 -- Communicator tag for regular (non-periodic) messages
	 local recvReq = {} -- list of recv requests
	 -- post a non-blocking recv request
	 for _, recvId in ipairs(neigIds) do
            local dataType = self._recvMPIDataType[recvId]
            local loc = self._recvMPILoc[recvId]
	    -- recv data: (its from recvId-1 as MPI ranks are zero indexed)
	    recvReq[recvId] = Mpi.Irecv(dataPtr+loc, 1, dataType, recvId-1, tag, comm)
	 end
	 
	 -- do a blocking send (does not really block as recv requests
	 -- are already posted)
	 for _, sendId in ipairs(neigIds) do
            local dataType = self._sendMPIDataType[sendId]
            local loc = self._sendMPILoc[sendId]
	    -- send data: (its to sendId-1 as MPI ranks are zero indexed)
	    Mpi.Send(dataPtr+loc, 1, dataType, sendId-1, tag, comm)
	 end

	 -- complete recv
         -- since MPI DataTypes eliminate the need for buffers, 
         -- all we have to do is wait for non-blocking receives to finish
	 for _, recvId in ipairs(neigIds) do
	    Mpi.Wait(recvReq[recvId], nil)
	 end
      end,
      _field_periodic_sync = function (self, dataPtr)
	 local comm = self._grid:commSet().nodeComm -- communicator to use
	 if not Mpi.Is_comm_valid(comm) then
	    return -- no need to do anything if communicator is not valid
	 end
	 
	 -- immediately return if nothing to sync
	 if self._lowerGhost == 0 and self._upperGhost == 0 then return end

	 local grid = self._grid

	 -- Steps: (1) Post non-blocking recv requests. (2) Do
	 -- blocking sends, (3) Complete recv and copy data into ghost
	 -- cells

	 local decomposedRange = self._grid:decomposedRange()
	 local myId = self._grid:subGridId() -- grid ID on this processor
	 local basePerTag = 53 -- tag for periodic BCs

	 -- Note on tags: Each MPI message (send/recv pair) must have
	 -- a unique tag. With periodic BCs it is possible that a rank
	 -- may send/recv more than one message and hence some way to
	 -- distinguishing various messages is needed. The
	 -- non-periodic messages are all tagged 42. The periodic
	 -- messages have a base tag of 53, with up->lo tags being
	 -- 53+dir+10, while lo->up tags being 53+dir. As at most we
	 -- will have 6 directions, this generates enough unique tags
	 -- to do periodic communication safely.

	 local recvUpperReq, recvLowerReq  = {}, {}
	 -- post non-blocking recv requests for periodic directions
	 for dir = 1, self._ndim do
	    if grid:isDirPeriodic(dir) then
	       local skelIds = decomposedRange:boundarySubDomainIds(dir)
	       for i = 1, #skelIds do
		  local loId, upId = skelIds[i].lower, skelIds[i].upper

		  if myId == loId then
		     local loTag = basePerTag+dir+10
		     local dataType = self._recvLowerPerMPIDataType[dir]
		     local loc = self._recvLowerPerMPILoc[dir]
		     recvLowerReq[dir] = Mpi.Irecv(
			dataPtr+loc, 1, dataType, upId-1, loTag, comm)
		  end
		  if myId == upId then
		     local upTag = basePerTag+dir
		     local dataType = self._recvUpperPerMPIDataType[dir]
		     local loc = self._recvUpperPerMPILoc[dir]
		     recvUpperReq[dir] = Mpi.Irecv(
			dataPtr+loc, 1, dataType, loId-1, upTag, comm)
		  end
	       end
	    end
	 end
	 
	 -- do a blocking send for periodic directions (does not
	 -- really block as recv requests are already posted)
	 for dir = 1, self._ndim do
	    if grid:isDirPeriodic(dir) then
	       local skelIds = decomposedRange:boundarySubDomainIds(dir)
	       for i = 1, #skelIds do
		  local loId, upId = skelIds[i].lower, skelIds[i].upper

		  if myId == loId then
		     local loTag = basePerTag+dir -- this must match recv tag posted above
		     local dataType = self._sendLowerPerMPIDataType[dir]
                     local loc = self._sendLowerPerMPILoc[dir]
		     Mpi.Send(dataPtr+loc, 1, dataType, upId-1, loTag, comm)
		  end
		  if myId == upId then
		     local upTag = basePerTag+dir+10 -- this must match recv tag posted above
		     local dataType = self._sendUpperPerMPIDataType[dir]
		     local loc = self._sendUpperPerMPILoc[dir]
		     Mpi.Send(dataPtr+loc, 1, dataType, loId-1, upTag, comm)
		  end
	       end
	    end
	 end

	 -- complete recv for periodic directions
	 -- since MPI DataTypes eliminate the need for buffers,
         -- all we have to do is wait for non-blocking receives to finish
	 for dir = 1, self._ndim do
	    if grid:isDirPeriodic(dir) then
	       local skelIds = decomposedRange:boundarySubDomainIds(dir)
	       for i = 1, #skelIds do
		  local loId, upId = skelIds[i].lower, skelIds[i].upper
		  if myId == loId then
		     Mpi.Wait(recvLowerReq[dir], nil)
		  end
		  if myId == upId then
		     Mpi.Wait(recvUpperReq[dir], nil)
		  end
	       end
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
