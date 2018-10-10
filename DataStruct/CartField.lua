-- Gkyl ------------------------------------------------------------------------
--
-- Multi-component fields on cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, sizeof, typeof, metatype")

-- Gkyl libraries
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Alloc = require "Lib.Alloc"
local AllocShared = require "Lib.AllocShared"
local CartDecompNeigh = require "Lib.CartDecompNeigh"
local LinearDecomp = require "Lib.LinearDecomp"
local Mpi = require "Comm.Mpi"
local Range = require "Lib.Range"

-- C interfaces
ffi.cdef [[
    // s: start index. nv: number of values to copy
    void gkylCartFieldAccumulate(unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldAssign(unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldScale(unsigned s, unsigned nv, double fact, double *out);
    void gkylCartFieldAbs(unsigned s, unsigned nv, double *out);
    void copyFromField(double *data, double *f, unsigned numComponents, unsigned c);
    void copyToField(double *f, double *data, unsigned numComponents, unsigned c);
]]

-- Local definitions
local rowMajLayout, colMajLayout = 1, 2 -- data layout
local indexerMakerFuncs = {} -- list of functions that make indexers
indexerMakerFuncs[rowMajLayout] = Range.makeRowMajorIndexer
indexerMakerFuncs[colMajLayout] = Range.makeColMajorIndexer
-- Default layout
local defaultLayout = colMajLayout

local genIndexerMakerFuncs = {} -- list of functions that make generic indexers
genIndexerMakerFuncs[rowMajLayout] = Range.makeRowMajorGenIndexer
genIndexerMakerFuncs[colMajLayout] = Range.makeColMajorGenIndexer

-- helper to check if two field are compatible
local function field_compatible(y, x)
   return y:localRange() == x:localRange() and y:numComponents() == x:numComponents()
end

-- copy field x into field y
local function field_memcpy(y, x)
   assert(field_compatible(y,x), "Can't copy incompatible fields")
   local sz = y:size()
   copy(y._data, x._data, sizeof(y:elemType())*sz)
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
   return metatype(typeof("struct { int32_t numComponents; $* _cdata; }", elct), field_comp_mt)
end

-- A function to create constructors for Field objects
local function Field_meta_ctor(elct)
   local fcompct = new_field_comp_ct(elct) -- Ctor for component data

   local isNumberType = false
   -- MPI data-types
   local elctCommType, elcCommSize = nil, 1
   if ffi.istype(new(elct), new("double")) then
      elctCommType = Mpi.DOUBLE
      isNumberType = true
   elseif ffi.istype(new(elct), new("float")) then
      elctCommType = Mpi.FLOAT
      isNumberType = true
   elseif ffi.istype(new(elct), new("int")) then
      elctCommType = Mpi.INT
      isNumberType = true
   else
      elctCommType = Mpi.BYTE -- by default, send stuff as byte array
      elcCommSize = sizeof(elct)
   end

   -- allocator for use in non-shared applications
   local function allocatorFunc(comm, numElem)
      local alloc = Alloc.Alloc_meta_ctor(elct)
      return alloc(numElem)
   end
   -- allocator for use in shared applications
   local function sharedAllocatorFunc(comm, numElem)
      local alloc = AllocShared.AllocShared_meta_ctor(elct)
      return alloc(comm, numElem)
   end

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
      
      -- setup object
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

      -- pre-allocate memory for send/recv calls when doing ghost-cell
      -- sync(). This prevents memory fragmentation in the C memory
      -- system as otherwise one would need to malloc/free every time
      -- sync() is called.
      self._sendData, self._recvData = {}, {}      
      local decomposedRange = self._grid:decomposedRange()
      local myId = self._grid:subGridId() -- grid ID on this processor
      local neigIds = self._decompNeigh:neighborData(myId) -- list of neighbors

      for _, sendId in ipairs(neigIds) do
	 local neighRgn = decomposedRange:subDomain(sendId)
	 local sendRgn = localRange:intersect(
	    neighRgn:extend(self._lowerGhost, self._upperGhost))
	 local sz = sendRgn:volume()*self._numComponents
	 self._sendData[sendId] = allocator(shmComm, sz)
      end

      local localExtRange = self:localExtRange()
      for _, recvId in ipairs(neigIds) do
	 local neighRgn = decomposedRange:subDomain(recvId)
	 local recvRgn = localExtRange:intersect(neighRgn)
	 local sz = recvRgn:volume()*self._numComponents
	 self._recvData[recvId] = allocator(shmComm, sz)
      end

      -- pre-allocate memory for send/recv calls for periodic BCs
      self._sendLowerPerData, self._recvLowerPerData = {}, {}
      self._sendUpperPerData, self._recvUpperPerData = {}, {}

      -- Following loop allocates memory for periodic directions. This
      -- is complicated as one needs to treat lower -> upper transfers
      -- differently than upper -> lower as the number of ghost cells
      -- may be different on each lower/upper side. (AHH)
      for dir = 1, self._ndim do
	 if grid:isDirPeriodic(dir) then
	    local skelIds = decomposedRange:boundarySubDomainIds(dir)
	    for i = 1, #skelIds do
	       local loId, upId = skelIds[i].lower, skelIds[i].upper

	       -- only allocate if we are on proper ranks
	       if myId == loId then
		  local rgnSend = decomposedRange:subDomain(loId):lowerSkin(dir, self._upperGhost)
		  local szSend = rgnSend:volume()*self._numComponents
		  self._sendLowerPerData[dir] = allocator(shmComm, szSend)
		  local rgnRecv = decomposedRange:subDomain(upId):upperSkin(dir, self._lowerGhost)
		  local szRecv = rgnRecv:volume()*self._numComponents
		  self._recvLowerPerData[dir] = allocator(shmComm, szRecv)
	       end
	       if myId == upId then
		  local rgnSend = decomposedRange:subDomain(upId):upperSkin(dir, self._lowerGhost)
		  local szSend = rgnSend:volume()*self._numComponents
		  self._sendUpperPerData[dir] = allocator(shmComm, szSend)
		  local rgnRecv = decomposedRange:subDomain(loId):lowerSkin(dir, self._upperGhost)
		  local szRecv = rgnRecv:volume()*self._numComponents
		  self._recvUpperPerData[dir] = allocator(shmComm, szRecv)
	       end	       
	    end
	 end
      end
      
      -- create IO object
      self._adiosIo = AdiosCartFieldIo { elemType = elct }
      -- tag to identify basis used to set this field
      self._basisId = "none"

      return self
   end
   setmetatable(Field, { __call = function (self, o) return self.new(self, o) end })

   -- set callable methods
   Field.__index = {
      elemType = function (self)
	 return elct
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
	 self:_assign(1.0, fIn) --field_memcpy(self, fIn)
      end,
      clear = function (self, val)
	 self._allocData:fill(val)
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

	 ffi.C.gkylCartFieldAssign(self:_localLower(), self:_localShape(), fact, fld._data, self._data)
      end,      
      _accumulateOneFld = function(self, fact, fld)
	 assert(field_compatible(self, fld),
		"CartField:accumulate/combine: Can only accumulate/combine compatible fields")
	 assert(type(fact) == "number",
		"CartField:accumulate/combine: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:accumulate/combine: Fields should have same layout for sums to make sense")

	 ffi.C.gkylCartFieldAccumulate(self:_localLower(), self:_localShape(), fact, fld._data, self._data)
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
      scale = isNumberType and
	 function (self, fact)
	    ffi.C.gkylCartFieldScale(self:_localLower(), self:_localShape(), fact, self._data)
	 end or
	 function (self, fact)
	    assert(false, "CartField:scale: Scale only works on numeric fields")
	 end,      
      abs = isNumberType and
         function (self)
            ffi.C.gkylCartFieldAbs(self:_localLower(), self:_localShape(), self._data)
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
      sync = function (self, syncPeriodicDirs)
         local syncPeriodicDirs = xsys.pickBool(syncPeriodicDirs, true)
	 self._field_sync(self)
	 if self._syncPeriodicDirs and syncPeriodicDirs then
	    self._field_periodic_sync(self)
	 end
	 -- this barrier is needed as when using MPI-SHM some
	 -- processors will not participate in sync()
	 Mpi.Barrier(self._grid:commSet().comm)
      end,
      setBasisId = function(self, basisId)
         self._basisId = basisId
      end,
      getBasisId = function(self)
         return self._basisId
      end,
      compatible = function(self, fld)
         return field_compatible(self, fld)
      end,
      _copy_from_field_region = function (self, rgn, data)
	 local indexer = self:genIndexer()
	 local c = 0
         local fitr = self:get(1)
	 for idx in rgn:rowMajorIter() do
	    self:fill(indexer(idx), fitr)
            ffi.C.copyFromField(data:data(), fitr:data(), self._numComponents, c)
            c = c + self._numComponents
	    --for k = 1, self._numComponents do
	    --   data[c] = fitr[k]; c = c+1
	    --end
	 end
      end,
      _copy_to_field_region = function (self, rgn, data)
	 local indexer = self:genIndexer()
	 local c = 0
         local fitr = self:get(1)
	 for idx in rgn:rowMajorIter() do
	    self:fill(indexer(idx), fitr)
            ffi.C.copyToField(fitr:data(), data:data(), self._numComponents, c)
            c = c + self._numComponents
	    --for k = 1, self._numComponents do
	    --   fitr[k] = data[c]; c = c+1
	    --end
	 end
      end,
      _field_sync = function (self)
	 local comm = self._grid:commSet().nodeComm -- communicator to use
	 if not Mpi.Is_comm_valid(comm) then
	    return -- no need to do anything if communicator is not valid
	 end
	 -- immediately return if nothing to sync
	 if self._lowerGhost == 0 and self._upperGhost == 0 then return end

	 -- Steps: (1) Post non-blocking recv requests. (2) Do
	 -- blocking sends, (3) Complete recv and copy data into ghost
	 -- cells

	 local decomposedRange = self._grid:decomposedRange()
	 local myId = self._grid:subGridId() -- grid ID on this processor
	 local neigIds = self._decompNeigh:neighborData(myId) -- list of neighbors
	 local tag = 42 -- Communicator tag for regular (non-periodic) messages
	 local localExtRange = self:localExtRange()

	 local recvReq = {} -- list of recv requests
	 -- post a non-blocking recv request
	 for _, recvId in ipairs(neigIds) do
	    local buff = self._recvData[recvId]
	    local sz = buff:size()
	    -- recv data: (its from recvId-1 as MPI ranks are zero indexed)
	    recvReq[recvId] = Mpi.Irecv(buff:data(), sz*elcCommSize, elctCommType, recvId-1, tag, comm)
	 end
	 
	 -- do a blocking send (does not really block as recv requests
	 -- are already posted)
	 for _, sendId in ipairs(neigIds) do
	    local neighRgn = decomposedRange:subDomain(sendId)
	    local sendRgn = self._localRange:intersect(
	       neighRgn:extend(self._lowerGhost, self._upperGhost))
	    local sz = sendRgn:volume()*self._numComponents

	    local buff = self._sendData[sendId]
	    self:_copy_from_field_region(sendRgn, buff) -- copy from skin cells
	    -- send data: (its to sendId-1 as MPI ranks are zero indexed)
	    Mpi.Send(buff:data(), sz*elcCommSize, elctCommType, sendId-1, tag, comm)
	 end

	 -- complete recv
	 for _, recvId in ipairs(neigIds) do
	    local neighRgn = decomposedRange:subDomain(recvId)
	    local recvRgn = localExtRange:intersect(neighRgn)
	    local buff = self._recvData[recvId]
	    Mpi.Wait(recvReq[recvId], nil)
	    self:_copy_to_field_region(recvRgn, buff) -- copy data into ghost cells
	 end
      end,
      _field_periodic_sync = function (self)
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
		     local loBuff = self._recvLowerPerData[dir]
		     local sz = loBuff:size()
		     --print(string.format("Recv request %d <- %d. Tag %d", myId, upId, loTag))
		     recvLowerReq[dir] = Mpi.Irecv(
			loBuff:data(), sz*elcCommSize, elctCommType, upId-1, loTag, comm)
		  end
		  if myId == upId then
		     local upTag = basePerTag+dir
		     local upBuff = self._recvUpperPerData[dir]
		     local sz = upBuff:size()
		     --print(string.format("Recv request %d <- %d. Tag %d", myId, loId, upTag))
		     recvUpperReq[dir] = Mpi.Irecv(
			upBuff:data(), sz*elcCommSize, elctCommType, loId-1, upTag, comm)
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
		     local loRgn = decomposedRange:subDomain(loId):lowerSkin(dir, self._upperGhost)
		     local loTag = basePerTag+dir -- this must match recv tag posted above
		     local loBuff = self._sendLowerPerData[dir]
		     self:_copy_from_field_region(loRgn, loBuff) -- copy from skin cells
		     local sz = loBuff:size()
		     --print(string.format("Sending %d -> %d. Tag %d", myId, upId, loTag))
		     Mpi.Send(loBuff:data(), sz*elcCommSize, elctCommType, upId-1, loTag, comm)
		  end
		  if myId == upId then
		     local upRgn = decomposedRange:subDomain(upId):upperSkin(dir, self._lowerGhost)
		     local upTag = basePerTag+dir+10 -- this must match recv tag posted above
		     local upBuff = self._sendUpperPerData[dir]
		     self:_copy_from_field_region(upRgn, upBuff) -- copy from skin cells
		     local sz = upBuff:size()
		     --print(string.format("Sending %d -> %d. Tag %d", myId, loId, upTag))
		     Mpi.Send(upBuff:data(), sz*elcCommSize, elctCommType, loId-1, upTag, comm)
		  end
	       end
	    end
	 end

	 -- complete recv for periodic directions
	 for dir = 1, self._ndim do
	    if grid:isDirPeriodic(dir) then
	       local skelIds = decomposedRange:boundarySubDomainIds(dir)
	       for i = 1, #skelIds do
	 	  local loId, upId = skelIds[i].lower, skelIds[i].upper

		  if myId == loId then
		     local loRgn = decomposedRange:subDomain(loId):lowerGhost(dir, self._lowerGhost)
		     local loBuff = self._recvLowerPerData[dir]
		     --print(string.format("Waiting for recv on %d", loId))
		     Mpi.Wait(recvLowerReq[dir], nil)
		     self:_copy_to_field_region(loRgn, loBuff) -- copy data into ghost cells
		  end
		  if myId == upId then
		     local upRgn = decomposedRange:subDomain(upId):upperGhost(dir, self._upperGhost)
		     local upBuff = self._recvUpperPerData[dir]
		     --print(string.format("Waiting for recv on %d", upId))
		     Mpi.Wait(recvUpperReq[dir], nil)
		     self:_copy_to_field_region(upRgn, upBuff) -- copy data into ghost cells
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
