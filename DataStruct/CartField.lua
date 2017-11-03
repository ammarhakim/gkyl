-- Gkyl ------------------------------------------------------------------------
--
-- Multi-component fields on cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries
local Adios = require "Io.Adios"
local Alloc = require "Lib.Alloc"
local AllocShared = require "Lib.AllocShared"
local CartDecompNeigh = require "Lib.CartDecompNeigh"
local LinearDecomp = require "Lib.LinearDecomp"
local Mpi = require "Comm.Mpi"
local Range = require "Lib.Range"

-- C interfaces
ffi.cdef [[
    // s: start index. sv: number of values to copy
    void gkylCartFieldAccumulate(unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldAssign(unsigned s, unsigned nv, double fact, const double *inp, double *out);
]]

-- Code from Lua wiki to convert table to comma-seperated-values
-- string.
-- Used to escape "'s by toCSV
local function escapeCSV (s)
  if string.find(s, '[,"]') then
    s = '"' .. string.gsub(s, '"', '""') .. '"'
  end
  return s
end
-- Convert from table to CSV string
local function toCSV (tt)
  local s = ""
  -- ChM 23.02.2014: changed pairs to ipairs assumption is that
  -- fromCSV and toCSV maintain data as ordered array
  for _,p in ipairs(tt) do  
    s = s .. "," .. escapeCSV(p)
  end
  return string.sub(s, 2)      -- remove first comma
end

-- Local definitions
local rowMajLayout, colMajLayout = 1, 2 -- data layout
local indexerMakerFuncs = {} -- list of functions that make indexers
indexerMakerFuncs[rowMajLayout] = Range.makeRowMajorIndexer
indexerMakerFuncs[colMajLayout] = Range.makeColMajorIndexer

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
   -- MPI and ADIOS data-types
   local elctCommType, elcCommSize = nil, 1
   if ffi.istype(new(elct), new("double")) then
      elctCommType = Mpi.DOUBLE
      elctIoType = Adios.double
      isNumberType = true
   elseif ffi.istype(new(elct), new("float")) then
      elctCommType = Mpi.FLOAT
      elctIoType = Adios.real
      isNumberType = true
   elseif ffi.istype(new(elct), new("int")) then
      elctCommType = Mpi.INT
      elctIoType = Adios.integer
      isNumberType = true
   else
      elctCommType = Mpi.BYTE -- by default, send stuff as byte array
      elcCommSize = sizeof(elct)
      elctIoType = Adios.byte
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
      local syncCorners = tbl.syncCorners and tbl.syncCorners or false -- Don't sync() corners by default
      local syncPeriodicDirs = tbl.syncPeriodicDirs and tbl.syncPeriodicDirs or true -- sync() periodic directions by default
      -- local and global ranges
      local globalRange = grid:globalRange()
      local localRange = grid:localRange()
      -- various communicators for use in shared allocator
      local comm = grid:commSet().comm     
      local shmComm = grid:commSet().sharedComm
      local nodeComm = grid:commSet().nodeComm

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
      self._localRange = localRange
      self._layout = colMajLayout -- default layout is column-major
      if tbl.layout then
	 if tbl.layout == "row-major" then
	    self._layout = rowMajLayout
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
	 local sendRgn = self._localRange:intersect(
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

      -- allocate space for IO (used to send data to ADIOS)
      self._outBuff = allocator(shmComm, localRange:volume()*self._numComponents)

      -- for use in ADIOS output
      local adLocalSz, adGlobalSz, adOffset = {}, {}, {}
      for d = 1, self._ndim do
	 adLocalSz[d] = localRange:shape(d)
	 adGlobalSz[d] = globalRange:shape(d)
	 adOffset[d] = localRange:lower(d)-1
      end
      adLocalSz[self._ndim+1] = self._numComponents
      adGlobalSz[self._ndim+1] = self._numComponents
      adOffset[self._ndim+1] = 0

      -- convert tables to comma-seperated-string. For some strange
      -- reasons, this is what ADIOS expects
      self._adLocalSz = toCSV(adLocalSz)
      self._adGlobalSz = toCSV(adGlobalSz)
      self._adOffset = toCSV(adOffset)
      
      return self
   end
   setmetatable(Field, { __call = function (self, o) return self.new(self, o) end })

   -- set callable methods
   Field.__index = {
      elemType = function(self)
	 return elct
      end,
      ndim = function (self)
	 return self._ndim
      end,
      numComponents = function (self)
	 return self._numComponents
      end,
      copy = function (self, fIn)
	 field_memcpy(self, fIn)
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
	 assert(field_compatible(self, fld), "CartField:accumulate/combine: Can only accumulate/combine compatible fields")
	 assert(type(fact) == "number", "CartField:accumulate/combine: Factor not a number")

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
	 return self:localRange():extend(self:lowerGhost(), self:upperGhost())
      end,      
      globalRange = function (self)
	 return self._globalRange
      end,
      globalExtRange = function (self) -- includes ghost cells
	 return self:globalRange():extend(self:lowerGhost(), self:upperGhost())
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
      sync = function (self)
	 return self._field_sync(self)
      end,
      write = function (self, outNm, tmStamp)
	 local comm = Mpi.getComm(self._grid:commSet().nodeComm)
	 -- (the extra getComm() is needed as Lua has no concept of
	 -- pointers and hence we don't know before hand if nodeComm
	 -- is a pointer or an object)
	 
	 if not Mpi.Is_comm_valid(comm) then
	    return -- no need to do anything if communicator is not valid
	 end

	 if not tmStamp then tmStamp = 0.0 end -- default time-stamp 
	 
	 local rank = Mpi.Comm_rank(comm)
	 -- setup ADIOS for IO
	 Adios.init_noxml(comm)
	 --Adios.set_max_buffer_size(16) -- 16 MB chunks	 

	 local ndim = self:ndim()
	 local localRange, globalRange = self:localRange(), self:globalRange()
	 -- setup group and set I/O method
	 local grpId = Adios.declare_group("CartField", "", Adios.flag_no)
	 Adios.select_method(grpId, "MPI", "", "")

	 -- field attributes
	 local cells = new("int[?]", ndim)
	 for d = 1, ndim do cells[d-1] = globalRange:shape(d) end
	 Adios.define_attribute_byvalue(grpId, "numCells", "", Adios.integer, ndim, cells)

	 local lower = new("double[?]", ndim)
	 for d = 1, ndim do lower[d-1] = self._grid:lower(d) end
	 Adios.define_attribute_byvalue(grpId, "lowerBounds", "", Adios.double, ndim, lower)

	 local upper = new("double[?]", ndim)
	 for d = 1, ndim do upper[d-1] = self._grid:upper(d) end
	 Adios.define_attribute_byvalue(grpId, "upperBounds", "", Adios.double, ndim, upper)

	 -- define data to write
	 Adios.define_var(
	    grpId, "time", "", Adios.double, "", "", "")
	 Adios.define_var(
	    grpId, "CartGridField", "", elctIoType, self._adLocalSz, self._adGlobalSz, self._adOffset)

	 -- copy field into output buffer (this copy is needed as
	 -- field also contains ghost-cell data, and, in addition,
	 -- ADIOS expects data to be laid out in row-major order)
	 self:_copy_from_field_region(self:localRange(), self._outBuff)

	 local fullNm = GKYL_OUT_PREFIX .. "_" .. outNm -- concatenate prefix
	 -- open file to write out group
	 local fd = Adios.open("CartField", fullNm, "w", comm)
	 local tmStampBuff = new("double[1]"); tmStampBuff[0] = tmStamp
	 Adios.write(fd, "time", tmStampBuff)
	 Adios.write(fd, "CartGridField", self._outBuff:data())
	 Adios.close(fd)
	 
	 Adios.finalize(rank)
      end,
      _copy_from_field_region = function (self, rgn, data)
	 local indexer = self:genIndexer()
	 local c = 1
	 for idx in rgn:rowMajorIter() do
	    local fitr = self:get(indexer(idx))
	    for k = 1, self._numComponents do
	       data[c] = fitr[k]; c = c+1
	    end
	 end
      end,
      _copy_to_field_region = function (self, rgn, data)
	 local indexer = self:genIndexer()
	 local c = 1
	 for idx in rgn:rowMajorIter() do
	    local fitr = self:get(indexer(idx))
	    for k = 1, self._numComponents do
	       fitr[k] = data[c]; c = c+1
	    end
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
	 local basePerTag = 50 -- Base tag for periodic messages
	 local localExtRange = self:localExtRange()

	 local recvReq = {} -- list of recv requests
	 -- post a non-blocking recv request
	 for _, recvId in ipairs(neigIds) do
	    local neighRgn = decomposedRange:subDomain(recvId)
	    local recvRgn = localExtRange:intersect(neighRgn)
	    local sz = recvRgn:volume()*self._numComponents
	    local buff = self._recvData[recvId]
	    -- recv data: (its from recvId-1 as MPI ranks are zero indexed)
	    recvReq[recvId] = Mpi.Irecv(buff:data(), sz*elcCommSize, elctCommType, recvId-1, tag, comm)
	 end
	 
	 -- do a blocking send (does not really block as the recv
	 -- requests are already posted)
	 for _, sendId in ipairs(neigIds) do
	    local neighRgn = decomposedRange:subDomain(sendId)
	    local sendRgn = self._localRange:intersect(
	       neighRgn:extend(self._lowerGhost, self._upperGhost))
	    local sz = sendRgn:volume()*self._numComponents

	    local buff = self._sendData[sendId]
	    self:_copy_from_field_region(sendRgn, buff)
	    -- send data: (its to sendId-1 as MPI ranks are zero indexed)
	    Mpi.Send(buff:data(), sz*elcCommSize, elctCommType, sendId-1, tag, comm)
	 end	 

	 -- complete recv
	 for _, recvId in ipairs(neigIds) do
	    local neighRgn = decomposedRange:subDomain(recvId)
	    local recvRgn = localExtRange:intersect(neighRgn)
	    local buff = self._recvData[recvId]
	    Mpi.Wait(recvReq[recvId], nil)
	    -- copy data into ghost cells
	    self:_copy_to_field_region(recvRgn, buff)	    
	 end
      end,
   }
   
   return Field
end

return {
   new_field_ct = Field_meta_ctor,
   Field = Field_meta_ctor(typeof("double")),
}
