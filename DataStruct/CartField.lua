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
local Alloc = require "Lib.Alloc"
local Range = require "Lib.Range"
local Grid = require "Grid"
local Mpi = require "Comm.Mpi"
local CartDecompNeigh = require "Lib.CartDecompNeigh"

-- Local definitions
local rowMajLayout, colMajLayout = 1, 2 -- data layout
local indexerMakerFuncs = {} -- list of functions that make indexers
indexerMakerFuncs[rowMajLayout] = Range.makeRowMajorIndexer
indexerMakerFuncs[colMajLayout] = Range.makeColMajorIndexer

local genIndexerMakerFuncs = {} -- list of functions that make generic indexers
genIndexerMakerFuncs[rowMajLayout] = Range.makeRowMajorGenIndexer
genIndexerMakerFuncs[colMajLayout] = Range.makeColMajorGenIndexer

-- copy field x into field y
local function field_memcpy(y, x)
   assert(y:localRange() == x:localRange() and y:numComponents() == x:numComponents(), "Can't copy incompatible fields")
   local sz = y:size()
   copy(y._data, x._data, sizeof(y:elemType())*sz)
end

-- Field accessor object: allows access to field values in cell
local function new_field_comp_ct(elct)
   local field_comp_mt = {
      __index = function(self, k)
	 return self._cdata[k-1]
      end,
      __newindex = function(self, k, v)
	 self._cdata[k-1] = v
      end,
   }
   return metatype(typeof("struct { int32_t numComponents; $* _cdata; }", elct), field_comp_mt)
end

-- A function to create constructors for Field objects
local function Field_meta_ctor(elct)
   local fcompct = new_field_comp_ct(elct)
   local allocator = Alloc.Alloc_meta_ctor(elct)
   local elctCommType, elcCommSize = nil, 1
   if ffi.istype(new(elct), new("double")) then
      elctCommType = Mpi.DOUBLE
   elseif ffi.istype(new(elct), new("float")) then
      elctCommType = Mpi.FLOAT
   elseif ffi.istype(new(elct), new("int")) then
      elctCommType = Mpi.INT
   else
      elctCommType = Mpi.BYTE -- by default, send stuff as byte array
      elcCommSize = sizeof(elct)
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

      -- local and global ranges
      local globalRange = grid:globalRange()
      local localRange = grid:localRange()

      -- allocate memory: this is NOT managed by the LuaJIT GC, allowing fields to be arbitrarly large
      local sz = localRange:extend(ghost[1], ghost[2]):volume()*nc -- amount of data in field
      self._allocData = allocator(sz) -- store this so it does not vanish under us
      self._data = self._allocData:data() -- pointer to data
      
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
      
      -- compute communication neighbors
      self._decompNeigh = CartDecompNeigh(grid:decomposedRange())
      if syncCorners then
	 self._decompNeigh:calcAllCommNeigh(ghost[1], ghost[2])
      else
	 self._decompNeigh:calcFaceCommNeigh(ghost[1], ghost[2])
      end
      
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
      fill = function (self, k, fc)
	 local loc = (k-1)*self._numComponents -- (k-1) as k is 1-based index	 
	 fc._cdata = self._data+loc
      end,
      sync = function (self)
	 return self._field_sync(self)
      end,
      _copy_from_field_region = function (self, rgn, data)
	 local indexer = self:genIndexer()	 
	 local c = 1
	 for idx in rgn:colMajorIter() do
	    local fitr = self:get(indexer(idx))
	    for k = 1, self._numComponents do
	       data[c] = fitr[k]; c = c+1
	    end
	 end
      end,
      _copy_to_field_region = function (self, rgn, data)
	 local indexer = self:genIndexer()
	 local c = 1
	 for idx in rgn:colMajorIter() do
	    local fitr = self:get(indexer(idx))
	    for k = 1, self._numComponents do
	       fitr[k] = data[c]; c = c+1
	    end
	 end
      end,
      _field_sync = function (self)
	 -- immediately return if nothing to sync
	 if self._lowerGhost == 0 and self._upperGhost == 0 then return end
	 
	 local comm = self._grid:comm() -- communicator to use
	 local decomposedRange = self._grid:decomposedRange()
	 local myId = self._grid:subGridId() -- grid ID on this processor
	 local neigIds = self._decompNeigh:neighborData(myId) -- list of neighbors
	 local tag = 11 -- Communicator tag for messages

	 -- send data of our skin cells to neighbor ghost cells
	 for _, sendId in ipairs(neigIds) do
	    local neighRgn = decomposedRange:subDomain(sendId)
	    local sendRgn = self._localRange:intersect(
	       neighRgn:extend(self._lowerGhost, self._upperGhost))
	    local sz = sendRgn:volume()*self._numComponents
	    local data = allocator(sz)
	    self:_copy_from_field_region(sendRgn, data)
	    -- send data: (its to sendId-1 as MPI ranks are zero indexed)
	    Mpi.Send(data:data(), sz*elcCommSize, elctCommType, sendId-1, tag, comm)
	 end

	 local localExtRange = self:localExtRange()
	 -- recv data from neighbor skin cells and copy into our ghost cells
	 for _, recvId in ipairs(neigIds) do
	    local neighRgn = decomposedRange:subDomain(recvId)
	    local recvRgn = localExtRange:intersect(neighRgn)
	    local sz = recvRgn:volume()*self._numComponents
	    local data = allocator(sz)
	    -- recv data: (its from recvId-1 as MPI ranks are zero indexed)
	    Mpi.Recv(data:data(), sz*elcCommSize, elctCommType, recvId-1, tag, comm, nil)
	    -- copy it into field
	    self:_copy_to_field_region(recvRgn, data)
	 end	 
      end,
   }
   return Field
end

return {
   new_field_ct = Field_meta_ctor,
   Field = Field_meta_ctor(typeof("double")),
}
