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

-- CartField -------------------------------------------------------------------
--
-- Multi-component field on cartesian grids
--------------------------------------------------------------------------------

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

   -- make constructor for Field 
   local Field = {}
   function Field:new(tbl)
      local self = setmetatable({}, Field)

      -- read data from input table
      local grid = tbl.onGrid
      local nc = tbl.numComponents and tbl.numComponents or 1 -- default numComponents=1
      local ghost = tbl.ghost and tbl.ghost or {0, 0} -- No ghost cells by default

      -- setup object: FIX ONCE CODE IS IN PARALLEL. (Need to get
      -- ranges from grid's decomp object)
      local l, u = {}, {}
      for dir = 1, grid:ndim() do
	 l[dir], u[dir] = 1, grid:numCells(dir)
      end
      local globalRange = grid:globalRange()
      local localRange = grid:localRange()

      -- allocate memory: this is not managed by the LuaJIT GC
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
   }
   return Field
end


return {
   new_field_ct = Field_meta_ctor,
   Field = Field_meta_ctor(typeof("double")),
}
