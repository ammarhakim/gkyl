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
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"
local Grid = require "Grid"

-- CartField -------------------------------------------------------------------
--
-- Multi-component field on cartesian grids
--------------------------------------------------------------------------------

-- Declare meta-data for field. Actual data to field is not stored in
-- this structure.
ffi.cdef [[
typedef struct {
    uint8_t _ndim; /* Dimension */
    uint16_t _lowerGhost, _upperGhost; /* Ghost cells in lower/upper faces */
    uint16_t _numComponents; /* Number of components */
    uint64_t _size; /* Size of field */
    Range_t _localRange, _globalRange; /* Local/global range */
    uint8_t _layout; /* Layout: 1 row-major, 2: column-major */
} FieldData_t;
]]

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

local function new_field_ct(elct)
   local fct = typeof("struct { FieldData_t _m; $ _data[?]; }", elct)
   local fcompct = new_field_comp_ct(elct)
   -- field object methods
   local field_mf = {
      elemType = function(self)
	 return elct
      end,
      ndim = function (self)
	 return self._m._ndim
      end,
      numComponents = function (self)
	 return self._m._numComponents
      end,
      copy = function (self, fIn)
	 field_memcpy(self, fIn)
      end,
      layout = function (self)
	 if self._m._layout == rowMajLayout then
	    return "row-major"
	 end
	 return "col-major"
      end,
      lowerGhost = function (self)
	 return self._m._lowerGhost
      end,
      upperGhost = function (self)
	 return self._m._upperGhost
      end,
      localRange = function (self)
	 return self._m._localRange
      end,
      localExtRange = function (self) -- includes ghost cells
	 return self:localRange():extend(self:lowerGhost(), self:upperGhost())
      end,      
      globalRange = function (self)
	 return self._m._globalRange
      end,
      globalExtRange = function (self) -- includes ghost cells
	 return self:globalRange():extend(self:lowerGhost(), self:upperGhost())
      end,
      size = function (self)
	 return self._m._size
      end,
      indexer = function (self) -- linear indexer taking (i,j,...)
	 return indexerMakerFuncs[self._m._layout](self:localExtRange())
      end,
      genIndexer = function (self) -- linear indexer taking indices as a vector
	 return genIndexerMakerFuncs[self._m._layout](self:localExtRange())
      end,
      get = function (self, k) -- k is an integer returned by a linear indexer
	 local loc = (k-1)*self._m._numComponents -- (k-1) as k is 1-based index
	 return fcompct(self._m._numComponents, self._data+loc)
      end,
      fill = function (self, k, fc)
	 local loc = (k-1)*self._m._numComponents -- (k-1) as k is 1-based index	 
	 fc._cdata = self._data+loc
      end,
   }
   -- field object meta-type
   local field_mt = {
      __new = function (self, tbl)
	 local grid = tbl.onGrid
	 local nc = tbl.numComponents and tbl.numComponents or 1 -- default numComponents=1
	 local ghost = tbl.ghost and tbl.ghost or {0, 0} -- No ghost cells by default
	 local l, u = {}, {}
	 for dir = 1, grid:ndim() do
	    l[dir], u[dir] = 1, grid:numCells(dir)
	 end
	 local globalRange = Range.Range(l, u)
	 local localRange = Range.Range(l, u) -- ADJUST WHEN DOING PARALLEL
	 local sz = localRange:extend(ghost[1], ghost[2]):volume()*nc -- amount of data in field
	 local f = new(self, sz) -- REPLACE WITH MALLOC-ED DATA TO OVERCOME 2GB LIMIT
	 f._m._ndim = grid:ndim()
	 f._m._lowerGhost, f._m._upperGhost = ghost[1], ghost[2]
	 f._m._numComponents = nc
	 f._m._size = sz

	 f._m._globalRange = globalRange
	 f._m._localRange = localRange
	 f._m._layout = colMajLayout -- default layout is column-major
	 if tbl.layout then
	    if tbl.layout == "row-major" then
	       f._m._layout = rowMajLayout
	    end
	 end
	 return f
      end,
      __index = field_mf,
   }
   return metatype(fct, field_mt)   
end

return {
   new_field_ct = new_field_ct,
   Field = new_field_ct(typeof("double")),
}
