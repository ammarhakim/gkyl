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
-- Multi-component on cartesian grids
--------------------------------------------------------------------------------

-- Declare meta-data for field. Actual data to field is not stored in
-- this structure.
 ffi.cdef [[
typedef struct {
    uint8_t _ndim; /* Dimension */
    uint16_t _lowerGhost, _upperGhost; /* Ghost cells in lower/upper faces */
    uint16_t _numComponents; /* Number of components */
    Range_t _localRange, _globalRange; /* Local/global range */
} FieldData_t;
]]

-- Field accessor object: allows access to field values in a given
-- cell.
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
      lowerGhost = function (self)
	 return self._m._lowerGhost
      end,
      upperGhost = function (self)
	 return self._m._upperGhost
      end,
      localRange = function (self)
	 return self._m._localRange
      end,
      localExtRange = function (self)
	 return self:localRange():extend(self:lowerGhost(), self:upperGhost())
      end,      
      globalRange = function (self)
	 return self._m._globalRange
      end,
      globalExtRange = function (self)
	 return self:globalRange():extend(self:lowerGhost(), self:upperGhost())
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
	 local sz = globalRange:extend(ghost[1], ghost[2]):volume()*nc -- amount of data in field
	 local f = new(self, nc)
	 f._m._ndim = grid:ndim()
	 f._m._lowerGhost, f._m._upperGhost = ghost[1], ghost[2]
	 f._m._globalRange = globalRange
	 f._m._localRange = globalRange -- ADJUST WHEN BRINGING IN PARALLEL
	 return f
      end,
      __index = field_mf,
   }
   return metatype(fct, field_mt)   
end

return {
   new_field_ct = new_field_ct,
   Field = new_field_ct(typeof("double")),
   FloatField = new_field_ct(typeof("float")),
}
