-- Gkyl ------------------------------------------------------------------------
--
-- Dynamically growing 1D vector
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
local Mpi = require "Comm.Mpi"
local Adios = require "Io.Adios"

-- Template to copy from table/vector
local copyTempl = xsys.template [[
return function (src, dest)
| for i = 1, NCOMP do
  dest[${i}] = src[${i}]
| end
end 	    
]]

-- DynVector -------------------------------------------------------------------
--
-- Dynamically growing 1D vector. Used to store diagnostics
--------------------------------------------------------------------------------

-- make constructor for DynVector
local DynVector = {}
function DynVector:new(tbl)
   local self = setmetatable({}, DynVector)
   
   self._numComponents = tbl.numComponents and tbl.numComponents or 1 -- default numComponents=1

   -- We store 1 extra element than requested to allow for 1-based indexing of returned values
   local allocator = Alloc.createAllocator(string.format("double[%d]", self._numComponents+1))
   self._timeMesh = Alloc.Double() -- time-mesh
   self._data = allocator() -- data
   self._tmpData = new(string.format("double[%d]", self._numComponents+1)) -- temp storage for single entry

   -- construct various functions from template representations
   self._copyToTempData = loadstring( copyTempl {NCOMP=self._numComponents} )()
   
   return self
end
setmetatable(DynVector, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
DynVector.__index = {
   numComponents = function(self)
      return self._numComponents
   end,
   appendData = function(self, t, v)
      self._timeMesh:push(t)
      self._copyToTempData(v, self._tmpData)
      self._data:push(self._tmpData)
   end,
   lastTime = function (self)
      return self._timeMesh:last()
   end,
   lastData = function(self)
      return self._data:last()
   end,
   timeMesh = function(self)
      return self._timeMesh
   end,
   data = function(self)
      return self._data
   end,
}

return { DynVector = DynVector }
