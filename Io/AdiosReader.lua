-- Gkyl ------------------------------------------------------------------------
--
-- Interface to ADIOS read API
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Mpi = require "Comm.Mpi"
local Adios = require "Io.Adios"
local ffi = require "ffi"
local Proto = require "Lib.Proto"
local Alloc = require "Lib.Alloc"
local Lin = require "Lib.Linalg"

-- table of type-strings and cast-strings for use in LuaJIT
local typeStringMap = { 
   [Adios.unsigned_byte] = {
      "unsigned_byte", "uint8_t*", Alloc.createAllocator("uint8_t")
   },
   [Adios.byte] = {
      "byte", "int8_t*", Alloc.createAllocator("int8_t")
   },
   [Adios.short] = {
      "short", "int16_t*"
   },
   [Adios.unsigned_short] = {
      "unsigned_short", "uint16_t*", Alloc.createAllocator("uint16_t")
   },
   [Adios.integer] = {
      "integer", "int32_t*", Alloc.createAllocator("int32_t")
   },
   [Adios.unsigned_integer] = {
      "unsigned_integer", "uint32_t*", Alloc.createAllocator("uint32_t")
   },
   [Adios.real] = {
      "real", "float*", Alloc.createAllocator("float")
   },
   [Adios.double] = {
      "double", "double*", Alloc.createAllocator("double")
   },
   [Adios.long_double] = {
      "long_double", "long double*", Alloc.createAllocator("long double")
   },
   [Adios.string] = {
      "string", "char*", Alloc.createAllocator("char *")
   },
}

-- wrapper function to convert enum to integer before indexing table
local function getTypeStringMap(ty) return typeStringMap[tonumber(ty)] end

-- method to read data and casting it to LuaJIT types
local function getValue(ty, data)
   local tmap = getTypeStringMap(ty)
   local v = ffi.cast(tmap[2], data)
   return v[0]
end

-- Object to represent a variable in an ADIOS file (used internally
-- and can not be instantiated by the user)
local AdiosVar = Proto()

function AdiosVar:init(fd, vName, vObj)
   self._fd, self.name, self._obj = fd, vName, vObj
   -- store shape of object and compute its size
   self.shape, self.size = { }, 1
   for j = 0, vObj.ndim-1 do
      self.shape[j+1] = vObj.dims[j]
      self.size = self.size*vObj.dims[j]
   end
   local t = getTypeStringMap(self._obj.type)
   self.type = t and t[1] or nil
end

-- Returned data is the value if var is a scalar or is the data
-- stored as a flat array otherwise.
function AdiosVar:read()
   if self._obj.ndim == 0 then
      return getValue(self._obj.type, self._obj.value)
   end
   
   local ndim = self._obj.ndim
   -- need to read stuff from file: allocate memory first
   local allocator = getTypeStringMap(self._obj.type)[3]
   local data = allocator(tonumber(self.size))
   
   -- create a bounding box for region of grid to read
   local start, count = Lin.UIntVec64(ndim), Lin.UIntVec64(ndim)
   for d = 1, ndim do
      start[d] = 0
      count[d] = self._obj.dims[d-1] -- dims is 0-indexed
   end
   -- perform read
   local sel = Adios.selection_boundingbox(ndim, start, count)
   Adios.schedule_read_byid(self._fd, sel, self._obj.varid, 0, 1, data)
   Adios.perform_reads(self._fd, 1)

   return data
end

-- Reader ----------------------------------------------------------------------
--
-- Reader for ADIOS BP file
--------------------------------------------------------------------------------
local Reader = Proto()

function Reader:init(fName, comm)
   if not comm then comm = Mpi.COMM_WORLD end

   self.fName = fName
   self.fd = Adios.read_open_file(fName, comm)
   assert(not (self.fd == nil), string.format(
	     "Unable to open ADIOS file %s for reading", fName))

   self.varList = { }
   -- read list of variables in file
   for i = 0, self.fd.nvars-1 do
      local nm = ffi.string(self.fd.var_namelist[i])
      self.varList[nm] = AdiosVar(self.fd, nm, Adios.inq_var_byid(self.fd, i))
   end
end

function Reader:hasVar(nm)
   if self.varList[nm] then return true end
   return false
end

function Reader:getVar(nm)
   return self.varList[nm]
end

return { Reader = Reader }
