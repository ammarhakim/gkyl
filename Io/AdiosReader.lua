-- Gkyl ------------------------------------------------------------------------
--
-- Interface to ADIOS read API.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Mpi   = require "Comm.Mpi"
local Adios = require "Io.Adios"
local ffi   = require "ffi"
local Proto = require "Lib.Proto"
local Alloc = require "Lib.Alloc"
local Lin   = require "Lib.Linalg"

-- Table of type-strings and cast-strings for use in LuaJIT.
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

-- Wrapper function to convert enum to integer before indexing table.
local function getTypeStringMap(ty) return typeStringMap[tonumber(ty)] end

-- AdiosVar --------------------------------------------------------------------
--
-- Object to represent a variable in an ADIOS file (used internally
-- and can not be instantiated outside file).
--------------------------------------------------------------------------------

local AdiosVar = Proto()

function AdiosVar:init(fd, vName, vid)
   self._fd, self.name  = fd, vName
   self._obj            = Adios.inq_var_byid(fd, vid)

   self.shape, self.size = { }, 1
   for j = 0, self._obj.ndim-1 do
      self.shape[j+1] = tonumber(self._obj.dims[j])
      self.size       = self.size*tonumber(self._obj.dims[j])
   end

   local t   = getTypeStringMap(self._obj.type)
   self.type = t and t[1] or nil
end

-- Returned data is the value if var is a scalar or is the data
-- stored as a flat array otherwise.
function AdiosVar:read()
   local tmap = getTypeStringMap(self._obj.type)

   local ndim = self._obj.ndim
   if ndim == 0 then
      return ffi.cast(tmap[2], self._obj.value)[0]
   end
   
   local allocator = tmap[3]
   local data      = allocator(self.size)
   
   -- (ADIOS expects input to be const uint64_t* objects, hence vector
   -- types below).
   local start, count = Lin.UInt64Vec(ndim), Lin.UInt64Vec(ndim)
   for d = 1, ndim do
      start[d] = 0
      count[d] = self._obj.dims[d-1] -- dims is 0-indexed
   end
   local sel = Adios.selection_boundingbox(ndim, start, count)

   Adios.schedule_read_byid(self._fd, sel, self._obj.varid, 0, 1, data)
   Adios.perform_reads(self._fd, 1)

   return data
end

-- AdiosAttr -------------------------------------------------------------------
--
-- Object to represent an attribute in an ADIOS file (used internally
-- and can not be instantiated outside this file).
--------------------------------------------------------------------------------

local AdiosAttr = Proto()

function AdiosAttr:init(fd, aName, aid)
   self.name                 = aName
   local atype, asize, adata = Adios.get_attr_byid(fd, aid)
   local type_size           = Adios.type_size(atype, adata)
   self._values              = { }

   local tmap = getTypeStringMap(atype)

   if tmap[1] ~= "string" then
      local v = ffi.cast(tmap[2], adata) -- Cast to proper type.
      for i = 1, asize/type_size do
	 self._values[i] = v[i-1]
      end
   else
      -- For string attribute we need to cast it to Lua string object.
      self._values[1] = ffi.string(adata)
   end
   
   ffi.C.free(adata)
end

function AdiosAttr:read() return self._values end

-- Reader ----------------------------------------------------------------------
--
-- Reader for ADIOS BP file.
--------------------------------------------------------------------------------
local Reader = Proto()

function Reader:init(fName, comm)
   if not comm then comm = Mpi.COMM_WORLD end

   self._fd = Adios.read_open_file(fName, comm)
   assert(not (self._fd == nil), string.format(
	     "Unable to open ADIOS file %s for reading", fName))

   self.varList = { }
   for i = 0, self._fd.nvars-1 do
      local nm         = ffi.string(self._fd.var_namelist[i])
      self.varList[nm] = AdiosVar(self._fd, nm, i)
   end

   self.attrList =  { }
   for i = 0, self._fd.nattrs-1 do
      local nm          = ffi.string(self._fd.attr_namelist[i])
      self.attrList[nm] = AdiosAttr(self._fd, nm, i)
   end
end

function Reader:hasVar(nm)
   if self.varList[nm] then return true end
   return false
end
function Reader:getVar(nm)
   return self.varList[nm]
end

function Reader:hasAttr(nm)
   if self.attrList[nm] then return true end
   return false
end
function Reader:getAttr(nm)
   return self.attrList[nm]
end

function Reader:close() Adios.read_close(self._fd) end

return { Reader = Reader }
