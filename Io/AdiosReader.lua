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

-- method to read data and casting it to LuaJIT types
local function getValue(ty, data)
   local tyStr, value
   
   if ty == Adios.unsigned_byte then
      tyStr = "unsigned_byte"
      value = ffi.cast("uint8_t*", data)
   elseif ty == Adios.byte then
      tyStr = "byte"
      value = ffi.cast("int8_t*", data)
   elseif ty == Adios.short then
      tyStr = "short"
      value = ffi.cast("int16_t*", data)
   elseif ty == Adios.unsigned_short then
      tyStr = unsigned_short
      value = ffi.cast("uint16_t*", data)
   elseif ty == Adios.integer then
      tyStr = "integer"
      value = ffi.cast("int32_t*", data)
   elseif ty == Adios.unsigned_integer then
      tyStr = "unsigned_integer"
      value = ffi.cast("uint32_t*", data)      
   elseif ty == Adios.real then
      tyStr = "real"
      value = ffi.cast("float*", data)
   elseif ty == Adios.double then
      tyStr = "double"
      value = ffi.cast("double*", data)
   elseif ty == Adios.long_double then
      tyStr = "long_double"
      value = ffi.cast("long double*", data)
   elseif ty == Adios.string then
      tyStr = "string"
      value = ffi.string(data)
   else
      assert(false, "Unknown/NYI ADIOS type")
   end

   return tyStr, value
end

-- Object to represent a variable in an ADIOS file (used internally
-- and can not be instantiated by the user)
local AdiosVar = Proto()

function AdiosVar:init(vName, vObj)
   self.name = vName
   self._obj = vObj
   -- store shape of object
   self.shape = { }
   for j = 0, vObj.ndim-1 do
      self.shape[j+1] = vObj.dims[j]
   end
   self.type = nil
end

-- Returned data is the value if var is a scalar or is a pointer to
-- the data otherwise,
function AdiosVar:read() -->  type, data
   local tyStr, value = getValue(self._obj.type, self._obj.value)
   self.type = tyStr -- store this for further use
   -- return scalar value if var is a scalar
   if self._obj.ndim == 0 then
      return value[0]
   end
   return value
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
      local v = Adios.inq_var_byid(self.fd, i)
      self.varList[nm] = AdiosVar(nm, v)
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
