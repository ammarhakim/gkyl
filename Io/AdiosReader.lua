-- Gkyl ------------------------------------------------------------------------
--
-- Interface to ADIOS read API.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Mpi   = require "Comm.Mpi"
local Adios = require "Io.Adios"
local Proto = require "Lib.Proto"
local Lin   = require "Lib.Linalg"

-- AdiosVar --------------------------------------------------------------------
--
-- Object to represent a variable in an ADIOS file (used internally
-- and can not be instantiated outside file).
--------------------------------------------------------------------------------

local AdiosVar = Proto()

function AdiosVar:init(ad_io, ad_engine, vName)
   self._ad_engine, self.name = ad_engine, vName
   self._ad_var = Adios.inquire_variable(ad_io, vName)

   self.ndim  = Adios.variable_ndims(self._ad_var)
   self.shape = Adios.variable_shape(self._ad_var)

   self.size = 1
   for j = 1, self.ndim do self.size = self.size*self.shape[j] end

   self.type = Adios.variable_type(self._ad_var)
end

-- Returned data is the value if var is a scalar or is the data
-- stored as a flat array otherwise.
function AdiosVar:read()
   local ndim = self.ndim
   local data
   if self.type == Adios.type_string then
      data = Lin.CharVec(Adios.name_char_num_max)
   elseif self.type == Adios.type_float then
      data = Lin.FloatVec(self.size)
   elseif self.type == Adios.type_double then
      data = Lin.Vec(self.size)
   elseif self.type == Adios.type_int32_t then
      data = Lin.IntVec(self.size)
   elseif self.type == Adios.type_uint64_t then
      data = Lin.UInt64Vec(self.size)
   end

   if ndim == 0 then
      local _ = Adios.get(self._ad_engine, self._ad_var, data:data(), Adios.mode_sync)
   else
      local start, count = {}, {}
      for d = 1, ndim do
         start[d] = 0
         count[d] = self.shape[d]
      end
      local _ = Adios.set_selection(self._ad_var, ndim, start, count)
      local _ = Adios.get(self._ad_engine, self._ad_var, data:data(), Adios.mode_sync)
   end

   if self.type == Adios.type_string then data = {ffi.string(data:data())} end
   if ndim == 0 then data = data[1] end
   return data
end

-- AdiosAttr -------------------------------------------------------------------
--
-- Object to represent an attribute in an ADIOS file (used internally
-- and can not be instantiated outside this file).
--------------------------------------------------------------------------------

local AdiosAttr = Proto()

function AdiosAttr:init(ad_io, aName)
   self.name    = aName
   ad_attr      = Adios.inquire_attribute(ad_io, aName)
   self._values = Adios.attribute_data(ad_attr)
end

function AdiosAttr:read() return self._values end

-- Reader ----------------------------------------------------------------------
--
-- Reader for ADIOS BP file.
--------------------------------------------------------------------------------
local Reader = Proto()

function Reader:init(fName, comm)
   if not comm then comm = Mpi.COMM_WORLD end

   self._ad = Adios.init_mpi(comm)
   self._ad_io = Adios.declare_io(self._ad, fName)

   self._ad_engine = Adios.open(self._ad_io, fName, Adios.mode_readRandomAccess)

   local var_names = Adios.available_variables(self._ad_io)
   self.varList = { }
   for i, nm in ipairs(var_names) do
      self.varList[nm] = AdiosVar(self._ad_io, self._ad_engine, nm)
   end

   local attr_names = Adios.available_attributes(self._ad_io)
   self.attrList =  { }
   for i, nm in ipairs(attr_names) do
      self.attrList[nm] = AdiosAttr(self._ad_io, nm)
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

function Reader:close()
   Adios.close(self._ad_engine)
   Adios.finalize(self._ad)
end

return { Reader = Reader }
