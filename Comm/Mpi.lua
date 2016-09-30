-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for MPI 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local _M = {}

ffi.cdef [[
  typedef struct MPI_Comm_type *MPI_Comm;
  MPI_Comm get_MPI_COMM_WORLD();

  typedef struct MPI_Datatype_type *MPI_Datatype;
  MPI_Datatype get_MPI_CHAR();
  MPI_Datatype get_MPI_BYTE();
  MPI_Datatype get_MPI_SHORT();
  MPI_Datatype get_MPI_INT();
  MPI_Datatype get_MPI_LONG();
  MPI_Datatype get_MPI_FLOAT();
  MPI_Datatype get_MPI_DOUBLE();
  MPI_Datatype get_MPI_UNSIGNED_CHAR();
  MPI_Datatype get_MPI_UNSIGNED_SHORT();
  MPI_Datatype get_MPI_UNSIGNED();
  MPI_Datatype get_MPI_UNSIGNED_LONG();
  MPI_Datatype get_MPI_LONG_DOUBLE();
  MPI_Datatype get_MPI_LONG_LONG_INT();
  MPI_Datatype get_MPI_FLOAT_INT();
  MPI_Datatype get_MPI_LONG_INT();
  MPI_Datatype get_MPI_DOUBLE_INT();
  MPI_Datatype get_MPI_SHORT_INT();
  MPI_Datatype get_MPI_2INT();
  MPI_Datatype get_MPI_LONG_DOUBLE_INT();
  MPI_Datatype get_MPI_PACKED();  

  int MPI_Comm_rank(MPI_Comm comm, int *rank);
  int MPI_Comm_size(MPI_Comm comm, int *size);  
]]

-- Predefined communicators
_M.COMM_WORLD = ffi.C.get_MPI_COMM_WORLD()

-- Datatypes
_M.MPI_CHAR = ffi.C.get_MPI_CHAR()
_M.MPI_BYTE = ffi.C.get_MPI_BYTE()
_M.MPI_SHORT = ffi.C.get_MPI_SHORT()
_M.MPI_INT = ffi.C.get_MPI_INT()
_M.MPI_LONG = ffi.C.get_MPI_LONG()
_M.MPI_FLOAT = ffi.C.get_MPI_FLOAT()
_M.MPI_DOUBLE = ffi.C.get_MPI_DOUBLE()
_M.MPI_UNSIGNED_CHAR = ffi.C.get_MPI_UNSIGNED_CHAR()
_M.MPI_UNSIGNED_SHORT = ffi.C.get_MPI_UNSIGNED_SHORT()
_M.MPI_UNSIGNED = ffi.C.get_MPI_UNSIGNED()
_M.MPI_UNSIGNED_LONG = ffi.C.get_MPI_UNSIGNED_LONG()
_M.MPI_LONG_DOUBLE = ffi.C.get_MPI_LONG_DOUBLE()
_M.MPI_LONG_LONG_INT = ffi.C.get_MPI_LONG_LONG_INT()
_M.MPI_FLOAT_INT = ffi.C.get_MPI_FLOAT_INT()
_M.MPI_LONG_INT = ffi.C.get_MPI_LONG_INT()
_M.MPI_DOUBLE_INT = ffi.C.get_MPI_DOUBLE_INT()
_M.MPI_SHORT_INT = ffi.C.get_MPI_SHORT_INT()
_M.MPI_2INT = ffi.C.get_MPI_2INT()
_M.MPI_LONG_DOUBLE_INT = ffi.C.get_MPI_LONG_DOUBLE_INT()
_M.MPI_PACKED = ffi.C.get_MPI_PACKED()

-- some types for use in MPI functions
local int_1 = typeof("int[1]")

function _M.Comm_rank(comm)
   local r = int_1()
   local err = ffi.C.MPI_Comm_rank(comm, r)
   return r[0]
end

function _M.Comm_size(comm)
   local r = int_1()
   local err = ffi.C.MPI_Comm_size(comm, r)
   return r[0]
end

return _M


