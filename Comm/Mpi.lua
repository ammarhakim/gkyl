-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrappers for MPI 
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

  int MPI_Comm_rank(MPI_Comm comm, int *rank);
  int MPI_Comm_size(MPI_Comm comm, int *size);  
]]

-- Predefined communicators and constants
_M.COMM_WORLD = ffi.C.get_MPI_COMM_WORLD()

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


