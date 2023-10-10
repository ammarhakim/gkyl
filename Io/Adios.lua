-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for ADIOS.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Don't bother if ADIOS is not built in.
assert(GKYL_HAVE_ADIOS, "Gkyl was not built with ADIOS!")

local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

require "Io._AdiosCdef"     -- Load FFI binding.

-- de-reference if object is a pointer
local function getObj(obj, ptyp)
   return ffi.istype(typeof(ptyp), obj) and obj[0] or obj
end

local _M = {}

-- adios2_type
_M.unknown = -1

_M.string = 0
_M.float = 1
_M.double = 2
_M.float_complex = 3
_M.double_complex = 4

_M.int8_t = 5
_M.int16_t = 6
_M.int32_t = 7
_M.int64_t = 8

_M.uint8_t = 9
_M.uint16_t = 10
_M.uint32_t = 11
_M.uint64_t = 12
_M.long_double = 13

-- adios2_init_mpi
function _M.init_mpi(comm)
   ad_handler = ffi.new("adios2_adios")
   ad_handler = ffiC.adios2_init_mpi(getObj(comm, "MPI_Comm[1]"))
   return ad_handler
end

-- adios2_declare_io
function _M.declare_io(ad_handler, io_name)
   io_handler = ffi.new("adios2_io")
   io_handler = ffiC.adios2_declare_io(handler, io_name)
   return io_handler
end

function _M.finalize(ad_handler)
   ad_error = ffiC.adios2_finalize(ad_handler)
end

return _M
