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

ffi.cdef [[

// Opaque types.
typedef struct adios2_adios *adios2_adios;
typedef struct adios2_io *adios2_io;
typedef struct adios2_variable *adios2_variable;
typedef struct adios2_attribute *adios2_attribute;
typedef struct adios2_engine *adios2_engine;
typedef struct adios2_operator *adios2_operator;

// Size of opaque objects.
int sizeof_adios2_adios();
int sizeof_adios2_io();
int sizeof_adios2_variable();
int sizeof_adios2_attribute();
int sizeof_adios2_engine();
int sizeof_adios2_operator();

adios2_adios get_adios2_adios();
adios2_io get_adios2_io();

/**
 * @brief adios2_error return types for all ADIOS2 C API functions
 * Based on the library C++ standardized exceptions
 * https://en.cppreference.com/w/cpp/error/exception
 * Each error will issue a more detailed description in the standard error
 * output, stderr
 */
typedef enum
{
    /** success */
    adios2_error_none = 0,

    /**
     * user input error
     */
    adios2_error_invalid_argument = 1,

    /** low-level system error, e.g. system IO error */
    adios2_error_system_error = 2,

    /** runtime errors other than system errors, e.g. memory overflow */
    adios2_error_runtime_error = 3,

    /** any other error exception */
    adios2_error_exception = 4

} adios2_error;

/**
 * Starting point for MPI apps. Creates an ADIOS handler.
 * MPI collective and it calls MPI_Comm_dup
 * @param comm defines domain scope from application
 * @return success: handler, failure: NULL
 */
adios2_adios *adios2_init_mpi(MPI_Comm comm);

/**
 * Declares a new io handler
 * @param adios owner the io handler
 * @param name unique io identifier within current adios handler
 * @return success: handler, failure: NULL
 */
adios2_io *adios2_declare_io(adios2_adios *adios, const char *name);

/**
 * Final point for adios handler. Deallocates adios pointer. Required to avoid
 * memory leaks.
 * MPI collective and it calls MPI_Comm_free
 * @param adios handler to be deallocated, must be initialized with
 * adios2_init or adios2_init_config
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_finalize(adios2_adios *adios);

]]
-- de-reference if object is a pointer
local function getObj(obj, ptyp)
   return ffi.istype(typeof(ptyp), obj) and obj[0] or obj
end

local _M = {}

--_M.adios = ffiC.get_adios2_adios()
--_M.io = ffiC.get_adios2_io()

local function new_adios2_adios()
   return new("adios2_adios[1]")
end

local function new_adios2_io()
   return new("adios2_io[1]")
end

-- Extract object from (potentially) a pointer
function _M.get_adios(obj)
   return getObj(obj, "adios2_adios[1]")
end
function _M.get_io(obj)
   return getObj(obj, "adios2_io[1]")
end

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
   local ad_handler = new_adios2_adios()
   ad_handler = ffiC.adios2_init_mpi(getObj(comm, "MPI_Comm[1]"))
   return ad_handler
end

-- adios2_declare_io
function _M.declare_io(ad_handler, io_name)
   local io_handler = new_adios2_io()
   io_handler = ffiC.adios2_declare_io(_M.get_io(ad_handler), io_name)
   return io_handler
end

function _M.finalize(ad_handler)
   ad_error = ffiC.adios2_finalize(_M.get_adios(ad_handler))
end

return _M
