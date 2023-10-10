-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for ADIOS: header CDEFs
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local Mpi = require "Comm.Mpi"

-- Copied from adios2 header files. Not complete, only pasted the functions I think we'll need.

ffi.cdef [[

typedef struct adios2_adios adios2_adios;
typedef struct adios2_io adios2_io;
typedef struct adios2_variable adios2_variable;
typedef struct adios2_attribute adios2_attribute;
typedef struct adios2_engine adios2_engine;
typedef struct adios2_operator adios2_operator;

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

typedef enum
{
    adios2_type_unknown = -1,

    adios2_type_string = 0,
    adios2_type_float = 1,
    adios2_type_double = 2,
    adios2_type_float_complex = 3,
    adios2_type_double_complex = 4,

    adios2_type_int8_t = 5,
    adios2_type_int16_t = 6,
    adios2_type_int32_t = 7,
    adios2_type_int64_t = 8,

    adios2_type_uint8_t = 9,
    adios2_type_uint16_t = 10,
    adios2_type_uint32_t = 11,
    adios2_type_uint64_t = 12,
    adios2_type_long_double = 13 // junmin added
} adios2_type;

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
