-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for ADIOS.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Don't bother if ADIOS is not built in.
assert(GKYL_HAVE_ADIOS, "Gkyl was not built with ADIOS!")

local Lin  = require "Lib.Linalg"
local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Wrap bits of adios2 we think we need. There are more features not wrapped.
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

typedef enum
{
    adios2_constant_dims_false = 0,
    adios2_constant_dims_true = 1,
} adios2_constant_dims;

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

typedef enum
{
    adios2_mode_undefined = 0,
    adios2_mode_write = 1,
    adios2_mode_read = 2,
    adios2_mode_append = 3,
    adios2_mode_readRandomAccess = 6,

    adios2_mode_deferred = 4,
    adios2_mode_sync = 5
} adios2_mode;

typedef enum
{
    adios2_step_mode_append = 0,
    adios2_step_mode_update = 1,
    adios2_step_mode_read = 2,
} adios2_step_mode;

typedef enum
{
    adios2_step_status_other_error = -1,
    adios2_step_status_ok = 0,
    adios2_step_status_not_ready = 1,
    adios2_step_status_end_of_stream = 2
} adios2_step_status;

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

/**
 * @brief Begin a logical adios2 step stream
 * Check each engine documentation for MPI collective/non-collective
 * behavior.
 * @param engine handler
 * @param mode see enum adios2_step_mode in adios2_c_types.h for options,
 * read is the common use case
 * @param timeout_seconds provide a time out in Engine opened in read mode
 * @param status output from enum adios2_step_status in adios2_c_types.h
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_begin_step(adios2_engine *engine,
                               const adios2_step_mode mode,
                               const float timeout_seconds,
                               adios2_step_status *status);

//***************** PUT *****************
/**
 * Put data associated with a Variable in an engine, used for engines with
 * adios2_mode_write at adios2_open
 * @param engine handler for a particular engine where data will be put
 * @param variable contains variable metadata information
 * @param data user data to be associated with a variable, must be the same type
 * passed to adios2_define_variable
 * @param launch mode launch policy
 * <pre>
 *              adios2_mode_deferred: lazy evaulation, do not use data until
 * first adios2_perform_puts, adios2_end_step, or adios2_close. This is the
 * preferred way.
 *              adios_mode_sync, data is consumed by the engine and can be
 * reused immediately. Special case, only use if necessary.
 * </pre>
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_put(adios2_engine *engine, adios2_variable *variable,
                        const void *data, const adios2_mode launch);

/**
 * Put data associated with a Variable in an engine, used for engines with
 * adios2_mode_write at adios2_open. This is the name string version
 * @param engine handler for a particular engine where data will be put
 * @param variable_name variable with this name must exists in adios2_io that
 * opened the engine handler (1st parameter)
 * @param data user data to be associated with a variable, must be the same type
 * passed to adios2_define_variable
 * @param launch mode launch policy
 * <pre>
 *              adios2_mode_deferred: lazy evaulation, do not use data until
 * first adios2_perform_puts, adios2_end_step, or adios2_close. This is the
 * preferred way.
 *              adios_mode_sync, data is consumed by the engine and can be
 * reused immediately. Special case, only use if necessary.
 * </pre>
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_put_by_name(adios2_engine *engine,
                                const char *variable_name, const void *data,
                                const adios2_mode launch);

/**
 * Performs all the adios2_put and adios2_put_by_name called with mode
 * adios2_mode_deferred, up to this point, by copying user data into internal
 * ADIOS buffers. User data can be reused after this point.
 * @param engine handler for a particular engine where data will be put
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_perform_puts(adios2_engine *engine);

/**
 * Write array data to disk.  This may relieve memory pressure by clearing ADIOS
 * buffers.  It is a collective call. User data can be reused after this point.
 * @param engine handler for a particular engine where data will be put
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_perform_data_write(adios2_engine *engine);

//***************** GET *****************
/**
 * Gets data associated with a Variable from an engine, used for engines with
 * adios2_mode_read at adios2_open.
 * This is the name string version
 * @param engine handler for a particular engine where data will be put
 * @param variable handler must exists in adios2_io that
 * opened the engine handler (1st parameter). Typically from
 * adios2_inquire_variable
 * @param data user data to be associated with a variable, must be the same type
 * passed to adios2_define_variable. Must be pre-allocated for the required
 * variable selection.
 * @param launch mode launch policy
 * <pre>
 *              adios2_mode_deferred: lazy evaluation, do not use data until
 * first adios2_perform_puts, adios2_end_step, or adios2_close. This is the
 * preferred way.
 *              adios_mode_sync: data is populated by the engine and can be
 * reused immediately. Special case, only use if necessary.
 * </pre>
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_get(adios2_engine *engine, adios2_variable *variable,
                        void *data, const adios2_mode launch);

/**
 * Gets data associated with a Variable from an engine, used for engines with
 * adios2_mode_read at adios2_open.
 * This is the name string version
 * @param engine handler for a particular engine where data will be put
 * @param variable_name variable with this name must exists in adios2_io that
 * opened the engine handler (1st parameter).
 * @param data user data to be associated with a variable, must be the same type
 * passed to adios2_define_variable. Must be pre-allocated for the required
 * variable selection.
 * @param launch mode launch policy
 * <pre>
 *              adios2_mode_deferred: lazy evaluation, do not use data until
 * first adios2_perform_puts, adios2_end_step, or adios2_close. This is the
 * preferred way.
 *              adios_mode_sync, data is populated by the engine and can be
 * reused
 * immediately. Special case, only use if necessary.
 * </pre>
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_get_by_name(adios2_engine *engine,
                                const char *variable_name, void *data,
                                const adios2_mode launch);

/**
 * Performs all the adios2_get and adios2_get_by_name called with mode
 * adios2_mode_deferred up to this point by getting the data from the Engine.
 * User data can be reused after this point.
 * @param engine handler for a particular engine where data will be obtained
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_perform_gets(adios2_engine *engine);

/**
 * Terminates interaction with current step. By default puts/gets data to/from
 * all transports
 * Check each engine documentation for MPI collective/non-collective behavior.
 * @param engine handler executing IO tasks
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_end_step(adios2_engine *engine);

/**
 * Close all transports in adios2_Engine. Call is required to close system
 * resources.
 * MPI Collective, calls MPI_Comm_free for duplicated communicator at Open
 * @param engine handler containing all transports to
 * be closed. NOTE: engines NEVER become NULL after this function is called.
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_close(adios2_engine *engine);

/**
 * @brief Define a variable within io
 * @param io handler that owns the variable
 * @param name unique variable identifier
 * @param type primitive type from enum adios2_type in adios2_c_types.h
 * @param ndims number of dimensions
 * @param shape global dimension
 * @param start local offset
 * @param count local dimension
 * @param constant_dims adios2_constant_dims_true:: shape, start, count
 * won't change; adios2_constant_dims_false: shape, start, count will change
 * after definition
 * @return success: handler, failure: NULL
 */
adios2_variable *
adios2_define_variable(adios2_io *io, const char *name, const adios2_type type,
                       const size_t ndims, const size_t *shape,
                       const size_t *start, const size_t *count,
                       const adios2_constant_dims constant_dims);

/**
 * @brief Define an attribute value inside io
 * @param io handler that owns the attribute
 * @param name unique attribute name inside IO handler
 * @param type primitive type from enum adios2_type in adios2_c_types.h
 * @param value attribute single value
 * @return success: handler, failure: NULL
 */
adios2_attribute *adios2_define_attribute(adios2_io *io, const char *name,
                                          const adios2_type type,
                                          const void *value);

/**
 * @brief Define an attribute array inside io
 * @param io handler that owns the attribute
 * @param name unique attribute name inside IO handler
 * @param type primitive type from enum adios2_type in adios2_c_types.h
 * @param data attribute data array
 * @param size number of elements of data array
 * @return success: handler, failure: NULL
 */
adios2_attribute *adios2_define_attribute_array(adios2_io *io, const char *name,
                                                const adios2_type type,
                                                const void *data,
                                                const size_t size);
/**
 * Open an Engine to start heavy-weight input/output operations.
 * In MPI version reuses the communicator from adios2_init or adios2_init_config
 * MPI Collective function as it calls MPI_Comm_dup
 * @param io engine owner
 * @param name unique engine identifier
 * @param mode adios2_mode_write, adios2_mode_read, adios2_mode_append, and
 * adios2_mode_readRandomAccess
 * @return success: handler, failure: NULL
 */
adios2_engine *adios2_open(adios2_io *io, const char *name,
                           const adios2_mode mode);
]]
-- de-reference if object is a pointer
local function getObj(obj, ptyp)
   return ffi.istype(typeof(ptyp), obj) and obj[0] or obj
end

local _M = {}

--_M.adios = ffiC.get_adios2_adios()
--_M.io = ffiC.get_adios2_io()

local function new_adios2_adios() return new("adios2_adios[1]") end
local function new_adios2_io() return new("adios2_io[1]") end
local function new_adios2_variable() return new("adios2_variable[1]") end
local function new_adios2_engine() return new("adios2_engine[1]") end

-- Extract object from (potentially) a pointer
function _M.get_adios(obj) return getObj(obj, "adios2_adios[1]") end
function _M.get_io(obj) return getObj(obj, "adios2_io[1]") end
function _M.get_variable(obj) return getObj(obj, "adios2_variable[1]") end
function _M.get_engine(obj) return getObj(obj, "adios2_engine[1]") end

-- adios2_type
_M.type_unknown = -1
_M.type_string = 0
_M.type_float = 1
_M.type_double = 2
_M.type_float_complex = 3
_M.type_double_complex = 4
_M.type_int8_t = 5
_M.type_int16_t = 6
_M.type_int32_t = 7
_M.type_int64_t = 8
_M.type_uint8_t = 9
_M.type_uint16_t = 10
_M.type_uint32_t = 11
_M.type_uint64_t = 12
_M.type_long_double = 13

-- adios2_mode
_M.mode_undefined = 0
_M.mode_write = 1
_M.mode_read = 2
_M.mode_append = 3
_M.mode_readRandomAccess = 6
_M.mode_deferred = 4
_M.mode_sync = 5

-- adios2_init_mpi
function _M.init_mpi(comm)
   local ad_h = new_adios2_adios()
   ad_h = ffiC.adios2_init_mpi(getObj(comm, "MPI_Comm[1]"))
   return ad_h
end

-- adios2_declare_io
function _M.declare_io(ad_h, io_name)
   local io_h = new_adios2_io()
   io_h = ffiC.adios2_declare_io(_M.get_adios(ad_h), io_name)
   return io_h
end

-- adios2_finalize
function _M.finalize(ad_h)
   ad_error = ffiC.adios2_finalize(_M.get_adios(ad_h))
end

-- adios2_define_variable
function _M.define_variable(io_h, varName, varType, ndims, shape, start, count, are_dims_constant)
   local vecsIn = {shape = shape, start = start, count = count}
   local vecs   = {shape = type(shape)=='table' and Lin.UInt64Vec(ndims) or shape,
                   start = type(start)=='table' and Lin.UInt64Vec(ndims) or start,
                   count = type(count)=='table' and Lin.UInt64Vec(ndims) or count,}
   for k, v in pairs(vecsIn) do 
      if type(v) == 'table' then
         for d = 1, ndims do vecs[k][d] = v[d] end
      end
   end
   local constDims = are_dims_constant and 1 or 0

   local var = new_adios2_variable()
   var = ffiC.adios2_define_variable(_M.get_io(io_h), varName, varType, ndims, vecs.shape:data(),
                       vecs.start:data(), vecs.count:data(), constDims)
   return var
end

-- adios2_open
function _M.open(io_h, fileName, access_mode) 
   local engine_h = new_adios2_engine()
   engine_h = ffiC.adios2_open(_M.get_io(io_h), fileName, access_mode)
   return engine_h
end

-- adios2_put
function _M.put(engine_h, variable_h, data, launch_mode)
   local err = ffiC.adios2_put(_M.get_engine(engine_h), _M.get_variable(variable_h),
                               data, launch_mode)
   return err
end

-- adios2_close
function _M.close(engine_h)
   local err = ffiC.adios2_close(_M.get_engine(engine_h))
   return err
end

return _M
