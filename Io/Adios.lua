-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for ADIOS.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Don't bother if ADIOS is not built in.
assert(GKYL_HAVE_ADIOS, "Gkyl was not built with ADIOS!")

local Mpi  = require "Comm.Mpi"
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
    adios2_false = 0,
    adios2_true = 1,
} adios2_bool;

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

typedef enum
{
    adios2_shapeid_unknown = -1,
    adios2_shapeid_global_value = 0,
    adios2_shapeid_global_array = 1,
    adios2_shapeid_joined_array = 2,
    adios2_shapeid_local_value = 3,
    adios2_shapeid_local_array = 4
} adios2_shapeid;

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
 * Retrieves a previously declared io handler by name
 * @param adios owner the io handler
 * @param name unique name for the previously declared io handler
 * @return success: handler, failure: NULL
 */
adios2_io *adios2_at_io(adios2_adios *adios, const char *name);

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
 * @brief Retrieve a variable handler within current io handler
 * @param io handler to variable io owner
 * @param name unique variable identifier within io handler
 * @return found: handler, not found: NULL
 */
adios2_variable *adios2_inquire_variable(adios2_io *io, const char *name);

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
 * Returns a handler to a previously defined attribute by name
 * @param io handler to attribute io owner
 * @param name unique attribute identifier within io handler
 * @return found: handler, not found: NULL
 */
adios2_attribute *adios2_inquire_attribute(adios2_io *io, const char *name);

/**
 * Returns an array of attribute handlers for all attribute present in the io
 * group
 * @param attributes output array of attribute pointers (pointer to an
 * adios2_attribute**)
 * @param size output number of attributes
 * @param io handler to attributes io owner
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_inquire_all_attributes(adios2_attribute ***attributes,
                                           size_t *size, adios2_io *io);

/**
 * @brief returns an array or c strings for names of available variables
 * Might create dangling pointers
 * @param io handler variables io owner
 * @param length of array of strings
 * @return names of variables as an array of strings
 */
char **adios2_available_variables(adios2_io *io, size_t *size);

/**
 * @brief returns an array or c strings for names of available attributes
 * Might create dangling pointers
 * @param io handler variables io owner
 * @param length of array of strings
 * @return names of variables as an array of strings
 */
char **adios2_available_attributes(adios2_io *io, size_t *size);

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

/**
 * Open an Engine to start heavy-weight input/output operations.
 * MPI Collective function as it calls MPI_Comm_dup
 * @param io engine owner
 * @param name unique engine identifier
 * @param mode adios2_mode_write, adios2_mode_read, adios2_mode_append, and
 * adios2_mode_readRandomAccess
 * @param comm communicator other than adios' handler comm. MPI only.
 * @return success: handler, failure: NULL
 */
adios2_engine *adios2_open_new_comm(adios2_io *io, const char *name,
                                    const adios2_mode mode, MPI_Comm comm);

/**
 * Set new start and count dimensions
 * @param variable handler for which new selection will be applied to
 * @param ndims number of dimensions for start and count
 * @param start new start dimensions array
 * @param count new count dimensions array
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_selection(adios2_variable *variable, const size_t ndims,
                                  const size_t *start, const size_t *count);

/**
 * Retrieve attribute name
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param name output, string without trailing '\0', NULL or preallocated buffer
 * @param size output, name size without '\0'
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_name(char *name, size_t *size,
                                   const adios2_attribute *attribute);

/**
 * Retrieve attribute type
 * @param type
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_type(adios2_type *type,
                                   const adios2_attribute *attribute);

/**
 * Retrieve attribute type in string form "char", "unsigned long", etc.
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param type output, string without trailing '\0', NULL or preallocated buffer
 * @param size output, type size without '\0'
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_type_string(char *type, size_t *size,
                                          const adios2_attribute *attribute);

/**
 * Checks if attribute is a single value or an array
 * @param result output, adios2_true: single value, adios2_false: array
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_is_value(adios2_bool *result,
                                       const adios2_attribute *attribute);

/**
 * Returns the number of elements (as in C++ STL size() function) if attribute
 * is a 1D array. If single value returns 1
 * @param size output, number of elements in attribute
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_size(size_t *size,
                                   const adios2_attribute *attribute);

/**
 * Retrieve attribute data pointer (read-only)
 * @param data output attribute values, must be pre-allocated
 * @param size data size
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_data(void *data, size_t *size,
                                   const adios2_attribute *attribute);

/**
 * Retrieve variable name
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param name output, string without trailing '\0', NULL or preallocated buffer
 * @param size output, name size without '\0'
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_name(char *name, size_t *size,
                                  const adios2_variable *variable);

/**
 * Retrieve variable type
 * @param type output, from enum adios2_type
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_type(adios2_type *type,
                                  const adios2_variable *variable);

/**
 * Retrieve variable type in string form "char", "unsigned long", etc.
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param type output, string without trailing '\0', NULL or preallocated buffer
 * @param size output, type size without '\0'
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_type_string(char *type, size_t *size,
                                         const adios2_variable *variable);

/**
 * Retrieve variable shapeid
 * @param shapeid output, from enum adios2_shapeid
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_shapeid(adios2_shapeid *shapeid,
                                     const adios2_variable *variable);

/**
 * Retrieve current variable number of dimensions
 * @param ndims output
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_ndims(size_t *ndims,
                                   const adios2_variable *variable);

/**
 * Retrieve current variable shape
 * @param shape output, must be pre-allocated with ndims
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_shape(size_t *shape,
                                   const adios2_variable *variable);
]]
-- de-reference if object is a pointer
local function getObj(obj, ptyp)
   return ffi.istype(typeof(ptyp), obj) and obj[0] or obj
end

local _M = {}

-- Maximum number of characters in attribute/variable name,
-- and maximum number of attributes/variables.
_M.name_char_num_max = 20
_M.varattr_num_max   = 20
_M.string_num_max    = 100000

local function new_adios2_adios() return new("adios2_adios[1]") end
local function new_adios2_io() return new("adios2_io[1]") end
local function new_adios2_variable() return new("adios2_variable[1]") end
local function new_adios2_attribute() return new("adios2_attribute[1]") end
local function new_adios2_engine() return new("adios2_engine[1]") end
local function new_adios2_bool() return new("adios2_bool[1]") end
local function new_adios2_type() return new("adios2_type[1]") end
local function new_adios2_step_status() return new("adios2_step_status[1]") end
local function new_adios2_shapeid() return new("adios2_shapeid[1]") end

-- Extract object from (potentially) a pointer
function _M.get_adios(obj) return getObj(obj, "adios2_adios[1]") end
function _M.get_io(obj) return getObj(obj, "adios2_io[1]") end
function _M.get_variable(obj) return getObj(obj, "adios2_variable[1]") end
function _M.get_attribute(obj) return getObj(obj, "adios2_attribute[1]") end
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

-- adios2_step_mode
_M.step_mode_append = 0
_M.step_mode_update = 1
_M.step_mode_read = 2

-- adios2_step_status
_M.step_status_other_error = -1
_M.step_status_ok = 0
_M.step_status_not_ready = 1
_M.step_status_end_of_stream = 2

-- Create a new step_status.
function _M.new_step_status() return new_adios2_step_status() end

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

-- adios2_at_io
function _M.at_io(ad_h, io_name)
   local io_h = new_adios2_io()
   io_h = ffiC.adios2_at_io(_M.get_adios(ad_h), io_name)
   return io_h
end

-- adios2_finalize
function _M.finalize(ad_h)
   local err = ffiC.adios2_finalize(_M.get_adios(ad_h))
   assert(err == 0, string.format("Io.Adios: error in finalize. Error code: %d.",tonumber(err)))
end

-- adios2_define_variable
function _M.define_variable(io_h, varName, varType, ndims, shape, start, count, are_dims_constant)
   -- Two options:
   --   1. pass only the first 3 arguments, or pass ndims=0 and shape=start=count=nil, to define a scalar.
   --   2. pass all arguments to define an array. In this case shape/start/count can be tables,
   --      a number (for ndim=1) or a Lin.UInt64Vec.
   local vecsIn = {shape = shape, start = start, count = count}
   local vecs   = {
      shape = type(shape)=='table' and Lin.UInt64Vec(ndims) or
             ((type(shape)=='number' or type(shape)=='nil') and {data=function(self) return shape end} or shape),
      start = type(start)=='table' and Lin.UInt64Vec(ndims) or
             ((type(start)=='number' or type(start)=='nil') and {data=function(self) return start end} or start),
      count = type(count)=='table' and Lin.UInt64Vec(ndims) or
             ((type(count)=='number' or type(count)=='nil') and {data=function(self) return count end} or count),
   }
   for k, v in pairs(vecsIn) do 
      if type(v) == 'table' then
         for d = 1, ndims do vecs[k][d] = v[d] end
      end
   end
   local constDims = are_dims_constant and 1 or 0

   local var = new_adios2_variable()
   var = ffiC.adios2_define_variable(_M.get_io(io_h), varName, varType, ndims or 0, vecs.shape:data(),
                       vecs.start:data(), vecs.count:data(), constDims)
   return var
end


-- adios2_inquire_variable
-- MF 2023/10/10: doesn't seem to work.
function _M.inquire_variable(io_h, varName)
   local var = new_adios2_variable()
   var = ffiC.adios2_inquire_variable(_M.get_io(io_h), varName)
   return var
end

-- adios2_define_attribute
function _M.define_attribute(io_h, attrName, attrType, dataPointer)
   local attr = new_adios2_attribute()
   attr = ffiC.adios2_define_attribute(_M.get_io(io_h), attrName, attrType, dataPointer)
   return attr
end

-- adios2_define_attribute_array
function _M.define_attribute_array(io_h, attrName, attrType, data, size)
   local dataVec
   if type(data) == 'table' then
      if     attrType == _M.type_float   then dataVec = Lin.Vec(size)
      elseif attrType == _M.type_double  then dataVec = Lin.Vec(size)
      elseif attrType == _M.type_int32_t then dataVec = Lin.IntVec(size)
      end

      for i = 1, size do dataVec[i] = data[i] end
   elseif type(data) == 'number' then
      dataVec = {data = function(self) return data end}
   else
      dataVec = data
   end
   local attr = new_adios2_attribute()
   attr = ffiC.adios2_define_attribute_array(_M.get_io(io_h), attrName, attrType, dataVec:data(), size)
   return attr
end

-- adios2_open
function _M.open(io_h, fileName, access_mode) 
   local engine_h = new_adios2_engine()
   engine_h = ffiC.adios2_open(_M.get_io(io_h), fileName, access_mode)
   return engine_h
end

-- adios2_open_new_comm
function _M.open_new_comm(io_h, fileName, access_mode, comm) 
   local engine_h = new_adios2_engine()
   engine_h = ffiC.adios2_open_new_comm(_M.get_io(io_h), fileName, access_mode, getObj(comm, "MPI_Comm[1]"))
   return engine_h
end

-- adios2_set_selection
function _M.set_selection(variable_h, ndims, start, count)
   local vecsIn = {shape = shape, start = start, count = count}
   local vecs   = {
      start = type(start)=='table' and Lin.UInt64Vec(ndims) or
             (type(start)=='number' and {data=function(self) return start end} or start),
      count = type(count)=='table' and Lin.UInt64Vec(ndims) or
             (type(count)=='number' and {data=function(self) return count end} or count),
   }
   for k, v in pairs(vecsIn) do 
      if type(v) == 'table' then
         for d = 1, ndims do vecs[k][d] = v[d] end
      end
   end
   local err = ffiC.adios2_set_selection(variable_h, ndims, vecs.start:data(), vecs.count:data())
   assert(err == 0, string.format("Io.Adios: error in set_selection. Error code: %d.",tonumber(err)))
end

-- adios2_put
function _M.put(engine_h, variable_h, dataPointer, launch_mode)
   local err = ffiC.adios2_put(_M.get_engine(engine_h), _M.get_variable(variable_h),
                               dataPointer, launch_mode)
   assert(err == 0, string.format("Io.Adios: error in put. Error code: %d.",tonumber(err)))
end

-- adios2_get
function _M.get(engine_h, variable_h, dataPointer, launch_mode)
   local err = ffiC.adios2_get(_M.get_engine(engine_h), _M.get_variable(variable_h),
                               dataPointer, launch_mode)
   assert(err == 0, string.format("Io.Adios: error in get. Error code: %d.",tonumber(err)))
end

-- adios2_close
function _M.close(engine_h)
   local err = ffiC.adios2_close(_M.get_engine(engine_h))
   assert(err == 0, string.format("Io.Adios: error in close. Error code: %d.",tonumber(err)))
end

-- adios2_begin_step
function _M.begin_step(engine_h, step_mode, timeout, step_stat)
   local err = ffiC.adios2_begin_step(_M.get_engine(engine_h), step_mode, timeout, step_stat)
   assert(err == 0, string.format("Io.Adios: error in begin_step. Error code: %d.",tonumber(err)))
end

-- adios2_end_step
function _M.end_step(engine_h)
   local err = ffiC.adios2_end_step(_M.get_engine(engine_h))
   assert(err == 0, string.format("Io.Adios: error in end_step. Error code: %d.",tonumber(err)))
end

-- adios2_inquire_attribute
function _M.inquire_attribute(io_h, attrName)
   local attr = new_adios2_attribute()
   attr = ffiC.adios2_inquire_attribute(_M.get_io(io_h), attrName)
   return attr
end

-- MF 2023/10/14: commenting this out because I don't know of an easy way to allocate adios2_attribute ***).
---- adios2_inquire_all_attributes
----adios2_error adios2_inquire_all_attributes(adios2_attribute ***attributes,
----                                           size_t *size, adios2_io *io);
--function _M.inquire_all_attributes(io_h)
--   local attrs, size = new("adios2_attribute[1]"), Lin.UInt64Vec(1)
--   local err = ffiC.adios2_inquire_all_attributes(attrs, size:data(), _M.get_io(io_h))
--   assert(err == 0, string.format("Io.Adios: error in inquire_all_attributes. Error code: %d.",tonumber(err)))
--   return attrs, size[1]
--end

-- adios2_available_variables
-- char **adios2_available_variables(adios2_io *io, size_t *size);
function _M.available_variables(io_h)
   local size = Lin.UInt64Vec(1)
   local strs = new("char[?]", _M.name_char_num_max*_M.varattr_num_max)
   strs = ffiC.adios2_available_variables(_M.get_io(io_h), size:data())
   local strsOut = {}
   for i = 1, tonumber(size[1]) do strsOut[i] = ffi.string(strs[i-1]) end
   return strsOut
end

-- adios2_available_attributes
function _M.available_attributes(io_h)
   local size = Lin.UInt64Vec(1)
   local strs = new("char[?]", _M.name_char_num_max*_M.varattr_num_max)
   strs = ffiC.adios2_available_attributes(_M.get_io(io_h), size:data())
   local strsOut = {}
   for i = 1, tonumber(size[1]) do strsOut[i] = ffi.string(strs[i-1]) end
   return strsOut
end

-- adios2_attribute_name
function _M.attribute_name(attr_h)
   local size = Lin.UInt64Vec(1)
   local nameC = new("char [?]", _M.name_char_num_max)
   local err = ffiC.adios2_attribute_name(nameC, size:data(), _M.get_attribute(attr_h))
   assert(err == 0, string.format("Io.Adios: error in attribute_name. Error code: %d.",tonumber(err)))
   local nameOut = ffi.string(nameC)
   return nameOut
end

-- adios2_attribute_type
function _M.attribute_type(attr_h)
   local attrType = new_adios2_type()
   local err = ffiC.adios2_attribute_type(attrType, _M.get_attribute(attr_h))
   assert(err == 0, string.format("Io.Adios: error in attribute_type. Error code: %d.",tonumber(err)))
   return attrType[0]
end

-- adios2_attribute_type_string
function _M.attribute_type_string(attr_h)
   local size = Lin.UInt64Vec(1)
   -- Call it once with ia null pointer to get the size.
   local err = ffiC.adios2_attribute_type_string(nil, size:data(), _M.get_attribute(attr_h))
   -- Call it again to fill the new string with type of attribute.
   local typeStrC = Lin.CharVec(tonumber(size[1])+1)
   local err = ffiC.adios2_attribute_type_string(typeStrC:data(), size:data(), _M.get_attribute(attr_h))
   assert(err == 0, string.format("Io.Adios: error in attribute_type_string. Error code: %d.",tonumber(err)))
   return ffi.string(typeStrC:data())
end

-- adios2_attribute_is_value
function _M.attribute_is_value(attr_h)
   local is_value = new_adios2_bool()
   local err = ffiC.adios2_attribute_is_value(is_value, _M.get_attribute(attr_h))
   assert(err == 0, string.format("Io.Adios: error in attribute_is_value. Error code: %d.",tonumber(err)))
   return is_value[0]
end

-- adios2_attribute_size
function _M.attribute_size(attr_h)
   local sizeOut = Lin.UInt64Vec(1)
   local err = ffiC.adios2_attribute_size(sizeOut:data(), _M.get_attribute(attr_h))
   assert(err == 0, string.format("Io.Adios: error in attribute_size. Error code: %d.",tonumber(err)))
   return tonumber(sizeOut[1])
end

-- adios2_attribute_data
function _M.attribute_data(attr_h)
   local typ, sz = _M.attribute_type(attr_h), _M.attribute_size(attr_h)
   -- Call it once with ia null pointer to get the size.
   local dataOut, sizeOut = nil, Lin.UInt64Vec(1)
   if typ == _M.type_string then
      dataOut = Lin.CharVec(_M.string_num_max)
   elseif typ == _M.type_float then
      dataOut = Lin.FloatVec(sz)
   elseif typ == _M.type_double then
      dataOut = Lin.Vec(sz)
   elseif typ == _M.type_int32_t then
      dataOut = Lin.IntVec(sz)
   elseif typ == _M.type_uint64_t then
      dataOut = Lin.UInt64Vec(sz)
   end
   local err = ffiC.adios2_attribute_data(dataOut:data(), sizeOut:data(), _M.get_attribute(attr_h))
   assert(err == 0, string.format("Io.Adios: error in attribute_data. Error code: %d.",tonumber(err)))

   if typ == _M.type_string then dataOut = {ffi.string(dataOut:data())} end
   return dataOut, tonumber(sizeOut[1])
end

-- adios2_variable_name
function _M.variable_name(var_h)
   local size = Lin.UInt64Vec(1)
   local nameC = new("char [?]", _M.name_char_num_max)
   local err = ffiC.adios2_variable_name(nameC, size:data(), _M.get_variable(var_h))
   assert(err == 0, string.format("Io.Adios: error in variable_name. Error code: %d.",tonumber(err)))
   local nameOut = ffi.string(nameC)
   return nameOut
end

-- adios2_variable_type
function _M.variable_type(var_h)
   local varType = new_adios2_type()
   local err = ffiC.adios2_variable_type(varType, _M.get_variable(var_h))
   assert(err == 0, string.format("Io.Adios: error in variable_type. Error code: %d.",tonumber(err)))
   return varType[0]
end

-- adios2_variable_type_string
function _M.variable_type_string(var_h)
   local size = Lin.UInt64Vec(1)
   -- Call it once with ia null pointer to get the size.
   local err = ffiC.adios2_variable_type_string(nil, size:data(), _M.get_variable(var_h))
   -- Call it again to fill the new string with type of attribute.
   local typeStrC = Lin.CharVec(tonumber(size[1])+1)
   local err = ffiC.adios2_variable_type_string(typeStrC:data(), size:data(), _M.get_variable(var_h))
   assert(err == 0, string.format("Io.Adios: error in variable_type_string. Error code: %d.",tonumber(err)))
   return ffi.string(typeStrC:data())
end

-- adios2_variable_shapeid
function _M.variable_shapeid(var_h)
   local shapeid = new_adios2_shapeid()
   local err = ffiC.adios2_variable_shapeid(shapeid, _M.get_variable(var_h))
   assert(err == 0, string.format("Io.Adios: error in variable_shapeid. Error code: %d.",tonumber(err)))
   return shapeid[0]
end

-- adios2_variable_ndims
function _M.variable_ndims(var_h)
   local ndims = Lin.UInt64Vec(1)
   local err = ffiC.adios2_variable_ndims(ndims:data(), _M.get_variable(var_h))
   assert(err == 0, string.format("Io.Adios: error in variable_ndims. Error code: %d.",tonumber(err)))
   return tonumber(ndims[1])
end

-- adios2_variable_shape
function _M.variable_shape(var_h)
   local ndims = _M.variable_ndims(var_h)
   local shape = Lin.UInt64Vec(ndims)
   local err = ffiC.adios2_variable_shape(shape:data(), _M.get_variable(var_h))
   assert(err == 0, string.format("Io.Adios: error in variable_shape. Error code: %d.",tonumber(err)))
   local shapeOut = Lin.IntVec(ndims)
   for d = 1, ndims do shapeOut[d] = tonumber(shape[d]) end
   return shapeOut
end


return _M
