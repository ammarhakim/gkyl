-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for ADIOS: header CDEFs
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local Mpi = require "Comm.Mpi"

-- Cut and paste from ADIOS header file, comments unchanged: Perhaps
-- this is not a good idea and this file could be auto-generated. Its
-- unlikely, though, that the ADIOS API will change.
ffi.cdef [[

/* global defines needed for the type creation/setup functions */
enum ADIOS_DATATYPES {adios_unknown = -1             /* (size) */

                     ,adios_byte = 0                 /* (1) */
                     ,adios_short = 1                /* (2) */
                     ,adios_integer = 2              /* (4) */
                     ,adios_long = 4                 /* (8) */

                     ,adios_unsigned_byte = 50       /* (1) */
                     ,adios_unsigned_short = 51      /* (2) */
                     ,adios_unsigned_integer = 52    /* (4) */
                     ,adios_unsigned_long = 54       /* (8) */

                     ,adios_real = 5                 /* (4) */
                     ,adios_double = 6               /* (8) */
                     ,adios_long_double = 7          /* (16) */

                     ,adios_string = 9               /* (?) */
                     ,adios_complex = 10             /* (8) */
                     ,adios_double_complex = 11      /* (16) */

                     /* Only for attributes: char** array of strings.
                        Number of elements must be known externally */
                     ,adios_string_array = 12        /* (sizeof(char*)) usually 4 */
                     };

enum ADIOS_FLAG {adios_flag_unknown = 0
                ,adios_flag_yes = 1
                ,adios_flag_no = 2
                };

enum ADIOS_BUFFER_ALLOC_WHEN {ADIOS_BUFFER_ALLOC_UNKNOWN
                             ,ADIOS_BUFFER_ALLOC_NOW
                             ,ADIOS_BUFFER_ALLOC_LATER
                             };
      
int adios_init (const char * config, MPI_Comm comm);
int adios_is_initialized(void);

int adios_finalize (int mype);

// end user calls for each I/O operation
// modes = "r" = "read", "w" = "write", "a" = "append", "u" = "update"
int adios_open (int64_t * fd, 
                const char * group_name, 
                const char * name,
                const char * mode, 
                MPI_Comm comm
               );

int adios_group_size (int64_t fd_p, 
                      uint64_t data_size,
                      uint64_t * total_size
                     ); 

int adios_write (int64_t fd_p, const char * name, const void * var);

int adios_get_write_buffer (int64_t fd_p, 
                            const char * name,
                            uint64_t * size,
                            void ** buffer
                           );

int adios_read (int64_t fd_p, 
                const char * name, 
                void * buffer,
                uint64_t buffer_size
               );

// OBSOLETE, kept only for backward compatibility of old codes
// Inconsistent behavior with new ADIOS variable naming
int adios_set_path (int64_t fd_p, const char * path);
// OBSOLETE, kept only for backward compatibility of old codes
// name should be the full path of the variable (as it was defined)
// in adios_write(), still the OLD full path should be used otherwise
// the variable is not found in the hash table
int adios_set_path_var (int64_t fd_p, const char * path, const char * name);

int adios_end_iteration (void);

int adios_start_calculation (void);

int adios_stop_calculation (void);

int adios_close (int64_t fd_p);

// ADIOS No-XML API's
int adios_init_noxml (MPI_Comm comm);

// To allocate ADIOS buffer OBSOLETE
int adios_allocate_buffer (
        enum ADIOS_BUFFER_ALLOC_WHEN adios_buffer_alloc_when,
        uint64_t buffer_size
        );

// To set maximum buffer size for each adios_open()...adios_close() operation.
void adios_set_max_buffer_size (uint64_t max_buffer_size_MB);

// To declare a ADIOS group
int adios_declare_group (int64_t * id, 
                         const char * name,
                         const char * time_index, 
                         enum ADIOS_FLAG stats
                        );

// To free a ADIOS group
int adios_free_group (int64_t id);

// To select a I/O method for a ADIOS group
int adios_select_method (int64_t group, 
                         const char * method,
                         const char * parameters,
                         const char * base_path
                        );

// To define a ADIOS variable
// Returns a variable ID, which can be used in adios_write_byid()
// 0 return value indicates an error
int64_t adios_define_var (int64_t group_id, 
                          const char * name,
                          const char * path,
                          enum ADIOS_DATATYPES type,
                          const char * dimensions,
                          const char * global_dimensions,
                          const char * local_offsets
                         );

// To remove all variable definitions from a group.
// Use it if you want to have a new set of variables defined
// for the next output step.
int adios_delete_vardefs (int64_t id);

// Return the expected size (in bytes) of a defined variable.
// It is simply the product of local dimensions and byte-size of type.
// It works only if the variable is defined with numeric dimensions or
// the dimension variables are already written with adios_write().
// Input is id returned by adios_define_var().
// Returns 0 on error, check adios_errno for error code
uint64_t adios_expected_var_size (int64_t var_id);

// To set the transform method for a variable just defined 
// var_id is the value returned by adios_define_var
// returns adios_errno (0=OK)
int adios_set_transform (int64_t var_id, const char *transform_type_str);

int adios_define_attribute (int64_t group, 
                            const char * name,
                            const char * path, 
                            enum ADIOS_DATATYPES type,
                            const char * value, 
                            const char * var
                           );

// define an attribute with values. 
// it can define an (1D) array of scalars, 'nelems' elements
// 'values' should point to the array of 'nelems' number of values of the specified type
int adios_define_attribute_byvalue (int64_t group, 
                            const char * name,
                            const char * path, 
                            enum ADIOS_DATATYPES type,
                            int  nelems,
                            const void * values
                           );

int adios_delete_attrdefs (int64_t id);

/** This function does similar function as adios_write. It is, however, used
 * in the following scenario that
 * 1. numbers, instead of a variable, are used to annotate array dimensions, and
 * 2. a variable is written mutiple times on a processor (e.g., AMR codes)
 */
int adios_write_byid (int64_t fd_p, int64_t id, const void * var);

/** Set the application's ID for adios_read_init()
 *  when using a staging method (DATASPACES, DIMES, NSSI or DATATAP).
 *  The ID should be unique for each application accessing the staging area
 *  IN:  id   a number unique for this application
 *  RETURN:       0 if accepted, <0 on error
 *  It is optional to use it before calling adios_init. Default is 1. 
 *  It has no effect for file based methods.
 *  Note: this function is defined both in adios.h and adios_read.h so that
 *  writing-only and reading-only applications can both use it.
 */ 
/*int adios_set_application_id (int id);*/

void adios_timing_write_xml (int64_t fd_p, const char* filename);

// no-xml schema API
// Define adios schema version
// The function implements the same as "schema version="1.1 ""in xml
int adios_define_schema_version (int64_t group_id, char * schema_version);

// Assign mesh to a variable
// The function implements the same as "var name="Var1" mesh="meshname" " in xml 
int adios_define_var_mesh(int64_t group_id, const char * varname, const char * meshname);

// Define centering of the variable value onto the mesh, centering is "cell" or "point"
int adios_define_var_centering(int64_t group_id, const char * varname, const char * centering);

// Define a external file where mesh variables are written 
int adios_define_mesh_file(int64_t group_id, char * name, char * file);

// The time-­steps points to time variables using steps, starting from step 0
int adios_define_var_timesteps (const char * timesteps, int64_t group_id, const char * name);

// The time-­steps points to time variables using real time, starting from time 
// Exactly like the time steps except with real numbers
int adios_define_var_timescale (const char * timescale, int64_t group_id, const char * name);

// Describe the padding pattern for output images
// If this number is 4, then the time-steps for images will be padded with 0 up to 4 digit numbers
int adios_define_var_timeseriesformat (const char * timeseries, int64_t group_id, const char * name);

// Use the concept of start, stride and count in all dimensions of a variable to identify a subset of a dataset
int adios_define_var_hyperslab (const char * hyperslab, int64_t group_id, const char * name);

// Describe if the mesh changes over time, and the option is "yes" or "no" 
int adios_define_mesh_timevarying (const char * timevarying, int64_t group_id, const char * name);

// The time-­steps points to time variables using steps, starting from step 0
int adios_define_mesh_timesteps (const char * timesteps, int64_t group_id, const char * name);

// The time-­steps points to time variables using real time, starting from time 0
int adios_define_mesh_timescale (const char * timescale, int64_t group_id, const char * name);

// Represent an integer for padding and formatting output image files
// If this number is 4, then the time-steps for images will be padded with 0 up to 4 digit number
int adios_define_mesh_timeseriesformat (const char * timeseries, int64_t group_id, const char * name);

// Indicates where (which ADIOS group) mesh variables are stored
int adios_define_mesh_group (const char * group, int64_t group_id, const char * name);

// Defines a uniform mesh
// For not requried attributes in this function, please use 0 instead
int adios_define_mesh_uniform (char * dimensions, 
                               char * origin, 
                               char * spacing, 
                               char * maximum, 
                               char * nspace,
                               int64_t group_id,
                               const char * name);

// Defines a rectilinear mesh
// For not requried attributes in this function, please use 0 instead
int adios_define_mesh_rectilinear (char * dimensions, 
                                   char * coordinates,
                                   char * nspace,
                                   int64_t group_id,
                                   const char * name);

// Defines a structured mesh
// For not requried attributes in this function, please use 0 instead
int adios_define_mesh_structured (char * dimensions, 
                                  char * points,
                                  char * nspace,
                                  int64_t group_id,
                                  const char * name);

// Define an unstructured mesh
// For not requried attributes in this function, please use 0 instead
int adios_define_mesh_unstructured (char * points,
                                    char * data, 
                                    char * count, 
                                    char * cell_type,
                                    char * npoints,
                                    char * nspace,
                                    int64_t group_id,
                                    const char * name);
]]
