-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for ADIOS: header CDEFs
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

-- Cut and paste from ADIOS header file, comments unchanged: Perhaps
-- this is not a good idea and this file could be auto-generated. Its
-- unlikely, though, that the ADIOS API will change.
ffi.cdef [[

/* adios_read_v2_fwd */

struct _ADIOS_FILE;
typedef struct _ADIOS_FILE ADIOS_FILE;

struct _ADIOS_VARSTAT;
typedef struct _ADIOS_VARSTAT ADIOS_VARSTAT;

struct _ADIOS_VARBLOCK;
typedef struct _ADIOS_VARBLOCK ADIOS_VARBLOCK;

struct _ADIOS_VARMESH;
typedef struct _ADIOS_VARMESH ADIOS_VARMESH;

struct _ADIOS_VARINFO;
typedef struct _ADIOS_VARINFO ADIOS_VARINFO;

struct _ADIOS_VARCHUNK;
typedef struct _ADIOS_VARCHUNK ADIOS_VARCHUNK;

/* adios_selection.h */

/*************************/
/* Types used in the API */
/*************************/

typedef struct ADIOS_SELECTION_STRUCT  ADIOS_SELECTION;

/* Type of selection */
enum ADIOS_SELECTION_TYPE {
    ADIOS_SELECTION_BOUNDINGBOX  = 0,   /* Contiguous block of data defined by offsets and counts in each dimension */
    ADIOS_SELECTION_POINTS       = 1,   /* List of individual points                                                */
    ADIOS_SELECTION_WRITEBLOCK   = 2,   /* Selection of an individual block written by a writer process             */
    ADIOS_SELECTION_AUTO         = 3    /* Let the method decide what to return                                     */
};


/* A Bounding Box */
typedef struct { 
    int       ndim;
    uint64_t *start;
    uint64_t *count;
} ADIOS_SELECTION_BOUNDINGBOX_STRUCT;

/* A list of points.
 * It is a 1D array of indices, linearized for all dimension
 *     (e.g.  [i1,j1,k1,i2,j2,k2,...,in,jn,kn] for n points in a 3D space.
 * If a container selection is given, points are relative coordinates/offsets in the
 * container box/writeblock.
 * 1D offsets in N-dimensional container are allowed. Such selections are returned by
 * FASTBIT and ALACRITY query method. File reading is supported for such selections.
 * adios_selection_points_1DtoND() can be used to convert 1D to N-D points.
 */
typedef struct { 
    int       ndim;
    int       _free_points_on_delete;     // user provided points are not copied, won't free either
    uint64_t  npoints;
    uint64_t *points;
    ADIOS_SELECTION *container_selection; // a writeblock, a bounding box, or NULL
} ADIOS_SELECTION_POINTS_STRUCT;

/* A selected block produced by a writer
 * Identified with an index in current versions.
 */
typedef struct { 
    int index;

    /* NCSU ALACRITY-ADIOS:
     *     Adding timestep-relative vs. absolute writeblock selections, as
     *     well as sub-PG selection support. Both of these are currently only
     *     used by the transform layer
     */
    int is_absolute_index;   /* 0 if 'index' is relative to the current timestep, != 0
                                otherwise (i.e., absolute index) */
    int is_sub_pg_selection; /* Whether this writeblock selection contains sub-PG bounds.
                                The following fields only matter if this is != 0 */

    /* Reads the linear range of elements in [element_offset, element_offset + nelements) */
    uint64_t element_offset;
    uint64_t nelements;
} ADIOS_SELECTION_WRITEBLOCK_STRUCT;

/* Let the read method decide what to return to each reading client. 
 * Hints are method-dependent parameters to influence what and how to 
 * return (e.g. the ordering of returned chunks, decomposition among
 * read processes, etc.)
 */
typedef struct { 
    char     *hints;
} ADIOS_SELECTION_AUTO_STRUCT;

/** Selection for reading a subset of a variable. 
 *   A selection is an additive list of bounding boxes and point-sets 
 */
struct ADIOS_SELECTION_STRUCT  {
       enum ADIOS_SELECTION_TYPE    type; /* Type of selection */
       union {
            ADIOS_SELECTION_BOUNDINGBOX_STRUCT bb;
            ADIOS_SELECTION_POINTS_STRUCT points;
            ADIOS_SELECTION_WRITEBLOCK_STRUCT block;
            ADIOS_SELECTION_AUTO_STRUCT autosel;
       } u;
       //ADIOS_SELECTION             *next;
};

/** Boundingbox selection to read a subset of a non-scalar variable.
 *  IN:
 *       ndim      Number of dimensions
 *       start     array of offsets to start reading in each dimension
 *       count     number of data elements to read in each dimension
 *  RETURN:        A selection which can be used to read variables
 */
ADIOS_SELECTION * adios_selection_boundingbox (int ndim, const uint64_t *start, const uint64_t *count);


/** Selection for a selection of an enumeration of positions.
 *  IN:
 *       ndim      Number of dimensions
 *       npoints   Number of points of the selection
 *       points    1D array of indices, compacted for all dimension
 *                 (e.g.  [i1,j1,k1,i2,j2,k2,...,in,jn,kn] for
 *                 n points in a 3D space.
 */
ADIOS_SELECTION* adios_selection_points (int ndim, uint64_t npoints, const uint64_t *points);

/** Selection for a block of data coming from a certain producer.
 *
 *  IN:
 *       index      Identifier of the written block, starting from 0 for the first block
 *                  written by producer rank 0. If each writer outputs one block
 *                  then the index equals to the rank of the write process. 
 *                  With multi-var writing and multiple steps in a file, index should be
 *                  calculated by the reading application using outside information beyond
 *                  what is provided by the ADIOS Read API.
 */
ADIOS_SELECTION* adios_selection_writeblock (int index);

/** Let the method decide what data gets to what reader process.
 *  This selection enables each reading method to provide an 'optimal'
 *  data transfer from writers to readers. It depends on the method and the 
 *  circumstances, what this selection actually means.
 *
 *  E.g. in situ processing: readers on a compute node will receive all data 
 *       from the writers on the same compute node.
 *
 *  IN:
 *       hints    Method dependent parameters to influence what and how to 
 *                return (e.g. decomposition; ordering of returned chunks)
 */
ADIOS_SELECTION* adios_selection_auto (char * hints);

/** Make a strided hyperslab selection the same way as in HDF5. 
 *  IN:
 *       ndim      Number of dimentsions
 *       start     array of offsets to start reading in each dimension
 *       strides   striding steps, NULL=1 in all dimensions (= boundingbox) 
 *       count     number of data elements to read in each dimension
 *       blocks    block size at each stride, NULL=1 in all dimensions
 */
/*
   No support: Klasky, Podhorszki
   Support: 

ADIOS_SELECTION* adios_selection_hyperslab (uint64_t ndim, uint64_t *start, uint64_t *strides, 
                                    uint64_t *count, uint64_t *blocks);
*/


/** Delete a selection and free up memory used by the selection.
  * IN: selection
  * RESULT: None
  * The ADIOS_SELECTION object can be simply freed by free(), too. 
  */
void adios_selection_delete (ADIOS_SELECTION *selection);

/* Convert one selection of 1D point offsets in a bounding box,
 * returned by FASTBIT and ALACRITY query methods, to N-dimensional points.
 * This function works only if there is a bounding box in pointsinbox1D->u.points.container!
 * It allocates memory for the result selection, after which the original can be freed.
 * Return:
 * If global==0, the points will be relative to the container, if not, the points will be
 * global coordinates (container's starting offsets added to each point) and result's container
 * will be NULL.
 */
ADIOS_SELECTION * adios_selection_points_1DtoND (ADIOS_SELECTION * pointsinbox1D, int global);

/* adios_schema.h */

enum ADIOS_MESH_TYPE
{
     ADIOS_MESH_UNIFORM      = 1
    ,ADIOS_MESH_STRUCTURED   = 2
    ,ADIOS_MESH_RECTILINEAR  = 3
    ,ADIOS_MESH_UNSTRUCTURED = 4
};

typedef struct
{
    int num_dimensions;
    uint64_t * dimensions;
    double * origins;
    double * spacings;
    double * maximums;
} MESH_UNIFORM;

typedef struct
{
    int use_single_var;        /* 1 means coordinates-single-var,0 means coordinates-multi-var */
    int num_dimensions;
    uint64_t * dimensions;
    char ** coordinates;       /* name of the variable(s) containing the rectilinear spacing values */
} MESH_RECTILINEAR;

typedef struct
{
    int use_single_var;        /* 1 means points-single-var, 0 mean points-multi-var */
    int num_dimensions;
    uint64_t * dimensions;
    int nspaces;
    char ** points;            /* name of the variable(s) containing the point coordinates  */
} MESH_STRUCTURED;

/* ADIOS Schema: supported cell types */
enum ADIOS_CELL_TYPE
{
     ADIOS_CELL_PT         = 1
    ,ADIOS_CELL_LINE       = 2
    ,ADIOS_CELL_TRI        = 3
    ,ADIOS_CELL_QUAD       = 4
    ,ADIOS_CELL_HEX        = 5
    ,ADIOS_CELL_PRI        = 6
    ,ADIOS_CELL_TET        = 7
    ,ADIOS_CELL_PYR        = 8
};

typedef struct
{
    int nspaces;
    uint64_t npoints;
    int nvar_points;           /* how much vars for points-multi-var, 1 for points-single-var */
    char ** points;
    int ncsets;
    uint64_t * ccounts;
    char ** cdata;
    enum ADIOS_CELL_TYPE * ctypes;
} MESH_UNSTRUCTURED;


typedef struct {   /*type returned by adios_inq_mesh for read method */
    int id;
    char * name;
    char * file_name; /* 0 means mesh struct from the same file, otherwise mesh struct from externel file  */
    int time_varying;           /*0 means not time-varying, 1 means time-varying */
    enum ADIOS_MESH_TYPE type;
    union
    {
        MESH_UNIFORM * uniform;
        MESH_RECTILINEAR * rectilinear;
        MESH_STRUCTURED * structured;
        MESH_UNSTRUCTURED * unstructured;
    };
} ADIOS_MESH;

/* adios_link.h */

enum ADIOS_LINK_TYPE {
    LINK_VAR = 1,
    LINK_IMAGE = 2
    /* expand supported link types here */
};

typedef struct
{
    int id;
    char * name;
    int nrefs;
    enum ADIOS_LINK_TYPE * type;
    char ** ref_names;    /* the linked variable name referred from this var */
    char ** ref_files;    /* full path, 0 means link from the same file, otherwise link from externel file */
} ADIOS_LINK;

/* adios_read_v2.h */

/*************************/
/* Types used in the API */
/*************************/

struct _ADIOS_FILE {
        uint64_t fh;                /* File handler                                                   */
        int      nvars;             /* Number of variables in all groups (with full path)             */
        char     ** var_namelist;   /* Variable names in a char* array                                */
        int      nattrs;            /* Number of attributes in all groups                             */
        char     ** attr_namelist;  /* Attribute names in a char* array                               */
        int      nmeshes;           /* Number of meshes in all groups                                 */
        char     ** mesh_namelist;  /* Mesh names in a char* array                                    */
        int      nlinks;            /* Number of links in all groups                                  */
        char     ** link_namelist;  /* link names in a char* array                                    */

        /* Stream step information */
        int      current_step;      /* The current step in a stream. For a file, it is always 0.      */
        int      last_step;         /* The currently available latest step in the stream/file.        */
        int      is_streaming;      /* Non-zero if in streaming mode, zero if in non-streaming mode   */

        /* Information about file/stream */
        char     *path;             /* Full path file name (as passed at open)                        */
        int      endianness;        /* 0: little endian, 1: big endian                                */
                                    /*   the read API takes care of conversion automatically          */
        int      version;           /* Version of ADIOS-BP format                                     */
        uint64_t file_size;         /* Size of file in bytes not including subfiles                   */

        /* Internals */
        void     * internal_data;   /* Data for internal use                                          */
};

struct _ADIOS_VARSTAT {
        void     * min;            /* minimum value in an array variable, = value for a scalar        */
        void     * max;            /* maximum value of an array variable (over all steps)             */
        double   * avg;            /* average value of an array variable (over all steps)             */
        double   * std_dev;        /* standard deviation value of an array variable (over all steps)  */

        struct ADIOS_STAT_STEP     /* per step statistics (if requested and recorded at writing)      */
        {
            void     ** mins;      /* minimum per each step (array of 'nsteps' elements)              */
            void     ** maxs;      /* maximum per each step (array of 'nsteps' elements)              */
            double   ** avgs;      /* average per each step (array of 'nsteps' elements)              */
            double   ** std_devs;  /* standard deviation per each step (array of 'nsteps' elements)   */
        } *steps;

        struct ADIOS_STAT_BLOCK    /* per block statistics (if requested and recorded at writing)     */
        {
            void     ** mins;      /* minimum per each block (array of 'sum_nblocks' elements)        */
            void     ** maxs;      /* maximum per each block (array of 'sum_nblocks' elements)        */
            double   ** avgs;      /* average per each block (array of 'sum_nblocks' elements)        */
            double   ** std_devs;  /* std deviation per each block (array of 'sum_nblocks' elements)  */
        } *blocks;

        struct ADIOS_HIST           /* Histogram if recorded at writing */
        {
            uint32_t    num_breaks;
            double      max;
            double      min;
            double *    breaks;
            uint32_t ** frequencies;
            uint32_t *  gfrequencies;
        } *histogram;
};

struct _ADIOS_VARBLOCK {
    uint64_t * start;      /* offset start point in global array ('ndim' elements)                    */
    uint64_t * count;      /* local sizes in global array ('ndim' elements)                           */
    uint32_t process_id;   /* a kind of ID of the writing process (likely MPI rank)                   */
    uint32_t time_index;   /* a kind of timestep info of the writing process >= step of variable      */
};

enum var_centering
{
    point = 1,            /* unstructured mesh point centering */
    cell = 2              /* unstructured mesh cell centering */
};

struct _ADIOS_VARMESH {
    int meshid;
    enum var_centering centering;
};

struct _ADIOS_VARINFO {
        int        varid;           /* variable index (0..ADIOS_FILE.nvars-1)                         */
        enum ADIOS_DATATYPES type;  /* type of variable                                               */
        int        ndim;            /* number of dimensions, 0 for scalars                            */
        uint64_t * dims;            /* size of each dimension.
                                       If variable has no global view 'dims' report the size of the 
                                       local array written by process rank 0. 
                                    */
        int        nsteps;          /* Number of steps of the variable in file. 
                                       There is always at least one step.                             */
                                    /* In streams it always equals 1.                                 */
        void     * value;           /* value of a scalar variable, NULL for array.                    */
        int        global;          /* 1: global view (was defined by writer), 
                                       0: pieces written by writers without defining a global array   */
        int      * nblocks;         /* Number of blocks that comprise this variable in a step
                                       It is an array of 'nsteps' integers                            */
        int        sum_nblocks;     /* Number of all blocks of all steps, the sum of the nblocks array*/
        int        nattrs;          /* Number of attributes with the name <variable_fullpath>/<name>  */
        int      * attr_ids;        /* Attribute ids in an array, use fp->attr_namelist[<id>] for names */
        ADIOS_VARSTAT  *statistics; /* Statistics, retrieved in separate call: adios_inq_var_stat()   */
        ADIOS_VARBLOCK *blockinfo;  /* Spatial arrangement of written blocks, 
                                       retrieved in separate call: adios_inq_var_blockinfo()       
                                       It is an array of 'sum_nblocks' elements                       */
        ADIOS_VARMESH *meshinfo;    /* structure in this file,
                                       retrieved in separate call: adios_inq_var_meshinfo()          */ 
};

struct _ADIOS_VARCHUNK {
        int                   varid;    /* variable index (0..ADIOS_FILE.nvars-1)              */
        enum ADIOS_DATATYPES  type;     /* type of variable                                    */
        /* NCSU ALACRITY-ADIOS - Added timestep information into varchunks */
        int                   from_steps; /* the first timestep in the returned data             */
        int                   nsteps;     /* the number of timesteps in the returned data        */
        ADIOS_SELECTION     * sel;      /* sub-selection of requested selection                */
        void                * data;     /* pointer to data, at next adios_read_check() memory 
                                           will likely be overwritten                          */
};


/* The list of the available read methods */
enum ADIOS_READ_METHOD {
        ADIOS_READ_METHOD_BP            = 0,  /* Read from ADIOS BP file (written by POSIX, MPI etc methods) */
        ADIOS_READ_METHOD_BP_AGGREGATE  = 1,  /* Read from ADIOS BP file (written by POSIX, MPI_AMR etc methods)  */
        ADIOS_READ_METHOD_DATASPACES    = 3,  /* Read from memory written by DATASPACES method               */
        ADIOS_READ_METHOD_DIMES         = 4,  /* Read from memory written by DIMES method                    */
        ADIOS_READ_METHOD_FLEXPATH      = 5,  /* Read from memory written by FLEXPATH method                 */
        ADIOS_READ_METHOD_ICEE          = 6,  /* Read from memory written by ICEE method                 */
};

/** Functions that return a pointer to some data structures (e.g. adios_fopen),
 *  return NULL on error and set adios_errno to a non-zero value and writes an
 *  error string.
 *  You have direct access to that string so you can print it.
 *  Do not write anything into it.
 *  The last error message is always available; it is not cleared until another error is detected.
 */
extern int adios_errno;
const char *adios_errmsg();

/** Initialize a reading method before opening a file/stream with using 
 *  the method.
 *  IN:  method     read method to use
 *       comm       MPI communicator of all processes participating
 *                  in a file/stream operation
 *       parameters A series of name=value pairs separated by ;
 *                  E.g. "max_memory=200; app_id = 1"
 *                  List of parameters is documented for each method 
 *                  separately. 
 *                  
 *  RETURN:       0 if accepted, !=0 on error (adios_errno value)
 *  Initialization is required for the staging methods, where init/finalize 
 *  perform the connection/disconnection to the staging server once.
 *  The ADIOS_READ_METHOD_BP does not need to be initialized/finalized.
 */
int adios_read_init_method (enum ADIOS_READ_METHOD method, 
                            MPI_Comm comm, 
                            const char * parameters);
 
/*  max_memory 
 *                  maximum size of memory to be used by the read method,
 *                  in MB, for caching. The methods use this size for 
 *                  caching/staging the data as a limit but use other 
 *                  memory for all internal data structures.
 *                  0 means, it can use as much memory as necessary
 *  app_id
 *      when using a staging method (DATASPACES, DIMES, NSSI or DATATAP).
 *      The ID should be unique for each application accessing the staging area
 */

/** Finalize the selected method. Required for all methods that are initialized. 
 */
int adios_read_finalize_method(enum ADIOS_READ_METHOD method);

/** Open an adios file/stream as a stream.
 *  Only one step at a time can be read. The list of variables will change when
 *  advancing the step if the writing application writes different variables at
 *  different times. 
 *
 *  IN:  fname    pathname of file/stream to be opened
 *       method   read method to use for this particular stream
 *       comm     the MPI communicator of all processes that want to read data from the stream
 *                If compiled with -D_NOMPI, pass any integer here.
 *       lock_mode   In case of real streams, a step may need to be locked in memory to be able
 *                      to read all data of the step completely.
 *       timeout_sec  >=0.0: block until the stream becomes available but 
 *                           for max 'timeout_sec' seconds.
 *                           0.0 means return immediately if stream is not available
 *                     <0.0: block possibly forever
 *                     Note: <0.0 does not ever return with err_file_not_found error, 
 *                        which is dangerous if the stream name is simply mistyped in the code.
 *  RETURN:       pointer to an ADIOS_FILE struct, NULL on error (sets adios_errno)
 *
 *  Note: the current_step field of the returned struct indicates which step the stream is at.
 *        Steps start from 0.
 *
 *  Possible errors (adios_errno values):
 *       err_file_not_found_error  File/stream does not exist / not yet available
 *       err_end_of_stream         Stream has ended, nothing is available and
 *                                 no more steps should be expected.
 */
ADIOS_FILE * adios_read_open (const char * fname, 
                              enum ADIOS_READ_METHOD method, 
                              MPI_Comm comm, 
                              enum ADIOS_LOCKMODE lock_mode,
                              float timeout_sec);


/** Open an adios file as a file.
 *  Each variable can have different number of steps, which is reported in adios_inq_var(). 
 *  Arbitrary steps of a variable can be read at any time. 
 *  In general, the 'last_step' field of the returned ADIOS_FILE struct indicates
 *  the number of open-write-close cycles of the writing application, while 
 *  'current_step' is always 0. If all variables have been written at each writing step, 
 *  'last_step' equals the 'nsteps' available for each variable reported by adios_inq_var().
 *  
 *
 *  IN:  fname    pathname of file to be opened
 *       method   read method to use for this particular file
 *       comm     the MPI communicator of all processes that want to read data from the file
 *                If compiled with -D_NOMPI, pass any integer here.
 *  RETURN:       pointer to an ADIOS_FILE struct, NULL on error (sets adios_errno)
 *
 *  Possible errors (adios_errno values):
 *       err_file_not_found_error  File does not exist 
 */
ADIOS_FILE * adios_read_open_file (const char * fname, 
                                   enum ADIOS_READ_METHOD method, 
                                   MPI_Comm comm);

/** Close an adios file.
 *  It will free the content of the underlying data structures and the fp pointer itself.
 *  IN:   fp       pointer to an ADIOS_FILE struct
 *  RETURN: 0 OK, !=0 on error (adios_errno value)
 */
int adios_read_close (ADIOS_FILE *fp);


/** Advance the current step of a stream.
 *  In case of streams, 
 *     - An error should be expected for any step, 
 *       since that might not yet be available 
 *     - Also, only the current step can be read.
 *     - Seeking to step N informs the read method that all steps 
 *       before N can be removed if space is needed.
 *  For files opened as file, stepping has no effect.
 *
 *  IN:   fp       pointer to an ADIOS_FILE struct
 *        last     0: next available step, !=0: newest available step 
 *        timeout_sec  >=0.0: block until the next step becomes available but 
 *                            for max 'timeout_sec' seconds.
 *                            0.0 means return immediately if step is not available
 *                 <0.0: block forever if necessary
 *  RETURN: 0 OK, !=0 on error (adios_errno value)
 *      
 *  Possible errors (adios_errno values):
 *       err_end_of_stream    Stream has ended, no more steps should be expected
 *       err_step_notready    The requested step is not yet available
 *       err_step_disappeared The requested step is not available anymore
 */
int adios_advance_step (ADIOS_FILE *fp, int last, float timeout_sec); 

/** Release a step in a stream without seeking to the next step.
  * This function is to inform the read method that the current step is
  * no longer needed, but the reader does not yet want to read another step.
  * This function releases the lock on the step only. The current step is not
  * changed in the ADIOS_FILE struct, but reading any variable may fail due to
  * removal. If no locking is requested at open time, this function has no
  * effect at all.
  *
  * adios_advance_step() also releases the step from which one advances 
  * forward so it is not necessary to call this function.
  */
void adios_release_step (ADIOS_FILE *fp);

/** Inquiry about a variable.
 *  This function does not read anything from the file but processes info
 *  already in memory after fopen.
 *  It allocates memory for the ADIOS_VARINFO struct and content, so
 *  you need to free resources later with adios_free_varinfo().
 *
 *  Note that you can get a scalar variable's value (including strings)
 *  with this operation without touching the file/stream.
 *  The 'stats' element will be NULL after this call. To get the statistics, 
 *  another call must be made after this: adios_inq_var_stat().
 *  The 'blocks' element will be NULL after this call. To get the decomposition
 *  of a variable in the file/stream, another call must be made after this: 
 *  adios_inq_var_blockinfo().
 *
 *  IN:  fp       pointer to an (opened) ADIOS_FILE struct
 *       varname  name of the variable
 *  RETURN:       pointer to and ADIOS_VARINFO struct, NULL on error (sets adios_errno)
 */
ADIOS_VARINFO * adios_inq_var (ADIOS_FILE *fp, const char * varname);  

/** Inquiry a variable by index
 *       varid    index of variable (0..fp->nvars-1)
 *                in fp->vars_namelist of ADIOS_FILE struct
 */
ADIOS_VARINFO * adios_inq_var_byid (ADIOS_FILE *fp, int varid);   

/** Free memory used by an ADIOS_VARINFO struct */
void adios_free_varinfo (ADIOS_VARINFO *cp);

/** Get statistics recorded about a variable.
 *  The information to calculate the statistics are recorded in the metadata,
 *  so no extra file access is necessary after adios_fopen() for this operation.
 *  The result is stored in the ADIOS_VARSTAT struct under varinfo.stats. 
 * 
 *  adios_free_varinfo() will free the extra memory allocated in this call.
 *
 *  IN:  fp             pointer to an (opened) ADIOS_FILE struct
 *       varinfo        result of adios_inq_var() 
 *       per_step_stat  !=0: return statistics also per step
 *       per_block_stat !=0: return statistics also per writer block 
 *  RETURN: 0 OK, !=0 on error (adios_errno value)
 */
int adios_inq_var_stat (ADIOS_FILE *fp, ADIOS_VARINFO * varinfo,
                        int per_step_stat, int per_block_stat);

/** Get the block-decomposition of the variable about how it is stored in 
 *  the file or stream. The decomposition information are recorded in the
 *  metadata, so no extra file access is necessary after adios_fopen() for 
 *  this operation. The result is stored in the array of 
 *  ADIOS_VARBLOCK structs under varinfo.blocks. 
 * 
 *  adios_free_varinfo() will free the extra memory allocated in this call.
 *
 *  IN:  fp       pointer to an (opened) ADIOS_FILE struct
 *       varinfo  result of adios_inq_var() 
 *  RETURN: 0 OK, !=0 on error (adios_errno value)
 */
int adios_inq_var_blockinfo (ADIOS_FILE *fp, ADIOS_VARINFO * varinfo);

/** Inquiry a link by index
*       linkid   index of link (0..fp->nlinks-1)
*                in fp->link_namelist of ADIOS_FILE struct
*/
ADIOS_LINK * adios_inq_link_byid (ADIOS_FILE *fp, int linkid);

/** Free memory used by an ADIOS_LINK struct */
void adios_free_linkinfo (ADIOS_LINK *linkinfo);

/** Inquiry a mesh by index
*       meshid   index of mesh (0..fp->nmeshes-1)
*                in fp->mesh_namelist of ADIOS_FILE struct
*/
ADIOS_MESH * adios_inq_mesh_byid (ADIOS_FILE *fp, int meshid);

/** fill in the complete mesh structure if mesh struct is stored in external file  
    call adios_inq_mesh_byid() first to determine if external file is NULL
*/
int adios_complete_meshinfo (ADIOS_FILE *datafile, ADIOS_FILE *meshfile, ADIOS_MESH *meshinfo);

/** Free memory used by an ADIOS_MESH struct */
void adios_free_meshinfo (ADIOS_MESH *meshinfo);

/** Get the mesh information of the variable about how it is stored in
 * ADIOS_MESH structure
 * IN:  fp       pointer to an (opened) ADIOS_FILE struct
 *      varinfo  result of adios_inq_var() 
 * RETURN: 0 OK, !=0 on error
 */
int adios_inq_var_meshinfo (ADIOS_FILE *fp, ADIOS_VARINFO * varinfo);

/** Schedule reading a variable (slice) from the file.
 *  You need to call adios_perform_reads() to do the reading of
 *  variables.
 *  IN:  fp         pointer to an (opened) ADIOS_FILE struct
 *       sel        selection created beforehand with adios_selection...().
 *                  sel=NULL means global selection (whole variable)
 *       varname    name of the variable
 *       from_step  File mode only: read the 'nsteps' consecutive steps
 *                  of a variable from 'from_step'.
 *                  It is not used in case of a stream.
 *       nsteps     Read 'nsteps' consecutive steps of a variable.
 *                  Must be 1 for a stream. 
 *  OUT: data       pointer to the memory to hold data of the variable
 *                  In blocking read mode, the memory should be 
 *                  pre-allocated. In non-blocking mode, memory can be
 *                  allocated or not, and that changes the behavior of
 *                  the chunked read. If memory is allocated, 
 *                  adios_check_read() returns a variable if it is completed.
 *                  If memory is not allocated, the check returns any chunk
 *                  already available of a variable (in ADIOS own memory)
 *                  and the application has to rearrange the data. The user
 *                  has to process/copy the data before getting new chunks.
 *  RETURN: 0 OK, !=0 on error (adios_errno value)
 */
int adios_schedule_read (const ADIOS_FILE * fp,
                         const ADIOS_SELECTION * sel,
                         const char            * varname,
                         int                     from_steps,
                         int                     nsteps,
                         void                  * data);

/** Read a variable by index
 *       varid    index of variable (0..fp->nvars-1)
 *                in fp->vars_namelist of ADIOS_FILE struct
 */
int adios_schedule_read_byid (const ADIOS_FILE * fp, 
                              const ADIOS_SELECTION * sel,
                              int                     varid,
                              int                     from_steps,
                              int                     nsteps,
                              void                  * data);

/* NCSU ALACRITY-ADIOS: Support for those transforms that can change reading behavior (e.g., level-of-detail) */
int adios_schedule_read_param (const ADIOS_FILE * fp,
                               const ADIOS_SELECTION * sel,
                               const char            * varname,
                               int                     from_steps,
                               int                     nsteps,
                               const char            * param,
                               void                  * data);

int adios_schedule_read_byid_param (const ADIOS_FILE * fp,
                                    const ADIOS_SELECTION * sel,
                                    int                     varid,
                                    int                     from_steps,
                                    int                     nsteps,
                                    const char            * param,
                                    void                  * data);


/** Let ADIOS perform the scheduled reads 
 *  IN:  blocking  If non-zero, return only when all reads are completed.
 *                 If zero, return immediately and report partial completions
 *                 through adios_check_reads()
 */
int adios_perform_reads (const ADIOS_FILE *fp, int blocking);

/** Get a chunk of completed read(s) in a chunking read scenario.
 *  Note that a selection of an array specified in a read may be completed in
 *  multiple chunks (usually when they come from multiple sources, like
 *  different disks or different application processes). 
 *  This function should be called in a loop until all chunks are processed. 
 *  That is indicated by a 0 return value. A NULL result for chunk only
 *  indicates that no chunk is available at the time of call. 
 *
 *  IN:  fp         Handler to file or stream
 *  OUT: chunk      A chunk completed by the time of calling this function.
 *                  It is NULL if no chunk is returned.
 *  RETURN:         0: all chunks have been returned previously, 
 *                     no need to call again (chunk is NULL, too)
 *                  1: some chunks are/will be available, call again
 *                  <0 on error (adios_errno value)
 */
int adios_check_reads (const ADIOS_FILE  * fp, 
                       ADIOS_VARCHUNK   ** chunk);

/** Free the memory of a chunk allocated inside adios_check_reads().
  * It only frees the ADIOS_VARCHUNK struct and the ADIOS_SELECTION struct
  * pointed by the chunk. The data pointer should never be freed since
  * that memory belongs to the reading method.
  */
void adios_free_chunk (ADIOS_VARCHUNK *chunk);


/** Get an attribute in a file.
 *  This function does not read anything from the file but processes info
 *  already in memory after fopen.
 *  The memory for the data is allocated within the library.
 *  You can use free() to free the memory after use.
 *
 *  IN:  fp       pointer to an (opened) ADIOS_FILE struct
 *       attrname name of the attribute
 *  OUT: type     adios type of attribute (see enum ADIOS_DATATYPES in adios_types.h)
 *       size     memory size of value (n+1 for a string of n characters)
 *       data     pointer to the value. You need to cast it afterward according to the type.
 *  RETURN: 0 OK, !=0 on error (adios_errno value)
 */
int adios_get_attr (ADIOS_FILE            * fp,
                    const char            * attrname,
                    enum ADIOS_DATATYPES  * type,
                    int                   * size,
                    void                 ** data);

/** Convenience function to get an attribute by name
 *       attrid   index of attribute (0..fp->nattrs-1)
 *                in fp->attr_namelist of ADIOS_FILE struct
 */
int adios_get_attr_byid (ADIOS_FILE  * fp, int attrid, 
                         enum ADIOS_DATATYPES * type,
                         int * size, void ** data);

/** Return the name of an adios type */
const char * adios_type_to_string (enum ADIOS_DATATYPES type);

/** Return the memory size of one data element of an adios type.
 *  If the type is adios_string, and the second argument is
 *  the string itself, it returns strlen(data)+1.
 *  For other types, it does not care about data and returns
 *  the size occupied by one element.
 */
int adios_type_size(enum ADIOS_DATATYPES type, const void *data);



/* ADIOS Read API v1 style grouping  */

/** Return the list of groups (names) that are written into
 *  the file. There is always at least one group there.
 *
 *  IN:   fp              pointer to an (opened) ADIOS_FILE struct
 *  OUT:  group_namelist  list of strings,
 *                        Note: one should pass a pointer to char** list
 *  RETURN:               number of groups, <0 on error (adios_errno value)
 *  
 */
int adios_get_grouplist (ADIOS_FILE  *fp, char ***group_namelist);

/** Restrict the view of variables/attributes to a certain group.
 *  The provided ADIOS_FILE structure is directly modified but
 *  another calls can change to a different group view, or reset
 *  back to full view.
 *
 *  IN:       groupid    id of the selected group (0..# of groups-1)
 *                       use -1 to reset to the complete list
 *  IN/OUT:   fp         pointer to an (opened) ADIOS_FILE struct
 *                       nvars, var_namelist, nattrs, and
 *                       attr_namelist will be modified.
 *  RETURN:   0 OK, !=0 on error (adios_errno value)
 *  Note: a stream does not have groups. Only a file can have
 *  multiple groups (from separate adios_open/adios_close operations)
 */
int adios_group_view (ADIOS_FILE  *fp, int groupid);

/** OTHER API FUNCTIONS **/

/* NCSU - Timeseries analysis functions */
double adios_stat_cor (ADIOS_VARINFO * vix, ADIOS_VARINFO * viy, char * characteristic, uint32_t time_start, uint32_t time_end, uint32_t lag);
double adios_stat_cov (ADIOS_VARINFO * vix, ADIOS_VARINFO * viy, char * characteristic, uint32_t time_start, uint32_t time_end, uint32_t lag);

/** This function can be called if user defines the sequences of 
 *  dimensions in reversed order for a variable in the ADIOS XML 
 *  definition to view the variable with the correct order of 
 *  dimensions. Call this before reading variables.
 */
void adios_reset_dimension_order (ADIOS_FILE *fp, int is_fortran);

/** Test function to print basic info about the file to stdout */
void adios_print_fileinfo (ADIOS_FILE *fp);


typedef struct
{
    int nmethods; 		// number of available read methods
    char ** name;     	// array of read method names
    enum ADIOS_READ_METHOD * methodID; // enum ID of the method in source code
} ADIOS_AVAILABLE_READ_METHODS;

/* Provide the names of read methods available in the running application
 */
ADIOS_AVAILABLE_READ_METHODS * adios_available_read_methods();

void adios_available_read_methods_free (ADIOS_AVAILABLE_READ_METHODS *);
]]
