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
]]
