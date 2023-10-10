-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for NCCL.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- don't bother if MPI if not built in
assert(GKYL_HAVE_MPI, "Gkyl was not built with MPI!")

local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, typeof = xsys.from(ffi, "new, typeof")

local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"
local cuda = require "Cuda.RunTime"

-- global count of barriers
local numMpiBarrier = 0
-- time spent in barrier function
local timeMpiBarrier = 0.0

local _M = {}

ffi.cdef [[

/* Opaque handle to communicator */
typedef struct ncclComm* ncclComm_t;

// #define NCCL_UNIQUE_ID_BYTES 128
typedef struct { char internal[128]; } ncclUniqueId;

/* Error type */
typedef enum { ncclSuccess                 =  0,
               ncclUnhandledCudaError      =  1,
               ncclSystemError             =  2,
               ncclInternalError           =  3,
               ncclInvalidArgument         =  4,
               ncclInvalidUsage            =  5,
               ncclRemoteError             =  6,
               ncclInProgress              =  7,
               ncclNumResults              =  8 } ncclResult_t;

/* Communicator configuration. Users can assign value to attributes to specify the
 * behavior of a communicator. */
typedef struct ncclConfig_v21700 {
  /* attributes that users should never touch. */
  size_t size;
  unsigned int magic;
  unsigned int version;
  /* attributes that users are able to customize. */
  int blocking;
  int cgaClusterSize;
  int minCTAs;
  int maxCTAs;
  const char *netName;
  int splitShare;
} ncclConfig_t;

/* Data types */
typedef enum { ncclInt8       = 0, ncclChar       = 0,
               ncclUint8      = 1,
               ncclInt32      = 2, ncclInt        = 2,
               ncclUint32     = 3,
               ncclInt64      = 4,
               ncclUint64     = 5,
               ncclFloat16    = 6, ncclHalf       = 6,
               ncclFloat32    = 7, ncclFloat      = 7,
               ncclFloat64    = 8, ncclDouble     = 8,
// MF Commented these out presumably because LuaJIT can't handle preprocessor instructions.
// #if defined(__CUDA_BF16_TYPES_EXIST__)
//                ncclBfloat16   = 9,
//                ncclNumTypes   = 10
// #else
//                ncclNumTypes   = 9
// #endif
               ncclNumTypes   = 9
} ncclDataType_t;

/* Reduction operation selector */
typedef enum { ncclNumOps_dummy = 5 } ncclRedOp_dummy_t;
typedef enum { ncclSum        = 0,
               ncclProd       = 1,
               ncclMax        = 2,
               ncclMin        = 3,
               ncclAvg        = 4,
               /* ncclNumOps: The number of built-in ncclRedOp_t values. Also
                * serves as the least possible value for dynamic ncclRedOp_t's
                * as constructed by ncclRedOpCreate*** functions. */
               ncclNumOps     = 5,
               /* ncclMaxRedOp: The largest valid value for ncclRedOp_t.
                * It is defined to be the largest signed value (since compilers
                * are permitted to use signed enums) that won't grow
                * sizeof(ncclRedOp_t) when compared to previous NCCL versions to
                * maintain ABI compatibility. */
               ncclMaxRedOp   = 0x7fffffff>>(32-8*sizeof(ncclRedOp_dummy_t))
             } ncclRedOp_t;


/* Generates an Id to be used in ncclCommInitRank. */
ncclResult_t  ncclGetUniqueId(ncclUniqueId* uniqueId);

/* Comm creation. */
ncclResult_t  ncclCommInitRank(ncclComm_t* comm, int nranks, ncclUniqueId commId, int rank);

/* Create a new communicator (multi thread/process version) with a configuration
 * set by users. */
ncclResult_t  ncclCommInitRankConfig(ncclComm_t* comm, int nranks, ncclUniqueId commId, int rank, ncclConfig_t* config);

/* Frees communicator. */
ncclResult_t  ncclCommDestroy(ncclComm_t comm);

/* Checks whether the comm has encountered any asynchronous errors */
ncclResult_t  ncclCommGetAsyncError(ncclComm_t comm, ncclResult_t *asyncError);

/**
* AllGather
* Gather sendcount values from all GPUs into recvbuff, receiving data from rank i at offset i*sendcount.
*
* Note: This assumes the receive count is equal to nranks*sendcount, which means that
* recvbuff should have a size of at least nranks*sendcount elements.
*
* In-place operation will happen if sendbuff == recvbuff + rank * sendcount.
*/
ncclResult_t ncclAllGather(const void* sendbuff, void* recvbuff, size_t count, ncclDataType_t datatype, ncclComm_t comm, cudaStream_t stream);


/**
 * Allreduce
 *
 * Reduce data arrays of length count in sendbuff using op operation
 * and leaves identical copies of the result on each recvbuff.
 *
 * In-place operation will happen if sendbuff == recvbuff.
 */
ncclResult_t ncclAllReduce(const void* sendbuff, void* recvbuff, size_t count, ncclDataType_t datatype, ncclRedOp_t op, ncclComm_t comm, cudaStream_t stream);

/*
 * Receive
 *
 * Receive data from rank peer into recvbuff.
 *
 * Rank peer needs to call ncclSend with the same datatype and the same count to this
 * rank.
 *
 * This operation is blocking for the GPU. If multiple ncclSend and ncclRecv operations
 * need to progress concurrently to complete, they must be fused within a ncclGroupStart/
 * ncclGroupEnd section.
 */
 ncclResult_t  ncclRecv(void* recvbuff, size_t count, ncclDataType_t datatype, int peer,
 ncclComm_t comm, cudaStream_t stream);

/*
 * Send
 *
 * Send data from sendbuff to rank peer.
 *
 * Rank peer needs to call ncclRecv with the same datatype and the same count from this
 * rank.
 *
 * This operation is blocking for the GPU. If multiple ncclSend and ncclRecv operations
 * need to progress concurrently to complete, they must be fused within a ncclGroupStart/
 * ncclGroupEnd section.
 */
ncclResult_t  ncclSend(const void* sendbuff, size_t count, ncclDataType_t datatype, int peer,
    ncclComm_t comm, cudaStream_t stream);

/*
 * Group Start
 *
 * Start a group call. All calls to NCCL until ncclGroupEnd will be fused into
 * a single NCCL operation. Nothing will be started on the CUDA stream until
 * ncclGroupEnd.
 */
ncclResult_t  ncclGroupStart();

/*
 * Group End
 *
 * End a group call. Start a fused NCCL operation consisting of all calls since
 * ncclGroupStart. Operations on the CUDA stream depending on the NCCL operations
 * need to be called after ncclGroupEnd.
 */
ncclResult_t  ncclGroupEnd();

/*
   Gkyl functions to wrap other NCCL objects/functions.
*/
void gkyl_NCCL_CONFIG_INITIALIZER(ncclConfig_t *nc);
]]

-- ncclResult_t enums. Have to match what's in nccl.h
_M.Success                 =  0
_M.UnhandledCudaError      =  1
_M.SystemError             =  2
_M.InternalError           =  3
_M.InvalidArgument         =  4
_M.InvalidUsage            =  5
_M.RemoteError             =  6
_M.InProgress              =  7
_M.NumResults              =  8

-- ncclDataType_t enums. Have to match what's in nccl.h
_M.Int    = 2
_M.Float  = 7
_M.Double = 8

-- ncclRedOp_t enum have to match the definition in nccl.h
_M.Sum  = 0
_M.Prod = 1
_M.Max  = 2
_M.Min  = 3
_M.Avg  = 4

local function new_ncclComm()
   local comm = new("ncclComm_t[1]")
   ffi.gc(comm, function(c) ffiC.ncclCommDestroy(c[0]) end)
   return comm
end

local function getObj(obj, ptyp)
   return ffi.istype(typeof(ptyp), obj) and obj[0] or obj
end

local function new_ncclConfig()
   return new("ncclConfig_t[1]")
end

-- ncclUniqueId object.
function _M.UniqueId() return new("ncclUniqueId") end

-- ncclComm_t object.
function _M.Comm() return new_ncclComm() end

-- ncclResult_t enum.
function _M.Result() return ffi.new("ncclResult_t[1]") end

-- ncclConfig_t object.
function _M.Config()
   local conf = new_ncclConfig()
   ffiC.gkyl_NCCL_CONFIG_INITIALIZER(conf)
   return conf
end

-- ncclGetUniqueId.
function _M.GetUniqueId(id)
   return ffiC.ncclGetUniqueId(id)
end

-- ncclCommInitRank.
function _M.CommInitRank(comm, nranks, commId, rank)
   return ffiC.ncclCommInitRank(comm, nranks, commId, rank)
end

-- ncclCommInitRankConfig.
function _M.CommInitRankConfig(comm, nranks, commId, rank, config)
   return ffiC.ncclCommInitRankConfig(comm, nranks, commId, rank, config)
end

-- ncclAllGather
function _M.AllGather(sendbuff, recvbuff, count, datatype, comm, stream)
   return ffiC.ncclAllGather(sendbuff, recvbuff, count, datatype,
      getObj(comm, "ncclComm_t[1]"), getObj(stream, "cudaStream_t[1]"))
end

-- ncclAllReduce.
function _M.AllReduce(sendbuff, recvbuff, count, datatype, op, comm, stream)
   return ffiC.ncclAllReduce(sendbuff, recvbuff, count, datatype, op,
      getObj(comm, "ncclComm_t[1]"), getObj(stream, "cudaStream_t[1]"))
end

-- ncclSend
function _M.Send(sendbuff, count, datatype, peer, comm, stream)
   return ffiC.ncclSend(sendbuff, count, datatype, peer,
      getObj(comm, "ncclComm_t[1]"), getObj(stream, "cudaStream_t[1]"))
end

-- ncclGroupStart
function _M.GroupStart()
   return ffiC.ncclGroupStart()
end

-- ncclGroupEnd
function _M.GroupEnd()
   return ffiC.ncclGroupEnd()
end

-- ncclRecv
function _M.Recv(recvbuff, count, datatype, peer, comm, stream)
   return ffiC.ncclRecv(recvbuff, count, datatype, peer,
      getObj(comm, "ncclComm_t[1]"), getObj(stream, "cudaStream_t[1]"))
end

-- ncclCommDestroy.
function _M.CommDestroy(comm)
   return ffiC.ncclCommDestroy(getObj(comm, "ncclComm_t[1]"))
end

-- ncclCommGetAsyncError.
function _M.CommGetAsyncError(comm, asyncError)
   return ffiC.ncclCommGetAsyncError(getObj(comm, "ncclComm_t[1]"), asyncError)
end

return _M
