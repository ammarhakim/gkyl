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
local new, typeof = xsys.from(ffi,
     "new, typeof")

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
               ncclNumResults              =  6 } ncclResult_t;

/* Communicator configuration. Users can assign value to attributes to specify the
 * behavior of a communicator. */
typedef struct ncclConfig_v21400 {
  /* attributes that users should never touch. */
  size_t size;
  unsigned int magic;
  unsigned int version;
  /* attributes that users are able to customize. */
  int blocking;
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

/* Frees communicator. */
ncclResult_t  ncclCommDestroy(ncclComm_t comm);

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
]]

-- ncclDataType_t enums. Have to match what's in nccl.h
_M.Int = 2
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

-- ncclConfig_t object.
function _M.Config() return new_ncclConfig() end

-- ncclGetUniqueId.
function _M.GetUniqueId(id)
   return ffiC.ncclGetUniqueId(id)
end

-- ncclCommInitRank.
function _M.CommInitRank(comm, nranks, commId, rank)
   return ffiC.ncclCommInitRank(comm, nranks, commId, rank)
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

-- ncclRecv
function _M.Recv(recvbuff, count, datatype, peer, comm, stream)
   return ffiC.ncclRecv(recvbuff, count, datatype, peer,
      getObj(comm, "ncclComm_t[1]"), getObj(stream, "cudaStream_t[1]"))
end

-- ncclCommDestroy.
function _M.CommDestroy(comm)
   return ffiC.ncclCommDestroy(getObj(comm, "ncclComm_t[1]"))
end

return _M
