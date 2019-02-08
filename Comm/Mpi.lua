-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for MPI
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

local _M = {}

ffi.cdef [[
  // Opaque types
  typedef struct MPI_Comm_type *MPI_Comm;
  typedef struct MPI_Datatype_type *MPI_Datatype;
  typedef struct MPI_Op_type *MPI_Op;
  typedef struct MPI_Status_type *MPI_Status;
  typedef struct MPI_Group_type *MPI_Group;
  typedef struct MPI_Request_type *MPI_Request;
  typedef struct MPI_Info_type *MPI_Info;
  typedef struct MPI_Win_type *MPI_Win;
  typedef unsigned MPI_Aint; /* Not sure if this is correct, but it seems to work */

  // size of various objects
  int sizeof_MPI_Status();
  int sizeof_MPI_Group();
  int sizeof_MPI_Comm();
  int sizeof_MPI_Request();

  // Pre-defined objects and constants
  MPI_Comm get_MPI_COMM_WORLD();
  MPI_Comm get_MPI_COMM_NULL();
  MPI_Comm get_MPI_COMM_SELF();
  MPI_Request get_MPI_REQUEST_NULL();
  MPI_Status *getPtr_MPI_STATUS_IGNORE();
  MPI_Info get_MPI_INFO_NULL();
  int get_MPI_COMM_TYPE_SHARED();
  int get_MPI_UNDEFINED();
  int get_MPI_ORDER_C();
  int get_MPI_ORDER_FORTRAN();

  // Datatypes
  MPI_Datatype get_MPI_CHAR();
  MPI_Datatype get_MPI_BYTE();
  MPI_Datatype get_MPI_SHORT();
  MPI_Datatype get_MPI_INT();
  MPI_Datatype get_MPI_LONG();
  MPI_Datatype get_MPI_FLOAT();
  MPI_Datatype get_MPI_DOUBLE();
  MPI_Datatype get_MPI_UNSIGNED_CHAR();
  MPI_Datatype get_MPI_UNSIGNED_SHORT();
  MPI_Datatype get_MPI_UNSIGNED();
  MPI_Datatype get_MPI_UNSIGNED_LONG();
  MPI_Datatype get_MPI_LONG_DOUBLE();
  MPI_Datatype get_MPI_LONG_LONG_INT();
  MPI_Datatype get_MPI_FLOAT_INT();
  MPI_Datatype get_MPI_LONG_INT();
  MPI_Datatype get_MPI_DOUBLE_INT();
  MPI_Datatype get_MPI_SHORT_INT();
  MPI_Datatype get_MPI_2INT();
  MPI_Datatype get_MPI_LONG_DOUBLE_INT();
  MPI_Datatype get_MPI_PACKED();

  // Operators
  MPI_Op get_MPI_MAX();
  MPI_Op get_MPI_MIN();
  MPI_Op get_MPI_SUM();
  MPI_Op get_MPI_PROD();
  MPI_Op get_MPI_LAND();
  MPI_Op get_MPI_BAND();
  MPI_Op get_MPI_LOR();
  MPI_Op get_MPI_BOR();
  MPI_Op get_MPI_LXOR();
  MPI_Op get_MPI_BXOR();
  MPI_Op get_MPI_MINLOC();
  MPI_Op get_MPI_MAXLOC();

  // Errors
  int get_MPI_SUCCESS();

  // Communicators
  int MPI_Comm_rank(MPI_Comm comm, int *rank);
  int MPI_Comm_size(MPI_Comm comm, int *size);
  int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
  int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);

  // MPI_Datatype handling
  int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);
  int MPI_Type_vector(int count,
    int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype * newtype);
  int MPI_Type_indexed(int count, const int array_of_blocklengths[],
    const int array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype);
  int MPI_Type_create_subarray(int ndims, const int array_of_sizes[], const
    int array_of_subsizes[], const int array_of_starts[], int order, MPI_Datatype
    oldtype, MPI_Datatype *newtype);
  int MPI_Type_commit(MPI_Datatype * datatype);
  int MPI_Type_free(MPI_Datatype *datatype);

  // Win & SHM calls
  int MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm);
  int MPI_Win_allocate_shared (MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
  int MPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr);
  int MPI_Win_lock_all(int assert, MPI_Win win);
  int MPI_Win_unlock_all(MPI_Win win);
  int MPI_Win_free(MPI_Win *win);

  // Point-to-point communication
  int MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count);
  int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
		    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
  int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
                  MPI_Comm comm );
  int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
	       MPI_Comm comm);
  int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
	       MPI_Comm comm, MPI_Status *status);
  int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
		int tag, MPI_Comm comm, MPI_Request *request);
  int MPI_Wait(MPI_Request *request, MPI_Status *status);

  // Groups
  int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
  int MPI_Group_rank(MPI_Group group, int *rank);
  int MPI_Group_size(MPI_Group group, int *size);
  int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup);
  int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
  int MPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[], MPI_Group group2, int ranks2[]);

  // Global operators
  int MPI_Barrier(MPI_Comm comm);
  int MPI_Abort(MPI_Comm comm, int errorcode);

  // Gkyl utility functions
  void GkMPI_fillStatus(const MPI_Status* inStatus, int *outStatus);
]]

-- Predefined objects and constants
_M.COMM_WORLD = ffiC.get_MPI_COMM_WORLD()
_M.COMM_NULL = ffiC.get_MPI_COMM_NULL()
_M.COMM_SELF = ffiC.get_MPI_COMM_SELF()
_M.REQUEST_NULL = ffiC.get_MPI_REQUEST_NULL()
_M.STATUS_IGNORE = ffiC.getPtr_MPI_STATUS_IGNORE()

_M.INFO_NULL = ffiC.get_MPI_INFO_NULL()
_M.COMM_TYPE_SHARED = ffiC.get_MPI_COMM_TYPE_SHARED()
_M.UNDEFINED = ffiC.get_MPI_UNDEFINED()
_M.ORDER_C = ffiC.get_MPI_ORDER_C()
_M.ORDER_FORTRAN = ffiC.get_MPI_ORDER_FORTRAN()

-- Object sizes
_M.SIZEOF_STATUS = ffiC.sizeof_MPI_Status()

-- Datatypes
_M.CHAR = ffiC.get_MPI_CHAR()
_M.BYTE = ffiC.get_MPI_BYTE()
_M.SHORT = ffiC.get_MPI_SHORT()
_M.INT = ffiC.get_MPI_INT()
_M.LONG = ffiC.get_MPI_LONG()
_M.FLOAT = ffiC.get_MPI_FLOAT()
_M.DOUBLE = ffiC.get_MPI_DOUBLE()
_M.UNSIGNED_CHAR = ffiC.get_MPI_UNSIGNED_CHAR()
_M.UNSIGNED_SHORT = ffiC.get_MPI_UNSIGNED_SHORT()
_M.UNSIGNED = ffiC.get_MPI_UNSIGNED()
_M.UNSIGNED_LONG = ffiC.get_MPI_UNSIGNED_LONG()
_M.LONG_DOUBLE = ffiC.get_MPI_LONG_DOUBLE()
_M.LONG_LONG_INT = ffiC.get_MPI_LONG_LONG_INT()
_M.FLOAT_INT = ffiC.get_MPI_FLOAT_INT()
_M.LONG_INT = ffiC.get_MPI_LONG_INT()
_M.DOUBLE_INT = ffiC.get_MPI_DOUBLE_INT()
_M.SHORT_INT = ffiC.get_MPI_SHORT_INT()
_M.TWOINT = ffiC.get_MPI_2INT()
_M.LONG_DOUBLE_INT = ffiC.get_MPI_LONG_DOUBLE_INT()
_M.PACKED = ffiC.get_MPI_PACKED()

-- Operators
_M.MAX = ffiC.get_MPI_MAX()
_M.MIN = ffiC.get_MPI_MIN()
_M.SUM = ffiC.get_MPI_SUM()
_M.PROD = ffiC.get_MPI_PROD()
_M.LAND = ffiC.get_MPI_LAND()
_M.BAND = ffiC.get_MPI_BAND()
_M.LOR = ffiC.get_MPI_LOR()
_M.BOR = ffiC.get_MPI_BOR()
_M.LXOR = ffiC.get_MPI_LXOR()
_M.BXOR = ffiC.get_MPI_BXOR()
_M.MINLOC = ffiC.get_MPI_MINLOC()
_M.MAXLOC = ffiC.get_MPI_MAXLOC()

-- Error codes
_M.SUCCESS = ffiC.get_MPI_SUCCESS()

-- some types for use in MPI functions
local int_1 = typeof("int[1]")
local uint_1 = typeof("unsigned[1]")
local voidp = typeof("void *[1]")

-- ctors for various MPI objects: MPI_Status object is not a opaque
-- pointer and so we need to allocate memory for it. Other object
-- resources are managed by MPI itself, so we can't allocate space for
-- them. Doing so will lead to program crash as the LJ GC will try to
-- free memory which is managed by MPI. A very bad thing.
local function new_MPI_Status()
   return ffi.cast("MPI_Status*", new("uint8_t[?]", _M.SIZEOF_STATUS))
end
local function new_MPI_Comm()
   return new("MPI_Comm[1]")
end
local function new_MPI_Group()
   return new("MPI_Group[1]")
end
local function new_MPI_Request()
   return new("MPI_Request[1]")
end
local function new_MPI_Win()
   return new("MPI_Win[1]")
end
local function new_MPI_Datatype()
   return new("MPI_Datatype[1]")
end

-- de-reference if object is a pointer
local function getObj(obj, ptyp)
   return ffi.istype(typeof(ptyp), obj) and obj[0] or obj
end

-- Extract object from (potentially) a pointer
function _M.getComm(obj)
   return getObj(obj, "MPI_Comm[1]")
end

-- MPI_Status object
_M.Status = {}
function _M.Status:new()
   local self = setmetatable({}, _M.Status)
   self.SOURCE, self.TAG, self.ERROR = 0, 0, 0
   self.mpiStatus = new_MPI_Status()
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(_M.Status, { __call = function (self, o) return self.new(self, o) end })

-- MPI_Comm object
function _M.Comm()  return new_MPI_Comm() end
-- MPI_Group object
function _M.Group()  return new_MPI_Group() end
-- MPI_Datatype object
function _M.MPI_Datatype() return new_MPI_Datatype() end

-- MPI_Comm_rank
function _M.Comm_rank(comm)
   local r = int_1()
   local _ = ffiC.MPI_Comm_rank(getObj(comm, "MPI_Comm[1]"), r)
   return r[0]
end
-- MPI_Comm_size
function _M.Comm_size(comm)
   local r = int_1()
   local _ = ffiC.MPI_Comm_size(getObj(comm, "MPI_Comm[1]"), r)
   return r[0]
end
-- MPI_Comm_dup
function _M.Comm_dup(comm)
   local c = new_MPI_Comm()
   local _ = ffiC.MPI_Comm_dup(getObj(comm, "MPI_Comm[1]"), c)
   return c
end
-- MPI_Comm_split
function _M.Comm_split(comm, color, key)
   local c = new_MPI_Comm()
   local _ = ffiC.MPI_Comm_split(getObj(comm, "MPI_Comm[1]"), color, key, c)
   return c
end
-- MPI_Comm_split_type
function _M.Comm_split_type(comm, split_type, key, info)
   local c = new_MPI_Comm()
   local _ = ffiC.MPI_Comm_split_type(getObj(comm, "MPI_Comm[1]"), split_type, key, info, c)
   return c
end
-- MPI_Win_allocate_shared
function _M.Win_allocate_shared(sz, disp_unit, info, comm)
   local w = new_MPI_Win()
   local baseptr = voidp()
   local _ = ffiC.MPI_Win_allocate_shared(sz, disp_unit, info, getObj(comm, "MPI_Comm[1]"), baseptr, w)
   return baseptr[0], w
end
-- MPI_Win_shared_query
function _M.Win_shared_query(win, rank)
   local sz, du = uint_1(), int_1()
   local baseptr = voidp()
   local _ = ffiC.MPI_Win_shared_query(getObj(win, "MPI_Win[1]"), rank, sz, du, baseptr)
   return sz[0], du[0], baseptr[0]
end
-- MPI_Win_lock_all
function _M.Win_lock_all(assert, win)
   local _ = ffiC.MPI_Win_lock_all(assert, getObj(win, "MPI_Win[1]"))
end
-- MPI_Win_unlock_all
function _M.Win_unlock_all(win)
   local _ = ffiC.MPI_Win_unlock_all(getObj(win, "MPI_Win[1]"))
end
-- MPI_Win_free
function _M.Win_free(win)
   local _ = ffiC.MPI_Win_free(win)
end

-- MPI_Type_contiguous
function _M.Type_contiguous(count, oldtype)
   local t = new_MPI_Datatype()
   local _ = ffiC.MPI_Type_contiguous(count, getObj(oldtype, "MPI_Datatype[1]"), t)
   return t
end
-- MPI_Type_vector
function _M.Type_vector(count, blocklength, stride, oldtype)
   local t = new_MPI_Datatype()
   local _ = ffiC.MPI_Type_vector(
      count, blocklength, stride,
      getObj(oldtype, "MPI_Datatype[1]"), t)
   return t
end
-- MPI_Type_indexed
function _M.Type_indexed(array_of_blocklengths, array_of_displacements, oldtype)
   local t = new_MPI_Datatype()
   local count = #array_of_blocklengths
   local _ = ffiC.MPI_Type_indexed(
      count, array_of_blocklengths:data(), array_of_displacements:data(),
      getObj(oldtype, "MPI_Datatype[1]"), t)
   return t
end
-- MPI_Type_create_subarray
function _M.Type_create_subarray(array_of_sizes, array_of_subsizes, array_of_starts, order, oldtype)
   local t = new_MPI_Datatype()
   local ndims = #array_of_sizes
   local _ = ffiC.MPI_Type_create_subarray(
      ndims, array_of_sizes,
      array_of_subsizes, array_of_starts, order,
      getObj(oldtype, "MPI_Datatype[1]"), t)
   return t
end
-- MPI_Type_commit
function _M.Type_commit(datatype)
   local _ = ffiC.MPI_Type_commit(datatype)
   return datatype
end
-- MPI_Type_free
function _M.Type_free(datatype)
   local _ = ffiC.MPI_Type_free(datatype)
end

-- MPI_Get_count
function _M.Get_count(status, datatype)
   local r = int_1()
   local _ = ffiC.MPI_Get_count(status.mpiStatus, datatype, r)
   return r[0]
end
-- MPI_Allreduce
function _M.Allreduce(sendbuf, recvbuf, count, datatype, op, comm)
   local _ = ffiC.MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Bcast
function _M.Bcast(buffer, count, datatype, root, comm)
   local _ = ffiC.MPI_Bcast(buffer, count, datatype, root, getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Send
function _M.Send(buf, count, datatype, dest, tag, comm)
   local _ = ffiC.MPI_Send(
      buf, count, getObj(datatype, "MPI_Datatype[1]"), dest, tag, getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Recv
function _M.Recv(buf, count, datatype, source, tag, comm, status)
   local st = status and status.mpiStatus or _M.STATUS_IGNORE
   local _ = ffiC.MPI_Recv(
      buf, count, getObj(datatype, "MPI_Datatype[1]"), source, tag, getObj(comm, "MPI_Comm[1]"), st)
   -- store MPI_Status
   if status ~= nil then
      local gks = new("int[3]")
      ffiC.GkMPI_fillStatus(st, gks)
      status.SOURCE, status.TAG, status.ERROR = gks[0], gks[1], gks[2]
   end
end
-- MPI_Irecv
function _M.Irecv(buf, count, datatype, source, tag, comm)
   local req = new_MPI_Request()
   local _ = ffiC.MPI_Irecv(
      buf, count, datatype, source, tag, getObj(comm, "MPI_Comm[1]"), req)
   return req
end
-- MPI_Wait
function _M.Wait(request, status)
   local st = status and status.mpiStatus or _M.STATUS_IGNORE
   local _ = ffiC.MPI_Wait(request, st)
   -- store MPI_Status
   if status ~= nil then
      local gks = new("int[3]")
      ffiC.GkMPI_fillStatus(st, gks)
      status.SOURCE, status.TAG, status.ERROR = gks[0], gks[1], gks[2]
   end
end

-- MPI_Barrier
function _M.Barrier(comm)
   local _ = ffiC.MPI_Barrier(getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Abort
function _M.Abort(comm, errCode)
   local _ = ffiC.MPI_Abort(getObj(comm, "MPI_Comm[1]"), errCode)
end

-- MPI_Comm_group
function _M.Comm_group(comm)
   local grp = new_MPI_Group()
   local _ = ffiC.MPI_Comm_group(getObj(comm, "MPI_Comm[1]"), grp)
   return grp
end
-- MPI_Group_rank
function _M.Group_rank(group)
   local r = int_1()
   local _ = ffiC.MPI_Group_rank(getObj(group, "MPI_Group[1]"), r)
   return r[0]
end
-- MPI_Group_size
function _M.Group_size(group)
   local r = int_1()
   local _ = ffiC.MPI_Group_size(getObj(group, "MPI_Group[1]"), r)
   return r[0]
end
-- MPI_Group_incl
function _M.Group_incl(group, n, ranks)
   local newgroup = new_MPI_Group()
   local _ = ffiC.MPI_Group_incl(getObj(group, "MPI_Group[1]"), n, ranks, newgroup)
   return newgroup
end
-- MPI_Comm_create
function _M.Comm_create(comm, group)
   local c = new_MPI_Comm()
   local _ = ffiC.MPI_Comm_create(
      getObj(comm, "MPI_Comm[1]"), getObj(group, "MPI_Group[1]"), c)
   return c
end
-- MPI_Group_translate_ranks
function _M.Group_translate_ranks(group1, ranks1, group2)
   local n = #ranks1
   local ranks2 = Lin.IntVec(n)
   local _ = ffiC.MPI_Group_translate_ranks(
      getObj(group1, "MPI_Group[1]"), n, ranks1:data(), getObj(group2, "MPI_Group[1]"), ranks2:data())
   return ranks2
end

-- Convenience functions (these are not wrappers over MPI but make
-- some things a little cleaner)

-- Collect 'ranks' from 'comm' and create a new communicator with just
-- those ranks. The 'ranks' parameter must be of a Linalg.IntVec object
function _M.Split_comm(comm, ranks)
   local grp = _M.Comm_group(comm)
   local newGrp = _M.Group_incl(grp, #ranks, ranks:data())
   return _M.Comm_create(comm, newGrp)
end

-- Check if the communicator is NULL
function _M.Is_comm_valid(comm)
   return getObj(comm, "MPI_Comm[1]") ~=  _M.COMM_NULL
end

-- Makes a MPI_Datatype object from block sizes and offsets
function _M.createDataTypeFromBlockSizeAndOffset(blockSize, blockOffset, oldtype)
   if #blockSize == 1 then
      return _M.Type_contiguous(blockSize[1], oldtype)
   end

   -- check if sizes of all blocks and relative offsets are the same
   local isVectorType = true
   local sz = blockSize[1]
   for i = 2, #blockSize do
      if blockSize[i] ~= sz then
	 isVectorType = false
      end
   end
   if #blockSize > 1 then
      local stride = blockOffset[2]-blockOffset[1]
      for i = 2, #blockSize do
	 if stride ~= blockOffset[i]-blockOffset[i-1] then
	    isVectorType = false
	 end
      end
   end

   if isVectorType then
      return _M.Type_vector(
	 #blockOffset,
	 blockSize[1],
	 #blockOffset>1 and blockOffset[2]-blockOffset[1] or 0,
	 oldtype)
   else
      -- must be contructed as an indexed type
      local array_of_blocklengths, array_of_displacements
	 = Lin.IntVec(#blockSize), Lin.IntVec(#blockSize)
      array_of_blocklengths:setValues(blockSize)
      array_of_displacements:setValues(blockOffset)
      return _M.Type_indexed(array_of_blocklengths, array_of_displacements, oldtype)
   end
end

-- Constructs block offsets and sizes from range object
function _M.createBlockInfoFromRange(dir, range, nlayer, numComponent, ordering)
   local indexer = range:genIndexer(ordering)
   local rFace = range:shorten(dir, nlayer)
   local lo = rFace:lowerAsVec()
   local up = rFace:upperAsVec()

   local blockSize, blockOffset = {}, {}

   local dist = indexer(up)-indexer(lo)+1
   if dist == rFace:volume() then
      blockSize[1] = rFace:volume()*numComponent
      blockOffset[1] = 0
   else
      -- determine size of each block and its zero-based index (MPI
      -- expects zero-based indices)
      local currBlockSize, lastLinIdx, lastOffset
      local count = 0
      for idx in rFace:iter(ordering) do
	 if count == 0 then
	    -- reset things if first time in this loop
	    count = 1
	    currBlockSize = 1
	    lastLinIdx = indexer(idx)
	    lastOffset = 0
	 else
	    local linIdx = indexer(idx)
	    if linIdx-lastLinIdx == 1 then -- indices are contiguous
	       currBlockSize = currBlockSize+1
	    else
	       -- store size and index
	       blockSize[count] = currBlockSize*numComponent
	       blockOffset[count] = lastOffset

	       -- prep for next round
	       currBlockSize = 1
	       lastOffset = (linIdx-1)*numComponent
	       count = count+1
	    end
	    lastLinIdx = linIdx -- for next round
	 end
	 blockSize[count] = currBlockSize*numComponent
	 blockOffset[count] = lastOffset
      end
   end
   return blockSize, blockOffset
end

-- Creates an MPI_Datatype object from a range object in a specified
-- direction and ordering. 'ordering' must be one of Range.rowMajor or
-- Range.colMajor.
--
-- 'nlayer' is number of layers in range to include in datatype and
-- 'numComponent' is the number of components in field (usually 1)
function _M.createDataTypeFromRange(dir, range, nlayer, numComponent, ordering, oldtype)
   local blockSize, blockOffset = _M.createBlockInfoFromRange(dir, range, nlayer, numComponent, ordering)
   return _M.Type_commit(
      _M.createDataTypeFromBlockSizeAndOffset(blockSize, blockOffset, oldtype)
   )
end

return _M
