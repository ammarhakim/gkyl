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
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

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

  // Win & SHM calls
  int MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm);
  int MPI_Win_allocate_shared (MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
  int MPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr);
  int MPI_Win_lock_all(int assert, MPI_Win win);
  int MPI_Win_unlock_all(MPI_Win win);
  int MPI_Win_free(MPI_Win *win);

  // point-to-point communication
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

-- MPI_Comm_rank
function _M.Comm_rank(comm)
   local r = int_1()
   local err = ffiC.MPI_Comm_rank(getObj(comm, "MPI_Comm[1]"), r)
   return r[0]
end
-- MPI_Comm_size
function _M.Comm_size(comm)
   local r = int_1()
   local err = ffiC.MPI_Comm_size(getObj(comm, "MPI_Comm[1]"), r)
   return r[0]
end
-- MPI_Comm_dup
function _M.Comm_dup(comm)
   local c = new_MPI_Comm()
   local err = ffiC.MPI_Comm_dup(getObj(comm, "MPI_Comm[1]"), c)
   return c
end
-- MPI_Comm_split
function _M.Comm_split(comm, color, key)
   local c = new_MPI_Comm()
   local err = ffiC.MPI_Comm_split(getObj(comm, "MPI_Comm[1]"), color, key, c)
   return c
end
-- MPI_Comm_split_type
function _M.Comm_split_type(comm, split_type, key, info)
   local c = new_MPI_Comm()
   local err = ffiC.MPI_Comm_split_type(getObj(comm, "MPI_Comm[1]"), split_type, key, info, c)
   return c
end
-- MPI_Win_allocate_shared
function _M.Win_allocate_shared(sz, disp_unit, info, comm)
   local w = new_MPI_Win()
   local baseptr = voidp()
   local err = ffiC.MPI_Win_allocate_shared(sz, disp_unit, info, getObj(comm, "MPI_Comm[1]"), baseptr, w)
   return baseptr[0], w
end
-- MPI_Win_shared_query
function _M.Win_shared_query(win, rank)
   local sz, du = uint_1(), int_1()
   local baseptr = voidp()
   local err = ffiC.MPI_Win_shared_query(getObj(win, "MPI_Win[1]"), rank, sz, du, baseptr)
   return sz[0], du[0], baseptr[0]
end
-- MPI_Win_lock_all
function _M.Win_lock_all(assert, win)
   local err = ffiC.MPI_Win_lock_all(assert, getObj(win, "MPI_Win[1]"))
end
-- MPI_Win_unlock_all
function _M.Win_unlock_all(win)
   local err = ffiC.MPI_Win_unlock_all(getObj(win, "MPI_Win[1]"))
end
-- MPI_Win_free
function _M.Win_free(win)
   local err = ffiC.MPI_Win_free(win)
end

-- MPI_Get_count
function _M.Get_count(status, datatype)
   local r = int_1()
   local err = ffiC.MPI_Get_count(status.mpiStatus, datatype, r)
   return r[0]
end
-- MPI_Allreduce
function _M.Allreduce(sendbuf, recvbuf, count, datatype, op, comm)
   local err = ffiC.MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Bcast
function _M.Bcast(buffer, count, datatype, root, comm)
   local err = ffiC.MPI_Bcast(buffer, count, datatype, root, getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Send
function _M.Send(buf, count, datatype, dest, tag, comm)
   local err = ffiC.MPI_Send(buf, count, datatype, dest, tag, getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Recv
function _M.Recv(buf, count, datatype, source, tag, comm, status)
   local st = status and status.mpiStatus or _M.STATUS_IGNORE
   local err = ffiC.MPI_Recv(buf, count, datatype, source, tag, getObj(comm, "MPI_Comm[1]"), st)
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
   local err = ffiC.MPI_Irecv(
      buf, count, datatype, source, tag, getObj(comm, "MPI_Comm[1]"), req)
   return req
end
-- MPI_Wait
function _M.Wait(request, status)
   local st = status and status.mpiStatus or _M.STATUS_IGNORE
   local err = ffiC.MPI_Wait(request, st)
   -- store MPI_Status
   if status ~= nil then
      local gks = new("int[3]")
      ffiC.GkMPI_fillStatus(st, gks)
      status.SOURCE, status.TAG, status.ERROR = gks[0], gks[1], gks[2]
   end
end

-- MPI_Barrier
function _M.Barrier(comm)
   local err = ffiC.MPI_Barrier(getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Abort
function _M.Abort(comm, errCode)
   local err = ffiC.MPI_Abort(getObj(comm, "MPI_Comm[1]"), errCode)
end

-- MPI_Comm_group
function _M.Comm_group(comm)
   local grp = new_MPI_Group()
   local err = ffiC.MPI_Comm_group(getObj(comm, "MPI_Comm[1]"), grp)
   return grp
end
function _M.Group_rank(group)
   local r = int_1()
   local err = ffiC.MPI_Group_rank(getObj(group, "MPI_Group[1]"), r)
   return r[0]
end
function _M.Group_size(group)
   local r = int_1()
   local err = ffiC.MPI_Group_size(getObj(group, "MPI_Group[1]"), r)
   return r[0]
end
function _M.Group_incl(group, n, ranks)
   local newgroup = new_MPI_Group()
   local err = ffiC.MPI_Group_incl(getObj(group, "MPI_Group[1]"), n, ranks, newgroup)
   return newgroup
end
function _M.Comm_create(comm, group)
   local c = new_MPI_Comm()
   local err = ffiC.MPI_Comm_create(
      getObj(comm, "MPI_Comm[1]"), getObj(group, "MPI_Group[1]"), c)
   return c
end
function _M.Group_translate_ranks(group1, ranks1, group2)
   local n = #ranks1
   local ranks2 = Lin.IntVec(n)
   local err = ffiC.MPI_Group_translate_ranks(
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

return _M
