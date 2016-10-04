-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for MPI 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local _M = {}

ffi.cdef [[
  // Opaque types
  typedef struct MPI_Comm_type *MPI_Comm;
  typedef struct MPI_Datatype_type *MPI_Datatype;        
  typedef struct MPI_Op_type *MPI_Op;
  typedef struct MPI_Status_type *MPI_Status;
  typedef struct MPI_Group_type *MPI_Group;

  // size of various objects
  int sizeof_MPI_Status();
  int sizeof_MPI_Group();
  int sizeof_MPI_Comm();

  // Pre-defined objects and constants
  MPI_Comm get_MPI_COMM_WORLD();
  MPI_Comm get_MPI_COMM_NULL();
  MPI_Comm get_MPI_COMM_SELF();
  MPI_Status *get_MPI_STATUS_IGNORE();

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

  // Communicators
  int MPI_Comm_rank(MPI_Comm comm, int *rank);
  int MPI_Comm_size(MPI_Comm comm, int *size);
  int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);  

  // point-to-point communication
  int MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count);
  int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
		    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
  int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
	       MPI_Comm comm);
  int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
	       MPI_Comm comm, MPI_Status *status);

  // Groups
  int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
  int MPI_Group_rank(MPI_Group group, int *rank);
  int MPI_Group_size(MPI_Group group, int *size);
  int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup);
  int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
  
  // Global operators
  int MPI_Barrier(MPI_Comm comm);
  int MPI_Abort(MPI_Comm comm, int errorcode);
  
  // Gkyl utility functions
  void GkMPI_fillStatus(const MPI_Status* inStatus, int *outStatus);
]]

-- Predefined objects and constants
_M.COMM_WORLD = ffi.C.get_MPI_COMM_WORLD()
_M.COMM_NULL = ffi.C.get_MPI_COMM_NULL()
_M.COMM_SELF = ffi.C.get_MPI_COMM_SELF()
_M.STATUS_IGNORE = ffi.C.get_MPI_STATUS_IGNORE()

-- Object sizes
_M.SIZEOF_STATUS = ffi.C.sizeof_MPI_Status()
_M.SIZEOF_GROUP = ffi.C.sizeof_MPI_Group()
_M.SIZEOF_COMM = ffi.C.sizeof_MPI_Comm()

-- Datatypes
_M.CHAR = ffi.C.get_MPI_CHAR()
_M.BYTE = ffi.C.get_MPI_BYTE()
_M.SHORT = ffi.C.get_MPI_SHORT()
_M.INT = ffi.C.get_MPI_INT()
_M.LONG = ffi.C.get_MPI_LONG()
_M.FLOAT = ffi.C.get_MPI_FLOAT()
_M.DOUBLE = ffi.C.get_MPI_DOUBLE()
_M.UNSIGNED_CHAR = ffi.C.get_MPI_UNSIGNED_CHAR()
_M.UNSIGNED_SHORT = ffi.C.get_MPI_UNSIGNED_SHORT()
_M.UNSIGNED = ffi.C.get_MPI_UNSIGNED()
_M.UNSIGNED_LONG = ffi.C.get_MPI_UNSIGNED_LONG()
_M.LONG_DOUBLE = ffi.C.get_MPI_LONG_DOUBLE()
_M.LONG_LONG_INT = ffi.C.get_MPI_LONG_LONG_INT()
_M.FLOAT_INT = ffi.C.get_MPI_FLOAT_INT()
_M.LONG_INT = ffi.C.get_MPI_LONG_INT()
_M.DOUBLE_INT = ffi.C.get_MPI_DOUBLE_INT()
_M.SHORT_INT = ffi.C.get_MPI_SHORT_INT()
_M.TWOINT = ffi.C.get_MPI_2INT()
_M.LONG_DOUBLE_INT = ffi.C.get_MPI_LONG_DOUBLE_INT()
_M.PACKED = ffi.C.get_MPI_PACKED()

-- Operators
_M.MAX = ffi.C.get_MPI_MAX();
_M.MIN = ffi.C.get_MPI_MIN();
_M.SUM = ffi.C.get_MPI_SUM();
_M.PROD = ffi.C.get_MPI_PROD();
_M.LAND = ffi.C.get_MPI_LAND();
_M.BAND = ffi.C.get_MPI_BAND();
_M.LOR = ffi.C.get_MPI_LOR();
_M.BOR = ffi.C.get_MPI_BOR();
_M.LXOR = ffi.C.get_MPI_LXOR();
_M.BXOR = ffi.C.get_MPI_BXOR();
_M.MINLOC = ffi.C.get_MPI_MINLOC();
_M.MAXLOC = ffi.C.get_MPI_MAXLOC();

-- some types for use in MPI functions
local int_1 = typeof("int[1]")

-- ctors for various MPI objects: these functions work by first
-- allocating space for the object and then casting the memory to the
-- appropriate type. One needs to do this as MPI defines these objects
-- as "implementation defined" and hence we can't assume anything
-- about their internal representation. (Even for MPI_Status which has
-- public members, an implementation is free to put other elements in
-- the struct and hence we can't use the size of the public members to
-- allocate memory).
local function new_MPI_Status()
   return ffi.cast("MPI_Status*", new("uint8_t[?]", _M.SIZEOF_STATUS))
end
local function new_MPI_Comm()
   return ffi.cast("MPI_Comm*", new("uint8_t[?]", _M.SIZEOF_COMM))
end
local function new_MPI_Group()
   return ffi.cast("MPI_Group*", new("uint8_t[?]", _M.SIZEOF_GROUP))
end

-- de-reference if object is a pointer
local function getObj(obj, ptyp)
   return ffi.istype(typeof(ptyp), obj) and obj[0] or obj
end

-- MPI_Status object
_M.Status = {}
function _M.Status:new()
   local self = setmetatable({}, Status)
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
   local err = ffi.C.MPI_Comm_rank(getObj(comm, "MPI_Comm*"), r)
   return r[0]
end
-- MPI_Comm_size
function _M.Comm_size(comm)
   local r = int_1()
   local err = ffi.C.MPI_Comm_size(getObj(comm, "MPI_Comm*"), r)
   return r[0]
end
-- MPI_Comm_dup
function _M.Comm_dup(comm)
   local c = new_MPI_Comm()
   local err = ffi.C.MPI_Comm_dup(getObj(comm, "MPI_Comm*"), c)
   return c
end

-- MPI_Allreduce
function _M.Allreduce(sendbuf, recvbuf, count, datatype, op, comm)
   local err = ffi.C.MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, getObj(comm, "MPI_Comm*"))
end
-- MPI_Get_count
function _M.Get_count(status, datatype)
   local r = int_1()
   local err = ffi.C.MPI_Get_count(status.mpiStatus, datatype, r)
   return r[0]
end
-- MPI_Send
function _M.Send(buf, count, datatype, dest, tag, comm)
   local err = ffi.C.MPI_Send(buf, count, datatype, dest, tag, getObj(comm, "MPI_Comm*"))
end
-- MPI_Recv
function _M.Recv(buf, count, datatype, source, tag, comm, status)
   local st = status and status.mpiStatus or _M.STATUS_IGNORE
   local err = ffi.C.MPI_Recv(buf, count, datatype, source, tag, getObj(comm, "MPI_Comm*"), st)
   -- store MPI_Status
   if status ~= nil then
      local gks = new("int[3]")
      ffi.C.GkMPI_fillStatus(st, gks)
      status.SOURCE, status.TAG, status.ERROR = gks[0], gks[1], gks[2]
   end
end

-- MPI_Barrier
function _M.Barrier(comm)
   local err = ffi.C.MPI_Barrier(getObj(comm, "MPI_Comm*"))
end
-- MPI_Abort
function _M.Abort(comm, err)
   local err = ffi.C.MPI_Abort(getObj(comm, "MPI_Comm*"), err)
end

-- MPI_Comm_group
function _M.Comm_group(comm)
   local grp = new_MPI_Group()
   local err = ffi.C.MPI_Comm_group(getObj(comm, "MPI_Comm*"), grp)
   return grp
end
function _M.Group_rank(group)
   local r = int_1()
   local err = ffi.C.MPI_Group_rank(getObj(group, "MPI_Group*"), r)
   return r[0]
end
function _M.Group_size(group)
   local r = int_1()
   local err = ffi.C.MPI_Group_size(getObj(group, "MPI_Group*"), r)
   return r[0]
end
function _M.Group_incl(group, n, ranks)
   local newgroup = new_MPI_Group()
   local err = ffi.C.MPI_Group_incl(getObj(group, "MPI_Group*"), n, ranks, newgroup)
   return newgroup
end
function _M.Comm_create(comm, group)
   local c = new_MPI_Comm()
   local err = ffi.C.MPI_Comm_create(
      getObj(comm, "MPI_Comm*"), getObj(group, "MPI_Group*"), c)
   return c
end

-- Convenience functions

-- Collect 'ranks' from 'comm' and create a new communicator with just
-- those ranks. The 'ranks' parameter must be of a Linalg.IntVec object
function _M.Split_comm(comm, ranks)
   local grp = _M.Comm_group(comm)
   local newGrp = _M.Group_incl(grp, #ranks, ranks:data())
   return _M.Comm_create(comm, newGrp)
end

return _M
