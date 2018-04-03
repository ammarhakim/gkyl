// Gkyl ------------------------------------------------------------------------
//
// Functions for use in MPI LuaJIT binding
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <iostream>

#include <mpi.h>
#include <GkMpiFuncs.h>

// Sizeof operators for various objects
GET_MPI_OBJ_SIZE(MPI_Status);

// Pre-defined objects and constants
GET_MPI_OBJECT(Comm, MPI_COMM_WORLD);
GET_MPI_OBJECT(Comm, MPI_COMM_NULL);
GET_MPI_OBJECT(Comm, MPI_COMM_SELF)
GET_MPI_OBJECT(Request, MPI_REQUEST_NULL);
GET_MPI_OBJECT_PTR(Status, MPI_STATUS_IGNORE);
GET_MPI_OBJECT(Info, MPI_INFO_NULL);

// Datatypes
GET_MPI_OBJECT(Datatype, MPI_CHAR);
GET_MPI_OBJECT(Datatype, MPI_BYTE);
GET_MPI_OBJECT(Datatype, MPI_SHORT);
GET_MPI_OBJECT(Datatype, MPI_INT);
GET_MPI_OBJECT(Datatype, MPI_LONG);
GET_MPI_OBJECT(Datatype, MPI_FLOAT);
GET_MPI_OBJECT(Datatype, MPI_DOUBLE);
GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_CHAR);
GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_SHORT);
GET_MPI_OBJECT(Datatype, MPI_UNSIGNED);
GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_LONG);
GET_MPI_OBJECT(Datatype, MPI_LONG_DOUBLE);
GET_MPI_OBJECT(Datatype, MPI_LONG_LONG_INT);
GET_MPI_OBJECT(Datatype, MPI_FLOAT_INT);
GET_MPI_OBJECT(Datatype, MPI_LONG_INT);
GET_MPI_OBJECT(Datatype, MPI_DOUBLE_INT);
GET_MPI_OBJECT(Datatype, MPI_SHORT_INT);
GET_MPI_OBJECT(Datatype, MPI_2INT);
GET_MPI_OBJECT(Datatype, MPI_LONG_DOUBLE_INT);
GET_MPI_OBJECT(Datatype, MPI_PACKED);

// Operators
GET_MPI_OBJECT(Op, MPI_MAX);
GET_MPI_OBJECT(Op, MPI_MIN);
GET_MPI_OBJECT(Op, MPI_SUM);
GET_MPI_OBJECT(Op, MPI_PROD);
GET_MPI_OBJECT(Op, MPI_LAND);
GET_MPI_OBJECT(Op, MPI_BAND);
GET_MPI_OBJECT(Op, MPI_LOR);
GET_MPI_OBJECT(Op, MPI_BOR);
GET_MPI_OBJECT(Op, MPI_LXOR);
GET_MPI_OBJECT(Op, MPI_BXOR);
GET_MPI_OBJECT(Op, MPI_MINLOC);
GET_MPI_OBJECT(Op, MPI_MAXLOC);

// Error codes
GET_INT_OBJECT(MPI_SUCCESS);

// Constants
GET_INT_OBJECT(MPI_COMM_TYPE_SHARED);
GET_INT_OBJECT(MPI_UNDEFINED);

// utility functions
void
GkMPI_fillStatus(const MPI_Status* inStatus, int *outStatus) {
  outStatus[0] = inStatus->MPI_SOURCE;
  outStatus[1] = inStatus->MPI_TAG;
  outStatus[2] = inStatus->MPI_ERROR;
}
