// Gkyl ------------------------------------------------------------------------
//
// Functions for use in MPI LuaJIT binding
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <mpi.h>
#include <GkMpiFuncs.h>

GET_MPI_OBJECT(Comm, MPI_COMM_WORLD);

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
