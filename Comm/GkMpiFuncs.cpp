// Gkyl ------------------------------------------------------------------------
//
// Functions for use in MPI LuaJIT binding
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <mpi.h>
#include <GkMpiFuncs.h>

// Communicators
GET_MPI_OBJECT(Comm, MPI_COMM_WORLD);

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
