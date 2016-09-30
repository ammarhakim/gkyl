// Gkyl ------------------------------------------------------------------------
//
// Functions for use in MPI LuaJIT binding
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_MPI_FUNCS_H
#define GK_MPI_FUNCS_H

#include <GkMpiMacros.h>

extern "C" {
    // Communicators    
    DECL_GET_MPI_OBJECT(Comm, MPI_COMM_WORLD);

    // Datatypes
    DECL_GET_MPI_OBJECT(Datatype, MPI_CHAR);
    DECL_GET_MPI_OBJECT(Datatype, MPI_BYTE);
    DECL_GET_MPI_OBJECT(Datatype, MPI_SHORT);
    DECL_GET_MPI_OBJECT(Datatype, MPI_INT);
    DECL_GET_MPI_OBJECT(Datatype, MPI_LONG);
    DECL_GET_MPI_OBJECT(Datatype, MPI_FLOAT);
    DECL_GET_MPI_OBJECT(Datatype, MPI_DOUBLE);
    DECL_GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_CHAR);
    DECL_GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_SHORT);
    DECL_GET_MPI_OBJECT(Datatype, MPI_UNSIGNED);
    DECL_GET_MPI_OBJECT(Datatype, MPI_UNSIGNED_LONG);
    DECL_GET_MPI_OBJECT(Datatype, MPI_LONG_DOUBLE);
    DECL_GET_MPI_OBJECT(Datatype, MPI_LONG_LONG_INT);
    DECL_GET_MPI_OBJECT(Datatype, MPI_FLOAT_INT);
    DECL_GET_MPI_OBJECT(Datatype, MPI_LONG_INT);
    DECL_GET_MPI_OBJECT(Datatype, MPI_DOUBLE_INT);
    DECL_GET_MPI_OBJECT(Datatype, MPI_SHORT_INT);
    DECL_GET_MPI_OBJECT(Datatype, MPI_2INT);
    DECL_GET_MPI_OBJECT(Datatype, MPI_LONG_DOUBLE_INT);
    DECL_GET_MPI_OBJECT(Datatype, MPI_PACKED);

    // Ops
    DECL_GET_MPI_OBJECT(Op, MPI_MAX);
    DECL_GET_MPI_OBJECT(Op, MPI_MIN);
    DECL_GET_MPI_OBJECT(Op, MPI_SUM);
    DECL_GET_MPI_OBJECT(Op, MPI_PROD);
    DECL_GET_MPI_OBJECT(Op, MPI_LAND);
    DECL_GET_MPI_OBJECT(Op, MPI_BAND);
    DECL_GET_MPI_OBJECT(Op, MPI_LOR);
    DECL_GET_MPI_OBJECT(Op, MPI_BOR);
    DECL_GET_MPI_OBJECT(Op, MPI_LXOR);
    DECL_GET_MPI_OBJECT(Op, MPI_BXOR);
    DECL_GET_MPI_OBJECT(Op, MPI_MINLOC);
    DECL_GET_MPI_OBJECT(Op, MPI_MAXLOC);
}

#endif // GK_MPI_FUNCS_H
