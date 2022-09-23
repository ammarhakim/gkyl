// Gkyl ------------------------------------------------------------------------
//
// Functions for use in MPI LuaJIT binding
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GKYL_MPI_FUNCS_H
#define GKYL_MPI_FUNCS_H

#include <GkylMpiMacros.h>

extern "C" {
    // Sizes of various objects
    DECL_GET_MPI_OBJ_SIZE(MPI_Status);
    DECL_GET_MPI_OBJ_SIZE(MPI_Request);
    DECL_GET_MPI_OBJ_PTR_SIZE(MPI_Status);
    DECL_GET_MPI_OBJ_PTR_SIZE(MPI_Request);
    
    // Pre-defined objects and constants
    DECL_GET_MPI_OBJECT(Comm, MPI_COMM_WORLD);
    DECL_GET_MPI_OBJECT(Comm, MPI_COMM_NULL);
    DECL_GET_MPI_OBJECT(Comm, MPI_COMM_SELF);
    DECL_GET_MPI_OBJECT(Info, MPI_INFO_NULL);
    
    DECL_GET_MPI_OBJECT(Request, MPI_REQUEST_NULL);
    DECL_GET_MPI_OBJECT_PTR(Status, MPI_STATUS_IGNORE);
    DECL_INT_OBJECT(MPI_PROC_NULL);

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

    // error codes
    DECL_INT_OBJECT(MPI_SUCCESS);

    // Constants
    DECL_INT_OBJECT(MPI_COMM_TYPE_SHARED);
    DECL_INT_OBJECT(MPI_UNDEFINED);
    DECL_INT_OBJECT(MPI_ORDER_C);
    DECL_INT_OBJECT(MPI_ORDER_FORTRAN);

    // Some utility functions to allow accessing non-opaque MPI types

    /* Fill the outStatus with MPI_Status public data */
    void GkMPI_fillStatus(const MPI_Status* inStatus, int *outStatus);
    void GkMPI_fillStatusArray(int count, const MPI_Status* inStatus, int *outStatus);
}

#endif // GK_MPI_FUNCS_H
