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
    DECL_GET_MPI_OBJECT(Comm, MPI_COMM_WORLD);
}

#endif // GK_MPI_FUNCS_H
