// Gkyl ------------------------------------------------------------------------
//
// Macros for internal use in MPI LuaJIT bindings
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_MPI_MACROS_H
#define GK_MPI_MACROS_H

#define GET_MPI_OBJECT(type, value) MPI_##type get_##value() { return value; }
#define DECL_GET_MPI_OBJECT(type, value) MPI_##type get_##value()

#define LUA_SET_MPI_OBJECT(type, value) _M.##value = ffi.C.get_##value()

#endif // GK_MPI_MACROS_H
