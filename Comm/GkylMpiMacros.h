// Gkyl ------------------------------------------------------------------------
//
// Macros for internal use in MPI LuaJIT bindings
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GKYL_MPI_MACROS_H
#define GKYL_MPI_MACROS_H

// Macros to declare/define functions to get MPI objects
#define GET_MPI_OBJECT(type, value) MPI_##type get_##value() { return value; }
#define DECL_GET_MPI_OBJECT(type, value) MPI_##type get_##value()
#define LUA_SET_MPI_OBJECT(type, value) _M.value = ffi.C.get_##value();

#define GET_INT_OBJECT(value) int get_##value() { return value; }
#define DECL_INT_OBJECT(value) int get_##value()

// Macros to declare/define functions to get MPI object pointers
#define GET_MPI_OBJECT_PTR(type, value) MPI_##type *getPtr_##value() { return value; }
#define DECL_GET_MPI_OBJECT_PTR(type, value) MPI_##type *getPtr_##value()

// Macros to declare/define functions to get sizes of various MPI objects
#define GET_MPI_OBJ_SIZE(type) int sizeof_##type() { return sizeof(type); }
#define GET_MPI_OBJ_PTR_SIZE(type) int sizeof_ptr_##type() { return sizeof(type *); }
#define DECL_GET_MPI_OBJ_SIZE(type) int sizeof_##type()
#define DECL_GET_MPI_OBJ_PTR_SIZE(type) int sizeof_ptr_##type()

#endif // GKYL_MPI_MACROS_H
