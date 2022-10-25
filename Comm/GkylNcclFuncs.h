// Gkyl ------------------------------------------------------------------------
//
// Functions for use in NCCL LuaJIT binding.
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GKYL_NCCL_FUNCS_H
#define GKYL_NCCL_FUNCS_H

#include <mpi.h>
#include <nccl.h>

// Initialize a ncclConfig_t.
extern "C" {
  void gkyl_NCCL_CONFIG_INITIALIZER(ncclConfig_t *nc);
}

#endif // GKYL_NCCL_FUNCS_H
