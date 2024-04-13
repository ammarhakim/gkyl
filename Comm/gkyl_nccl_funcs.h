// Gkyl ------------------------------------------------------------------------
//
// Functions for use in NCCL LuaJIT binding.
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

#ifdef GKYL_HAVE_NCCL

#include <mpi.h>
#include <nccl.h>

// Initialize a ncclConfig_t.
void gkyl_NCCL_CONFIG_INITIALIZER(ncclConfig_t *nc);

#endif
