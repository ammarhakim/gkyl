// Gkyl ------------------------------------------------------------------------
//
// Functions for use in NCCL LuaJIT binding.
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifdef GKYL_HAVE_NCCL

#include <gkyl_nccl_funcs.h>

// Initialize a ncclConfig_t.
void
gkyl_NCCL_CONFIG_INITIALIZER(ncclConfig_t *nc)
{
  nc[0] = NCCL_CONFIG_INITIALIZER;
}

#endif
