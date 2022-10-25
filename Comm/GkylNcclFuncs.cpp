// Gkyl ------------------------------------------------------------------------
//
// Functions for use in NCCL LuaJIT binding.
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <GkylNcclFuncs.h>

// Initialize a ncclConfig_t.
void
gkyl_NCCL_CONFIG_INITIALIZER(ncclConfig_t *nc) {

  nc[0] = NCCL_CONFIG_INITIALIZER;
}
