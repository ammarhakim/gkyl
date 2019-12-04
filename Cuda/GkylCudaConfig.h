// Gkyl ------------------------------------------------------------------------
//
// Configuration header for CUDA (or lack of CUDA)
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_CUDA_CONFIG_H
#define GK_CUDA_CONFIG_H

#include <gkylconfig.h>

#ifndef HAVE_CUDA_H
# define __global__
# define __device__
# define __host__
#endif

#endif // GK_CUDA_CONFIG_H
