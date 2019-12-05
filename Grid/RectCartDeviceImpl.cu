#include "RectCartDeviceImpl.h"

__device__ void cellCenter(RectCart_t* grid, int* idx, double* xc)
{
  for(unsigned int d=0; d<grid->ndim; d++) {
    xc[d] = grid->lower[d] + (idx[d]-0.5)*grid->dx[d];
  }
}

