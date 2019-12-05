#include "RectCartDeviceImpl.h"

__device__ void RectCart_t::cellCenter(int* idx, double* xc)
{
  for(unsigned int d=0; d<this->ndim; d++) {
    xc[d] = this->lower[d] + (idx[d]-0.5)*this->dx[d];
  }
}

