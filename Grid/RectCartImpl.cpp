#include <RectCartImpl.h>

void cellCenter(int ndim, double* lower, int* currIdx, double* dx, double* xc) {
  for(unsigned int d=0; d<ndim; d++) {
    xc[d] = lower[d] + (currIdx[d]-0.5)*dx[d];
  }
}
