/* -*- c++ -*- */

#include <cstdio>
#include "RectCartDeviceImpl.h"

extern "C" 
{
    void unit_sumArray(int numBlocks, int numThreads, int n, double a, double *x, double *y);
    void unit_sayHello();
    void unit_showRange(GkylRange_t *range);
    void unit_showGrid(RectCart_t *grid);
}

__global__ void ker_unit_sumArray(int n, double a, double *x, double *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i<n) y[i] = a*x[i]+y[i];
}

__global__ void ker_unit_sayHello()
{
  printf("Hello!\n");
}

__global__ void ker_unit_showRange(GkylRange_t *range)
{
  printf("Range ndim: %d\n", range->ndim);
  for (unsigned i=0; i<range->ndim; ++i)
    printf(" %d, %d\n", range->lower[i], range->upper[i]);
}

__global__ void ker_unit_showGrid(RectCart_t *grid)
{
  printf("Grid ndim: %d\n", grid->ndim);
  for (unsigned i=0; i<grid->ndim; ++i)
    printf(" %g, %g\n", grid->lower[i], grid->upper[i]);

  double xc[6];
  int idx[6];
  cellCenter(grid, idx, xc);
}

void unit_sumArray(int numBlocks, int numThreads, int n, double a, double *x, double *y)
{
   ker_unit_sumArray<<<numBlocks, numThreads>>>(n, a, x, y);
}

void unit_sayHello()
{
  ker_unit_sayHello<<<1, 1>>>();
}

void unit_showRange(GkylRange_t *devRange)
{
  ker_unit_showRange<<<1, 1>>>(devRange);
}

void unit_showGrid(RectCart_t *devGrid)
{
  ker_unit_showGrid<<<1, 1>>>(devGrid);
}
