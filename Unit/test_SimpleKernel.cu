/* -*- c++ -*- */

#include <cstdio>
#include <RectCartDeviceImpl.h>
#include <RangeDeviceImpl.h>

extern "C" 
{
    void unit_sumArray(int numBlocks, int numThreads, int n, double a, double *x, double *y);
    void unit_sayHello();
    void unit_showRange(Range_t *range);
    void unit_showGrid(RectCart_t *grid);
    void unit_getCellCenter(RectCart_t *grid, int *idx, double *xc);
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

__global__ void ker_unit_showRange(Range_t *range)
{
  printf("Range ndim: %d\n", range->ndim);
  for (unsigned i=0; i<range->ndim; ++i)
    printf(" %d, %d\n", range->lower[i], range->upper[i]);

  Gkyl::Indexer<1> idxr1;
  Gkyl::Indexer<2> idxr2;
  int lin1 = idxr1.index(range, 0);
  int lin2 = idxr2.index(range, 0, 1);
  printf("  --- Calls to indexers %d, %d\n", lin1, lin2);
}

__global__ void ker_unit_showGrid(RectCart_t *grid)
{
  printf("Grid ndim: %d\n", grid->ndim);
  for (unsigned i=0; i<grid->ndim; ++i)
    printf(" %g, %g\n", grid->lower[i], grid->upper[i]);
}

__global__ void ker_unit_getCellCenter(RectCart_t *grid, int *idx, double *xc)
{
  grid->cellCenter(idx, xc);
}

void unit_sumArray(int numBlocks, int numThreads, int n, double a, double *x, double *y)
{
   ker_unit_sumArray<<<numBlocks, numThreads>>>(n, a, x, y);
}

void unit_sayHello()
{
  ker_unit_sayHello<<<1, 1>>>();
}

void unit_showRange(Range_t *devRange)
{
  ker_unit_showRange<<<1, 1>>>(devRange);
}

void unit_showGrid(RectCart_t *devGrid)
{
  ker_unit_showGrid<<<1, 1>>>(devGrid);
}

void unit_getCellCenter(RectCart_t *devGrid, int *idx, double *xc)
{
  ker_unit_getCellCenter<<<1, 1>>>(devGrid, idx, xc);
}
