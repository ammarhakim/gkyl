/* -*- c++ -*- */

#include <cstdio>
#include "RectCartDeviceImpl.h"

template <unsigned NDIM>
class SimpleIndexer
{
  public:
    __device__ inline int indexer(const GkylRange_t *range, const int indx[NDIM]) {
      return 0;
    }
};

template <>
class SimpleIndexer<1>
{
  public:
    __device__ inline int indexer(const GkylRange_t *range, int i1) {
      return 1;
    }
};

template <>
class SimpleIndexer<2>
{
  public:
    __device__ inline int indexer(const GkylRange_t *range, int i1, int i2) {
      return 2;
    }
};

extern "C" 
{
    void unit_sumArray(int numBlocks, int numThreads, int n, double a, double *x, double *y);
    void unit_sayHello();
    void unit_showRange(GkylRange_t *range);
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

__global__ void ker_unit_showRange(GkylRange_t *range)
{
  printf("Range ndim: %d\n", range->ndim);
  for (unsigned i=0; i<range->ndim; ++i)
    printf(" %d, %d\n", range->lower[i], range->upper[i]);

  SimpleIndexer<1> idxr1;
  SimpleIndexer<2> idxr2;
  int lin1 = idxr1.indexer(range, 0);
  int lin2 = idxr2.indexer(range, 0, 1);
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

void unit_showRange(GkylRange_t *devRange)
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
