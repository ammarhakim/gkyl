/* -*- c++ -*- */

#include <cstdio>

extern "C" 
{

    typedef struct { int32_t ndim; int32_t lower[6]; int32_t upper[6]; } GkylRange;
    typedef struct {
        int32_t ndim;
        int32_t cells[6];
        double lower[6], upper[6];
        double vol, dx[6];
    } GkylRectCart;

    void unit_sumArray(int numBlocks, int numThreads, int n, double a, double *x, double *y);
    void unit_sayHello();
    void unit_showRange(GkylRange *range);
    void unit_showGrid(GkylRectCart *grid);
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

__global__ void ker_unit_showRange(GkylRange *range)
{
  printf("Range ndim: %d\n", range->ndim);
  for (unsigned i=0; i<range->ndim; ++i)
    printf(" %d, %d\n", range->lower[i], range->upper[i]);
}

__global__ void ker_unit_showGrid(GkylRectCart *grid)
{
  printf("Grid ndim: %d\n", grid->ndim);
  for (unsigned i=0; i<grid->ndim; ++i)
    printf(" %g, %g\n", grid->lower[i], grid->upper[i]);
}

void unit_sumArray(int numBlocks, int numThreads, int n, double a, double *x, double *y)
{
   ker_unit_sumArray<<<numBlocks, numThreads>>>(n, a, x, y);
}

void unit_sayHello()
{
  ker_unit_sayHello<<<1, 1>>>();
}

void unit_showRange(GkylRange *devRange)
{
  ker_unit_showRange<<<1, 1>>>(devRange);
}

void unit_showGrid(GkylRectCart *devGrid)
{
  ker_unit_showGrid<<<1, 1>>>(devGrid);
}
