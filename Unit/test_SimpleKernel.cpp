/* -*- c++ -*- */

#include <iostream>
#include <cstdio>
#include <GkylBasisTypes.h>
#include <GkylCartField.h>
#include <GkylRange.h>
#include <GkylRectCart.h>

extern "C" 
{
    void unit_sumArray(int numBlocks, int numThreads, int n, double a, double *x, double *y);
    void unit_sayHello();
    void unit_showRange(GkylRange_t *range);
    void unit_showGrid(GkylRectCart_t *grid);
    void unit_getCellCenter(GkylRectCart_t *grid, int *idx, double *xc);

    void unit_showFieldRange(GkylCartField_t *f, double *g);
    void unit_showFieldGrid(GkylCartField_t *f);

    void unit_test_BasisTypes_1xp1_ser();
    void unit_test_BasisTypes_2xp2_ser();
    void unit_test_BasisTypes_5xp2_ser();

    // Type of function that "adds" two numbers
    typedef double (*sumFunc_t)(void *self, double a, double b);
    // Base equation object
    struct SimpleEquation {
        void *child { }; // pointer to child object
        sumFunc_t sumFunc { }; // provided by "children"
        
        __host__ __device__ double sum(double a, double b) {
          return sumFunc(child, a, b);
        }
    };
    // Euler equation object
    struct SimpleEulerEquation {
        double gamma { 1.4 };
    };

    double unit_test_SimpleEquation(SimpleEquation *eqn, double* ab);
    sumFunc_t getEulerSumFuncOnDevice();
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

  Gkyl::Indexer<1> idxr(range);
  int lin1 = idxr.index(range->lower[0]);
  int lin2 = idxr.index(range->upper[0]);
  printf("  --- Calls to indexers %d, %d (%d)\n", lin1, lin2, lin2-lin1);

  int idx[1];
  for (int i=0; i<range->upper[0]-range->lower[0]+1; ++i) {
    idxr.invIndex(i, idx);
    printf("Lin index %d -> index %d\n", i, idx[0]);
  }

  // can also instantiate a C++ range object using the GkylRange_t range pointer (from Lua)
  Gkyl::Range range_cpp(range);
  printf("Range ndim: %d\n", range_cpp.ndim());
  for (unsigned i=0; i<range_cpp.ndim(); ++i)
    printf(" %d, %d\n", range_cpp.lower(i), range_cpp.upper(i));
  printf("Volume = %d\n", range_cpp.volume());
}

__global__ void ker_unit_showGrid(GkylRectCart_t *grid)
{
  printf("Grid ndim: %d\n", grid->ndim);
  for (unsigned i=0; i<grid->ndim; ++i)
    printf(" %g, %g\n", grid->lower[i], grid->upper[i]);
}

__global__ void ker_unit_getCellCenter(GkylRectCart_t *grid, int *idx, double *xc)
{
  grid->cellCenter(idx, xc);
}

__global__ void ker_unit_showFieldRange(GkylCartField_t *f, double *g)
{
  GkylRange_t *range = f->localRange;
  GkylRange_t *rangeExt = f->localExtRange;
  printf("Field ndim: %d\n", f->ndim);
  printf("localRange ndim: %d\n", range->ndim);
  for (unsigned i=0; i<range->ndim; ++i)
    printf("localRange, d=%d: %d, %d\n", i, range->lower[i], range->upper[i]);
  for (unsigned i=0; i<rangeExt->ndim; ++i)
    printf("localExtRange, d=%d: %d, %d\n", i, rangeExt->lower[i], rangeExt->upper[i]);

  Gkyl::GenIndexer fIdxr = f->genIndexer();
  Gkyl::GenIndexer idxr = Gkyl::GenIndexer(range);
  int lin1 = fIdxr.index(rangeExt->lower);
  int lin2 = fIdxr.index(rangeExt->upper);
  printf("  --- Calls to indexers %d, %d (%d)\n", lin1, lin2, lin2-lin1);

  int idx[1];
  for (int i=0; i<range->volume(); ++i) {
    idxr.invIndex(i, idx);
    printf("Lin index %d -> index %d, f = %f, g = %f\n", i, idx[0], f->getDataPtrAt(i)[0], g[i*3]);
  }
}

__global__ void ker_unit_showFieldGrid(GkylCartField_t *f)
{
  GkylRectCart_t *grid = f->grid;
  printf("Grid ndim: %d\n", grid->ndim);
  for (unsigned i=0; i<grid->ndim; ++i)
    printf(" %g, %g\n", grid->lower[i], grid->upper[i]);
}

// example of templated kernel function
template <int NDIM, int POLYORDER, int BASISTYPE>
__global__ void ker_test_BasisTypes()
{
  int numBasis = Gkyl::BasisCount<NDIM, POLYORDER, BASISTYPE>::numBasis();
  printf("NDIM = %d, polyOrder = %d : numBasis = %d\n", NDIM, POLYORDER, numBasis);
}


// "add" function specific to Euler equation
__device__ double simpleEulerSumFunc(SimpleEulerEquation *self, double a, double b) {
  printf("simpleEulerSumFunc: %g*(%g + %g)\n", self->gamma, a, b);
  return self->gamma*(a+b);
}

sumFunc_t getEulerSumFuncOnDevice() {
  sumFunc_t eulerSumFunc_p;
  auto err = cudaMemcpyFromSymbol(&eulerSumFunc_p, simpleEulerSumFunc, sizeof(sumFunc_t));
  std::cout << "getEulerSumFuncOnDevice " << err << std::endl;
  return eulerSumFunc_p;
}

__global__ void ker_updateEquation(SimpleEquation *eqn, double *ab, double *out) {
  double v = eqn->sum(ab[0], ab[1]);
  (*out) = v;
}

// Functions callable from LuaJIT
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

void unit_showGrid(GkylRectCart_t *devGrid)
{
  ker_unit_showGrid<<<1, 1>>>(devGrid);
}

void unit_getCellCenter(GkylRectCart_t *devGrid, int *idx, double *xc)
{
  ker_unit_getCellCenter<<<1, 1>>>(devGrid, idx, xc);
}

void unit_showFieldRange(GkylCartField_t *f, double *g)
{
  ker_unit_showFieldRange<<<1, 1>>>(f, g);
}

void unit_showFieldGrid(GkylCartField_t *f)
{
  ker_unit_showFieldGrid<<<1, 1>>>(f);
}

void unit_test_BasisTypes_1xp1_ser()
{
  ker_test_BasisTypes<1,1,Gkyl::G_SERENDIPITY_C><<<1,1>>>();
}

void unit_test_BasisTypes_2xp2_ser()
{
  ker_test_BasisTypes<2,2,Gkyl::G_SERENDIPITY_C><<<1,1>>>();
}

void unit_test_BasisTypes_5xp2_ser()
{
  ker_test_BasisTypes<5,2,Gkyl::G_SERENDIPITY_C><<<1,1>>>();
}

double unit_test_SimpleEquation(SimpleEquation *eqn, double* ab) 
{
  double *out;
  auto err = cudaMalloc(&out, sizeof(double));
  ker_updateEquation<<<1,1>>>(eqn, ab, out);

  double retVal = 0;
  cudaMemcpy(&retVal, out, sizeof(double), cudaMemcpyDeviceToHost);
  return retVal;
}
