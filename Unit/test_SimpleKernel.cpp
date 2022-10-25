/* -*- c++ -*- */

#include <iostream>
#include <cstdio>

extern "C" 
{
    void unit_sumArray(int numBlocks, int numThreads, int n, double a, double *x, double *y);
    void unit_sayHello();

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

    typedef struct {
        double x, y, z;
    } GkylTestParticle_t;

    void unit_showParticle(GkylTestParticle_t *ptcl);
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

__global__ void ker_unit_showParticle(GkylTestParticle_t *ptcl)
{
  printf("Particle: %g %g %g\n", ptcl->x, ptcl->y, ptcl->z);
}

// "add" function specific to Euler equation
__device__ double simpleEulerSumFunc(void *obj, double a, double b) {
  SimpleEulerEquation *self = (SimpleEulerEquation*) obj;
  return  self->gamma*(a+b);
}

__device__ sumFunc_t p_simpleEulerSumFunc = &simpleEulerSumFunc;

sumFunc_t getEulerSumFuncOnDevice() {
  sumFunc_t eulerSumFunc_p;
  auto err = cudaMemcpyFromSymbol(&eulerSumFunc_p, p_simpleEulerSumFunc, sizeof(sumFunc_t));
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

void unit_showParticle(GkylTestParticle_t *ptcl)
{
  ker_unit_showParticle<<<1, 1>>>(ptcl);
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

