/* -*- c++ -*- */
// Gkyl ------------------------------------------------------------------------
//
// CUDA back-end for CartField
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include "CartFieldDeviceImpl.h"

__global__ void ker_gkylCartFieldAccumulate(unsigned s, unsigned nv, double fact, const double *inp, double *out) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x + s; n < nv; n += blockDim.x * gridDim.x)
     out[n] += fact*inp[n];
}

__global__  void ker_gkylCartFieldAssign(unsigned s, unsigned nv, double fact, const double *inp, double *out) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x + s; n < nv; n += blockDim.x * gridDim.x)
    out[n] = fact*inp[n];
}

__global__  void ker_gkylCartFieldScale(unsigned s, unsigned nv, double fact, double *out) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x + s; n < nv; n += blockDim.x * gridDim.x)
    out[n] *= fact;
}

__global__  void ker_gkylCartFieldAbs(unsigned s, unsigned nv, double *out) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x + s; n < nv; n += blockDim.x * gridDim.x)
    out[n] = fabs(out[n]);
}

__global__  void ker_gkylCopyFromField(double *data, double *f, unsigned numComponents, unsigned offset) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x; n < numComponents; n += blockDim.x * gridDim.x)
    data[n+offset] = f[n];
}

__global__  void ker_gkylCopyToField(double *f, double *data, unsigned numComponents, unsigned offset) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x; n < numComponents; n += blockDim.x * gridDim.x)
    f[n] = data[n+offset];
}

__global__  void ker_gkylCartFieldAssignAll(unsigned s, unsigned nv, double val, double *out) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x + s; n < nv; n += blockDim.x * gridDim.x)
    out[n] = val;
}

void gkylCartFieldDeviceAccumulate(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out)
{
  ker_gkylCartFieldAccumulate<<<numBlocks, numThreads>>>(s, nv, fact, inp, out);
}

void gkylCartFieldDeviceAssign(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out)
{
  ker_gkylCartFieldAssign<<<numBlocks, numThreads>>>(s, nv, fact, inp, out);
}

void gkylCartFieldDeviceScale(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, double *out)
{
  ker_gkylCartFieldScale<<<numBlocks, numThreads>>>(s, nv, fact, out);
}

void gkylCartFieldDeviceAbs(int numBlocks, int numThreads, unsigned s, unsigned nv, double *out)
{
  ker_gkylCartFieldAbs<<<numBlocks, numThreads>>>(s, nv, out);
}

void gkylCopyFromFieldDevice(int numBlocks, int numThreads, double *data, double *f, unsigned numComponents, unsigned c)
{
  ker_gkylCopyFromField<<<numBlocks, numThreads>>>(data, f, numComponents, c);
}

void gkylCopyToFieldDevice(int numBlocks, int numThreads, double *f, double *data, unsigned numComponents, unsigned c)
{
  ker_gkylCopyToField<<<numBlocks, numThreads>>>(f, data, numComponents, c);
}

void gkylCartFieldDeviceAssignAll(int numBlocks, int numThreads, unsigned s, unsigned nv, double val, double *out)
{
  ker_gkylCartFieldAssignAll<<<numBlocks, numThreads>>>(s, nv, val, out);
}
