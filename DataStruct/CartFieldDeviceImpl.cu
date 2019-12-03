/* -*- c++ -*- */
// Gkyl ------------------------------------------------------------------------
//
// CUDA back-end for CartField
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <CartFieldDeviceImpl.h>

__global__ void ker_gkylCartFieldAccumulate(unsigned s, unsigned nv, double fact, const double *inp, double *out) {
  int n = blockIdx.x*blockDim.x + threadIdx.x;
  if (n >= s && n < (s+nv))
    out[n] += fact*inp[n];
}

__global__  void ker_gkylCartFieldAssign(unsigned s, unsigned nv, double fact, const double *inp, double *out) {
  int n = blockIdx.x*blockDim.x + threadIdx.x;
  if (n >= s && n < (s+nv))
    out[n] = fact*inp[n];
}

__global__  void ker_gkylCartFieldScale(unsigned s, unsigned nv, double fact, double *out) {
  int n = blockIdx.x*blockDim.x + threadIdx.x;
  if (n >= s && n < (s+nv))
    out[n] *= fact;
}

__global__  void ker_gkylCartFieldAbs(unsigned s, unsigned nv, double *out) {
  int n = blockIdx.x*blockDim.x + threadIdx.x;
  if (n >= s && n < (s+nv))
    out[n] = fabs(out[n]);
}

__global__  void ker_gkylCopyFromField(double *data, double *f, unsigned numComponents, unsigned offset) {
  int k = blockIdx.x*blockDim.x + threadIdx.x;
  if (k < numComponents)
    data[k+offset] = f[k];
}

__global__  void ker_gkylCopyToField(double *f, double *data, unsigned numComponents, unsigned offset) {
  int k = blockIdx.x*blockDim.x + threadIdx.x;
  if (k < numComponents)
    f[k] = data[k+offset];
}

__global__  void ker_gkylCartFieldAssignAll(unsigned s, unsigned nv, double val, double *out) {
  int n = blockIdx.x*blockDim.x + threadIdx.x;
  if (n >= s && n < (s+nv))
    out[n] = val;
}

void gkylCartFieldDeviceAccumulate(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out)
{
  ker_gkylCartFieldAccumulate<<<numBlocks, numThreads>>>(s, nv, fact, inp, out)
}

void gkylCartFieldDeviceAssign(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out)
{
  ker_gkylCartFieldAssign<<<numBlocks, numThreads>>>(s, nv, fact, inp, out)
}

void gkylCartFieldDeviceScale(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, double *out)
{
  ker_gkylCartFieldScale<<<numBlocks, numThreads>>>(s, nv, fact, out)
}

void gkylCartFieldDeviceAbs(int numBlocks, int numThreads, unsigned s, unsigned nv, double *out)
{
  ker_gkylCartFieldAbs<<<numBlocks, numThreads>>>(s, nv, out)
}

void gkylCopyFromFieldDevice(int numBlocks, int numThreads, double *data, double *f, unsigned numComponents, unsigned c)
{
  ker_gkylCopyFromField<<<numBlocks, numThreads>>>(data, f, numComponents, offset)
}

void gkylCopyToFieldDevice(int numBlocks, int numThreads, double *f, double *data, unsigned numComponents, unsigned c)
{
  ker_gkylCopyToField<<<numBlocks, numThreads>>>(f, data, numComponents, offset)
}

void gkylCartFieldDeviceAssignAll(int numBlocks, int numThreads, unsigned s, unsigned nv, double val, double *out)
{
  ker_gkylCartFieldAssignAll<<<numBlocks, numThreads>>>(s, nv, val, out)
}
