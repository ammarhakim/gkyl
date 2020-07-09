/* -*- c++ -*- */
// Gkyl ------------------------------------------------------------------------
//
// CUDA back-end for CartField
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <CartFieldDeviceImpl.h>

__global__ void ker_gkylCartFieldAccumulate(unsigned s, unsigned nv, double fact, const double *inp, double *out) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x + s; n < s + nv; n += blockDim.x * gridDim.x)
     out[n] += fact*inp[n];
}

__global__ void ker_gkylCartFieldAccumulateOffset(unsigned sInp, unsigned sOut, unsigned nCells, unsigned compStart, unsigned nCompInp, unsigned nCompOut, double fact, const double *inp, double *out) {
   if (nCompInp < nCompOut) {
      for (unsigned i=blockIdx.x*blockDim.x + threadIdx.x; i<nCells; i += blockDim.x * gridDim.x) {
         for (unsigned c=0; c<nCompInp; ++c) {
            out[sOut + i*nCompOut + compStart + c] += fact*inp[sInp + i*nCompInp + c];
         }
      }
   }
   else {
      for (unsigned i=blockIdx.x*blockDim.x + threadIdx.x; i<nCells; i += blockDim.x * gridDim.x) {
         for (unsigned c=0; c<nCompOut; ++c) {
            out[sOut + i*nCompOut + c] += fact*inp[sInp + i*nCompInp + compStart + c];
         }
      }
   }
}

__global__ void ker_gkylCartFieldAssign(unsigned s, unsigned nv, double fact, const double *inp, double *out) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x + s; n < s + nv; n += blockDim.x * gridDim.x)
    out[n] = fact*inp[n];
}

__global__ void ker_gkylCartFieldScale(unsigned s, unsigned nv, double fact, double *out) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x + s; n < s + nv; n += blockDim.x * gridDim.x)
    out[n] *= fact;
}

__global__ void ker_gkylCartFieldAbs(unsigned s, unsigned nv, double *out) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x + s; n < s + nv; n += blockDim.x * gridDim.x)
    out[n] = fabs(out[n]);
}

__global__ void ker_gkylPeriodicCopy(GkylRange_t *rangeSkin, GkylRange_t *rangeGhost, GkylCartField_t *f, unsigned numComponents) 
{
  // set up indexers for prescribed rangeSkin, rangeGhost, and f
  Gkyl::GenIndexer localIdxrSkin(rangeSkin);
  Gkyl::GenIndexer localIdxrGhost(rangeGhost);
  Gkyl::GenIndexer fIdxr = f->genIndexer();

  unsigned linearIdx = threadIdx.x + blockIdx.x*blockDim.x;
  int idxSkin[6];
  int idxGhost[6];

  // get i,j,k... index idxSkin and idxGhost from linear index linearIdx using range invIndexer
  localIdxrSkin.invIndex(linearIdx, idxSkin);
  localIdxrGhost.invIndex(linearIdx, idxGhost);

  // convert i,j,k... index idxSkin and idxGhost into a linear index linearIdxSkin/Ghost
  const int linearIdxSkin = fIdxr.index(idxSkin);
  const int linearIdxGhost = fIdxr.index(idxGhost);

  const double *fSkin = f->getDataPtrAt(linearIdxSkin);
  double *fGhost = f->getDataPtrAt(linearIdxGhost);

  // Copy the appropriate skin cell data into the corresponding ghost cells.
  for (unsigned k=0; k<numComponents; ++k) {
    fGhost[k] = fSkin[k];
  }
}

__global__ void ker_gkylCartFieldAssignAll(unsigned s, unsigned nv, double val, double *out) 
{
  for (int n = blockIdx.x*blockDim.x + threadIdx.x + s; n < s + nv; n += blockDim.x * gridDim.x)
    out[n] = val;
}

void gkylCartFieldDeviceAccumulate(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out)
{
  ker_gkylCartFieldAccumulate<<<numBlocks, numThreads>>>(s, nv, fact, inp, out);
}

void gkylCartFieldDeviceAccumulateOffset(int numBlocks, int numThreads, unsigned sInp, unsigned sOut, unsigned nCells, unsigned compStart, unsigned nCompInp, unsigned nCompOut, double fact, const double *inp, double *out)
{
  ker_gkylCartFieldAccumulateOffset<<<numBlocks, numThreads>>>(sInp, sOut, nCells, compStart, nCompInp, nCompOut, fact, inp, out);
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

void gkylDevicePeriodicCopy(int numBlocks, int numThreads, GkylRange_t *rangeSkin, GkylRange_t *rangeGhost, GkylCartField_t *f, unsigned numComponents)
{
  ker_gkylPeriodicCopy<<<numBlocks, numThreads>>>(rangeSkin, rangeGhost, f, numComponents);
}

void gkylCartFieldDeviceAssignAll(int numBlocks, int numThreads, unsigned s, unsigned nv, double val, double *out)
{
  ker_gkylCartFieldAssignAll<<<numBlocks, numThreads>>>(s, nv, val, out);
}
