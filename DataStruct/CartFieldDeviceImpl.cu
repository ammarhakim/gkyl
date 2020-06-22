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

__global__ void ker_gkylPeriodicCopy(int dir, GkylCartField_t *f) 
{
  // Get numComponents and number of ghost cells
  unsigned numComponents = f->numComponents;
  int lowerGhost = f->lowerGhost;
  int upperGhost = f->upperGhost;
  // Get local range and compute skin cell and ghost cell ranges
  GkylRange_t *localRange = f->localRange;
  GkylRange_t lowerSkinRange = localRange.lowerSkin(dir, upperGhost);
  GkylRange_t upperSkinRange = localRange.upperSkin(dir, lowerGhost);
  GkylRange_t lowerGhostRange = localRange.lowerGhost(dir, lowerGhost);
  GkylRange_t upperGhostRange = localRange.upperGhost(dir, upperGhost);
  // Set up indexers for prescribed skin and ghost ranges, and f
  Gkyl::GenIndexer localIdxrLowerSkin(lowerSkinRange);
  Gkyl::GenIndexer localIdxrUpperSkin(upperSkinRange);
  Gkyl::GenIndexer localIdxrLowerGhost(lowerGhostRange);
  Gkyl::GenIndexer localIdxrUpperGhost(upperGhostRange);
  Gkyl::GenIndexer fIdxr = f->genIndexer();

  unsigned linearIdx = threadIdx.x + blockIdx.x*blockDim.x;
  int idxLowerSkin[6];
  int idxUpperSkin[6];
  int idxLowerGhost[6];
  int idxUpperGhost[6];

  // Get i,j,k... index idxLower/UpperSkin and idxLower/UpperGhost from linear index linearIdx using range invIndexer
  localIdxrLowerSkin.invIndex(linearIdx, idxLowerSkin);
  localIdxrUpperSkin.invIndex(linearIdx, idxUpperSkin);
  localIdxrLowerGhost.invIndex(linearIdx, idxLowerGhost);
  localIdxrUpperGhost.invIndex(linearIdx, idxUpperGhost);

  // convert i,j,k... index idxLower/UpperSkin and idxLower/UpperGhost into a linear index linearIdxLower/UpperSkin/Ghost
  const int linearIdxLowerSkin = fIdxr.index(idxLowerSkin);
  const int linearIdxUpperSkin = fIdxr.index(idxUpperSkin);
  const int linearIdxLowerGhost = fIdxr.index(idxLowerGhost);
  const int linearIdxUpperGhost = fIdxr.index(idxUpperGhost);

  const double *fLowerSkin = f->getDataPtrAt(linearIdxLowerSkin);
  const double *fUpperSkin = f->getDataPtrAt(linearIdxUpperSkin);
  double *fLowerGhost = f->getDataPtrAt(linearIdxLowerGhost);
  double *fUpperGhost = f->getDataPtrAt(linearIdxUpperGhost);

  // Copy the appropriate skin cell data into the corresponding ghost cells.
  for (unsigned k=0; k<numComponents; ++k) {
    fUpperGhost[k] = fLowerSkin[k];
    fLowerGhost[k] = fUpperSkin[k];
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

void gkylDevicePeriodicCopy(int numBlocks, int numThreads, int dir, GkylCartField_t *f)
{
  ker_gkylPeriodicCopy<<<numBlocks, numThreads>>>(dir, f);
}

void gkylCartFieldDeviceAssignAll(int numBlocks, int numThreads, unsigned s, unsigned nv, double val, double *out)
{
  ker_gkylCartFieldAssignAll<<<numBlocks, numThreads>>>(s, nv, val, out);
}
