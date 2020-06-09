#include <cstdio>
#include <GkylWavePropagation.h>
#include <GkylEuler.h>

__global__ void cuda_WavePropagation(GkylWavePropagation_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut) {

  GkylRange_t *localRange = fIn->localRange;
  unsigned int ndim = localRange->ndim;

  // set up indexers for localRange and fIn (localExtRange)
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = fIn->genIndexer();

  // get setup data from GkylWavePropagation_t structure
  GkylRectCart_t *grid = fIn->grid;
  int *updateDirs = hyper->updateDirs;
  int numUpdateDirs = hyper->numUpdateDirs;
  Gkyl::Euler *eq = hyper->equation;
  GkylCartField_t *cflRateByCell = hyper->cflRateByCell;

  // declaring this dummy array shared seems to alleviate register pressure and improve performance a bit
  extern __shared__ double dummy[];
  unsigned linearIdx = threadIdx.x + blockIdx.x*blockDim.x;

  int idxC[6];
  int idxL[6];
  int idxR[6];

  double xcC[6];
  double xcL[6];
  double xcR[6];

  // get i,j,k... index idxC from linear index linearIdx using localRange invIndexer
  localIdxr.invIndex(linearIdx, idxC);
  // convert i,j,k... index idxC into a linear index linearIdxC
  // note that linearIdxC != linearIdx.
  // this is because linearIdxC will have jumps because of ghost cells
  const int linearIdxC = fIdxr.index(idxC);

  grid->cellCenter(idxC, xcC);
  const double *dx = grid->dx;

  const double *fInC = fIn->getDataPtrAt(linearIdxC);
  double *fRhsOutC = fRhsOut->getDataPtrAt(linearIdxC);
  // cflRateByCell->getDataPtrAt(linearIdxC)[0] += cflRate;

  for(int i=0; i<numUpdateDirs; i++) {
    int dir = updateDirs[i] - 1;

    for(int d=0; d<ndim; d++) {
      idxL[d] = idxC[d];
      idxR[d] = idxC[d];
    }
    idxL[dir] = idxC[dir] - 1;
    idxR[dir] = idxC[dir] + 1;

    const int linearIdxL = fIdxr.index(idxL);
    const int linearIdxR = fIdxr.index(idxR);
    grid->cellCenter(idxL, xcL);
    grid->cellCenter(idxR, xcR);
    const double *fInL = fIn->getDataPtrAt(linearIdxL);
    const double *fInR = fIn->getDataPtrAt(linearIdxR);

    // hyper->maxsByCell->getDataPtrAt(linearIdxC)[i] = max(maxsL, maxsR);
  }

}

void wavePropagationAdvanceOnDevice(int numBlocks, int numThreads, GkylWavePropagation_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut) {
  cudaFuncSetAttribute(cuda_WavePropagation, cudaFuncAttributeMaxDynamicSharedMemorySize, 32*sizeof(double));
  cuda_WavePropagation<<<numBlocks, numThreads, 32*sizeof(double)>>>(hyper, fIn, fRhsOut);
}

__global__ void wavePropagationSetDtAndCflRateOnDevice(GkylWavePropagation_t *hyper, double dt, GkylCartField_t *cflRate) {
  hyper->dt = dt;
  hyper->cflRateByCell = cflRate;
}

void wavePropagationSetDtAndCflRate(GkylWavePropagation_t *hyper, double dt, GkylCartField_t *cflRate) {
  wavePropagationSetDtAndCflRateOnDevice<<<1, 1>>>(hyper, dt, cflRate);
}

