/* -*- c++ -*- */

#include <cstdio>
#include <GkylHyperDisCont.h>
#include <GkylEquation.h>
#include <VlasovModDecl.h>

__global__ void cuda_HyperDisCont(GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut) {

  GkylRange_t *localRange = fIn->localRange;
  unsigned int ndim = localRange->ndim;
  
  // set up indexers for localRange and fIn (localExtRange)
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = fIn->genIndexer();

  // get setup data from GkylHyperDisCont_t structure
  GkylRectCart_t *grid = fIn->grid;
  int *updateDirs = hyper->updateDirs;
  int numUpdateDirs = hyper->numUpdateDirs;
  bool *zeroFluxFlags = hyper->zeroFluxFlags;
  GkylEquation_t *eq = hyper->equation;
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
    const double cflRate = eq->volTerm(xcC, dx, idxC, fInC, fRhsOutC);
    cflRateByCell->getDataPtrAt(linearIdxC)[0] += cflRate;

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
      
      // left (of C) surface update. use dummy in place of fRhsOutL (cell to left of surface) so that only current cell (C) is updated.
      double maxsL, maxsR;
      if(!(zeroFluxFlags[dir] && idxC[dir] == localRange->lower[dir])) {
        maxsL = eq->surfTerm(dir, xcL, xcC, dx, dx, hyper->maxs[i], idxL, idxC, fInL, fInC, dummy, fRhsOutC);
      } else if( zeroFluxFlags[dir]) {
        eq->boundarySurfTerm(dir, xcL, xcC, dx, dx, hyper->maxs[i], idxL, idxC, fInL, fInC, dummy, fRhsOutC);
      }

      // right (of C) surface update. use dummy in place of fRhsOutR (cell to left of surface) so that only current cell (C) is updated.
      if(!(zeroFluxFlags[dir] && idxC[dir] == localRange->upper[dir])) {
        maxsR = eq->surfTerm(dir, xcC, xcR, dx, dx, hyper->maxs[i], idxC, idxR, fInC, fInR, fRhsOutC, dummy);
      } else if( zeroFluxFlags[dir]) {
        eq->boundarySurfTerm(dir, xcC, xcR, dx, dx, hyper->maxs[i], idxC, idxR, fInC, fInR, fRhsOutC, dummy);
      }
      hyper->maxsByCell->getDataPtrAt(linearIdxC)[i] = max(maxsL, maxsR);
    }
  
} 

void advanceOnDevice(int numBlocks, int numThreads, int numComponents, GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut) {
  cudaFuncSetAttribute(cuda_HyperDisCont, cudaFuncAttributeMaxDynamicSharedMemorySize, numComponents*sizeof(double));
  cuda_HyperDisCont<<<numBlocks, numThreads, numComponents*sizeof(double)>>>(hyper, fIn, fRhsOut);
}

__global__ void setDtAndCflRateOnDevice(GkylHyperDisCont_t *hyper, double dt, GkylCartField_t *cflRate) {
  hyper->dt = dt;
  hyper->cflRateByCell = cflRate;
}

void setDtAndCflRate(GkylHyperDisCont_t *hyper, double dt, GkylCartField_t *cflRate) {
  setDtAndCflRateOnDevice<<<1, 1>>>(hyper, dt, cflRate);
}
