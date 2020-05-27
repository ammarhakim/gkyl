#include <cstdio>
#include <GkylRectCart.h>
#include <GkylCartField.h>
#include <GkylRange.h>

__global__ void cuda_HyperDisCont(GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut) {

  unsigned int linIdx = threadIdx.x + blockIdx.x*blockDim.x;
  GkylRange_t *localRange = fIn->localRange;
  GkylRange_t *localExtRange = fIn->localExtRange;
  unsigned int nComp = fIn->numComponents;
  unsigned int ndim = localRange->ndim;
  
  Gkyl::GenIndexer locIdxr(localRange);
  Gkyl::GenIndexer fIdxr = fIn->genIndexer();

  // get setup data from GkylHyperDisCont_t structure
  GkylRectCart_t *grid = fIn->grid;
  int *updateDirs = hyper->updateDirs;
  int numUpdateDirs = hyper->numUpdateDirs;
  bool *zeroFluxFlags = hyper->zeroFluxFlags;
  GkylEquation_t *eq = hyper->equation;
  
  // "loop" over cells in local range
  if(linIdx < localRange->volume()) {
    int *idxC = new int[ndim];
    int *idxL = new int[ndim];
    int *idxR = new int[ndim];
    
    double *xcC = new double[ndim];
    double *xcL = new double[ndim];
    double *xcR = new double[ndim];

    double *dx = grid->dx;
    locIdxr.invIndexer(linIdx, idxC);
    int linIdxC = fIdxr.indexer(idxC);
    grid->cellCenter(idxC, xcC);
    
    double *fInC = fIn->getDevicePtr(linIdxC);
    double *fRhsOutC = fRhsOut->getDevicePtr(linIdxC);
    double cflRate = eq->volumeTerm(xcC, dx, idxC, fInC, fRhsOutC);

    double *dummy = new double[nComp];
 
    for(int i=0; i<numUpdateDirs; i++) {
      int dir = updateDirs[i];

      for(int d=0; d<ndim; d++) {
        if(d!=dir) {
          idxL[d] = idxC[d];
          idxR[d] = idxC[d];
        } else {
          idxL[d] = idxC[d] - 1;
          idxR[d] = idxC[d] + 1;
        }
      }

      int linIdxL = fIdxr.indexer(idxL);
      int linIdxR = fIdxr.indexer(idxR);
      grid->cellCenter(idxL, xcL);
      grid->cellCenter(idxR, xcR);
      double *fInL = fIn->getDevicePtr(linIdxL);
      double *fRhsOutL = fRhsOut->getDevicePtr(linIdxL);
      double *fInR = fIn->getDevicePtr(linIdxR);
      double *fRhsOutR = fRhsOut->getDevicePtr(linIdxR);
      
      // left (of C) surface update. use dummy in place of fRhsOutL (cell to left of surface) so that only current cell (C) is updated.
      if(!(zeroFluxFlags[dir] && idxC[dir] == localRange->lower[dir])) {
        eq->surfTerm(dir, xcL, xcC, dx, dx, idxL, idxC, fInL, fInC, dummy, fRhsOutC);
      } else if zeroFluxFlags[dir] {
        eq->boundarySurfTerm(dir, xcL, xcC, dx, dx, idxL, idxC, fInL, fInC, dummy, fRhsOutC);
      }

      // right (of C) surface update. use dummy in place of fRhsOutR (cell to left of surface) so that only current cell (C) is updated.
      if(!(zeroFluxFlags[dir] && idxC[dir] == localRange->lower[dir])) {
        eq->surfTerm(dir, xcC, xcR, dx, dx, idxC, idxR, fInC, fInR, fRhsOutC, dummy);
      } else if zeroFluxFlags[dir] {
        eq->boundarySurfTerm(dir, xcC, xcR, dx, dx, idxC, idxR, fInC, fInR, fRhsOutC, dummy);
      }
    }
  }
} 
