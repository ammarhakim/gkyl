#include <cstdio>
#include <GkylHyperDisCont.h>
#include <GkylVlasov.h>
#include <VlasovModDecl.h>

__global__ void cuda_HyperDisCont(GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut) {

  GkylRange_t *localRange = fIn->localRange;
  unsigned int ndim = localRange->ndim;
  unsigned int numComponents = fRhsOut->numComponents;
  
  // set up indexers for localRange and fIn (localExtRange)
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = fIn->genIndexer();

  // get setup data from GkylHyperDisCont_t structure
  GkylRectCart_t *grid = fIn->grid;
  int *updateDirs = hyper->updateDirs;
  int numUpdateDirs = hyper->numUpdateDirs;
  bool *zeroFluxFlags = hyper->zeroFluxFlags;
  bool clearOut = hyper->clearOut;
  Gkyl::Vlasov *eq = hyper->equation;
  
  // CUDA thread "loop" over (non-ghost) cells in local range
  for(unsigned int linearIdx = threadIdx.x + blockIdx.x*blockDim.x; linearIdx < localRange->volume(); linearIdx += blockDim.x*gridDim.x) {
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
    int linearIdxC = fIdxr.index(idxC);

    grid->cellCenter(idxC, xcC);
    double *dx = grid->dx;
    
    double *fInC = fIn->getDataPtrAt(linearIdxC);
    double *fRhsOutC = fRhsOut->getDataPtrAt(linearIdxC);
    if(clearOut) {
      memset(fRhsOutC, 0., sizeof(double)*numComponents);
    }
    double cflRate = eq->volTerm(xcC, dx, idxC, fInC, fRhsOutC);

    // hard code this size for now. 
    // should be numComponents, but want to avoid dynamic memory alloc
    double dummy[200];
 
    for(int i=0; i<numUpdateDirs; i++) {
      int dir = updateDirs[i] - 1;

      for(int d=0; d<ndim; d++) {
        if(d!=dir) {
          idxL[d] = idxC[d];
          idxR[d] = idxC[d];
        } else {
          idxL[d] = idxC[d] - 1;
          idxR[d] = idxC[d] + 1;
        }
      }

      int linearIdxL = fIdxr.index(idxL);
      int linearIdxR = fIdxr.index(idxR);
      grid->cellCenter(idxL, xcL);
      grid->cellCenter(idxR, xcR);
      double *fInL = fIn->getDataPtrAt(linearIdxL);
      double *fInR = fIn->getDataPtrAt(linearIdxR);
      
      // left (of C) surface update. use dummy in place of fRhsOutL (cell to left of surface) so that only current cell (C) is updated.
      if(!(zeroFluxFlags[dir] && idxC[dir] == localRange->lower[dir])) {
        eq->surfTerm(dir, dummy, dummy, xcL, xcC, dx, dx, 0., idxL, idxC, fInL, fInC, dummy, fRhsOutC);
      } else if( zeroFluxFlags[dir]) {
        eq->boundarySurfTerm(dir, dummy, dummy, xcL, xcC, dx, dx, 0., idxL, idxC, fInL, fInC, dummy, fRhsOutC);
      }

      // right (of C) surface update. use dummy in place of fRhsOutR (cell to left of surface) so that only current cell (C) is updated.
      if(!(zeroFluxFlags[dir] && idxC[dir] == localRange->upper[dir])) {
        eq->surfTerm(dir, dummy, dummy, xcC, xcR, dx, dx, 0., idxC, idxR, fInC, fInR, fRhsOutC, dummy);
      } else if( zeroFluxFlags[dir]) {
        eq->boundarySurfTerm(dir, dummy, dummy, xcC, xcR, dx, dx, 0., idxC, idxR, fInC, fInR, fRhsOutC, dummy);
      }
    }
  }
} 

void advanceOnDevice(int numThreads, int numBlocks, GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut) {
  cuda_HyperDisCont<<<numThreads, numBlocks>>>(hyper, fIn, fRhsOut);
}
