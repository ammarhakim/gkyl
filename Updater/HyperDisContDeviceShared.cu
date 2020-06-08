#include <cstdio>
#include <GkylHyperDisCont.h>
#include <GkylVlasov.h>
#include <VlasovModDecl.h>

__global__ void cuda_HyperDisCont_shared(GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut) {

  GkylRange_t *localRange = fIn->localRange;
  unsigned int ndim = localRange->ndim;
  unsigned int numComponents = fIn->numComponents;
  
  // set up indexers for localRange and fIn (localExtRange)
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = fIn->genIndexer();

  // get setup data from GkylHyperDisCont_t structure
  GkylRectCart_t *grid = fIn->grid;
  int *updateDirs = hyper->updateDirs;
  int numUpdateDirs = hyper->numUpdateDirs;
  bool *zeroFluxFlags = hyper->zeroFluxFlags;
  Gkyl::Vlasov *eq = hyper->equation;
 
  // shared memeory blocks of size blockDim.x * blockDim.y * numComponents
  __shared__ double dummy[32];
  extern __shared__ double fIn_shared[];
  unsigned linearIdx = threadIdx.x + blockIdx.x*blockDim.x;

  int idxC[6];
  int idxL[6];
  int idxR[6];
  
  double xcC[6];
  double xcL[6];
  double xcR[6];

  // read f into shared memory
  const int jump = blockDim.x/numComponents;
  for(int j=0; j<numComponents; j++) {
    const int sharedIdx = jump*j + blockDim.x*blockIdx.x;
    localIdxr.invIndex(sharedIdx, idxC);
    const int linearIdxC = fIdxr.index(idxC);
    if(sharedIdx < localRange->volume()) {
      fIn_shared[threadIdx.x + j*blockDim.x] = fIn->_data[threadIdx.x + linearIdxC*numComponents]; 
    }
  }

  // sync to make sure all data has been loaded
  __syncthreads();

  if(linearIdx < localRange->volume()) {
    // get i,j,k... index idxC from linear index linearIdx using localRange invIndexer
    localIdxr.invIndex(linearIdx, idxC);
    grid->cellCenter(idxC, xcC);
    const double *dx = grid->dx;
    // convert i,j,k... index idxC into a linear index linearIdxC
    // note that linearIdxC != linearIdx.
    // this is because linearIdxC will have jumps because of ghost cells
    const int linearIdxC = fIdxr.index(idxC);
    
    const double *fInC = &fIn_shared[numComponents*threadIdx.x];
    double *fRhsOutC = fRhsOut->getDataPtrAt(linearIdxC);
    double cflRate = eq->volTerm_shared(xcC, dx, idxC, fInC, fRhsOutC);

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
      if(!(zeroFluxFlags[dir] && idxC[dir] == localRange->lower[dir])) {
        eq->surfTerm(dir, xcL, xcC, dx, dx, 0., idxL, idxC, fInL, fInC, dummy, fRhsOutC);
      } else if( zeroFluxFlags[dir]) {
        eq->boundarySurfTerm(dir, xcL, xcC, dx, dx, 0., idxL, idxC, fInL, fInC, dummy, fRhsOutC);
      }

      // right (of C) surface update. use dummy in place of fRhsOutR (cell to left of surface) so that only current cell (C) is updated.
      if(!(zeroFluxFlags[dir] && idxC[dir] == localRange->upper[dir])) {
        eq->surfTerm(dir, xcC, xcR, dx, dx, 0., idxC, idxR, fInC, fInR, fRhsOutC, dummy);
      } else if( zeroFluxFlags[dir]) {
        eq->boundarySurfTerm(dir, xcC, xcR, dx, dx, 0., idxC, idxR, fInC, fInR, fRhsOutC, dummy);
      }
    }
  } 
} 

void advanceOnDevice_shared(int numBlocks, int numThreads, int numComponents, GkylHyperDisCont_t *hyper, GkylCartField_t *fIn, GkylCartField_t *fRhsOut) {
  cudaFuncSetAttribute(cuda_HyperDisCont_shared, cudaFuncAttributeMaxDynamicSharedMemorySize, ((numComponents)*numThreads+numComponents)*sizeof(double));
  cuda_HyperDisCont_shared<<<numBlocks, numThreads, ((numComponents)*numThreads)*sizeof(double)>>>(hyper, fIn, fRhsOut);
}
