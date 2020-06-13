#include <cstdio>
#include <GkylWavePropagation.h>
#include <GkylEuler.h>

__device__ static void calcDelta(
  const double *ql, const double *qr, double *delta, const int meqn)
{
  for (int i = 0; i < meqn; i++) {
    delta[i] = qr[i] - ql[i];
  }
}

__device__ static void calcFirstOrderGud(
  const double dtdx, double *ql, double *qr, const double *amdq,
  const double *apdq, const int meqn)
{
  for (int i = 0; i< meqn; i++) {
    /* qr[i] -= dtdx * apdq[i]; */
    /* ql[i] -= dtdx * amdq[i]; */
    // XXX calling __threadfence_system() between two calcFirstOrderGud dummy
    // calls fails occasionally with small numThreads
    atomicAdd(qr+i, -dtdx * apdq[i]);
    atomicAdd(ql+i, -dtdx * amdq[i]);
  }
}

__device__ static double calcCfla(
  const double cfla, const double dtdx, const double *s, const int mwave)
{
  double c = cfla;
  for (int i = 0; i < mwave; i ++) {
    c = max(c, dtdx * abs(s[i]));
  }
  return c;
}

__global__ void cuda_WavePropagation(
  GkylWavePropagation_t *hyper, GkylCartField_t *qIn, GkylCartField_t *qOut)
{

  GkylRange_t *localRange = qIn->localRange;
  GkylRange_t *localExtRange = qIn->localExtRange;
  GkylRange_t *localExtEdgeRange = qIn->localExtEdgeRange;
  int ndim = localRange->ndim;

  // set up indexers for localRange and qIn (localExtRange)
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = qIn->genIndexer();

  // get setup data from GkylWavePropagation_t structure
  GkylRectCart_t *grid = qIn->grid;
  int *updateDirs = hyper->updateDirs;
  int numUpdateDirs = hyper->numUpdateDirs;
  Gkyl::Euler *eq = hyper->equation;
  GkylCartField_t *dtByCell = hyper->dtByCell;

  const int meqn = eq->numEquations();
  const int mwave = eq->numWaves();

  // XXX use meqn and mwave
  double delta[10];
  double waves[30];
  double s[3];
  double amdq[10];
  double apdq[10];

  // declaring this dummy array shared seems to alleviate register pressure and
  // improve performance a bit
  extern __shared__ double dummy[];
  int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;

  int idxC[3];
  int idxL[3];
  int idxR[3];

  // get i,j,k... index idxC
  localIdxr.invIndex(linearIdx, idxC);
  // if ndim>1, linearIdxC!=linearIdx since linearIdxC jumps due to ghost cells
  const int linearIdxC = fIdxr.index(idxC);

  const double *dx = grid->dx;

  double cfl = hyper->_cfl;
  double cflm = hyper->_cflm;
  double cfla = 0; // actual CFL number used

  const double *qInC = qIn->getDataPtrAt(linearIdxC);
  double *qOutC = qOut->getDataPtrAt(linearIdxC);

  if(linearIdx < localExtEdgeRange->volume()) {
    for(int i = 0; i < meqn; i++) {
      qOutC[i] = qInC[i];
    }
  }

  for(int i=0; i<numUpdateDirs; i++) {
    int dir = updateDirs[i] - 1;
    const double dtdx = hyper->dt / dx[dir];

    for(int d=0; d<ndim; d++) {
      idxL[d] = idxC[d];
      idxR[d] = idxC[d];
    }
    idxL[dir] = idxC[dir] - 1;
    idxR[dir] = idxC[dir];

    const int linearIdxL = fIdxr.index(idxL);
    const int linearIdxR = fIdxr.index(idxR);
    const double *qInL = qIn->getDataPtrAt(linearIdxL);
    const double *qInR = qIn->getDataPtrAt(linearIdxR);

    double *qOutL = qOut->getDataPtrAt(linearIdxL);
    double *qOutR = qOut->getDataPtrAt(linearIdxR);

    if(linearIdx < localExtEdgeRange->volume()) {
      calcDelta(qInL, qInR, delta, meqn);

      eq->rp(dir, delta, qInL, qInR, waves, s);
      eq->qFluctuations(dir, qInL, qInR, waves, s, amdq, apdq);

      calcFirstOrderGud(dtdx, qOutL, qOutR, amdq, apdq, meqn);
      // XXX following fails with small numThreads
      /* calcFirstOrderGud(dtdx, qOutL, dummy, amdq, apdq, meqn); */
      /* __threadfence_system(); */
      /* calcFirstOrderGud(dtdx, dummy, qOutR, amdq, apdq, meqn); */

      cfla = calcCfla(cfla, dtdx, s, mwave);
    }
  }

  dtByCell->getDataPtrAt(linearIdxC)[0] = hyper->dt * cfl/cfla;
}

void wavePropagationAdvanceOnDevice(
  int numBlocks, int numThreads, GkylWavePropagation_t *hyper,
  GkylCartField_t *qIn, GkylCartField_t *qOut)
{
  Gkyl::Euler *eq = hyper->equation;
  // XXX
  const int meqn = 5; // eq->numEquations();
  const int mwave = 1; // eq->numWaves();
  const int nComponents = (mwave + mwave * meqn);
  const int sharedMemSize = numThreads * nComponents;

  cudaFuncSetAttribute(
    cuda_WavePropagation, cudaFuncAttributeMaxDynamicSharedMemorySize,
    sharedMemSize*sizeof(double));
  cuda_WavePropagation<<<numBlocks, numThreads, sharedMemSize*sizeof(double)>>>(
    hyper, qIn, qOut);
}

__global__ void setDtOnDevice(GkylWavePropagation_t *hyper, double dt) {
  hyper->dt = dt;
}

void setDt(GkylWavePropagation_t *hyper, double dt) {
  setDtOnDevice<<<1, 1>>>(hyper, dt);
}
