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
    /* if (threadIdx.x ==0 && blockIdx.x==0) */
    /*   printf("ql[0] before %13g\n", ql[0]); */
  for (int i = 0; i< meqn; i++) {
    /* qr[i] -= dtdx * apdq[i]; */
    /* ql[i] -= dtdx * amdq[i]; */
    // XXX calling __threadfence_system() between two calcFirstOrderGud dummy
    // calls fails occasionally with small numThreads
    atomicAdd(qr+i, -dtdx * apdq[i]);
    atomicAdd(ql+i, -dtdx * amdq[i]);
  }
    /* if (threadIdx.x ==0 && blockIdx.x==0) */
    /*   printf("ql[0] after %13g\n", ql[0]); */
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

__device__ static double waveDotProd(
    const double *waves, const double *waves1, const int mw,
    const int meqn) {
  double result = 0.;
  for (int i = 0; i < meqn; i++) {
    result += waves[meqn*mw+i] * waves[meqn*mw+i];
  }
  return result;
}

__device__ static inline double limiter_minMod(const double r) {
   return max(0., min(1., r));
}

__device__ static void limitWaves(
    double *waveSlice, const double *speedSlice, const int mwave,
    const int meqn) {
  int i = threadIdx.x + 1;  // FIXME
  int jump = meqn * mwave;
  for (int mw = 0; mw < mwave; mw++ ){
    const double wnorm2 = waveDotProd(
        waveSlice+i*jump, waveSlice+i*jump, mw, meqn);
    double wlimitr = 1.;
    if (wnorm2 > 0) {
      double r;
        const double dotl = waveDotProd(
            waveSlice+(i-1)*jump, waveSlice+i*jump, mw, meqn);
        const double dotr = waveDotProd(
            waveSlice+(i+1)*jump, waveSlice+i*jump, mw, meqn);
        r = speedSlice[i] > 0 ? dotl/wnorm2 : dotr/wnorm2;
      wlimitr = limiter_minMod(r);
        printf("[%2d] dotl %13g dotr %13g wnorm2 %13g\n",
          i, dotl, dotr, wnorm2);
    }
    __threadfence_block();
    for (int me = 0; me < meqn; me++) {
      waveSlice[i*jump+mw*meqn] *= wlimitr;
    }
  }
}

__device__ static void secondOrderFlux(
  const double dtdx, const double s, const double *wave, double *fs,
  const int meqn) {
  double sfact = 0.5 * abs(s) * (1 - abs(s) * dtdx);
  for (int i = 0; i < meqn; i++) {
    fs[i] += sfact * wave[i];
  }
}

__device__ static void secondOrderUpdate(
    const double dtdx, const double *fs, const double *fs1, double *q,
    const int meqn) {
  for (int i = 0; i < meqn; i++) {
    q[i] -= dtdx * (fs1[i] - fs[i]);
  }
}

__device__ static void copyComponents(
    const double *ptrFrom, double *ptrTo, const int nComponents) {
  for (int i = 0; i < nComponents; i++) {
    ptrTo[i] = ptrFrom[i];
  }
}

__global__ void cuda_WavePropagation(
  GkylWavePropagation_t *hyper, GkylCartField_t *qIn, GkylCartField_t *qOut)
{

  GkylRange_t *localRange = qIn->localRange;
  GkylRange_t *localExtRange = qIn->localExtRange;
  GkylRange_t *localEdgeRange = qIn->localEdgeRange;
  GkylRange_t *localExtEdgeRange = qIn->localExtEdgeRange;
  int ndim = localRange->ndim;

  // set up indexers for localRange and qIn (localExtRange)
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer extIdxr(localExtRange);
  Gkyl::GenIndexer edgeIdxr(localExtEdgeRange);
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
  double delta[5];
  double amdq[5];
  double apdq[5];

  // declaring this dummy array shared seems to alleviate register pressure and
  // improve performance a bit
  extern __shared__ double dummy[];
  int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;

  // assign buffer space for different usages
  int base = 0;

  const int baseWaveSlice = base;
  base += (meqn * mwave) * blockDim.x;
  double *waveSlice = dummy + baseWaveSlice;

  const int baseSpeedSlice = base;
  base += (mwave) * blockDim.x;
  double *speedSlice = dummy + baseSpeedSlice;

  const int baseFluxSlice = base;
  base += (meqn) * blockDim.x;
  double *fluxSlice = dummy + baseFluxSlice;

  // find buffer addresses for each thread
  // FIXME shall waves and s be created on the fly and then copied into slices
  double *waves = waveSlice + (meqn * mwave) * threadIdx.x;
  double *s = speedSlice + (mwave) * threadIdx.x;
  double *flux = fluxSlice + (meqn) * threadIdx.x;

  int idxC[3];
  int idxL[3];
  int idxR[3];

  // get i,j,k... index idxC
  localIdxr.invIndex(linearIdx, idxC);
  //edgeIdxr.invIndex(linearIdx, idxC);
  // if ndim>1, linearIdxC!=linearIdx since linearIdxC jumps due to ghost cells
  const int linearIdxC = fIdxr.index(idxC);

  const double *dx = grid->dx;

  double cfl = hyper->_cfl;
  double cflm = hyper->_cflm;
  double cfla = 0; // actual CFL number used

  const double *qInC = qIn->getDataPtrAt(linearIdxC);
  double *qOutC = qOut->getDataPtrAt(linearIdxC);

  /* if(linearIdx < localRange->volume()) { */
  /*   for(int i = 0; i < meqn; i++) { */
  /*     qOutC[i] = qInC[i]; */
  /*   } */
  /*   if (threadIdx.x == 0 && blockIdx.x == 0) { */
  /*     for(int i = 0; i < meqn; i++) { */
  /*       qOutC[i] = qInC[i]; */
  /*     } */
  /*   } */
  /* } */
    /* __syncthreads(); */
    /* if(linearIdx < localExtEdgeRange->volume()) */
    /*   printf("[%2d] qOut [%2d] %13g\n", linearIdx, idxC[0], qOutC[0]); */

  for(int i=0; i<numUpdateDirs; i++) {
    int dir = updateDirs[i] - 1;
    const double dtdx = hyper->dt / dx[dir];

    for(int d=0; d<ndim; d++) {
      idxL[d] = idxC[d];
      idxR[d] = idxC[d];
    }
    /* idxL[dir] = idxC[dir] - 1; */
    /* idxR[dir] = idxC[dir]; */
    idxL[dir] = idxC[dir] - 2;
    idxR[dir] = idxC[dir] - 1;

    const int linearIdxL = fIdxr.index(idxL);
    const int linearIdxR = fIdxr.index(idxR);
    const double *qInL = qIn->getDataPtrAt(linearIdxL);
    const double *qInR = qIn->getDataPtrAt(linearIdxR);

    double *qOutL = qOut->getDataPtrAt(linearIdxL);
    double *qOutR = qOut->getDataPtrAt(linearIdxR);

    /* __syncthreads(); */
    /* if(linearIdx < localExtEdgeRange->volume()) */
    /*   printf("[%2d] before L [%2d] %13g, R [%2d] %13g; wave %13g; amdq %13g; qOutL %13g\n", */
    /*       linearIdx, idxL[dir], qInL[0], idxR[dir], qInR[0], waves[0], amdq[0], qOutL[0]); */

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

    /* __syncthreads(); */
    /* if(linearIdx < localExtEdgeRange->volume()) */
    /*   printf("[%2d] L [%2d] %13g, R [%2d] %13g; wave %13g; amdq %13g; qOutL %13g\n", */
    /*       linearIdx, idxL[dir], qInL[0], idxR[dir], qInR[0], waves[0], amdq[0], qOutL[0]); */

    if(linearIdx < localEdgeRange->volume()) {
      limitWaves(waveSlice, speedSlice, mwave, meqn);
    }

    __syncthreads();
    if(linearIdx < localExtEdgeRange->volume()) {
      printf("[%2d] limited wave %13g\n",
          linearIdx, waves[0]);
    }

    /* if(linearIdx < localEdgeRange->volume()) { */
    /*   int idx = threadIdx.x + 1; */
    /*   double *waves = waveSlice + (meqn * mwave) * i; */
    /*   double *s = speedSlice + (mwave) * i; */
    /*   double *flux = fluxSlice + (meqn) * i; */
    /*  */
    /*   for (int c = 0; c < meqn; c++) { */
    /*     flux[c] = 0; */
    /*   } */
    /*  */
    /*   for (int mw = 0; mw < mwave; mw++) { */
    /*     secondOrderFlux(dtdx, s[mw], waves+mw*meqn, flux, meqn); */
    /*   } */
    /* } */
    /*  */
    /* if(linearIdx < localRange->volume()) { */
    /*   int idx = threadIdx.x + 1; */
    /*   double *flux = fluxSlice + (meqn) * i; */
    /*   double *flux1 = fluxSlice + (meqn) * (i + 1); */
    /*   idxC[dir] += 1; */
    /*   const int linearIdxCC = fIdxr.index(idxC); */
    /*   double *qOutCC = qOut->getDataPtrAt(linearIdxCC); */
    /*   secondOrderUpdate(dtdx, flux, flux1, qOutCC, meqn); */
    /* } */
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
  const int nComponents = (mwave + mwave * meqn + meqn);
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
