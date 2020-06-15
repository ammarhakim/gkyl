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

__device__ static double waveDotProd(
    const double *waves, const double *waves1, const int mw, const int meqn) {
  double result = 0.;
  for (int i = 0; i < meqn; i++) {
    result += waves[meqn*mw+i] * waves1[meqn*mw+i];
  }
  return result;
}

__device__ static inline double limiter_minMod(const double r) {
   return max(0., min(1., r));
}

__device__ static void limitWaves(
    const double *waves, const double *speeds, double *limitedWaves,
    const int mwave, const int meqn) {
  int jump = meqn * mwave;
  for (int mw = 0; mw < mwave; mw++ ){
    const double wnorm2 = waveDotProd(waves, waves, mw, meqn);
    if (wnorm2 > 0) {
      double wlimitr = 1.;
      const double dotl = waveDotProd(waves-jump, waves, mw, meqn);
      const double dotr = waveDotProd(waves+jump, waves, mw, meqn);
      const double r = speeds[mw] > 0 ? dotl/wnorm2 : dotr/wnorm2;
      wlimitr = limiter_minMod(r);
      for (int me = 0; me < meqn; me++) {
        limitedWaves[mw*meqn+me] *= wlimitr;
      }
    }
  }
}

__device__ static void secondOrderFluxOneWave(
  const double dtdx, const double s, const double *wave, double *fs,
  const int meqn) {
  double sfact = 0.5 * abs(s) * (1 - abs(s) * dtdx);
  for (int i = 0; i < meqn; i++) {
    fs[i] += sfact * wave[i];
  }
}

__device__ static void secondOrderFlux(
  const double dtdx, const double *s, const double *waves, double *fs,
  const int meqn, const int mwave) {
    for (int mw = 0; mw < mwave; mw++) {
      secondOrderFluxOneWave(dtdx, s[mw], waves+mw*meqn, fs, meqn);
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
  double delta[5];
  double amdq[5];
  double apdq[5];

  // declaring this dummy array shared seems to alleviate register pressure and
  // improve performance a bit
  extern __shared__ double dummy[];
  int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;

  // assign buffer space for different usages
  int base = 0;

  // numThreads == numRealCells
  // waveSlice and speedSlice are defined on ghost-ghost, ghost-real, and
  // real-real faces, thus the +3
  const int baseWaveSlice = base;
  base += (meqn * mwave) * (blockDim.x + 3);
  double *waveSlice = dummy + baseWaveSlice;

  const int baseSpeedSlice = base;
  base += (mwave) * (blockDim.x + 3);
  double *speedSlice = dummy + baseSpeedSlice;

  // limitedWaves and second-order fluxSlice are defiend on ghost-real and
  // real-real faces, thus the +1
  const int baseLimitedWaveSlice = base;
  base += (meqn * mwave) * (blockDim.x + 1);
  double *limitedWaveSlice = dummy + baseLimitedWaveSlice;

  const int baseFluxSlice = base;
  base += (meqn) * (blockDim.x + 1);
  double *fluxSlice = dummy + baseFluxSlice;

  // find buffer addresses for each thread
  double *waves = waveSlice + (meqn * mwave) * (threadIdx.x+1);
  double *s = speedSlice + (mwave) * (threadIdx.x+1);
  double *limitedWaves = limitedWaveSlice + (meqn * mwave) * (threadIdx.x);
  double *flux = fluxSlice + (meqn) * (threadIdx.x);

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

  for(int i=0; i<numUpdateDirs; i++) {
    int dir = updateDirs[i] - 1;
    const double dtdx = hyper->dt / dx[dir];

    for(int d=0; d<ndim; d++) {
      idxL[d] = idxC[d];
      idxR[d] = idxC[d];
    }
    // XXX firstOrder stuff is over extended edges, but idxC was calculated
    // from localIdxr.invIndex
    idxL[dir] = idxC[dir] - 1;
    idxR[dir] = idxC[dir] - 0;

    const int linearIdxL = fIdxr.index(idxL);
    const int linearIdxR = fIdxr.index(idxR);
    const double *qInL = qIn->getDataPtrAt(linearIdxL);
    const double *qInR = qIn->getDataPtrAt(linearIdxR);

    double *qOutL = qOut->getDataPtrAt(linearIdxL);
    double *qOutR = qOut->getDataPtrAt(linearIdxR);

    if(linearIdx < localRange->volume()) {
      calcDelta(qInL, qInR, delta, meqn);

      eq->rp(dir, delta, qInL, qInR, waves, s);
      eq->qFluctuations(dir, qInL, qInR, waves, s, amdq, apdq);

      calcFirstOrderGud(dtdx, qOutL, qOutR, amdq, apdq, meqn);
      // XXX following fails with small numThreads
      /* calcFirstOrderGud(dtdx, qOutL, dummy, amdq, apdq, meqn); */
      /* __threadfence_system(); */
      /* calcFirstOrderGud(dtdx, dummy, qOutR, amdq, apdq, meqn); */

      cfla = calcCfla(cfla, dtdx, s, mwave);

      copyComponents(waves, limitedWaves, meqn * mwave);

      // can we avoid branching?
      // solve one additional Riemann problem on the lower side
      if (threadIdx.x==0) {
        int inc = -1;
        idxL[dir] += inc;
        idxR[dir] += inc;
        const int linearIdxL = fIdxr.index(idxL);
        const int linearIdxR = fIdxr.index(idxR);
        const double *qInL = qIn->getDataPtrAt(linearIdxL);
        const double *qInR = qIn->getDataPtrAt(linearIdxR);
        calcDelta(qInL, qInR, delta, meqn);
        eq->rp(dir, delta, qInL, qInR, waves+inc*meqn*mwave, s+inc*mwave);
        idxL[dir] -= inc;
        idxR[dir] -= inc;
      }

      // solve two additional Riemann problems on the higher side and update the
      // last real cell
      if (threadIdx.x==blockDim.x-1 || linearIdx==localRange->volume()-1) {
        for (int inc = 1; inc < 3; inc++) {
          idxL[dir] += inc;
          idxR[dir] += inc;
          const int linearIdxL = fIdxr.index(idxL);
          const int linearIdxR = fIdxr.index(idxR);
          const double *qInL = qIn->getDataPtrAt(linearIdxL);
          const double *qInR = qIn->getDataPtrAt(linearIdxR);
          calcDelta(qInL, qInR, delta, meqn);
          eq->rp(dir, delta, qInL, qInR, waves+inc*meqn*mwave, s+inc*mwave);
          copyComponents(waves+inc*meqn*mwave, limitedWaves+inc*meqn*mwave, meqn * mwave);

          if (linearIdx==localRange->volume()-1 && inc==1) {
            double *qOutL = qOut->getDataPtrAt(linearIdxL);
            double *qOutR = qOut->getDataPtrAt(linearIdxR);
            eq->qFluctuations(dir, qInL, qInR, waves+inc*meqn*mwave, s+inc*mwave, amdq, apdq);
            calcFirstOrderGud(dtdx, qOutL, qOutR, amdq, apdq, meqn);
          }
          idxL[dir] -= inc;
          idxR[dir] -= inc;
        }
      }
    }

    __syncthreads();
    if(linearIdx < localRange->volume()) {
      limitWaves(waves, s, limitedWaves, mwave, meqn);

      if (threadIdx.x==blockDim.x-1 || linearIdx==localRange->volume()-1) {
        limitWaves(
            waves+meqn*mwave, s+mwave, limitedWaves+meqn*mwave, mwave, meqn);
      }
    }

    __syncthreads();
    if(linearIdx < localRange->volume()) {
      for (int c = 0; c < meqn; c++) {
        flux[c] = 0;
      }
      secondOrderFlux(dtdx, s, limitedWaves, flux, meqn, mwave);

      if (threadIdx.x==blockDim.x-1 || linearIdx==localRange->volume()-1) {
        for (int c = 0; c < meqn; c++) {
          (flux+meqn)[c] = 0;
        }
        secondOrderFlux(
            dtdx, s+mwave, limitedWaves+meqn*mwave, flux+meqn, meqn, mwave);
      }
    }

    if(linearIdx < localRange->volume()) {
      //idxC[dir] -= 1;
      const int linearIdxCC = fIdxr.index(idxC);
      double *qOutCC = qOut->getDataPtrAt(linearIdxCC);
      secondOrderUpdate(dtdx, flux, flux+meqn, qOutCC, meqn);
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
  int sharedMemSize = 0;
  // numThreads == numRealCellsPerBlock
  // s & waves are needed on all real-real, real-ghost, and ghost-ghost faces
  sharedMemSize += (numThreads+3) * (mwave+mwave*meqn);
  // limitedWaves and 2nd-order flux are needed on real-real & real-ghost faces
  sharedMemSize += (numThreads+1) * (meqn+meqn);
  sharedMemSize *= sizeof(double);

  cudaFuncSetAttribute(
      cuda_WavePropagation, cudaFuncAttributeMaxDynamicSharedMemorySize,
      sharedMemSize);
  cuda_WavePropagation<<<numBlocks, numThreads, sharedMemSize>>>(
      hyper, qIn, qOut);
}

__global__ void setDtOnDevice(GkylWavePropagation_t *hyper, double dt) {
  hyper->dt = dt;
}

void setDt(GkylWavePropagation_t *hyper, double dt) {
  setDtOnDevice<<<1, 1>>>(hyper, dt);
}
