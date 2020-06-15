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
   // return 1;
}

__device__ static void limitWaves(
    const double *waves, const double *speeds, double *limitedWaves,
    const int mwave, const int meqn) {
  int jump = meqn * mwave;
  for (int mw = 0; mw < mwave; mw++ ){
    const double wnorm2 = waveDotProd(waves, waves, mw, meqn);
    double wlimitr = 1.;
    if (wnorm2 > 0) {
      const double dotl = waveDotProd(waves-jump, waves, mw, meqn);
      const double dotr = waveDotProd(waves+jump, waves, mw, meqn);
      const double r = speeds[mw] > 0 ? dotl/wnorm2 : dotr/wnorm2;
      wlimitr = limiter_minMod(r);

/* printf("[%2d, %2d] r %13g limitr %13g dotl %13g dotr %13g wnorm2 %13g waves %13g %13g; wavesr %13g %13g\n", */
/*     blockIdx.x, threadIdx.x, r, wlimitr, dotl, dotr, wnorm2, */
/*     waves[0], waves[4], */
/*     (waves+jump)[0], (waves+jump)[4] */
/*     ); */
    }
    for (int me = 0; me < meqn; me++) {
      limitedWaves[mw*meqn+me] *= wlimitr;
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
  GkylRange_t *localEdgeRange = qIn->localEdgeRange;
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
  double delta[5];
  double amdq[5];
  double apdq[5];

  // declaring this dummy array shared seems to alleviate register pressure and
  // improve performance a bit
  extern __shared__ double dummy[];
  int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;

  // assign buffer space for different usages
  int base = 0;

  // numThreads == numRealEdges (real-real and real-ghost faces)
  // waveSlice and speedSlice contain ghost-ghost faces, thus the +2
  const int baseWaveSlice = base;
  base += (meqn * mwave) * (blockDim.x + 3);
  double *waveSlice = dummy + baseWaveSlice;

  const int baseSpeedSlice = base;
  base += (mwave) * (blockDim.x + 3);
  double *speedSlice = dummy + baseSpeedSlice;

  const int baseLimitedWaveSlice = base;
  base += (meqn * mwave) * (blockDim.x + 1);
  double *limitedWaveSlice = dummy + baseLimitedWaveSlice;

  const int baseFluxSlice = base;
  base += (meqn) * (blockDim.x + 1);
  double *fluxSlice = dummy + baseFluxSlice;

  // find buffer addresses for each thread
  // FIXME shall waves and s be created on the fly and then copied into slices
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

  // XXX still missing ghost cells
  /* if(linearIdx < localEdgeRange->volume()) { */
  /*   for(int i = 0; i < meqn; i++) { */
  /*     qOutC[i] = qInC[i]; */
  /*   } */
  /* } */

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
        // printf("[%2d, %2d] inc %d\n", blockIdx.x, threadIdx.x, inc);
      }
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

          if (inc==1 && linearIdx==localRange->volume()-1) {
            double *qOutL = qOut->getDataPtrAt(linearIdxL);
            double *qOutR = qOut->getDataPtrAt(linearIdxR);
            eq->qFluctuations(dir, qInL, qInR, waves+inc*meqn*mwave, s+inc*mwave, amdq, apdq);
            calcFirstOrderGud(dtdx, qOutL, qOutR, amdq, apdq, meqn);
          }
          idxL[dir] -= inc;
          idxR[dir] -= inc;
          // printf("[%2d, %2d] inc %d\n", blockIdx.x, threadIdx.x, inc);
        }
      }
    }

    /* __syncthreads(); */
    /* if(linearIdx < localEdgeRange->volume()) */
    /* printf("[%2d; %2d, %2d] after  L [%2d] %13g, R [%2d] %13g; wave %13g %13g; amdq %13g; qOutL %13g\n", */
    /*     linearIdxC, blockIdx.x, threadIdx.x, idxL[dir], qInL[0], idxR[dir], qInR[0], waves[0], waves[4], amdq[0], qOutL[0]); */

    if(linearIdx < localRange->volume()) {
      limitWaves(waves, s, limitedWaves, mwave, meqn);
      /* if (threadIdx.x==0) { */
      /*   int inc = -1; */
      /*   limitWaves(waves+inc*meqn*mwave, s+inc*mwave, limitedWaves+inc*meqn*mwave, mwave, meqn); */
      /* } */
      if (threadIdx.x==blockDim.x-1 || linearIdx==localRange->volume()-1) {
        int inc = +1;
        limitWaves(waves+inc*meqn*mwave, s+inc*mwave, limitedWaves+inc*meqn*mwave, mwave, meqn);
      }
    }

    __syncthreads();
    /* if(linearIdx < localEdgeRange->volume()) */
    /* printf("[%2d] limited wave %13g %13g %13g %13g %13g\n", */
    /*     linearIdx+1, limitedWaves[0], */
    /*     limitedWaves[1], */
    /*     limitedWaves[2], */
    /*     limitedWaves[3], */
    /*     limitedWaves[4]); */
    /* if (threadIdx.x==blockDim.x-1 || linearIdx==localRange->volume()-1) { */
    /* printf("[%2d] limited wave %13g %13g %13g %13g %13g\n", */
    /*     linearIdx+2, (limitedWaves+meqn*mwave)[0], */
    /*     (limitedWaves+meqn*mwave)[1], */
    /*     (limitedWaves+meqn*mwave)[2], */
    /*     (limitedWaves+meqn*mwave)[3], */
    /*     (limitedWaves+meqn*mwave)[4]); */
    /*  */
    /* } */

    if(linearIdx < localRange->volume()) {
      for (int c = 0; c < meqn; c++) {
        flux[c] = 0;
      }

      for (int mw = 0; mw < mwave; mw++) {
        secondOrderFlux(dtdx, s[mw], limitedWaves+mw*meqn, flux, meqn);
      }

      if (threadIdx.x==blockDim.x-1 || linearIdx==localRange->volume()-1) {
        for (int c = 0; c < meqn; c++) {
          (flux+meqn)[c] = 0;
        }

        for (int mw = 0; mw < mwave; mw++) {
          secondOrderFlux(dtdx, (s+mwave)[mw], limitedWaves+meqn*mwave+mw*meqn, (flux+meqn), meqn);
        }
      }
    }

    /* __syncthreads(); */
    /* if(linearIdx < localEdgeRange->volume()) */
    /* printf("[%2d] fluxes %13g %13g %13g %13g %13g\n", */
    /*     linearIdx+1, flux[0], */
    /*     flux[1], */
    /*     flux[2], */
    /*     flux[3], */
    /*     flux[4]); */

    if(linearIdx < localRange->volume()) {
      //idxC[dir] -= 1;
      const int linearIdxCC = fIdxr.index(idxC);
      double *qOutCC = qOut->getDataPtrAt(linearIdxCC);
    /* printf("[%2d] cpp before fl %13g %13g fr %13g %13g; qOut %13g %13g %13g %13g %13g\n", */
    /*     linearIdx+1, */
    /*     flux[0], */
    /*     flux[4], */
    /*     (flux+meqn)[0], */
    /*     (flux+meqn)[4], */
    /*     qOutCC[0], */
    /*     qOutCC[1], */
    /*     qOutCC[2], */
    /*     qOutCC[3], */
    /*     qOutCC[4] */
    /*     ); */
      secondOrderUpdate(dtdx, flux, flux+meqn, qOutCC, meqn);
    /* printf("[%2d] after fl %13g %13g fr %13g %13g; qOut %13g %13g %13g %13g %13g\n", */
    /*     linearIdx+1, */
    /*     flux[0], */
    /*     flux[4], */
    /*     (flux+meqn)[0], */
    /*     (flux+meqn)[4], */
    /*     qOutCC[0], */
    /*     qOutCC[1], */
    /*     qOutCC[2], */
    /*     qOutCC[3], */
    /*     qOutCC[4] */
    /*     ); */

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
  // numThreads == numRealEdgesPerBlock == numRealCellsPerBlock + 1
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
