#include <cstdio>
#include <GkylWavePropagation.h>

__device__ inline  static void calcFirstOrderGud(
    const double dtdx, double *ql, double *qr, const double *amdq,
    const double *apdq, const int meqn)
{
  for (int i = 0; i < meqn; i++)
  {
    /* qr[i] -= dtdx * apdq[i]; */
    /* ql[i] -= dtdx * amdq[i]; */
    // XXX
    atomicAdd(qr+i, -dtdx * apdq[i]);
    atomicAdd(ql+i, -dtdx * amdq[i]);
  }
}

__device__ static inline double calcCfla(
    const double cfla, const double dtdx, const double *speeds, const int mwave)
{
  double c = cfla;
  for (int i = 0; i < mwave; i ++)
  {
    c = max(c, dtdx * abs(speeds[i]));
  }
  return c;
}

__device__ static inline double waveDotProd(
    const double *waves, const double *waves1, const int mw, const int meqn)
{
  double result = 0.;
  for (int i = 0; i < meqn; i++)
  {
    result += waves[meqn*mw+i] * waves1[meqn*mw+i];
  }
  return result;
}

__device__ static inline double limiter_minMod(const double r)
{
   return max(0., min(1., r));
}

__device__ static void limitWaves(
    const double *waves, const double *speeds, double *limitedWaves,
    const int mwave, const int meqn)
{
  int jump = meqn * mwave;
  for (int mw = 0; mw < mwave; mw++ )
  {
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

__device__ static inline void secondOrderFluxOneWave(
    const double dtdx, const double speed, const double *wave, double *fs,
    const int meqn) {
  double sfact = 0.5 * abs(speed) * (1 - abs(speed) * dtdx);
  for (int i = 0; i < meqn; i++)
  {
    fs[i] += sfact * wave[i];
  }
}

__device__ static void secondOrderFlux(
    const double dtdx, const double *speeds, const double *waves, double *fs,
    const int meqn, const int mwave)
{
    for (int mw = 0; mw < mwave; mw++)
    {
      secondOrderFluxOneWave(dtdx, speeds[mw], waves+mw*meqn, fs, meqn);
    }
}

__device__ static inline void secondOrderUpdate(
    const double dtdx, const double *fs, const double *fs1, double *q,
    const int meqn)
{
  for (int i = 0; i < meqn; i++)
  {
    // q[i] -= dtdx * (fs1[i] - fs[i]);
    atomicAdd(q+i, -dtdx * (fs1[i] - fs[i]));
  }
}

__device__ static inline void copyComponents(
    const double *ptrFrom, double *ptrTo, const int nComponents)
{
  for (int i = 0; i < nComponents; i++)
  {
    ptrTo[i] = ptrFrom[i];
  }
}

__global__ void cuda_WavePropagation(
    GkylWavePropagation_t *hyper, GkylCartField_t *qIn, GkylCartField_t *qOut)
{
  GkylRange_t *localRange = qIn->localRange;
  Gkyl::GenIndexer fIdxr = qIn->genIndexer();
  Gkyl::GenIndexer localIdxr(localRange);
  GkylEquationFv_t *eq = hyper->equation;
  GkylCartField_t *dtByCell = hyper->dtByCell;
  const int meqn = eq->numEquations;
  const int mwave = eq->numWaves;

  const int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;

  // assign buffer space for different usages
  extern __shared__ double dummy[];
  int base = 0;

  double *amdq = dummy + base + (meqn)*(threadIdx.x);
  base += (meqn) * (blockDim.x);

  double *apdq = dummy + base + (meqn)*(threadIdx.x);
  base += (meqn) * (blockDim.x);

  // numThreads == numRealCells
  // waves and speeds are defined on ghost-ghost, ghost-real, and
  // real-real faces, thus the jump is blockDim.x+3;
  // also, the 0th thread will work on one more Rp on its left, thus the address
  // for its own waves/speeds is threadIdx.x+1 (not threadIdx.x)
  double *waves = dummy + base + (meqn*mwave)*(threadIdx.x+1);
  base += (meqn*mwave) * (blockDim.x+3);

  double *speeds = dummy + base + (mwave)*(threadIdx.x+1);
  base += (mwave) * (blockDim.x+3);

  // limitedWaves and second-order fluxs are defiend on ghost-real and
  // real-real faces, thus the jump is blockDim.x+1
  double *limitedWaves = dummy + base + (meqn*mwave) * (threadIdx.x);
  base += (meqn*mwave) * (blockDim.x+1);

  base = 0; // reusing the shared memory since amdq/apdq are not longer needed

  double *flux = dummy + base + (meqn)*(threadIdx.x);
  base += (meqn) * (blockDim.x+1);

  int idxC[3];
  int idxL[3];
  int idxR[3];

  // get i,j,k... index idxC; real cells have idxC=1,2,...,shape
  localIdxr.invIndex(linearIdx, idxC);
  // linearIdxC knows ghots cells, the 1st real cell has linearIdxC=2
  // if ndim>1, linearIdxC jumps (i.e., not consecutive) due to ghost cells
  const int linearIdxC = fIdxr.index(idxC);

  double cfl = hyper->_cfl;
  double cflm = hyper->_cflm;
  double cfla = 0; // actual CFL number used

  const double *qInC = qIn->getDataPtrAt(linearIdxC);
  double *qOutC = qOut->getDataPtrAt(linearIdxC);

  // ghost cells are not copied, but this is OK because the waves and speeds are
  // computed using qIn anyway
  if(linearIdx < localRange->volume())
  {
    copyComponents(qInC, qOutC, meqn);
  }

  for(int i=0; i<hyper->numUpdateDirs; i++)
  {
    const int dir = hyper->updateDirs[i] - 1;
    const double dtdx = hyper->dt / qIn->grid->dx[dir];

#pragma unroll
    for(int d=0; d<3; d++)
    {
      idxL[d] = idxC[d];
      idxR[d] = idxC[d];
    }
    idxL[dir] = idxC[dir]-1;

    const bool isLoBlockOrDevice = threadIdx.x==0 || idxC[dir]==1;
    const bool isHiDevice = idxC[dir]==localRange->shape(dir);
    const bool isHiBlockOrDevice = threadIdx.x==blockDim.x-1 || isHiDevice;

    const int linearIdxL = fIdxr.index(idxL);
    const int linearIdxR = fIdxr.index(idxR);
    const double *qInL = qIn->getDataPtrAt(linearIdxL);
    const double *qInR = qIn->getDataPtrAt(linearIdxR);

    double *qOutL = qOut->getDataPtrAt(linearIdxL);
    double *qOutR = qOut->getDataPtrAt(linearIdxR);

    if(linearIdx < localRange->volume())
    {
      eq->rp(dir, qInL, qInR, waves, speeds);
      eq->qFluctuations(dir, qInL, qInR, waves, speeds, amdq, apdq);
      calcFirstOrderGud(dtdx, qOutL, qOutR, amdq, apdq, meqn);
      cfla = calcCfla(cfla, dtdx, speeds, mwave);
      copyComponents(waves, limitedWaves, meqn*mwave);

      // can we avoid branching?
      // solve one additional Riemann problem on the lower side
      if (isLoBlockOrDevice)
      {
        int inc = -1;
        idxL[dir] += inc;
        idxR[dir] += inc;
        const int linearIdxL = fIdxr.index(idxL);
        const int linearIdxR = fIdxr.index(idxR);
        const double *qInL = qIn->getDataPtrAt(linearIdxL);
        const double *qInR = qIn->getDataPtrAt(linearIdxR);

        eq->rp(dir, qInL, qInR, waves+inc*meqn*mwave, speeds+inc*mwave);
        cfla = calcCfla(cfla, dtdx, speeds+inc*mwave, mwave);

        idxL[dir] -= inc;
        idxR[dir] -= inc;
      }

      // solve two additional Riemann problems on the higher side and update the
      // last real cell
      if (isHiBlockOrDevice)
      {
        for (int inc = 1; inc < 3; inc++)
        {
          idxL[dir] += inc;
          idxR[dir] += inc;
          const int linearIdxL = fIdxr.index(idxL);
          const int linearIdxR = fIdxr.index(idxR);
          const double *qInL = qIn->getDataPtrAt(linearIdxL);
          const double *qInR = qIn->getDataPtrAt(linearIdxR);

          eq->rp(dir, qInL, qInR, waves+inc*meqn*mwave, speeds+inc*mwave);
          cfla = calcCfla(cfla, dtdx, speeds+inc*mwave, mwave);

          if (inc==1)
          {
            copyComponents(
                waves+inc*meqn*mwave, limitedWaves+inc*meqn*mwave, meqn*mwave);

            if (isHiDevice)
            {
              double *qOutL = qOut->getDataPtrAt(linearIdxL);
              double *qOutR = qOut->getDataPtrAt(linearIdxR);
              eq->qFluctuations(
                  dir, qInL, qInR, waves+inc*meqn*mwave, speeds+inc*mwave, amdq,
                  apdq);
              calcFirstOrderGud(dtdx, qOutL, qOutR, amdq, apdq, meqn);
            }
          }

          idxL[dir] -= inc;
          idxR[dir] -= inc;
        }
      }
    }

    __syncthreads();
    if(linearIdx < localRange->volume())
    {
      limitWaves(waves, speeds, limitedWaves, mwave, meqn);

      if (isHiBlockOrDevice)
      {
        limitWaves(
            waves+meqn*mwave, speeds+mwave, limitedWaves+meqn*mwave, mwave,
            meqn);
      }
    }

    __syncthreads();
    if(linearIdx < localRange->volume())
    {
      for (int c=0; c<meqn; c++)
      {
        flux[c] = 0;
      }
      secondOrderFlux(dtdx, speeds, limitedWaves, flux, meqn, mwave);

      if (isHiBlockOrDevice) {
        for (int c=0; c<meqn; c++)
        {
          (flux+meqn)[c] = 0;
        }
        secondOrderFlux(
            dtdx, speeds+mwave, limitedWaves+meqn*mwave, flux+meqn, meqn,
            mwave);
      }
    }

    __syncthreads();
    if(linearIdx < localRange->volume())
    {
      secondOrderUpdate(dtdx, flux, flux+meqn, qOutC, meqn);
    }
  }

  dtByCell->getDataPtrAt(linearIdxC)[0] = hyper->dt * cfl/cfla;
}

static int calcSharedMemSize(
    const int meqn, const int mwave, const int numThreads)
{
  int sharedMemSize = 0;
  // numThreads == numRealCellsPerBlock
  // amdq and apdq
  sharedMemSize += (numThreads) * (meqn+meqn);
  // waves & speeds are needed on all real-real, real-ghost, and ghost-ghost
  // cell faces
  sharedMemSize += (numThreads+3) * (mwave*meqn+mwave);
  // limitedWaves and 2nd-order flux are needed on real-real & real-ghost faces
  sharedMemSize += (numThreads+1) * (mwave*meqn+meqn*0);
  sharedMemSize *= sizeof(double);
  return sharedMemSize;
}

void wavePropagationAdvanceOnDevice(
    const int meqn, const int mwave, const int numBlocks, const int numThreads,
    GkylWavePropagation_t *hyper, GkylCartField_t *qIn, GkylCartField_t *qOut)
{
  const int sharedMemSize = calcSharedMemSize(
      meqn, mwave, numThreads);

  cudaFuncSetAttribute(
      cuda_WavePropagation, cudaFuncAttributeMaxDynamicSharedMemorySize,
      sharedMemSize);
  cuda_WavePropagation<<<numBlocks, numThreads, sharedMemSize>>>(
      hyper, qIn, qOut);
}

void wavePropagationInitOnDevice(
    const int meqn, const int mwave, const int numBlocks, const int numThreads,
    GkylWavePropagation_t *hyper)
{
  const int sharedMemSize = calcSharedMemSize(meqn, mwave, numThreads);
  cudacall(cudaMalloc(&hyper->buf, sharedMemSize*numBlocks));
  hyper->sharedMemSize = sharedMemSize;
}

void wavePropagationDestroyOnDevice(GkylWavePropagation_t *hyper) {
    cudacall(cudaFree(hyper->buf));
}

__global__ void setDtOnDevice(GkylWavePropagation_t *hyper, double dt)
{
  hyper->dt = dt;
}

void setDt(GkylWavePropagation_t *hyper, double dt)
{
  setDtOnDevice<<<1, 1>>>(hyper, dt);
}
