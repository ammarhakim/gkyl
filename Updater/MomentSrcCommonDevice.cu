#include <GkylMomentSrc.h>
#include <string.h>
#include <cstdio>

// Makes indexing cleaner
static const int X = 0;
static const int Y = 1;
static const int Z = 2;

static const int RHO = 0;
static const int MX = 1;
static const int MY = 2;
static const int MZ = 3;
static const int ER = 4;

static const int EX = 0;
static const int EY = 1;
static const int EZ = 2;
static const int BX = 3;
static const int BY = 4;
static const int BZ = 5;
static const int PHIE = 6;

#define fidx(n, c) (3 * (n) + (c))
#define eidx(c) (3 * nFluids + (c))

#define sq(x) ((x) * (x))

#define F2(base,i,j) (base)[(j)*(matSize)+(i)]

#define MOMENTSRC_SHM_LAYOUT SOA
#define SOA (1)
#define AOS (2)

#if MOMENTSRC_SHM_LAYOUT == AOS
#define SHAREDPAD (1)
#else
#define SHAREDPAD (0)
#endif

__global__ static void cuda_gkylMomentSrcTimeCenteredCublasSetPtrs(
   const int matSize,  double *d_lhs, double *d_rhs, double **d_lhs_ptr,
   double **d_rhs_ptr) {
  // numThreads*numBlocks == numRealCells
  const int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;
  d_lhs_ptr[linearIdx] = d_lhs + (sq(matSize))*linearIdx;
  d_rhs_ptr[linearIdx] = d_rhs + (matSize)*linearIdx;
}


// FIXME simplify by returning by reference?
// XXX don't init for direct solver!!!
GkylMomentSrcDeviceData_t *cuda_gkylMomentSrcInit(
    const char *scheme, const int nFluids, const int numBlocks,
    const int numThreads) {
  GkylMomentSrcDeviceData_t *context = new GkylMomentSrcDeviceData_t;
  context->scheme = scheme;
  context->nFluids = nFluids;
  context->numThreads = numThreads;
  context->numBlocks = numBlocks;

  if (strcmp(scheme, "time-centered")==0) {
    cublascall(cublasCreate(&(context->handle)));

    const int matSize = 3 * nFluids + 3;
    const int batchSize = numThreads*numBlocks;
    // device memory for actuall arrays and vectors
    cudacall(cudaMalloc(&context->d_lhs, batchSize*sq(matSize)*sizeof(double)));
    cudacall(cudaMalloc(&context->d_rhs, batchSize*matSize*sizeof(double)));
    cudacall(cudaMalloc(&context->d_info, batchSize*sizeof(int)));
    // device memory for pointers to the actuall arrays and vectors
    cudacall(cudaMalloc(&context->d_lhs_ptr, batchSize*sizeof(double*)));
    cudacall(cudaMalloc(&context->d_rhs_ptr, batchSize*sizeof(double*)));

    cuda_gkylMomentSrcTimeCenteredCublasSetPtrs<<<numBlocks, numThreads>>>(
        matSize, context->d_lhs, context->d_rhs, context->d_lhs_ptr,
        context->d_rhs_ptr);
  }

  return context;
}


void cuda_gkylMomentSrcDestroy(const GkylMomentSrcDeviceData_t *context) {
  if (strcmp(context->scheme, "time-centered")==0) {
    cudacall(cudaFree(context->d_lhs_ptr));
    cudacall(cudaFree(context->d_rhs_ptr));
    cudacall(cudaFree(context->d_lhs));
    cudacall(cudaFree(context->d_rhs));
    cudacall(cudaFree(context->d_info));
    cublascall(cublasDestroy(context->handle));

    delete(context);
  }
}


__device__ static void cuda_gkylMomentSrcTimeCenteredUpdateRhovE(
    const int linearIdx, const int linearIdxC, const MomentSrcData_t *sd,
    const FluidData_t *fd, const double dt, GkylCartField_t **fluidFlds,
    GkylCartField_t *emFld, double **d_rhs_ptr) {
  const double *sol = d_rhs_ptr[linearIdx];
  const int nFluids = sd->nFluids;

  //------------> update solution for fluids
  double chargeDens = 0.0;
  for (int n=0; n<nFluids; ++n)
  {
    double *f = fluidFlds[n]->getDataPtrAt(linearIdxC);
    double qbym = fd[n].charge/fd[n].mass;

    chargeDens += qbym*f[RHO];

    f[MX] = 2*sol[fidx(n,X)]/qbym - f[MX];
    f[MY] = 2*sol[fidx(n,Y)]/qbym - f[MY];
    f[MZ] = 2*sol[fidx(n,Z)]/qbym - f[MZ];
  }

  //------------> update electric field
  double *em = emFld->getDataPtrAt(linearIdxC);
  em[EX] = 2*sol[eidx(X)] - em[EX];
  em[EY] = 2*sol[eidx(Y)] - em[EY];
  em[EZ] = 2*sol[eidx(Z)] - em[EZ];

  //------------> update correction potential
  const double crhoc = sd->chi_e*chargeDens/sd->epsilon0;
  em[PHIE] += dt*crhoc;
}


__global__ static void cuda_gkylMomentSrcTimeCenteredCublasUpdate(
    const MomentSrcData_t *sd, const FluidData_t *fd, const double dt,
    GkylCartField_t **fluidFlds, GkylCartField_t *emFld, double **d_rhs_ptr) {
  GkylRange_t *localRange = emFld->localRange;
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = emFld->genIndexer();

  // numThreads*numBlocks == numRealCells
  const int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;
  int idxC[3];
  localIdxr.invIndex(linearIdx, idxC);
  const int linearIdxC = fIdxr.index(idxC);
  
  const int nFluids = sd->nFluids;

  double keOld[2]; // XXX
  for (int n = 0; n < nFluids; ++n)
  {
    double *f = fluidFlds[n]->getDataPtrAt(linearIdxC);
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  cuda_gkylMomentSrcTimeCenteredUpdateRhovE(
      linearIdx, linearIdxC, sd, fd, dt, fluidFlds, emFld, d_rhs_ptr);

  if (sd->hasPressure)
  {
    for (int n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = fluidFlds[n]->getDataPtrAt(linearIdxC);
      const double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 
}


__global__ static void cuda_gkylMomentSrcTimeCenteredCublasSetMat(
    const MomentSrcData_t *sd, const FluidData_t *fd, const double dt,
    GkylCartField_t **fluidFlds, GkylCartField_t *emFld, double **d_lhs_ptr,
    double **d_rhs_ptr) {
  GkylRange_t *localRange = emFld->localRange;
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = emFld->genIndexer();

  const int nFluids = sd->nFluids;
  const int matSize = 3 * nFluids + 3;

  // numThreads*numBlocks == numRealCells
  const int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;
  int idxC[3];
  localIdxr.invIndex(linearIdx, idxC);
  const int linearIdxC = fIdxr.index(idxC);
  const double *em = emFld->getDataPtrAt(linearIdxC);

  double *lhs = d_lhs_ptr[linearIdx];
  double *rhs = d_rhs_ptr[linearIdx];

  for (int c=0; c<sq(matSize); c++)
    lhs[c] = 0;

  const double dt1 = 0.5 * dt;
  const double dt2 = 0.5 * dt / sd->epsilon0;

  for (int n=0; n<nFluids; ++n)
  {
    double qbym = fd[n].charge/fd[n].mass;
    double qbym2 = sq(qbym);

    const double *f = fluidFlds[n]->getDataPtrAt(linearIdxC);
    if (fd[n].evolve) {
      // off-diagonal elements of lhs
      // eqn. for X-component of current
      F2(lhs, fidx(n,X), fidx(n,Y)) = -dt1*qbym*(em[BZ]);
      F2(lhs, fidx(n,X), fidx(n,Z)) = dt1*qbym*(em[BY]);
      F2(lhs, fidx(n,X), eidx(X)) = -dt1*qbym2*f[RHO];

      // eqn. for Y-component of current
      F2(lhs, fidx(n,Y), fidx(n,X)) = dt1*qbym*(em[BZ]);
      F2(lhs, fidx(n,Y), fidx(n,Z)) = -dt1*qbym*(em[BX]);
      F2(lhs, fidx(n,Y), eidx(Y)) = -dt1*qbym2*f[RHO];

      // eqn. for Z-component of current
      F2(lhs, fidx(n,Z), fidx(n,X)) = -dt1*qbym*(em[BY]);
      F2(lhs, fidx(n,Z), fidx(n,Y)) = dt1*qbym*(em[BX]);
      F2(lhs, fidx(n,Z), eidx(Z)) = -dt1*qbym2*f[RHO];
    }
    // diagonal elements of lhs
    F2(lhs, fidx(n,X), fidx(n,X)) = 1.0;
    F2(lhs, fidx(n,Y), fidx(n,Y)) = 1.0;
    F2(lhs, fidx(n,Z), fidx(n,Z)) = 1.0;

    // fill corresponding RHS elements
    rhs[fidx(n,X)] = qbym*f[MX];
    rhs[fidx(n,Y)] = qbym*f[MY];
    rhs[fidx(n,Z)] = qbym*f[MZ];

    // set current contribution to electric field equation
    F2(lhs, eidx(X), fidx(n,X)) = dt2;
    F2(lhs, eidx(Y), fidx(n,Y)) = dt2;
    F2(lhs, eidx(Z), fidx(n,Z)) = dt2;
  }

  // fill in elements for electric field equations
  F2(lhs, eidx(EX), eidx(EX)) = 1.0;
  F2(lhs, eidx(EY), eidx(EY)) = 1.0;
  F2(lhs, eidx(EZ), eidx(EZ)) = 1.0;

  rhs[eidx(EX)] = em[EX];
  rhs[eidx(EY)] = em[EY];
  rhs[eidx(EZ)] = em[EZ];
}


static void cuda_gkylMomentSrcTimeCenteredCublas(
    const int nFluids, int numBlocks, int numThreads, const MomentSrcData_t *sd, const FluidData_t *fd,
    const double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld,
    const GkylMomentSrcDeviceData_t *context) {
  const int matSize = 3 * nFluids + 3;

  double **d_lhs_ptr = context->d_lhs_ptr;
  double **d_rhs_ptr = context->d_rhs_ptr;
  int *d_info = context->d_info;
  const cublasHandle_t &handle = context->handle;

  int batchSize = numThreads*numBlocks;

  cuda_gkylMomentSrcTimeCenteredCublasSetMat<<<numBlocks, numThreads>>>(
      sd, fd, dt, fluidFlds, emFld, d_lhs_ptr, d_rhs_ptr);

  cublascall(cublasDgetrfBatched(
      handle,
      matSize,  // n
      d_lhs_ptr,  // A
      matSize,  // lda
      NULL,  // int *PivotArray
      d_info,  // int *infoArray
      batchSize // number of pointers contained in A
      ));

  int info;
  cublascall(cublasDgetrsBatched(
      handle,
      CUBLAS_OP_N,  // trans
      matSize,  // n
      1,  // nrhs
      d_lhs_ptr,  // matrix A
      matSize,  // lda
      NULL,  // const int *devIpiv
      d_rhs_ptr,  // double *Barray[]
      matSize,  // ldb
      &info,  // int *info
      batchSize // number of pointers contained in A
      ));

  // update solution
  cuda_gkylMomentSrcTimeCenteredCublasUpdate<<<numBlocks, numThreads>>>(
      sd, fd, dt, fluidFlds, emFld, d_rhs_ptr);
}


__device__ static void cuda_gkylMomentSrcTimeCenteredDirectUpdateRhovE(
    const MomentSrcData_t *sd, const FluidData_t *fd, const double dt,
    double **ff, double *em) {
  const int nFluids = sd->nFluids;
  const double epsilon0 = sd->epsilon0;

  const double Bx = (em[BX]);
  const double By = (em[BY]);
  const double Bz = (em[BZ]);
  const double Bmag = sqrt(Bx*Bx + By*By + Bz*Bz);
  double b[] = {0, 0, 0};
  if (Bmag > 0)
  {
    b[0] = Bx / Bmag;
    b[1] = By / Bmag;
    b[2] = Bz / Bmag;
  }

  extern __shared__ double dummy[];

  int base = 0;
  double *qbym = dummy + base;
  base += nFluids;

  if (threadIdx.x==0)
  {
    for (int n=0; n<nFluids; ++n)
    {
      qbym[n] = fd[n].qbym;
    }
  }
  __syncthreads();

#if MOMENTSRC_SHM_LAYOUT == AOS
  double *ptr = dummy + base + threadIdx.x * (nFluids * (3 + 1 + 1) + SHAREDPAD);

  double *JJ = ptr;
  ptr += nFluids*3;

  double *Wc_dt = ptr;
  ptr += nFluids;

  double *wp_dt2 = ptr;
  ptr += nFluids;
#else
  double *JJ = dummy + base + threadIdx.x * (nFluids * 3);
  base += blockDim.x * nFluids * 3;

  double *Wc_dt = dummy + base + threadIdx.x * nFluids;
  base += blockDim.x * nFluids;

  double *wp_dt2 = dummy + base + threadIdx.x * nFluids;
  base += blockDim.x * nFluids;
#endif

  double K[] = {0, 0, 0};
  double w02 = 0.;
  double gam2 = 0.;
  double delta = 0.;

  for (int n=0; n < nFluids; ++n)
  {
    const double *f = ff[n];
    double *J = JJ+n*3;
    J[0] = f[MX] * qbym[n];
    J[1] = f[MY] * qbym[n];
    J[2] = f[MZ] * qbym[n];
    if (!fd[n].evolve)
      continue;
    Wc_dt[n] = qbym[n] * Bmag * dt;
    wp_dt2[n] = f[RHO] * sq(qbym[n]) / epsilon0 * sq(dt);
    const double tmp = 1. / (1. + sq(Wc_dt[n]) / 4.);
    w02 += wp_dt2[n] * tmp;
    gam2 += wp_dt2[n] * sq(Wc_dt[n]) * tmp;
    delta += wp_dt2[n] * Wc_dt[n] * tmp;

    const double bDotJ = b[0]*J[0] + b[1]*J[1] + b[2]*J[2];
    const double bCrossJ[] = {
      b[1]*J[2]-b[2]*J[1], // by*Jz-bz*Jy
      b[2]*J[0]-b[0]*J[2], // bz*Jx-bx*Jz
      b[0]*J[1]-b[1]*J[0], // bx*Jy-by*Jx
    };

#pragma unroll
    for(int c=0; c<3; c++) {
      K[c] -= dt * tmp * (J[c] + sq(Wc_dt[n] / 2.) * b[c] * bDotJ
              - (Wc_dt[n] / 2.) * bCrossJ[c]);
    }
  }
  const double Delta2 = sq(delta) / (1. + w02 / 4.);

  const double F[] = {em[EX] * epsilon0, em[EY] * epsilon0, em[EZ] * epsilon0};
  double Fbar[3];
  {
    double F_halfK[3];
#pragma unroll
    for (int c=0; c<3; c++) {
      F_halfK[c] = F[c] + 0.5 * K[c];
      for (int n=0; n < nFluids; ++n)
      {
        if (fd[n].evolve)
          continue;
        F_halfK[c] -= (0.5 * dt) * JJ[n*3+c]; // from the un-evolved fluid
      }
    }

    const double tmp = 1. / (1. + w02 / 4. + Delta2 / 64.);
    const double bDotF_halfK = b[0]*F_halfK[0] + b[1]*F_halfK[1] + b[2]*F_halfK[2];
    const double bCrossF_halfK[] = {
      b[1]*F_halfK[2]-b[2]*F_halfK[1], // by*Fz-bz*Fy
      b[2]*F_halfK[0]-b[0]*F_halfK[2], // bz*Fx-bx*Fz
      b[0]*F_halfK[1]-b[1]*F_halfK[0], // bx*Fy-by*Fx
    };

#pragma unroll
    for (int c=0; c<3; c++) {
      Fbar[c] = tmp * (
        F_halfK[c]
        + ((Delta2 / 64. - gam2 / 16.) / (1. + w02 / 4. + gam2 / 16.))
           * b[c] * bDotF_halfK
        + (delta / 8. / (1. + w02 / 4.)) * bCrossF_halfK[c]
        );
    } 
  }

  //------------> update electric field
  em[EX] = (2*Fbar[0]-F[0]) / epsilon0;
  em[EY] = (2*Fbar[1]-F[1]) / epsilon0;
  em[EZ] = (2*Fbar[2]-F[2]) / epsilon0;

  //------------> update solution for fluids
  double chargeDens = 0.0;
  for (int n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    chargeDens += qbym[n] * f[RHO];
    if (!fd[n].evolve)
      continue;

    const double *J = JJ+n*3;
    double Jstar[3];
    double J_new[3];

#pragma unroll
    for (int c=0; c<3; c++) {
      Jstar[c] = J[c] + Fbar[c] * (wp_dt2[n] / dt / 2.);
    }
    const double bDotJstar = b[0]*Jstar[0] + b[1]*Jstar[1] + b[2]*Jstar[2];
    const double bCrossJstar[] = {
      b[1]*Jstar[2]-b[2]*Jstar[1], // by*Jz-bz*Jy
      b[2]*Jstar[0]-b[0]*Jstar[2], // bz*Jx-bx*Jz
      b[0]*Jstar[1]-b[1]*Jstar[0], // bx*Jy-by*Jx
    };

#pragma unroll
    for (int c=0; c<3; c++) {
      J_new[c] = 2. * (Jstar[c] + sq(Wc_dt[n] / 2.) * b[c] * bDotJstar
                 - (Wc_dt[n] / 2.) * bCrossJstar[c]) / (1. + sq(Wc_dt[n] / 2.))
                 - J[c];
    }

    f[MX] = J_new[0] / qbym[n];
    f[MY] = J_new[1] / qbym[n];
    f[MZ] = J_new[2] / qbym[n];
  } 

  //------------> update correction potential
  em[PHIE] += dt * sd->chi_e * chargeDens/epsilon0;
}


__global__ static void cuda_gkylMomentSrcTimeCenteredDirect(
    int numBlocks, int numThreads, const MomentSrcData_t *sd, const FluidData_t *fd,
    const double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld,
    const GkylMomentSrcDeviceData_t *context) {
  GkylRange_t *localRange = emFld->localRange;
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = emFld->genIndexer();

  // numThreads*numBlocks == numRealCells
  const int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;
  int idxC[3];
  localIdxr.invIndex(linearIdx, idxC);
  const int linearIdxC = fIdxr.index(idxC);

  const int nFluids = sd->nFluids;

  double *em = emFld->getDataPtrAt(linearIdxC);
  double *ff[2]; // XXX

  double keOld[2]; // XXX
  for (int n = 0; n < nFluids; ++n)
  {
    double *f = fluidFlds[n]->getDataPtrAt(linearIdxC);
    if (!fd[n].evolve)
      continue;
    keOld[n] = (sq(f[MX]) + sq(f[MY]) + sq(f[MZ]));
    ff[n] = f;
  }

  cuda_gkylMomentSrcTimeCenteredDirectUpdateRhovE(sd, fd, dt, ff, em);

  if (sd->hasPressure)
  {
    for (int n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = ff[n];
      const double keNew = (sq(f[MX]) + sq(f[MY]) + sq(f[MZ]));
      f[ER] += 0.5 * (keNew - keOld[n]) / f[RHO];
    }
  } 
}


void momentSrcAdvanceOnDevice(
    const MomentSrcData_t *sd, const FluidData_t *fd, const double dt,
    GkylCartField_t **fluidFlds, GkylCartField_t *emFld,
    const GkylMomentSrcDeviceData_t *context)
{
  const int nFluids = context->nFluids;
  const int numThreads = context->numThreads;
  const int numBlocks = context->numBlocks;

  if (strcmp(context->scheme, "time-centered")==0) {
    cuda_gkylMomentSrcTimeCenteredCublas(
        nFluids, numBlocks, numThreads, sd, fd, dt, fluidFlds, emFld, context);
  } else if (strcmp(context->scheme, "time-centered-direct")==0
             || strcmp(context->scheme, "direct")==0) {
    int sharedMemSize = 0;

#if MOMENTSRC_SHM_LAYOUT == AOS
    // qbym, J, Wc_dt, wp_dt2
    sharedMemSize += numThreads * (nFluids * (3 + 1 + 1) + SHAREDPAD) + nFluids;
#else
    // J, Wc_dt, wp_dt2 and shared qbym
    sharedMemSize += numThreads * (nFluids * (3 + 1 + 1) + SHAREDPAD) + nFluids;
#endif
    sharedMemSize *= sizeof(double);

    cudaFuncSetSharedMemConfig(cuda_gkylMomentSrcTimeCenteredDirect,
        cudaSharedMemBankSizeEightByte);

    cuda_gkylMomentSrcTimeCenteredDirect
      <<<numBlocks, numThreads, sharedMemSize>>>(
      numBlocks, numThreads, sd, fd, dt, fluidFlds, emFld, context);
  }
}

#undef SHAREDPAD
#undef MOMENTSRC_SHM_LAYOUT
#undef SOA
#undef AOS
