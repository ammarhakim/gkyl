#include <GkylMomentSrc.h>
#include <cstdio>

// Makes indexing cleaner
static const unsigned X = 0;
static const unsigned Y = 1;
static const unsigned Z = 2;

static const unsigned RHO = 0;
static const unsigned MX = 1;
static const unsigned MY = 2;
static const unsigned MZ = 3;
static const unsigned ER = 4;

static const unsigned EX = 0;
static const unsigned EY = 1;
static const unsigned EZ = 2;
static const unsigned BX = 3;
static const unsigned BY = 4;
static const unsigned BZ = 5;
static const unsigned PHIE = 6;
static const unsigned PHIM = 7;

#define fidx(n, c) (3 * (n) + (c))
#define eidx(c) (3 * nFluids + (c))

#define sq(x) ((x) * (x))

#define N (9)
#define F2(base,i,j) (base)[(j)*N+(i)]


__global__ static void cuda_gkylMomentSrcSetPtrs(
    double *d_lhs, double *d_rhs, double **d_lhs_ptr, double **d_rhs_ptr) {
  // numThreads*numBlocks == numRealCells
  int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;
  d_lhs_ptr[linearIdx] = d_lhs + (N*N)*linearIdx;
  d_rhs_ptr[linearIdx] = d_rhs + (N)*linearIdx;
}


GkylMomentSrcDeviceCUBLAS_t *cuda_gkylMomentSrcTimeCenteredInit(
    int numBlocks, int numThreads) {
  GkylMomentSrcDeviceCUBLAS_t *context = new GkylMomentSrcDeviceCUBLAS_t[1];
  cublascall(cublasCreate(&(context->handle)));

  int batchSize = numThreads*numBlocks;
  // device memory for actuall arrays and vectors
  cudacall(cudaMalloc(&context->d_lhs, batchSize*N*N*sizeof(double)));
  cudacall(cudaMalloc(&context->d_rhs, batchSize*N*sizeof(double)));
  cudacall(cudaMalloc(&context->d_info, batchSize*sizeof(int)));
  // device memory for pointers to the actuall arrays and vectors
  cudacall(cudaMalloc(&context->d_lhs_ptr, batchSize*sizeof(double*)));
  cudacall(cudaMalloc(&context->d_rhs_ptr, batchSize*sizeof(double*)));

  cuda_gkylMomentSrcSetPtrs<<<numBlocks, numThreads>>>(
      context->d_lhs, context->d_rhs, context->d_lhs_ptr, context->d_rhs_ptr);

  return context;
}


void cuda_gkylMomentSrcTimeCenteredDestroy(
    GkylMomentSrcDeviceCUBLAS_t *context) {
  cudacall(cudaFree(context->d_lhs_ptr));
  cudacall(cudaFree(context->d_rhs_ptr));
  cudacall(cudaFree(context->d_lhs));
  cudacall(cudaFree(context->d_rhs));
  cudacall(cudaFree(context->d_info));
  cublascall(cublasDestroy(context->handle));
}


// FIXME simplify; separate out pressure part
__global__ static void cuda_gkylMomentSrcUpdateRhovE(
    MomentSrcData_t *sd, FluidData_t *fd, double dt, GkylCartField_t **fluidFlds,
    GkylCartField_t *emFld, double **d_rhs_ptr) {
  GkylRange_t *localRange = emFld->localRange;
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = emFld->genIndexer();

  // numThreads*numBlocks == numRealCells
  int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;
  int idxC[3];
  localIdxr.invIndex(linearIdx, idxC);
  const int linearIdxC = fIdxr.index(idxC);
  double *em = emFld->getDataPtrAt(linearIdxC);

  double *sol = d_rhs_ptr[linearIdx];
  unsigned nFluids = sd->nFluids;

  double keOld[2]; // XXX
  for (int n = 0; n < nFluids; ++n)
  {
    double *f = fluidFlds[n]->getDataPtrAt(linearIdxC);
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

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
  em[EX] = 2*sol[eidx(X)] - em[EX];
  em[EY] = 2*sol[eidx(Y)] - em[EY];
  em[EZ] = 2*sol[eidx(Z)] - em[EZ];

  //------------> update correction potential
  double crhoc = sd->chi_e*chargeDens/sd->epsilon0;
  em[PHIE] += dt*crhoc;

  if (sd->hasPressure)
  {
    for (int n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = fluidFlds[n]->getDataPtrAt(linearIdxC);
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 

}


__global__ static void cuda_gkylMomentSrcTimeCenteredSetMat(
    MomentSrcData_t *sd, FluidData_t *fd, double dt, GkylCartField_t **fluidFlds,
    GkylCartField_t *emFld, double **d_lhs_ptr, double **d_rhs_ptr) {
  GkylRange_t *localRange = emFld->localRange;
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = emFld->genIndexer();

  // numThreads*numBlocks == numRealCells
  int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;
  int idxC[3];
  localIdxr.invIndex(linearIdx, idxC);
  const int linearIdxC = fIdxr.index(idxC);
  const double *em = emFld->getDataPtrAt(linearIdxC);

  double *lhs = d_lhs_ptr[linearIdx];
  double *rhs = d_rhs_ptr[linearIdx];

  for (int c=0; c<N*N; c++)
    lhs[c] = 0;

  unsigned nFluids = sd->nFluids;
  double dt1 = 0.5 * dt;
  double dt2 = 0.5 * dt / sd->epsilon0;

  for (unsigned n=0; n<nFluids; ++n)
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


static void cuda_gkylMomentSrcTimeCenteredPreAlloc(
    int numBlocks, int numThreads, MomentSrcData_t *sd, FluidData_t *fd,
    double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld,
    GkylMomentSrcDeviceCUBLAS_t *context) {
  double **d_lhs_ptr = context->d_lhs_ptr;
  double **d_rhs_ptr = context->d_rhs_ptr;
  int *d_info = context->d_info;
  cublasHandle_t &handle = context->handle;

  int batchSize = numThreads*numBlocks;

  cuda_gkylMomentSrcTimeCenteredSetMat<<<numBlocks, numThreads>>>(
      sd, fd, dt, fluidFlds, emFld, d_lhs_ptr, d_rhs_ptr);

  cublascall(cublasDgetrfBatched(
      handle,
      N,  // n
      d_lhs_ptr,  // A
      N,  // lda
      NULL,  // int *PivotArray
      d_info,  // int *infoArray
      batchSize // number of pointers contained in A
      ));

  int info;
  cublascall(cublasDgetrsBatched(
      handle,
      CUBLAS_OP_N,  // trans
      N,  // n
      1,  // nrhs
      d_lhs_ptr,  // matrix A
      N,  // lda
      NULL,  // const int *devIpiv
      d_rhs_ptr,  // double *Barray[]
      N,  // ldb
      &info,  // int *info
      batchSize // number of pointers contained in A
      ));

  // update solution
  cuda_gkylMomentSrcUpdateRhovE<<<numBlocks, numThreads>>>(
      sd, fd, dt, fluidFlds, emFld, d_rhs_ptr);
}


void momentSrcAdvanceOnDevicePreAlloc(
    int numBlocks, int numThreads, MomentSrcData_t *sd, FluidData_t *fd,
    double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld,
    GkylMomentSrcDeviceCUBLAS_t *context)
{
  cuda_gkylMomentSrcTimeCenteredPreAlloc(
      numBlocks, numThreads, sd, fd, dt, fluidFlds, emFld, context);
}

