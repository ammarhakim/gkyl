#include <GkylMomentSrc.h>
#include <cstdio>

// Makes indexing cleaner
#define X (0)
#define Y (1)
#define Z (2)

#define RHO (0)
#define MX (1)
#define MY (2)
#define MZ (3)
#define ER (4)

#define EX (0)
#define EY (1)
#define EZ (2)
#define BX (3)
#define BY (4)
#define BZ (5)
#define PHIE (6)
#define PHIM (7)

#define fidx(n, c) (3 * (n) + (c))
#define eidx(c) (3 * nFluids + (c))

#define N (9)
#define sq(x) ((x) * (x))

#define F2(base,i,j) (base)[(j)*N+(i)]


// FIXME simplify; separate out pressure part
__global__ void cuda_gkylMomentSrcUpdateRhovE(
    MomentSrcData_t *sd, FluidData_t *fd, double dt, GkylCartField_t **fluidFlds,
    GkylCartField_t *emFld, double *d_lhs, double *d_rhs, double **d_lhs_ptr,
    double **d_rhs_ptr) {
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

__global__ void cuda_gkylMomentSrcSetMat(
    MomentSrcData_t *sd, FluidData_t *fd, double dt, GkylCartField_t **fluidFlds,
    GkylCartField_t *emFld, double *d_lhs, double *d_rhs, double **d_lhs_ptr,
    double **d_rhs_ptr) {
  GkylRange_t *localRange = emFld->localRange;
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = emFld->genIndexer();

  // numThreads*numBlocks == numRealCells
  int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;
  int idxC[3];
  localIdxr.invIndex(linearIdx, idxC);
  const int linearIdxC = fIdxr.index(idxC);
  const double *em = emFld->getDataPtrAt(linearIdxC);

  double *lhs = d_lhs + (N*N)*linearIdx;
  double *rhs = d_rhs + (N)*linearIdx;
  d_lhs_ptr[linearIdx] = lhs;
  d_rhs_ptr[linearIdx] = rhs;

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

static cudaError_t cuda_gkylMomentSrcTimeCentered(
    int numBlocks, int numThreads, MomentSrcData_t *sd, FluidData_t *fd,
    double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld) {
  // FIXME save d_lhs and avoid reallocating?
  double *d_lhs = 0;
  double *d_rhs = 0;
  double **d_lhs_ptr = 0;
  double **d_rhs_ptr = 0;
  int *d_info = 0;
  cublasHandle_t handle;

  int batchSize = numThreads*numBlocks;

  cublascall(cublasCreate(&handle));

  // memory for actuall arrays and vectors
  cudacall(cudaMalloc(&d_lhs, batchSize*N*N*sizeof(double)));
  cudacall(cudaMalloc(&d_rhs, batchSize*N*sizeof(double)));
  cudacall(cudaMalloc(&d_info, batchSize*sizeof(int)));

  // memory for pointers to the actuall arrays and vectors
  cudacall(cudaMalloc(&d_lhs_ptr, batchSize*sizeof(double*)));
  cudacall(cudaMalloc(&d_rhs_ptr, batchSize*sizeof(double*)));

  cuda_gkylMomentSrcSetMat<<<numBlocks, numThreads>>>(
      sd, fd, dt, fluidFlds, emFld, d_lhs, d_rhs, d_lhs_ptr, d_rhs_ptr);

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
      sd, fd, dt, fluidFlds, emFld, d_lhs, d_rhs, d_lhs_ptr, d_rhs_ptr);

  cudacall(cudaFree(d_lhs_ptr));
  cudacall(cudaFree(d_rhs_ptr));
  cudacall(cudaFree(d_lhs));
  cudacall(cudaFree(d_rhs));
  cudacall(cudaFree(d_info));

  cublascall(cublasDestroy(handle));

  return cudaSuccess;
}

void momentSrcAdvanceOnDevice(
    int numBlocks, int numThreads, MomentSrcData_t *sd, FluidData_t *fd,
    double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld)
{
  cuda_gkylMomentSrcTimeCentered(
      numBlocks, numThreads, sd, fd, dt, fluidFlds, emFld);
}
