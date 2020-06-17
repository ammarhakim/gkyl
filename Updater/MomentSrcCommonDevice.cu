#include <cstdio>
#include <GkylMomentSrc.h>
#include <GkylCudaConfig.h>

#include <cublas_v2.h>
#include <cuda_runtime.h>

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

__global__ void cuda_gkylMomentSrcSetMat(
    MomentSrcData_t *sd, FluidData_t *fd, double dt, GkylCartField_t **fluidFlds,
    GkylCartField_t *emFld, double *d_lhs, double *d_rhs) {
  GkylRange_t *localRange = emFld->localRange;
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr = emFld->genIndexer();

  // numThreads*numBlocks == numRealCells
  int linearIdx = threadIdx.x + blockIdx.x*blockDim.x;
  int idxC[3];
  localIdxr.invIndex(linearIdx, idxC);
  const int linearIdxC = fIdxr.index(idxC);
  const double *em = emFld->getDataPtrAt(linearIdxC);

  double *lhs = d_lhs + (N*N)*sizeof(double)*linearIdx;
  double *rhs = d_rhs + (N)*sizeof(double)*linearIdx;

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

static int cuda_gkylMomentSrcTimeCentered(
    int numBlocks, int numThreads, MomentSrcData_t *sd, FluidData_t *fd,
    double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld) {
  // FIXME save d_lhs and avoid reallocating?
  cublasStatus_t status;
  double *d_lhs = 0;
  double *d_rhs = 0;
  cublasHandle_t handle;

  int batchSize = numThreads*numBlocks;

  status = cublasCreate(&handle);
  if (status != CUBLAS_STATUS_SUCCESS) {
    fprintf(stderr, "!!!! CUBLAS initialization error\n");
    return EXIT_FAILURE;
  }

  if (cudaMalloc(reinterpret_cast<void **>(&d_lhs), N*N * sizeof(d_lhs[0])) !=
      cudaSuccess) {
    fprintf(stderr, "!!!! device memory allocation error (allocate A)\n");
    return EXIT_FAILURE;
  }

  if (cudaMalloc(reinterpret_cast<void **>(&d_rhs), N * sizeof(d_rhs[0])) !=
      cudaSuccess) {
    fprintf(stderr, "!!!! device memory allocation error (allocate A)\n");
    return EXIT_FAILURE;
  }

  cuda_gkylMomentSrcSetMat<<<numBlocks, numThreads>>>(
      sd, fd, dt, fluidFlds, emFld, d_lhs, d_rhs);

  /* status = cublasDgetrfBatched( */
  /*     handle, */
  /*     N,  // n */
  /*     &d_lhs,  // A */
  /*     N,  // lda */
  /*     NULL,  // int *PivotArray */
  /*     NULL,  // int *infoArray */
  /*     batchSize // number of pointers contained in A */
  /*     ); */
  /* if (status != CUBLAS_STATUS_SUCCESS) { */
  /*   fprintf(stderr, "!!!! LU decomposition kernel execution error.\n"); */
  /*   return EXIT_FAILURE; */
  /* } */
  /*  */
  /* status = cublasDgetrsBatched( */
  /*     handle, */
  /*     CUBLAS_OP_N,  // trans */
  /*     N,  // n */
  /*     1,  // nrhs */
  /*     &d_lhs,  // matrix A */
  /*     N,  // lda */
  /*     NULL,  // const int *devIpiv */
  /*     &d_rhs,  // double *Barray[] */
  /*     N,  // ldb */
  /*     NULL,  // int *info */
  /*     batchSize // number of pointers contained in A */
  /*     ); */
  /* if (status != CUBLAS_STATUS_SUCCESS) { */
  /*   fprintf(stderr, "!!!! solve kernel execution error.\n"); */
  /*   return EXIT_FAILURE; */
  /* } */

  // update solution

  if (cudaFree(d_lhs) != cudaSuccess) {
    fprintf(stderr, "!!!! memory free error (d_lhs)\n");
    return EXIT_FAILURE;
  }

  if (cudaFree(d_rhs) != cudaSuccess) {
    fprintf(stderr, "!!!! memory free error (d_rhs)\n");
    return EXIT_FAILURE;
  }
 
  status = cublasDestroy(handle);
  if (status != CUBLAS_STATUS_SUCCESS) {
    fprintf(stderr, "!!!! shutdown error (A)\n");
    return EXIT_FAILURE;
  }
 
  exit(EXIT_SUCCESS);
}

void momentSrcAdvanceOnDevice(
    int numBlocks, int numThreads, MomentSrcData_t *sd, FluidData_t *fd,
    double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld)
{
  cuda_gkylMomentSrcTimeCentered(
      numBlocks, numThreads, sd, fd, dt, fluidFlds, emFld);
}
