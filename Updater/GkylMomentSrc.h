
#ifndef GKYL_MOMENT_SRC_H
#define GKYL_MOMENT_SRC_H

#include <MomentSrcCommon.h>
#include <GkylCartField.h>
#include <GkylCudaConfig.h>
#include <cublas_v2.h>

extern "C" {
  typedef struct {
    // TODO domain decomposition info
    const char *scheme;
    int nFluids;
    int numThreads;
    int numBlocks;

    // for time-centered scheme using cublas batched getrf and getrs
    double *d_lhs;
    double *d_rhs;
    double **d_lhs_ptr;
    double **d_rhs_ptr;
    int *d_info;
    cublasHandle_t handle;
  } GkylMomentSrcDeviceData_t;

  GkylMomentSrcDeviceData_t *cuda_gkylMomentSrcInit(
    const char *scheme, const int nFluids, const int numBlocks, const int numThreads);
  void cuda_gkylMomentSrcDestroy(const GkylMomentSrcDeviceData_t *context);
  void momentSrcAdvanceOnDevice(
      const MomentSrcData_t *sd, const FluidData_t *fd, const double dt,
      GkylCartField_t **fluidFlds, GkylCartField_t *emFld,
      const GkylMomentSrcDeviceData_t *context);
  
}

#endif // GKYL_MOMENT_SRC_H
