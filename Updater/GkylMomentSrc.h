
#ifndef GKYL_MOMENT_SRC_H
#define GKYL_MOMENT_SRC_H

#include <MomentSrcCommon.h>
#include <GkylCartField.h>
#include <GkylCudaConfig.h>
#include <cublas_v2.h>

extern "C" {
  typedef struct {
    // TODO domain decomposition info

    // for time-centered scheme using cublas batched getrf and getrs
    double *d_lhs;
    double *d_rhs;
    double **d_lhs_ptr;
    double **d_rhs_ptr;
    int *d_info;
    cublasHandle_t handle;
  } GkylMomentSrcDeviceData_t;

  GkylMomentSrcDeviceData_t *cuda_gkylMomentSrcInit(
      int numBlocks, int numThreads);
  void cuda_gkylMomentSrcDestroy(
      GkylMomentSrcDeviceData_t *context);
  void momentSrcAdvanceOnDevice(
    int numBlocks, int numThreads, MomentSrcData_t *sd, FluidData_t *fd,
    double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld,
    GkylMomentSrcDeviceData_t *context);
}

#endif // GKYL_MOMENT_SRC_H
