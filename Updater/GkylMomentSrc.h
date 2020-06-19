
#ifndef GKYL_MOMENT_SRC_H
#define GKYL_MOMENT_SRC_H

#include <MomentSrcCommon.h>
#include <GkylCartField.h>
#include <GkylCudaConfig.h>
#include <cublas_v2.h>

extern "C" {
  typedef struct {
    // TODO domain decomposition info
    double *d_lhs;
    double *d_rhs;
    double **d_lhs_ptr;
    double **d_rhs_ptr;
    int *d_info;
    cublasHandle_t handle;
  } GkylMomentSrcDeviceCUBLAS_t;

  GkylMomentSrcDeviceCUBLAS_t *cuda_gkylMomentSrcTimeCenteredInit(
      int numBlocks, int numThreads);
  void momentSrcAdvanceOnDevicePreAlloc(
    int numBlocks, int numThreads, MomentSrcData_t *sd, FluidData_t *fd,
    double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld,
    GkylMomentSrcDeviceCUBLAS_t *context);
  void momentSrcAdvanceOnDevice(
      int numBlocks, int numThreads, MomentSrcData_t *sd, FluidData_t *fd,
      double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld);
}

#endif // GKYL_MOMENT_SRC_H
