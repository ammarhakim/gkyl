
#ifndef GKYL_MOMENT_SRC_H
#define GKYL_MOMENT_SRC_H

#include <MomentSrcCommon.h>
#include <GkylCartField.h>
#include <GkylCudaConfig.h>
#include <cublas_v2.h>

extern "C" {
  void momentSrcAdvanceOnDevice(
      int numBlocks, int numThreads, MomentSrcData_t *sd, FluidData_t *fd,
      double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld);
}

#endif // GKYL_MOMENT_SRC_H
