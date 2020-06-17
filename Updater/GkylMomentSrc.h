
#ifndef GKYL_MOMENT_SRC_H
#define GKYL_MOMENT_SRC_H

#include <GkylCudaConfig.h>
#include <MomentSrcCommon.h>
#include <GkylCartField.h>

extern "C" {
    /* CUDA GPU */
    void momentSrcAdvanceOnDevice(
        int numBlocks, int numThreads, MomentSrcData_t *sd, FluidData_t *fd,
        double dt, GkylCartField_t **fluidFlds, GkylCartField_t *emFld);
}

#endif // GKYL_MOMENT_SRC_H
