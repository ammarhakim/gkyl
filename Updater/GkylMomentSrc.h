
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

#define cudacall(call)                                                                                                          \
    do                                                                                                                          \
    {                                                                                                                           \
        cudaError_t err = (call);                                                                                               \
        if(cudaSuccess != err)                                                                                                  \
        {                                                                                                                       \
            fprintf(stderr,"CUDA Error:\nFile = %s\nLine = %d\nReason = %s\n", __FILE__, __LINE__, cudaGetErrorString(err));    \
            cudaDeviceReset();                                                                                                  \
            exit(EXIT_FAILURE);                                                                                                 \
        }                                                                                                                       \
    }                                                                                                                           \
    while (0)

#define cublascall(call)                                                                                        \
    do                                                                                                          \
    {                                                                                                           \
        cublasStatus_t status = (call);                                                                         \
        if(CUBLAS_STATUS_SUCCESS != status)                                                                     \
        {                                                                                                       \
            fprintf(stderr,"CUBLAS Error:\nFile = %s\nLine = %d\nCode = %d\n", __FILE__, __LINE__, status);     \
            cudaDeviceReset();                                                                                  \
            exit(EXIT_FAILURE);                                                                                 \
        }                                                                                                       \
                                                                                                                \
    }                                                                                                           \
    while(0)

#endif // GKYL_MOMENT_SRC_H
