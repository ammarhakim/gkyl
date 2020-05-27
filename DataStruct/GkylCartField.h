#pragma once

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylRange.h>
#include <GkylRectCart.h>

extern "C" 
{
    typedef struct {
        int32_t numComponents;
        int32_t ndim;
        GkylRange_t *localRange, *localExtRange;
        GkylRange_t *globalRange, *globalExtRange;
        GkylRectCart_t *grid;
        double *_data;
        __host__ __device__ __inline__ double* getDataPtrAt(int linIdx)
        {
           return _data + linIdx*numComponents;
        }
        __host__ __device__ __inline__ Gkyl::GenIndexer genIndexer()
        {
           Gkyl::GenIndexer idxr(localExtRange);
           return idxr;
        }
    } GkylCartField_t;
}

