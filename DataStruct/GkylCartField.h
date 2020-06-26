#ifndef GKYL_CART_FIELD_H
#define GKYL_CART_FIELD_H

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylRange.h>
#include <GkylRectCart.h>

extern "C" 
{
    typedef struct {
        int ndim;
        int elemSize;
        int numComponents;
        GkylRange_t *localRange, *localExtRange;
        GkylRange_t *localExtEdgeRange;
        GkylRange_t *globalRange, *globalExtRange;
        GkylRectCart_t *grid;
        double *_data;
        
        __host__ __device__ __inline__ double* getDataPtrAt(int linIdx)
        {
           return _data + linIdx*numComponents;
        }
        __host__ __device__ __inline__ Gkyl::GenIndexer genIndexer()
        {
           return Gkyl::GenIndexer(localExtRange);
        }
    } GkylCartField_t;
}

#endif // GKYL_CART_FIELD_H
