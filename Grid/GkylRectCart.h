#pragma once

extern "C" 
{
    typedef struct {
        int ndim;
        int cells[6];
        double lower[6], upper[6];
        double vol, dx[6];
        __host__ __device__ __inline__ void cellCenter(int* idx, double* xc)
        {
          for(unsigned int d=0; d<ndim; d++) {
            xc[d] = lower[d] + (idx[d]-0.5)*dx[d];
          }
        }
    } GkylRectCart_t;
}

