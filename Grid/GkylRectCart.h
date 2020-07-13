// Gkyl ------------------------------------------------------------------------
//
// Rectangular grid objects
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GKYL_RECT_CART_H
#define GKYL_RECT_CART_H

extern "C" 
{
    typedef struct {
        int ndim;
        int cells[6];
        double lower[6], upper[6];
        double vol, dx[6];
        __host__ __device__ __inline__ void cellCenter(const int* __restrict__ idx, double* xc)
        {
          #pragma unroll
          for(unsigned int d=0; d<6; d++) {
            xc[d] = lower[d] + (idx[d]-0.5)*dx[d];
          }
        }
    } GkylRectCart_t;
}

#endif // GKYL_RECT_CART_H
