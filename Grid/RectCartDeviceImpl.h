#ifndef RECT_CART_DEVICE_IMPL_H
#define RECT_CART_DEVICE_IMPL_H

extern "C" 
{
    typedef struct {
        int32_t ndim;
        int32_t cells[6];
        double lower[6], upper[6];
        double vol, dx[6];
        __device__ __inline__ void cellCenter(int* idx, double* xc)
        {
          for(unsigned int d=0; d<this->ndim; d++) {
            xc[d] = this->lower[d] + (idx[d]-0.5)*this->dx[d];
          }
        }
    } RectCart_t;
}


#endif // RECT_CART_DEVICE_IMPL_H
