#ifndef RECT_CART_DEVICE_IMPL_H
#define RECT_CART_DEVICE_IMPL_H

extern "C" 
{
    typedef struct { int32_t ndim; int32_t lower[6]; int32_t upper[6]; } Range_t;
    typedef struct {
        int32_t ndim;
        int32_t cells[6];
        double lower[6], upper[6];
        double vol, dx[6];
    } RectCart_t;
}

__device__ void cellCenter(RectCart_t* grid, int* idx, double* xc);

#endif // RECT_CART_DEVICE_IMPL_H
