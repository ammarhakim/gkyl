#ifndef MAXWELL_MOD_DELC_H 
#define MAXWELL_MOD_DELC_H
#include <GkylCudaConfig.h> 
#include <cmath> 
extern "C" { 
typedef struct { double c, chi, gamma; } MaxwellEq_t; 
 
__host__ __device__ double MaxwellVol1xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
__host__ __device__ double MaxwellSurf1xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf1xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

__host__ __device__ double MaxwellVol1xSerP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
__host__ __device__ double MaxwellSurf1xSer_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf1xSer_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

__host__ __device__ double MaxwellVol1xSerP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
__host__ __device__ double MaxwellSurf1xSer_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf1xSer_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 


__host__ __device__ double MaxwellVol2xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
__host__ __device__ double MaxwellSurf2xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf2xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellSurf2xSer_Y_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf2xSer_Y_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

__host__ __device__ double MaxwellVol2xSerP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
__host__ __device__ double MaxwellSurf2xSer_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf2xSer_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellSurf2xSer_Y_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf2xSer_Y_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

__host__ __device__ double MaxwellVol2xSerP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
__host__ __device__ double MaxwellSurf2xSer_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf2xSer_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellSurf2xSer_Y_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf2xSer_Y_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 


__host__ __device__ double MaxwellVol3xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
__host__ __device__ double MaxwellSurf3xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf3xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellSurf3xSer_Y_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf3xSer_Y_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellSurf3xSer_Z_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf3xSer_Z_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 


__host__ __device__ double MaxwellVol1xTensorP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
__host__ __device__ double MaxwellSurf1xTensor_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf1xTensor_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

__host__ __device__ double MaxwellVol1xTensorP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
__host__ __device__ double MaxwellSurf1xTensor_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf1xTensor_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 


__host__ __device__ double MaxwellVol2xTensorP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
__host__ __device__ double MaxwellSurf2xTensor_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf2xTensor_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellSurf2xTensor_Y_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf2xTensor_Y_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

__host__ __device__ double MaxwellVol2xTensorP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
__host__ __device__ double MaxwellSurf2xTensor_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf2xTensor_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellSurf2xTensor_Y_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
__host__ __device__ double MaxwellCentralSurf2xTensor_Y_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 


} 
#endif 
