#ifndef MAXWELL_MOD_DELC_H 
#define MAXWELL_MOD_DELC_H 
#include <cmath> 
extern "C" { 
typedef struct { double c, chi, gamma; } MaxwellEq_t; 
 
double MaxwellVol1xMaxP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf1xMax_X_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf1xMax_X_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol1xMaxP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf1xMax_X_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf1xMax_X_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol1xMaxP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf1xMax_X_P3(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf1xMax_X_P3(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 


double MaxwellVol2xMaxP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf2xMax_X_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xMax_X_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf2xMax_Y_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xMax_Y_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol2xMaxP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf2xMax_X_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xMax_X_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf2xMax_Y_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xMax_Y_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol2xMaxP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf2xMax_X_P3(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xMax_X_P3(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf2xMax_Y_P3(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xMax_Y_P3(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 


double MaxwellVol3xMaxP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf3xMax_X_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf3xMax_X_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf3xMax_Y_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf3xMax_Y_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf3xMax_Z_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf3xMax_Z_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol3xMaxP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf3xMax_X_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf3xMax_X_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf3xMax_Y_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf3xMax_Y_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf3xMax_Z_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf3xMax_Z_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 



 
double MaxwellVol1xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf1xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf1xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol1xSerP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf1xSer_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf1xSer_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol1xSerP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf1xSer_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf1xSer_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 


double MaxwellVol2xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf2xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf2xSer_Y_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xSer_Y_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol2xSerP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf2xSer_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xSer_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf2xSer_Y_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xSer_Y_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol2xSerP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf2xSer_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xSer_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf2xSer_Y_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xSer_Y_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 


double MaxwellVol3xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf3xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf3xSer_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf3xSer_Y_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf3xSer_Y_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf3xSer_Z_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf3xSer_Z_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 



 
double MaxwellVol1xTensorP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf1xTensor_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf1xTensor_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol1xTensorP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf1xTensor_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf1xTensor_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol1xTensorP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf1xTensor_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf1xTensor_X_P3(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 


double MaxwellVol2xTensorP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf2xTensor_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xTensor_X_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf2xTensor_Y_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xTensor_Y_P1(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 

double MaxwellVol2xTensorP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
double MaxwellSurf2xTensor_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xTensor_X_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellSurf2xTensor_Y_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 
double MaxwellCentralSurf2xTensor_Y_P2(const MaxwellEq_t * meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr); 


} 
#endif 
