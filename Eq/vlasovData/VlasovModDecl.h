#ifndef VLASOV_MOD_DELC_H 
#define VLASOV_MOD_DELC_H 
#include <cmath> 
extern "C" { 
double VlasovVolStream1x1vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vMaxP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vMaxP5(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vMaxP6(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P6(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vMaxP7(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P7(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream1x2vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vMaxP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vMaxP5(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vMaxP6(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P6(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream1x3vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vMaxP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vMax_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vMaxP5(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vMax_X_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream2x2vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vMaxP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vMax_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vMax_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vMaxP5(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vMax_X_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vMax_Y_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream2x3vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x3vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x3vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x3vMaxP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vMax_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vMax_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream3x3vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream3x3vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream3x3vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Z_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 



 
double VlasovVolStream1x1vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vSerP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vSerP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vSerP5(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vSerP6(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P6(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream1x2vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vSerP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vSerP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vSerP5(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vSer_X_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream1x3vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vSerP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vSerP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream2x2vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vSerP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vSer_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vSerP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vSer_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream2x3vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x3vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x3vSerP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vSer_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream3x3vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream3x3vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 



 
double VlasovVolStream1x1vTensorP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vTensor_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vTensorP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vTensor_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vTensorP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vTensor_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vTensorP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vTensor_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream1x2vTensorP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vTensor_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vTensorP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vTensor_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vTensorP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vTensor_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream1x3vTensorP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vTensor_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vTensorP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vTensor_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream2x2vTensorP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vTensor_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vTensor_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vTensorP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vTensor_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vTensor_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr); 


 
double VlasovVol1x1vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vMax_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vMaxP5(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vMax_VX_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vMaxP6(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vMax_VX_P6(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vMaxP7(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vMax_VX_P7(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol1x2vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x2vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x2vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x2vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vMax_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vMax_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x2vMaxP5(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vMax_VX_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vMax_VY_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x2vMaxP6(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vMax_VX_P6(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vMax_VY_P6(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol1x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x3vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x3vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vMax_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VZ_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x3vMaxP5(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vMax_VX_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VY_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VZ_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol2x2vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x2vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x2vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x2vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vMax_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vMax_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x2vMaxP5(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vMax_VX_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vMax_VY_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol2x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x3vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x3vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vMax_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VZ_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol3x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol3x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol3x3vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 



 
double VlasovVol1x1vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vSer_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vSerP5(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vSer_VX_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vSerP6(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vSer_VX_P6(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol1x2vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x2vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x2vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vSer_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x2vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vSer_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vSer_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x2vSerP5(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vSer_VX_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vSer_VY_P5(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol1x3vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x3vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x3vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x3vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vSer_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VZ_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol2x2vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x2vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x2vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vSer_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x2vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vSer_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vSer_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol2x3vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x3vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x3vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol3x3vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol3x3vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 



 
double VlasovVol1x1vTensorP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vTensor_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vTensorP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vTensor_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vTensorP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vTensor_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x1vTensorP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vTensor_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol1x2vTensorP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vTensor_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vTensor_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x2vTensorP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vTensor_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vTensor_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x2vTensorP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vTensor_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vTensor_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol1x3vTensorP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vTensor_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vTensor_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vTensor_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol1x3vTensorP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vTensor_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vTensor_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vTensor_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVol2x2vTensorP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vTensor_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vTensor_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVol2x2vTensorP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vTensor_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vTensor_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


} 
#endif 
