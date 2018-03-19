#ifndef VLASOV_MOD_DELC_H 
#define VLASOV_MOD_DELC_H 
#include <cmath> 
extern "C" { 
double VlasovVolStream1x1vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vMaxP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream1x2vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vMaxP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream1x3vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vMax_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vMaxP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vMax_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream2x2vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vMax_Y_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vMax_Y_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vMax_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vMax_Y_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vMaxP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vMax_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vMax_Y_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream2x3vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vMax_Y_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x3vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vMax_Y_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x3vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vMax_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vMax_Y_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x3vMaxP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vMax_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vMax_Y_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream3x3vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Y_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Z_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream3x3vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Y_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Z_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream3x3vMaxP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vMax_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Y_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Z_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream3x3vMaxP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vMax_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Y_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vMax_Z_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 



 
double VlasovVolStream1x1vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vSerP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x1vSerP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream1x2vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vSerP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vSer_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x2vSerP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vSer_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream1x3vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vSerP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vSer_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream1x3vSerP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vSer_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream2x2vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vSer_Y_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vSer_Y_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vSerP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vSer_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vSer_Y_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x2vSerP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vSer_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vSer_Y_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream2x3vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vSer_Y_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x3vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vSer_Y_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x3vSerP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vSer_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vSer_Y_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream2x3vSerP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vSer_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vSer_Y_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolStream3x3vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Y_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Z_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream3x3vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Y_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Z_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream3x3vSerP3(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vSer_X_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Y_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Z_P3(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolStream3x3vSerP4(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream3x3vSer_X_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Y_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream3x3vSer_Z_P4(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 



 
void VlasovVolElc1x1vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vMax_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x1vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vMax_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x1vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vMax_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x1vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vMax_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc1x2vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vMax_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x2vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vMax_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x2vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vMax_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x2vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vMax_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc1x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vMax_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vMax_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x3vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vMax_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x3vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vMax_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc2x2vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vMax_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vMax_VY_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x2vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vMax_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vMax_VY_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x2vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vMax_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vMax_VY_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x2vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vMax_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vMax_VY_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc2x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vMax_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vMax_VY_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vMax_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vMax_VY_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x3vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vMax_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vMax_VY_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x3vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vMax_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vMax_VY_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc3x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc3x3vMax_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vMax_VY_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vMax_VZ_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc3x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc3x3vMax_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vMax_VY_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vMax_VZ_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc3x3vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc3x3vMax_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vMax_VY_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vMax_VZ_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc3x3vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc3x3vMax_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vMax_VY_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vMax_VZ_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 



 
void VlasovVolElc1x1vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vSer_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x1vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vSer_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x1vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vSer_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x1vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vSer_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc1x2vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vSer_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x2vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vSer_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x2vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vSer_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x2vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vSer_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc1x3vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vSer_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x3vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vSer_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x3vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vSer_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x3vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vSer_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc2x2vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vSer_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vSer_VY_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x2vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vSer_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vSer_VY_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x2vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vSer_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vSer_VY_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x2vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vSer_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vSer_VY_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc2x3vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vSer_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vSer_VY_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x3vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vSer_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vSer_VY_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x3vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vSer_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vSer_VY_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x3vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vSer_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vSer_VY_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc3x3vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc3x3vSer_VX_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vSer_VY_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vSer_VZ_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc3x3vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc3x3vSer_VX_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vSer_VY_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vSer_VZ_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc3x3vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc3x3vSer_VX_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vSer_VY_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vSer_VZ_P3(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc3x3vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc3x3vSer_VX_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vSer_VY_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc3x3vSer_VZ_P4(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolElcMag1x1vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x1vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x1vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vMax_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x1vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vMax_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolElcMag1x2vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vMax_VY_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x2vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vMax_VY_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x2vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vMax_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vMax_VY_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x2vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vMax_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vMax_VY_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolElcMag1x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VY_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VZ_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VY_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VZ_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x3vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vMax_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VY_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VZ_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x3vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vMax_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VY_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vMax_VZ_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolElcMag2x2vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vMax_VY_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x2vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vMax_VY_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x2vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vMax_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vMax_VY_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x2vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vMax_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vMax_VY_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolElcMag2x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VY_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VZ_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VY_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VZ_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x3vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vMax_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VY_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VZ_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x3vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vMax_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VY_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vMax_VZ_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolElcMag3x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VY_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VZ_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag3x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VY_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VZ_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag3x3vMaxP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vMax_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VY_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VZ_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag3x3vMaxP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vMax_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VY_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vMax_VZ_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 



 
double VlasovVolElcMag1x1vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x1vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vSer_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x1vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vSer_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x1vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x1vSer_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolElcMag1x2vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vSer_VY_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x2vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vSer_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vSer_VY_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x2vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vSer_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vSer_VY_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x2vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x2vSer_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x2vSer_VY_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolElcMag1x3vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VY_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VZ_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x3vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vSer_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VY_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VZ_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x3vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vSer_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VY_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VZ_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag1x3vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag1x3vSer_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VY_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag1x3vSer_VZ_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolElcMag2x2vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vSer_VY_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x2vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vSer_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vSer_VY_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x2vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vSer_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vSer_VY_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x2vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x2vSer_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x2vSer_VY_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolElcMag2x3vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VY_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VZ_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x3vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vSer_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VY_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VZ_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x3vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vSer_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VY_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VZ_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag2x3vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag2x3vSer_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VY_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag2x3vSer_VZ_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


double VlasovVolElcMag3x3vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VY_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VZ_P1(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag3x3vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vSer_VX_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VY_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VZ_P2(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag3x3vSerP3(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vSer_VX_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VY_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VZ_P3(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

double VlasovVolElcMag3x3vSerP4(const double *w, const double *dxv, const double *E, const double *f, double *out); 
double VlasovSurfElcMag3x3vSer_VX_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VY_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
double VlasovSurfElcMag3x3vSer_VZ_P4(const double *w, const double *dxv, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr); 



 
} 
#endif 
