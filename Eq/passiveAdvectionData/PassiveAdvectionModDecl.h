#ifndef PASSIVE_ADVECTION_MOD_DECL_H 
#define PASSIVE_ADVECTION_MOD_DECL_H 
#include <cmath> 
#include <algorithm> 
#include <Positivity.h> 
#define SURFAVG 1 
#define QUAD 2 
#define cflType QUAD 
#define upwindType QUAD 
template <typename T> int sgn(T val) { 
    return (T(0) < val) - (val < T(0)); 
}
extern "C" { 
double PassiveAdvectionVol1xSerP1(const double *w, const double *dxv, double *cfl, const double *f, double *out); 
double PassiveAdvectionSurf1xSer_X1_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity1xSer_X1_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol1xSerP2(const double *w, const double *dxv, double *cfl, const double *f, double *out); 
double PassiveAdvectionSurf1xSer_X1_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity1xSer_X1_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol2xSerP1(const double *w, const double *dxv, double *cfl, const double *f, double *out); 
double PassiveAdvectionSurf2xSer_X1_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity2xSer_X1_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf2xSer_X2_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity2xSer_X2_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol2xSerP2(const double *w, const double *dxv, double *cfl, const double *f, double *out); 
double PassiveAdvectionSurf2xSer_X1_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity2xSer_X1_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf2xSer_X2_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity2xSer_X2_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol3xSerP1(const double *w, const double *dxv, double *cfl, const double *f, double *out); 
double PassiveAdvectionSurf3xSer_X1_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity3xSer_X1_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf3xSer_X2_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity3xSer_X2_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf3xSer_X3_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity3xSer_X3_P1(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol3xSerP2(const double *w, const double *dxv, double *cfl, const double *f, double *out); 
double PassiveAdvectionSurf3xSer_X1_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity3xSer_X1_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf3xSer_X2_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity3xSer_X2_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf3xSer_X3_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurfPositivity3xSer_X3_P2(const double *cflL, const double *cflR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr); 

} 
#endif 
