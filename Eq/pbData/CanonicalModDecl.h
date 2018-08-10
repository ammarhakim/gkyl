#ifndef CANONICAL_MOD_DECL_H 
#define CANONICAL_MOD_DECL_H 
#include <cmath> 
#include <Positivity.h> 
extern "C" { 
double CanonicalVol1x1vSerP1(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf1x1vSer_X_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity1x1vSer_X_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf1x1vSer_VX_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity1x1vSer_VX_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 

double CanonicalVol1x1vSerP2(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf1x1vSer_X_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity1x1vSer_X_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf1x1vSer_VX_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity1x1vSer_VX_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 


double CanonicalVol2x2vSerP1(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf2x2vSer_X_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity2x2vSer_X_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_Y_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity2x2vSer_Y_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_VX_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity2x2vSer_VX_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_VY_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity2x2vSer_VY_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 

double CanonicalVol2x2vSerP2(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf2x2vSer_X_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity2x2vSer_X_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_Y_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity2x2vSer_Y_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_VX_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity2x2vSer_VX_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_VY_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurfPositivity2x2vSer_VY_P2(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 


} 
#endif 
