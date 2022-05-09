#ifndef EPAR_GYROKINETIC_MOD_DECL_H 
#define EPAR_GYROKINETIC_MOD_DECL_H 

#include <cmath>

#define SURFAVG 1 
#define QUAD 2 
#define cflType QUAD 
#define upwindType QUAD 

template <typename T> int sgn(T val) { 
  return (T(0) < val) - (val < T(0)); 
}

extern "C" { 

  double EparGyrokineticVol1x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *Epar, const double *f, double *out); 
  double EparGyrokineticSurf1x2vSer_xL_P1(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *Epar, const double *fL, const double *fR, double *outL, double *outR); 
  double EparGyrokineticSurf1x2vSer_xR_P1(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *Epar, const double *fL, const double *fR, double *outL, double *outR); 
  double EparGyrokineticSurf1x2vSer_vpar_P1(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *Epar, const double *fL, const double *fR, double *outL, double *outR); 

  double EparGyrokineticVol1x2vSerP2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *Epar, const double *f, double *out); 
  double EparGyrokineticSurf1x2vSer_xL_P2(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *Epar, const double *fL, const double *fR, double *outL, double *outR); 
  double EparGyrokineticSurf1x2vSer_xR_P2(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *Epar, const double *fL, const double *fR, double *outL, double *outR); 
  double EparGyrokineticSurf1x2vSer_vpar_P2(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *Epar, const double *fL, const double *fR, double *outL, double *outR); 


} 
#endif 
