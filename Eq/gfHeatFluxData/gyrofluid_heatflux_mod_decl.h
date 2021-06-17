#ifndef GYROFLUID_HEATFLUX_MOD_DECL_H 
#define GYROFLUID_HEATFLUX_MOD_DECL_H 

#include <math.h> 

extern "C" { 

  double gyrofluid_heatflux_vol_1x_p0_ser(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *w, const double *dx, const double *jac, const double *jacDbmag, const double *sMom, const double *phi, double *primMom, double *out); 
  double gyrofluid_heatflux_surf_1x_p0_ser_x(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *wL1, const double *dxL1, const double *wR1, const double *dxR1, const double uMaxIn, const double *jacL, const double *jacDbmagL, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, double *outL, double *outR); 

  double gyrofluid_heatflux_vol_1x_p1_ser(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *w, const double *dx, const double *jac, const double *jacDbmag, const double *sMom, const double *phi, double *primMom, double *out); 
  double gyrofluid_heatflux_surf_1x_p1_ser_x(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *wL1, const double *dxL1, const double *wR1, const double *dxR1, const double uMaxIn, const double *jacL, const double *jacDbmagL, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, double *outL, double *outR); 

  double gyrofluid_heatflux_vol_1x_p2_ser(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *w, const double *dx, const double *jac, const double *jacDbmag, const double *sMom, const double *phi, double *primMom, double *out); 
  double gyrofluid_heatflux_surf_1x_p2_ser_x(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *wL1, const double *dxL1, const double *wR1, const double *dxR1, const double uMaxIn, const double *jacL, const double *jacDbmagL, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, double *outL, double *outR); 

} 
#endif 
