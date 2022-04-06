#ifndef HASEGAWA_WAKATANI_MOD_DECL_H 
#define HASEGAWA_WAKATANI_MOD_DECL_H 

#include <cmath>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

extern "C" { 

  double hasegawa_wakatani_vol_2x_p1_Ser(const double C_, const double kappa_, const double *xc, const double *dx, const double *phi, const double *fIn, double *out); 
  double hasegawa_wakatani_surf_2x_p1_Ser_x(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR); 
  double hasegawa_wakatani_surf_2x_p1_Ser_y(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR); 

  double hasegawa_wakatani_vol_2x_p2_Ser(const double C_, const double kappa_, const double *xc, const double *dx, const double *phi, const double *fIn, double *out); 
  double hasegawa_wakatani_surf_2x_p2_Ser_x(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR); 
  double hasegawa_wakatani_surf_2x_p2_Ser_y(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR); 

  double hasegawa_wakatani_vol_2x_p1_Tensor(const double C_, const double kappa_, const double *xc, const double *dx, const double *phi, const double *fIn, double *out); 
  double hasegawa_wakatani_surf_2x_p1_Tensor_x(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR); 
  double hasegawa_wakatani_surf_2x_p1_Tensor_y(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR); 

  double hasegawa_wakatani_vol_2x_p2_Tensor(const double C_, const double kappa_, const double *xc, const double *dx, const double *phi, const double *fIn, double *out); 
  double hasegawa_wakatani_surf_2x_p2_Tensor_x(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR); 
  double hasegawa_wakatani_surf_2x_p2_Tensor_y(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR); 

} 
#endif 
