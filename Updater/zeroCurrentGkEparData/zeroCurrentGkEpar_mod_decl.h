#ifndef ZERO_CURRENT_GK_EPAR_MOD_DECL_H
#define ZERO_CURRENT_GK_EPAR_MOD_DECL_H

extern "C" { 

#include <math.h>
  void zeroCurrentGkEpar_divP_1x_p1_ser(const double dz, const double charge, const double *delparFac, const double *dlnbmagdz, const double *m2par, const double *m2perp, double *divP); 
  void zeroCurrentGkEpar_M1scale_1x2v_p1_ser(const double *xvc, const double *dxv, const double vparMax, const double *alpha, double *distf); 


  void zeroCurrentGkEpar_divP_1x_p2_ser(const double dz, const double charge, const double *delparFac, const double *dlnbmagdz, const double *m2par, const double *m2perp, double *divP); 
  void zeroCurrentGkEpar_M1scale_1x2v_p2_ser(const double *xvc, const double *dxv, const double vparMax, const double *alpha, double *distf); 


  void zeroCurrentGkEpar_divP_1x_p1_tensor(const double dz, const double charge, const double *delparFac, const double *dlnbmagdz, const double *m2par, const double *m2perp, double *divP); 
  void zeroCurrentGkEpar_M1scale_1x2v_p1_tensor(const double *xvc, const double *dxv, const double vparMax, const double *alpha, double *distf); 


  void zeroCurrentGkEpar_divP_1x_p2_tensor(const double dz, const double charge, const double *delparFac, const double *dlnbmagdz, const double *m2par, const double *m2perp, double *divP); 
  void zeroCurrentGkEpar_M1scale_1x2v_p2_tensor(const double *xvc, const double *dxv, const double vparMax, const double *alpha, double *distf); 



 
}

#endif 
