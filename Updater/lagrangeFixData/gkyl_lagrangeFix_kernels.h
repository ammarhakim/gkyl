#ifndef LAGRANGE_FIX_KERNELS_H
#define LAGRANGE_FIX_KERNELS_H


#include <math.h>
extern "C" {

  void lagrangeFix_vlasov_1x1v_ser_p1(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_1x1v_ser_p1(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x1v_ser_p2(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_1x1v_ser_p2(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x1v_ser_p3(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_1x1v_ser_p3(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x2v_ser_p1(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_1x2v_ser_p1(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x2v_ser_p2(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_1x2v_ser_p2(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x2v_ser_p3(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_1x2v_ser_p3(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x3v_ser_p1(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x3v_ser_p2(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x3v_ser_p3(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_2x2v_ser_p1(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_2x2v_ser_p1(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_2x2v_ser_p2(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_2x2v_ser_p2(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_2x2v_ser_p3(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_2x2v_ser_p3(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_2x3v_ser_p1(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_2x3v_ser_p2(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_3x3v_ser_p1(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_3x2v_ser_p1(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x1v_tensor_p1(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_1x1v_tensor_p1(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x1v_tensor_p2(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_1x1v_tensor_p2(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x2v_tensor_p1(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_1x2v_tensor_p1(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_1x2v_tensor_p2(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_1x2v_tensor_p2(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_2x2v_tensor_p1(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_2x2v_tensor_p1(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_vlasov_2x2v_tensor_p2(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);
  void lagrangeFix_gk_2x2v_tensor_p2(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f);

}
#endif