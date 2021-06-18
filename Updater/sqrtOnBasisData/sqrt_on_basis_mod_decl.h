#ifndef SQRT_ON_BASIS_MOD_DECL_H 
#define SQRT_ON_BASIS_MOD_DECL_H 

#include <cmath> 

extern "C" { 

  void sqrt_on_basis_gauss_1x_p1_ser(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_1x_p2_ser(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_1x_p3_ser(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_2x_p1_ser(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_2x_p2_ser(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_2x_p3_ser(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_3x_p1_ser(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_3x_p2_ser(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_3x_p3_ser(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_4x_p1_ser(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_4x_p2_ser(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_5x_p1_ser(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_5x_p2_ser(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_6x_p1_ser(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_1x_p1_max(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_1x_p2_max(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_1x_p3_max(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_2x_p1_max(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_2x_p2_max(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_2x_p3_max(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_3x_p1_max(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_3x_p2_max(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_3x_p3_max(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_4x_p1_max(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_4x_p2_max(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_1x_p1_tensor(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_1x_p2_tensor(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_1x_p3_tensor(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_2x_p1_tensor(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_2x_p2_tensor(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_2x_p3_tensor(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_3x_p1_tensor(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_3x_p2_tensor(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_3x_p3_tensor(const double qExp, const double *f, double *out); 

  void sqrt_on_basis_gauss_4x_p1_tensor(const double qExp, const double *f, double *out); 
  void sqrt_on_basis_gauss_4x_p2_tensor(const double qExp, const double *f, double *out); 

} 

#endif 
