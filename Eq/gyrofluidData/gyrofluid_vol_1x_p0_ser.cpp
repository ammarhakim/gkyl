#include <gyrofluid_mod_decl.h>

double gyrofluid_vol_1x_p0_ser(const double q_, const double m_, const double *w, const double *dx, const double *jac, const double *rBmag, const double *jacDbmag, const double *dBdz, const double *sMom, const double *phi, const double *primMom, const double *cs, double *out) 
{ 
  // q_,m_:   species charge and mass.
  // w:       cell-center.
  // dx:      cell length.
  // uMaxIn:  maximum speed.
  // jac:     jacobian.
  // rBmag:   reciprocal of magnetic field magnitude (1/B).
  // jacDbmag:jacobian divided by B (J/B).
  // sMom:    stepped moments (times Jacobian).
  // phi:     electrostatic potential.
  // primMom: primitive moments (upar, Tpar, Tperp).
  // cs:      sound speed.
  // out:     output increment.

  double wx = w[0];
  double rdx2 = 2.0/dx[0];
  double rdxSq4 = rdx2*rdx2;


  out[1] += -0.7071067811865475*dBdz[0]*sMom[3]; 



  return rdx2*(fabs(0.7071067811865475*primMom[0]) + 0.7071067811865475*cs[0]); 
}
