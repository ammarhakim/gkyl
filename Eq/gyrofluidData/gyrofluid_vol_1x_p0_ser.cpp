#include <gyrofluid_mod_decl.h>

double gyrofluid_vol_1x_p0_ser(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *w, const double *dx, const double *jac, const double *rBmag, const double *jacDbmag, const double *dBdz, const double *sMom, const double *phi, double *primMom, double *out) 
{ 
  // q_,m_:   species charge and mass.
  // kappa:   heat conductivity coefficients.
  // kperpSq: k_perp^2.
  // w:       cell-center.
  // dx:      cell length.
  // uMaxIn:  maximum speed.
  // jac:     jacobian.
  // rBmag:   reciprocal of magnetic field magnitude (1/B).
  // jacDbmag:jacobian divided by B (J/B).
  // sMom:    stepped moments (times Jacobian).
  // phi:     electrostatic potential.
  // primMom: primitive moments (upar, Tpar, Tperp).
  // out:     output increment.

  double wx = w[0];
  double rdx2 = 2.0/dx[0];
  double rdxSq4 = rdx2*rdx2;


  out[1] += -0.7071067811865475*dBdz[0]*sMom[3]; 



  return fabs((0.4854917717073234*sqrt((2.0*primMom[2])/m_+primMom[1]/m_)+0.7071067811865475*primMom[0])*rdx2) + 1.0*fmax(kappaPar*rdxSq4,kappaPerp*rdxSq4); 
}
