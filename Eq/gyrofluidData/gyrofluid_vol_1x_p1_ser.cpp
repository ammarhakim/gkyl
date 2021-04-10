#include <gyrofluid_mod_decl.h>

double gyrofluid_vol_1x_p1_ser(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *w, const double *dx, const double *jac, const double *rBmag, const double *jacDbmag, const double *dBdz, const double *sMom, const double *phi, double *primMom, double *out) 
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

  out[1] += 1.732050807568877*sMom[2]*rdx2; 

  out[2] += (1.224744871391589*phi[0]*sMom[1]*q_*rdx2)/m_-0.7071067811865475*(dBdz[1]*sMom[7]+dBdz[0]*sMom[6]); 
  out[3] += (((2.449489742783178*phi[1]*sMom[1]+1.224744871391589*phi[0]*sMom[0])*q_)/m_+3.464101615137754*sMom[4])*rdx2-0.7071067811865475*(dBdz[0]*sMom[7]+dBdz[1]*sMom[6]); 

  out[4] += (1.224744871391589*phi[0]*sMom[3]*q_*rdx2)/m_; 
  out[5] += (((2.449489742783178*phi[1]*sMom[3]+1.224744871391589*phi[0]*sMom[2])*q_+0.8660254037844386*((primMom[0]*sMom[1]+sMom[0]*primMom[1])*primMom[3]+(primMom[1]*sMom[1]+primMom[0]*sMom[0])*primMom[2]))/m_+1.224744871391589*(primMom[1]*sMom[5]+primMom[0]*sMom[4]))*rdx2; 

  out[7] += 1.224744871391589*(primMom[1]*sMom[7]+primMom[0]*sMom[6])*rdx2; 

  return fabs(0.7071067811865475*primMom[0]*rdx2) + 1.333333333333333*fmax(kappaPar*rdxSq4,kappaPerp*rdxSq4); 
}
