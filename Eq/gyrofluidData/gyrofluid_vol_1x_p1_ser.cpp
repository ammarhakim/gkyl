#include <gyrofluid_mod_decl.h>

double gyrofluid_vol_1x_p1_ser(const double q_, const double m_, const double *w, const double *dx, const double *jac, const double *rBmag, const double *jacDbmag, const double *dBdz, const double *sMom, const double *phi, const double *primMom, const double *cs, double *out) 
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

  out[1] += 1.732050807568877*sMom[2]*rdx2; 

  out[2] += (-(1.224744871391589*sMom[0]*phi[1]*q_*rdx2)/m_)-0.7071067811865475*(dBdz[1]*sMom[7]+dBdz[0]*sMom[6]); 
  out[3] += (((-1.224744871391589*phi[1]*sMom[1]*q_)-2.449489742783178*(sMom[1]*primMom[5]+sMom[0]*primMom[4]))/m_+3.464101615137754*sMom[4])*rdx2-0.7071067811865475*(dBdz[0]*sMom[7]+dBdz[1]*sMom[6]); 

  out[4] += -(1.224744871391589*phi[1]*sMom[2]*q_*rdx2)/m_; 
  out[5] += ((0.8660254037844386*((primMom[0]*sMom[1]+sMom[0]*primMom[1])*primMom[3]+(primMom[1]*sMom[1]+primMom[0]*sMom[0])*primMom[2])-1.224744871391589*phi[1]*sMom[3]*q_)/m_+1.224744871391589*(primMom[1]*sMom[5]+primMom[0]*sMom[4]))*rdx2; 

  out[7] += 1.224744871391589*(primMom[1]*sMom[7]+primMom[0]*sMom[6])*rdx2; 

  return rdx2*(fabs(0.7071067811865475*primMom[0]) + 0.7071067811865475*cs[0]); 
}
