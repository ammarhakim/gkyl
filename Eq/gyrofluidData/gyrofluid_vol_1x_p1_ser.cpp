#include <gyrofluid_mod_decl.h>

double gyrofluid_vol_1x_p1_ser(const double q_, const double m_, const double *w, const double *dx, const double *rJac, const double *rBmag, const double *dBdz, const double *sMom, const double *phi, const double *primMom, const double *cRus, double *out) 
{ 
  // q_,m_:   species charge and mass.
  // w:       cell-center.
  // dx:      cell length.
  // uMaxIn:  maximum speed.
  // rJac:    reciprocal of jacobian (1/B).
  // rBmag:   reciprocal of magnetic field magnitude (1/B).
  // sMom:    stepped moments (times Jacobian).
  // phi:     electrostatic potential.
  // primMom: primitive moments (upar, Tpar, Tperp).
  // cRus:    phase speed in Rusanov numerical flux.
  // out:     output increment.

  double wx = w[0];
  double rdx2 = 2.0/dx[0];
  double rdxSq4 = rdx2*rdx2;

  out[1] += 0.8660254037844386*((rBmag[0]*rJac[1]+rJac[0]*rBmag[1])*sMom[3]+(rBmag[1]*rJac[1]+rBmag[0]*rJac[0])*sMom[2])*rdx2; 

  out[2] += (-(1.224744871391589*sMom[0]*phi[1]*q_*rdx2)/m_)-0.7071067811865475*(dBdz[1]*sMom[7]+dBdz[0]*sMom[6]); 
  out[3] += (1.732050807568877*((rBmag[0]*rJac[1]+rJac[0]*rBmag[1])*sMom[5]+(rBmag[1]*rJac[1]+rBmag[0]*rJac[0])*sMom[4])-(1.224744871391589*(phi[1]*sMom[1]*q_+((rBmag[1]*rJac[1]+rBmag[0]*rJac[0])*sMom[1]+sMom[0]*(rBmag[0]*rJac[1]+rJac[0]*rBmag[1]))*primMom[5]+((rBmag[0]*rJac[1]+rJac[0]*rBmag[1])*sMom[1]+sMom[0]*(rBmag[1]*rJac[1]+rBmag[0]*rJac[0]))*primMom[4]))/m_)*rdx2-0.7071067811865475*(dBdz[0]*sMom[7]+dBdz[1]*sMom[6]); 

  out[4] += -(1.224744871391589*phi[1]*sMom[2]*q_*rdx2)/m_; 
  out[5] += (((-1.224744871391589*phi[1]*sMom[3]*q_)+(0.4330127018922193*((primMom[0]*rBmag[1]+rBmag[0]*primMom[1])*rJac[1]+rJac[0]*(primMom[1]*rBmag[1]+primMom[0]*rBmag[0]))*sMom[1]+sMom[0]*(0.4330127018922193*(primMom[0]*rBmag[0]*rJac[1]+rJac[0]*(primMom[0]*rBmag[1]+rBmag[0]*primMom[1]))+0.7794228634059945*primMom[1]*rBmag[1]*rJac[1]))*primMom[3]+(0.4330127018922193*((primMom[0]*rBmag[0]*rJac[1]+rJac[0]*(primMom[0]*rBmag[1]+rBmag[0]*primMom[1]))*sMom[1]+sMom[0]*((primMom[0]*rBmag[1]+rBmag[0]*primMom[1])*rJac[1]+rJac[0]*(primMom[1]*rBmag[1]+primMom[0]*rBmag[0])))+0.7794228634059945*primMom[1]*rBmag[1]*rJac[1]*sMom[1])*primMom[2])/m_+1.10227038425243*primMom[1]*rBmag[1]*rJac[1]*sMom[5]+0.6123724356957944*((primMom[0]*rBmag[0]*rJac[1]+rJac[0]*(primMom[0]*rBmag[1]+rBmag[0]*primMom[1]))*sMom[5]+((primMom[0]*rBmag[1]+rBmag[0]*primMom[1])*rJac[1]+rJac[0]*(primMom[1]*rBmag[1]+primMom[0]*rBmag[0]))*sMom[4]))*rdx2; 

  out[7] += (1.10227038425243*primMom[1]*rBmag[1]*rJac[1]*sMom[7]+0.6123724356957944*((primMom[0]*rBmag[0]*rJac[1]+rJac[0]*(primMom[0]*rBmag[1]+rBmag[0]*primMom[1]))*sMom[7]+((primMom[0]*rBmag[1]+rBmag[0]*primMom[1])*rJac[1]+rJac[0]*(primMom[1]*rBmag[1]+primMom[0]*rBmag[0]))*sMom[6]))*rdx2; 

  return rdx2*(fabs(0.7071067811865475*primMom[0]) + 0.7071067811865475*cRus[0]); 
}
