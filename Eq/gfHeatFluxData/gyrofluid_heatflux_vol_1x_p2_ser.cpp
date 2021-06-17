#include <gyrofluid_heatflux_mod_decl.h>

double gyrofluid_heatflux_vol_1x_p2_ser(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *w, const double *dx, const double *jac, const double *jacDbmag, const double *sMom, const double *phi, double *primMom, double *out) 
{ 
  // q_,m_:   species charge and mass.
  // kappa:   heat conductivity coefficients.
  // kperpSq: k_perp^2.
  // w:       cell-center.
  // dx:      cell length.
  // jac:     jacobian.
  // jacDbmag:jacobian divided by B (J/B).
  // sMom:    stepped moments (times Jacobian).
  // phi:     electrostatic potential.
  // primMom: primitive moments (upar, Tpar, Tperp).
  // out:     output increment.

  double rdx2 = 2.0/dx[0];
  double rdxSq4 = rdx2*rdx2;

  out[7] += (kappaPerp*(((-2.371708245126284*phi[1]*jac[2])-1.060660171779821*phi[0]*jac[1])*kperpSq*q_+4.743416490252569*jac[2]*primMom[7]+2.121320343559642*jac[1]*primMom[6])+(2.371708245126284*jac[2]*primMom[4]+1.060660171779821*jac[1]*primMom[3])*kappaPar)*rdxSq4; 
  out[8] += (kappaPerp*((jac[2]*((-7.115124735378852*phi[2])-5.303300858899105*phi[0])-4.743416490252569*jac[1]*phi[1]-2.371708245126284*jac[0]*phi[0])*kperpSq*q_+14.23024947075771*jac[2]*primMom[8]+9.48683298050514*jac[1]*primMom[7]+(10.60660171779821*jac[2]+4.743416490252569*jac[0])*primMom[6])+(7.115124735378852*jac[2]*primMom[5]+4.743416490252569*jac[1]*primMom[4]+(5.303300858899105*jac[2]+2.371708245126284*jac[0])*primMom[3])*kappaPar)*rdxSq4; 

  out[10] += kappaPerp*(((-2.371708245126284*phi[1]*jacDbmag[2])-1.060660171779821*phi[0]*jacDbmag[1])*kperpSq*q_+4.743416490252569*jacDbmag[2]*primMom[7]+2.121320343559642*jacDbmag[1]*primMom[6])*rdxSq4; 
  out[11] += kappaPerp*((jacDbmag[2]*((-7.115124735378852*phi[2])-5.303300858899105*phi[0])-4.743416490252569*jacDbmag[1]*phi[1]-2.371708245126284*jacDbmag[0]*phi[0])*kperpSq*q_+14.23024947075771*jacDbmag[2]*primMom[8]+9.48683298050514*jacDbmag[1]*primMom[7]+(10.60660171779821*jacDbmag[2]+4.743416490252569*jacDbmag[0])*primMom[6])*rdxSq4; 

  return 1.8*fmax(((0.7071067811865475*jac[0]-0.7905694150420947*jac[2])*kappaPar*m_*rdxSq4)/(0.7071067811865475*sMom[0]-0.7905694150420947*sMom[2]),((0.7071067811865475*jac[0]-0.7905694150420947*jac[2])*kappaPerp*m_*rdxSq4)/(0.7071067811865475*sMom[0]-0.7905694150420947*sMom[2])); 
}
