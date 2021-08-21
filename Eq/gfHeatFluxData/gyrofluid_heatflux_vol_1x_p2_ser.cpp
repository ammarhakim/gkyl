#include <gyrofluid_heatflux_mod_decl.h>

double gyrofluid_heatflux_vol_1x_p2_ser(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *w, const double *dx, const double *jac, const double *rBmag, const double *rBmagSq, const double *sMom, const double *phi, double *primMom, double *out) 
{ 
  // q_,m_:   species charge and mass.
  // kappa:   heat conductivity coefficients.
  // kperpSq: k_perp^2.
  // w:       cell-center.
  // dx:      cell length.
  // jac:     jacobian.
  // rBmag:   reciprocal of magnetic field magnitude (1/B).
  // rBmagSq: rBmag^2.
  // sMom:    stepped moments (times Jacobian).
  // phi:     electrostatic potential.
  // primMom: primitive moments (upar, Tpar, Tperp).
  // out:     output increment.

  double rdx2 = 2.0/dx[0];
  double rdxSq4 = rdx2*rdx2;

  out[7] += (kappaPerp*(((-2.371708245126284*phi[1]*rBmag[2])-1.060660171779821*phi[0]*rBmag[1])*kperpSq*q_+4.743416490252569*rBmag[2]*primMom[7]+2.121320343559642*rBmag[1]*primMom[6])+(2.371708245126284*rBmag[2]*primMom[4]+1.060660171779821*rBmag[1]*primMom[3])*kappaPar)*rdxSq4; 
  out[8] += (kappaPerp*(((-5.303300858899105*phi[0]*rBmag[2])-4.743416490252569*phi[1]*rBmag[1]-2.371708245126284*phi[0]*rBmag[0])*kperpSq*q_+rBmag[2]*(14.23024947075771*primMom[8]-7.115124735378852*phi[2]*kperpSq*q_)+9.48683298050514*rBmag[1]*primMom[7]+(10.60660171779821*rBmag[2]+4.743416490252569*rBmag[0])*primMom[6])+(7.115124735378852*rBmag[2]*primMom[5]+4.743416490252569*rBmag[1]*primMom[4]+(5.303300858899105*rBmag[2]+2.371708245126284*rBmag[0])*primMom[3])*kappaPar)*rdxSq4; 

  out[10] += kappaPerp*(((-2.371708245126284*phi[1]*rBmagSq[2])-1.060660171779821*phi[0]*rBmagSq[1])*kperpSq*q_+4.743416490252569*rBmagSq[2]*primMom[7]+2.121320343559642*rBmagSq[1]*primMom[6])*rdxSq4; 
  out[11] += kappaPerp*((((-7.115124735378852*phi[2])-5.303300858899105*phi[0])*rBmagSq[2]-4.743416490252569*phi[1]*rBmagSq[1]-2.371708245126284*phi[0]*rBmagSq[0])*kperpSq*q_+14.23024947075771*rBmagSq[2]*primMom[8]+9.48683298050514*rBmagSq[1]*primMom[7]+(10.60660171779821*rBmagSq[2]+4.743416490252569*rBmagSq[0])*primMom[6])*rdxSq4; 

  return 1.8*fmax(((0.7071067811865475*rBmag[0]-0.7905694150420947*rBmag[2])*kappaPar*m_*rdxSq4)/(0.7071067811865475*sMom[0]-0.7905694150420947*sMom[2]),((0.7071067811865475*rBmag[0]-0.7905694150420947*rBmag[2])*kappaPerp*m_*rdxSq4)/(0.7071067811865475*sMom[0]-0.7905694150420947*sMom[2])); 
}
