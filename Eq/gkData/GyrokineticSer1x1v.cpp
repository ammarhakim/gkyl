#include <GyrokineticModDecl.h> 
double GyrokineticVol1x1vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  out[1] += 1.732050807568877*f[0]*dfac_x*wv; 
  out[2] += -(2.121320343559642*f[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  out[3] += (0.7071067811865475*(2.449489742783178*f[2]*dfac_x*m_*wv-3.0*Phi[1]*f[1]*dfac_v*dfac_x*q_))/m_; 
  double cflFreq = 0.0; 
  cflFreq += fabs(0.7071067811865475*Bmag[0]*wv)*dxInv; 
  cflFreq += fabs(-(0.8660254037844386*Bmag[0]*Phi[1]*dfac_x*q_)/m_)*dvInv; 
  return cflFreq; 
} 
