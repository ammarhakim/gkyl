#include <GyrokineticModDecl.h> 
double GyrokineticVol2x0vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaM = 0.0; 
  double alphaP = 0.0; 
  double alphax[4]; 
  alphax[0] = -0.8660254037844386*BmagInv[0]*Phi[2]*dfac_y; 
  alphax[1] = -0.8660254037844386*BmagInv[0]*Phi[3]*dfac_y; 
  out[1] += 0.8660254037844386*(alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[3] += 0.8660254037844386*(alphax[1]*f[3]+alphax[0]*f[2])*dfac_x; 
  alphaM = 0.5*(0.5*alphax[0]-0.8660254037844386*alphax[1])*dfac_x; 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.5*(0.8660254037844386*alphax[1]+0.5*alphax[0])*dfac_x; 
  if(alphaP>0) cflFreq += alphaP; 
  double alphay[4]; 
  alphay[0] = 0.8660254037844386*BmagInv[0]*Phi[1]*dfac_x; 
  alphay[2] = 0.8660254037844386*BmagInv[0]*Phi[3]*dfac_x; 
  out[2] += 0.8660254037844386*(alphay[2]*f[2]+alphay[0]*f[0])*dfac_y; 
  out[3] += 0.8660254037844386*(alphay[2]*f[3]+alphay[0]*f[1])*dfac_y; 
  alphaM = 0.5*(0.5*alphay[0]-0.8660254037844386*alphay[2])*dfac_y; 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.5*(0.8660254037844386*alphay[2]+0.5*alphay[0])*dfac_y; 
  if(alphaP>0) cflFreq += alphaP; 
  return cflFreq; 
} 
double GyrokineticVol2x0vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaM = 0.0; 
  double alphaP = 0.0; 
  double alphax[8]; 
  alphax[0] = -0.8660254037844386*BmagInv[0]*Phi[2]*dfac_y; 
  alphax[1] = -0.8660254037844386*BmagInv[0]*Phi[3]*dfac_y; 
  alphax[2] = -1.936491673103709*BmagInv[0]*Phi[5]*dfac_y; 
  alphax[3] = -1.936491673103709*BmagInv[0]*Phi[7]*dfac_y; 
  alphax[4] = -0.8660254037844386*BmagInv[0]*Phi[6]*dfac_y; 
  out[1] += 0.8660254037844386*(alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[3] += 0.1*(7.745966692414834*alphax[3]*f[7]+8.660254037844387*alphax[4]*f[6]+7.745966692414834*alphax[2]*f[5]+8.660254037844386*alphax[1]*f[3]+8.660254037844386*f[1]*alphax[3]+8.660254037844386*alphax[0]*f[2]+8.660254037844386*f[0]*alphax[2])*dfac_x; 
  out[4] += 0.1*(17.32050807568877*alphax[3]*f[6]+17.32050807568877*alphax[1]*f[4]+17.32050807568877*f[1]*alphax[4]+19.36491673103709*alphax[2]*f[3]+19.36491673103709*f[2]*alphax[3]+19.36491673103709*alphax[0]*f[1]+19.36491673103709*f[0]*alphax[1])*dfac_x; 
  out[6] += 0.1*(17.32050807568877*alphax[2]*f[7]+17.32050807568877*alphax[1]*f[6]+17.32050807568877*alphax[3]*f[5]+17.32050807568877*alphax[3]*f[4]+17.32050807568877*f[3]*alphax[4]+19.36491673103708*alphax[0]*f[3]+19.36491673103708*f[0]*alphax[3]+19.36491673103708*alphax[1]*f[2]+19.36491673103708*f[1]*alphax[2])*dfac_x; 
  out[7] += 0.1*(8.660254037844386*alphax[1]*f[7]+8.660254037844387*alphax[0]*f[5]+7.745966692414834*alphax[3]*f[3]+7.745966692414834*alphax[2]*f[2])*dfac_x; 
  alphaM = 0.25*(0.25*(4.47213595499958*alphax[4]+3.0*alphax[3]-1.732050807568877*alphax[2]-3.464101615137754*alphax[1]+2.0*alphax[0])+0.25*(4.47213595499958*alphax[4]-3.0*alphax[3]+1.732050807568877*alphax[2]-3.464101615137754*alphax[1]+2.0*alphax[0]))*dfac_x; 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.25*(0.25*(4.47213595499958*alphax[4]+3.0*alphax[3]+1.732050807568877*alphax[2]+3.464101615137754*alphax[1]+2.0*alphax[0])+0.25*(4.47213595499958*alphax[4]-3.0*alphax[3]-1.732050807568877*alphax[2]+3.464101615137754*alphax[1]+2.0*alphax[0]))*dfac_x; 
  if(alphaP>0) cflFreq += alphaP; 
  double alphay[8]; 
  alphay[0] = 0.8660254037844386*BmagInv[0]*Phi[1]*dfac_x; 
  alphay[1] = 1.936491673103709*BmagInv[0]*Phi[4]*dfac_x; 
  alphay[2] = 0.8660254037844386*BmagInv[0]*Phi[3]*dfac_x; 
  alphay[3] = 1.936491673103709*BmagInv[0]*Phi[6]*dfac_x; 
  alphay[5] = 0.8660254037844387*BmagInv[0]*Phi[7]*dfac_x; 
  out[2] += 0.8660254037844386*(alphay[5]*f[5]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0])*dfac_y; 
  out[3] += 0.1*(8.660254037844387*alphay[5]*f[7]+7.745966692414834*alphay[3]*f[6]+7.745966692414834*alphay[1]*f[4]+8.660254037844386*alphay[2]*f[3]+8.660254037844386*f[2]*alphay[3]+8.660254037844386*alphay[0]*f[1]+8.660254037844386*f[0]*alphay[1])*dfac_y; 
  out[5] += 0.1*(17.32050807568877*alphay[3]*f[7]+17.32050807568877*alphay[2]*f[5]+17.32050807568877*f[2]*alphay[5]+19.36491673103709*alphay[1]*f[3]+19.36491673103709*f[1]*alphay[3]+19.36491673103709*alphay[0]*f[2]+19.36491673103709*f[0]*alphay[2])*dfac_y; 
  out[6] += 0.1*(8.660254037844386*alphay[2]*f[6]+8.660254037844387*alphay[0]*f[4]+7.745966692414834*alphay[3]*f[3]+7.745966692414834*alphay[1]*f[1])*dfac_y; 
  out[7] += 0.1*(17.32050807568877*alphay[2]*f[7]+17.32050807568877*alphay[1]*f[6]+17.32050807568877*alphay[3]*f[5]+17.32050807568877*f[3]*alphay[5]+17.32050807568877*alphay[3]*f[4]+19.36491673103708*alphay[0]*f[3]+19.36491673103708*f[0]*alphay[3]+19.36491673103708*alphay[1]*f[2]+19.36491673103708*f[1]*alphay[2])*dfac_y; 
  alphaM = 0.25*(0.25*(4.47213595499958*alphay[5]+3.0*alphay[3]-3.464101615137754*alphay[2]-1.732050807568877*alphay[1]+2.0*alphay[0])+0.25*(4.47213595499958*alphay[5]-3.0*alphay[3]-3.464101615137754*alphay[2]+1.732050807568877*alphay[1]+2.0*alphay[0]))*dfac_y; 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.25*(0.25*(4.47213595499958*alphay[5]+3.0*alphay[3]+3.464101615137754*alphay[2]+1.732050807568877*alphay[1]+2.0*alphay[0])+0.25*(4.47213595499958*alphay[5]-3.0*alphay[3]+3.464101615137754*alphay[2]-1.732050807568877*alphay[1]+2.0*alphay[0]))*dfac_y; 
  if(alphaP>0) cflFreq += alphaP; 
  return cflFreq; 
} 
