#include <GyrokineticModDecl.h> 
double GyrokineticVol1x1vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaM = 0.0; 
  double alphaP = 0.0; 
  double alphax[8]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*dfac_x*wv; 
  alphax[2] = (0.8164965809277261*Gradpar[0]*dfac_x)/dfac_v; 
  alphaM = 0.25*alphax[0]; 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.25*alphax[0]; 
  if(alphaP>0) cflFreq += alphaP; 
  double alphav[8]; 
  alphav[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  alphav[1] = -(3.872983346207417*Gradpar[0]*Phi[2]*dfac_v*dfac_x*q_)/m_; 
  alphaM = 0.25*alphav[0]; 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.25*alphav[0]; 
  if(alphaP>0) cflFreq += alphaP; 
  out[1] += 0.8660254037844386*(alphax[2]*f[2]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphav[1]*f[1]+alphav[0]*f[0]); 
  out[3] += 0.1732050807568877*(4.47213595499958*(alphax[2]*f[5]+alphav[1]*f[4])+5.0*(alphax[0]*f[2]+f[0]*alphax[2]+alphav[0]*f[1]+f[0]*alphav[1])); 
  out[4] += 1.936491673103709*(alphax[2]*f[3]+alphax[0]*f[1]); 
  out[5] += 1.936491673103709*(alphav[1]*f[3]+alphav[0]*f[2]); 
  out[6] += 0.1*(17.32050807568877*alphax[2]*f[7]+8.660254037844387*alphav[0]*f[4]+19.36491673103708*alphax[0]*f[3]+f[1]*(19.36491673103708*alphax[2]+7.745966692414834*alphav[1])); 
  out[7] += 0.1*(17.32050807568877*alphav[1]*f[6]+8.660254037844387*alphax[0]*f[5]+19.36491673103708*alphav[0]*f[3]+(7.745966692414834*alphax[2]+19.36491673103708*alphav[1])*f[2]); 
  return cflFreq; 
} 
double GyrokineticVol1x1vSerP2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaM = 0.0; 
  double alphaP = 0.0; 
  double alphax[8]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*dfac_x*wv; 
  alphax[1] = 1.414213562373095*Gradpar[1]*dfac_x*wv; 
  alphax[2] = (0.8164965809277261*Gradpar[0]*dfac_x)/dfac_v; 
  alphax[3] = (0.8164965809277261*Gradpar[1]*dfac_x)/dfac_v; 
  alphax[4] = 1.414213562373095*Gradpar[2]*dfac_x*wv; 
  alphax[6] = (0.816496580927726*Gradpar[2]*dfac_x)/dfac_v; 
  alphaM = 0.25*(2.23606797749979*alphax[4]-1.732050807568877*alphax[1]+alphax[0]); 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.25*(2.23606797749979*alphax[4]+1.732050807568877*alphax[1]+alphax[0]); 
  if(alphaP>0) cflFreq += alphaP; 
  double alphav[8]; 
  alphav[0] = -(1.732050807568877*(2.23606797749979*Gradpar[1]*Phi[2]+Gradpar[0]*Phi[1])*dfac_v*dfac_x*q_)/m_; 
  alphav[1] = -(1.732050807568877*((2.0*Gradpar[2]+2.23606797749979*Gradpar[0])*Phi[2]+Gradpar[1]*Phi[1])*dfac_v*dfac_x*q_)/m_; 
  alphav[4] = -(1.732050807568877*(2.0*Gradpar[1]*Phi[2]+Phi[1]*Gradpar[2])*dfac_v*dfac_x*q_)/m_; 
  alphaM = 0.25*alphav[0]; 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.25*alphav[0]; 
  if(alphaP>0) cflFreq += alphaP; 
  out[1] += 0.8660254037844386*(alphax[6]*f[6]+alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphav[4]*f[4]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[3] += 0.1*(7.745966692414834*alphax[3]*f[7]+8.660254037844387*(alphax[4]*f[6]+f[4]*alphax[6])+7.745966692414834*(alphax[2]*f[5]+alphav[1]*f[4]+f[1]*alphav[4])+8.660254037844386*(alphax[1]*f[3]+f[1]*alphax[3]+alphax[0]*f[2]+f[0]*alphax[2]+alphav[0]*f[1]+f[0]*alphav[1])); 
  out[4] += 0.1*(17.32050807568877*(alphax[3]*f[6]+f[3]*alphax[6])+17.32050807568877*(alphax[1]*f[4]+f[1]*alphax[4])+19.36491673103709*(alphax[2]*f[3]+f[2]*alphax[3]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[5] += 1.936491673103709*(alphav[4]*f[6]+alphav[1]*f[3]+alphav[0]*f[2]); 
  out[6] += 0.01428571428571429*(108.4435336938077*alphax[6]*f[7]+121.2435565298214*(alphax[2]*f[7]+alphax[1]*f[6]+f[1]*alphax[6])+121.2435565298214*alphax[3]*f[5]+(38.72983346207417*alphav[4]+121.2435565298214*alphax[3]+60.62177826491071*alphav[0])*f[4]+121.2435565298214*f[3]*alphax[4]+60.62177826491071*f[0]*alphav[4]+135.5544171172596*(alphax[0]*f[3]+f[0]*alphax[3]+alphax[1]*f[2])+f[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphav[1])); 
  out[7] += 0.1*(8.660254037844386*alphax[1]*f[7]+(7.745966692414834*alphax[6]+17.32050807568877*alphav[1])*f[6]+8.660254037844387*alphax[0]*f[5]+f[3]*(17.32050807568877*alphav[4]+7.745966692414834*alphax[3]+19.36491673103708*alphav[0])+(7.745966692414834*alphax[2]+19.36491673103708*alphav[1])*f[2]); 
  return cflFreq; 
} 
