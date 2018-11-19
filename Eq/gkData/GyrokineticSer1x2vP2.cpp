#include <GyrokineticModDecl.h> 
double GyrokineticVol1x2vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[20]; 
  alphax[0] = 2.0*Gradpar[0]*dfac_x*wv; 
  alphax[2] = (1.154700538379252*Gradpar[0]*dfac_x)/dfac_v; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.1767766952966368*alphax[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.1767766952966368*alphax[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0441941738241592*alphax[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0441941738241592*alphax[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0441941738241592*alphax[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0441941738241592*alphax[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphav[20]; 
  alphav[0] = -(2.449489742783178*Gradpar[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  alphav[1] = -(5.477225575051662*Gradpar[0]*Phi[2]*dfac_v*dfac_x*q_)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.1767766952966368*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.1767766952966368*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0441941738241592*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0441941738241592*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0441941738241592*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0441941738241592*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  out[1] += 0.6123724356957944*(alphax[2]*f[2]+alphax[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphav[1]*f[1]+alphav[0]*f[0]); 
  out[4] += 0.1224744871391589*(4.47213595499958*(alphax[2]*f[8]+alphav[1]*f[7])+5.0*(alphax[0]*f[2]+f[0]*alphax[2]+alphav[0]*f[1]+f[0]*alphav[1])); 
  out[5] += 0.6123724356957944*(alphax[2]*f[6]+alphax[0]*f[3]); 
  out[6] += 0.6123724356957944*(alphav[1]*f[5]+alphav[0]*f[3]); 
  out[7] += 1.369306393762915*(alphax[2]*f[4]+alphax[0]*f[1]); 
  out[8] += 1.369306393762915*(alphav[1]*f[4]+alphav[0]*f[2]); 
  out[10] += 0.07071067811865474*(7.745966692414834*(alphax[2]*f[14]+alphav[1]*f[13])+8.660254037844386*(alphax[0]*f[6]+alphav[0]*f[5]+(alphax[2]+alphav[1])*f[3])); 
  out[11] += 0.07071067811865474*(17.32050807568877*alphax[2]*f[12]+8.660254037844387*alphav[0]*f[7]+19.36491673103708*alphax[0]*f[4]+f[1]*(19.36491673103708*alphax[2]+7.745966692414834*alphav[1])); 
  out[12] += 0.07071067811865474*(17.32050807568877*alphav[1]*f[11]+8.660254037844387*alphax[0]*f[8]+19.36491673103708*alphav[0]*f[4]+(7.745966692414834*alphax[2]+19.36491673103708*alphav[1])*f[2]); 
  out[13] += 1.369306393762915*(alphax[2]*f[10]+alphax[0]*f[5]); 
  out[14] += 1.369306393762915*(alphav[1]*f[10]+alphav[0]*f[6]); 
  out[15] += 0.07071067811865474*(8.660254037844386*alphax[2]*f[16]+8.660254037844387*alphax[0]*f[9]); 
  out[16] += 0.07071067811865474*(8.660254037844386*alphav[1]*f[15]+8.660254037844387*alphav[0]*f[9]); 
  out[17] += 0.07071067811865474*(17.32050807568877*alphax[2]*f[18]+8.660254037844387*alphav[0]*f[13]+19.36491673103709*alphax[0]*f[10]+(19.36491673103709*alphax[2]+7.745966692414834*alphav[1])*f[5]); 
  out[18] += 0.07071067811865474*(17.32050807568877*alphav[1]*f[17]+8.660254037844387*alphax[0]*f[14]+19.36491673103709*alphav[0]*f[10]+(7.745966692414834*alphax[2]+19.36491673103709*alphav[1])*f[6]); 
  out[19] += 0.07071067811865474*(8.660254037844387*(alphax[0]*f[16]+alphav[0]*f[15])+8.660254037844386*(alphax[2]+alphav[1])*f[9]); 
  return cflFreq; 
} 
double GyrokineticVol1x2vSerP2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[20]; 
  alphax[0] = -1.0*((1.732050807568877*(Bmag[2]*(2.0*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2])+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[2]*geoY[2]+BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))*dfac_x2*m_*wv2)/q_-2.0*Gradpar[0]*dfac_x*wv+(0.5773502691896258*(Bmag[2]*(2.0*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2])+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[2]*geoY[2]+BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))*dfac_x2*m_)/(dfac_v2*q_)); 
  alphax[1] = (1.732050807568877*((-1.0*(2.0*(0.4472135954999579*Bmag[1]*BmagInv[1]*geoY[2]+(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2])+Bmag[2]*(2.0*BmagInv[0]*geoY[2]+4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))-3.51382110749967*Bmag[2]*BmagInv[2]*geoY[2])*dfac_x2*m_*wv2)/q_+2.0*Gradpar[1]*dfac_x*wv+(((-1.0*(0.5773502691896258*(0.8944271909999159*Bmag[1]*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2])+2.23606797749979*BmagInv[0]*geoY[0]*Bmag[2]+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[2]*(1.154700538379252*(BmagInv[0]*geoY[2]+geoY[0]*BmagInv[2])+2.32379000772445*BmagInv[1]*geoY[1])))-2.028705562299123*Bmag[2]*BmagInv[2]*geoY[2])*dfac_x2*m_)/(dfac_v2*q_); 
  alphax[2] = -(2.0*(((Bmag[2]*(2.0*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2])+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[2]*geoY[2]+BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))*dfac_x2*m_*wv)/q_-0.5773502691896258*Gradpar[0]*dfac_x))/dfac_v; 
  alphax[4] = (2.0*(0.5773502691896258*Gradpar[1]*dfac_x-(1.0*(2.0*(0.4472135954999579*Bmag[1]*BmagInv[1]*geoY[2]+(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2])+Bmag[2]*(2.0*BmagInv[0]*geoY[2]+4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*dfac_x2*m_*wv)/q_)-(7.027642214999339*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x2*m_*wv)/q_)/dfac_v; 
  alphax[7] = 2.0*(Gradpar[2]*dfac_x*wv-(1.732050807568877*(BmagInv[1]*(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])+BmagInv[0]*geoY[1]*Bmag[2])*dfac_x2*m_*wv2)/q_)+(dfac_x2*m_*(1.732050807568877*((2.23606797749979*((-0.2857142857142857*Bmag[1]*BmagInv[2])-1.571428571428571*BmagInv[1]*Bmag[2])-1.0*BmagInv[0]*Bmag[1])*geoY[2]+((-3.51382110749967*geoY[1]*Bmag[2])-1.0*geoY[0]*Bmag[1])*BmagInv[2])*wv2+(0.5773502691896258*(2.23606797749979*((-1.571428571428571*Bmag[2]*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2]))-0.2857142857142857*Bmag[1]*BmagInv[2]*geoY[2])-1.0*(Bmag[1]*(BmagInv[0]*geoY[2]+geoY[0]*BmagInv[2])+2.0*(BmagInv[1]*(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])+BmagInv[0]*geoY[1]*Bmag[2]))))/dfac_v2))/q_; 
  alphax[8] = -(1.154700538379252*(Bmag[2]*(0.8944271909999159*BmagInv[1]*geoY[2]+geoY[1]*(0.8944271909999159*BmagInv[2]+BmagInv[0])+geoY[0]*BmagInv[1])+0.4472135954999579*Bmag[1]*(BmagInv[2]*geoY[2]+BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))*dfac_x2*m_)/(dfac_v2*q_); 
  alphax[11] = (3.872983346207417*((0.5773502691896258*(((-0.5714285714285714*Bmag[1]*BmagInv[2])-3.142857142857143*BmagInv[1]*Bmag[2]-0.8944271909999159*BmagInv[0]*Bmag[1])*geoY[2]+((-3.142857142857143*geoY[1]*Bmag[2])-0.8944271909999159*geoY[0]*Bmag[1])*BmagInv[2]+BmagInv[1]*((-1.788854381999832*geoY[0]*Bmag[2])-0.8*Bmag[1]*geoY[1])-1.788854381999832*BmagInv[0]*geoY[1]*Bmag[2])*dfac_x2*m_*wv)/q_+0.2981423969999719*Gradpar[2]*dfac_x))/dfac_v; 
  alphax[12] = (0.2581988897471611*((-2.0*(2.0*(0.4472135954999579*Bmag[1]*BmagInv[1]*geoY[2]+(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2])+Bmag[2]*(2.0*BmagInv[0]*geoY[2]+4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))-7.027642214999339*Bmag[2]*BmagInv[2]*geoY[2])*dfac_x2*m_)/(dfac_v2*q_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.125*(3.16227766016838*alphax[7]-2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(3.16227766016838*alphax[7]+2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*((-0.5477225575051661*alphax[12])-1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.6846531968814574*alphax[12]-0.3952847075210473*alphax[8]+0.7905694150420947*alphax[7]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.5477225575051661*alphax[12])+1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.5477225575051661*alphax[12])-1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.5477225575051661*alphax[12])+1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.5477225575051661*alphax[12])-1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.6846531968814574*alphax[12]-0.3952847075210473*alphax[8]+0.7905694150420947*alphax[7]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.5477225575051661*alphax[12])+1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(0.5477225575051661*alphax[12]-1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*((-0.6846531968814574*alphax[12])-0.3952847075210473*alphax[8]+0.7905694150420947*alphax[7]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.5477225575051661*alphax[12]+1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.5477225575051661*alphax[12]-1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.5477225575051661*alphax[12]+1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.5477225575051661*alphax[12]-1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*((-0.6846531968814574*alphax[12])-0.3952847075210473*alphax[8]+0.7905694150420947*alphax[7]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.5477225575051661*alphax[12]+1.060660171779821*alphax[11]+0.3162277660168379*alphax[8]+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphav[20]; 
  alphav[0] = dfac_v*(dfac_x2*(3.0*((geoY[0]*Bmag[1]*(3.16227766016838*BmagInv[1]*Bmag[2]+0.7071067811865475*BmagInv[0]*Bmag[1])*wm)/q_+Phi[2]*(1.414213562373095*(2.23606797749979*BmagInv[0]*Bmag[2]*geoY[2]+(2.23606797749979*geoY[0]*Bmag[2]+Bmag[1]*geoY[1])*BmagInv[2])+0.7071067811865475*((9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+2.23606797749979*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))+(1.414213562373095*Bmag[1]*BmagInv[1]*Phi[2]+Phi[1]*(0.7071067811865475*Bmag[1]*BmagInv[2]+1.414213562373095*BmagInv[1]*Bmag[2]))*geoY[2]+Phi[1]*(1.414213562373095*geoY[1]*Bmag[2]*BmagInv[2]+0.7071067811865475*(2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))))+0.7071067811865475*((3.0*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*wm)/q_+23.57142857142857*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2])+((3.0*(Bmag[2]*(BmagInv[0]*(Bmag[2]*(3.16227766016838*geoY[2]+3.535533905932737*geoY[0])+3.16227766016838*Bmag[1]*geoY[1])+1.414213562373095*(2.0*Bmag[1]*BmagInv[1]*geoY[2]+(2.23606797749979*geoY[0]*Bmag[2]+2.0*Bmag[1]*geoY[1])*BmagInv[2])+6.363961030678928*BmagInv[1]*geoY[1]*Bmag[2])+0.7071067811865475*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2])+16.66751698511148*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2])*wm)/q_)*wv-(2.449489742783178*dfac_x*((2.23606797749979*Gradpar[1]*Bmag[2]+Gradpar[0]*Bmag[1])*wm+((2.23606797749979*Gradpar[1]*Phi[2]+Gradpar[0]*Phi[1])*q2)/q_))/m_); 
  alphav[1] = dfac_v*(dfac_x2*(Bmag[2]*(1.414213562373095*(3.0*((4.024922359499621*Bmag[1]*BmagInv[1]*geoY[1]*wm)/q_+BmagInv[0]*Phi[1]*geoY[2])+17.24966725499838*geoY[1]*BmagInv[2]*Phi[2])+7.453940198968323*Phi[1]*BmagInv[2]*geoY[2])+3.0*(Bmag[1]*(((3.16227766016838*BmagInv[0]*geoY[0]*Bmag[2]+0.7071067811865475*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*wm)/q_+1.414213562373095*(BmagInv[0]*Phi[2]+0.4472135954999579*BmagInv[1]*Phi[1])*geoY[2])+((1.414213562373095*Bmag[1]*(2.0*geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2]+6.363961030678928*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]*Bmag[2])*wm)/q_+(1.414213562373095*geoY[0]*Bmag[1]*BmagInv[2]+0.7071067811865475*(9.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])))*Phi[2]+Phi[1]*(1.414213562373095*(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2]+0.7071067811865475*((4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])*Bmag[2]+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))))+2.23606797749979*((10.90964748116388*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*wm)/q_+(3.333503397022295*Bmag[1]*BmagInv[2]+10.90964748116388*BmagInv[1]*Bmag[2])*Phi[2]*geoY[2])+(1.414213562373095*(2.23606797749979*Bmag[2]*(4.714285714285714*Bmag[1]*BmagInv[2]+7.714285714285714*BmagInv[1]*Bmag[2])+3.0*Bmag[1]*(2.0*BmagInv[0]*Bmag[2]+0.4472135954999579*Bmag[1]*BmagInv[1]))*geoY[2]*wm)/q_)*wv-(2.449489742783178*dfac_x*((Bmag[2]*(2.0*Gradpar[2]+2.23606797749979*Gradpar[0])+Bmag[1]*Gradpar[1])*wm+(((2.0*Gradpar[2]+2.23606797749979*Gradpar[0])*Phi[2]+Gradpar[1]*Phi[1])*q2)/q_))/m_); 
  alphav[2] = 1.732050807568877*dfac_x2*(((Bmag[2]*(1.414213562373095*(2.0*Bmag[1]*BmagInv[1]*geoY[2]+(2.23606797749979*geoY[0]*Bmag[2]+2.0*Bmag[1]*geoY[1])*BmagInv[2])+Bmag[2]*(3.16227766016838*BmagInv[0]*geoY[2]+0.7071067811865475*(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])))+Bmag[1]*(0.7071067811865475*Bmag[1]*BmagInv[2]*geoY[2]+3.16227766016838*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+0.7071067811865475*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))+5.555838995037158*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2])*wm)/q_+((5.555838995037158*Bmag[2]*BmagInv[2]+1.414213562373095*(2.23606797749979*BmagInv[0]*Bmag[2]+Bmag[1]*BmagInv[1]))*Phi[2]+Phi[1]*(0.7071067811865475*Bmag[1]*BmagInv[2]+1.414213562373095*BmagInv[1]*Bmag[2]))*geoY[2]+(1.414213562373095*(2.23606797749979*geoY[0]*Bmag[2]+Bmag[1]*geoY[1])*BmagInv[2]+0.7071067811865475*((9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+2.23606797749979*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*Phi[2]+Phi[1]*(1.414213562373095*geoY[1]*Bmag[2]*BmagInv[2]+0.7071067811865475*(2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))); 
  alphav[3] = (dfac_v*((1.732050807568877*(Bmag[2]*(1.414213562373095*(2.0*Bmag[1]*BmagInv[1]*geoY[2]+(2.23606797749979*geoY[0]*Bmag[2]+2.0*Bmag[1]*geoY[1])*BmagInv[2])+Bmag[2]*(3.16227766016838*BmagInv[0]*geoY[2]+0.7071067811865475*(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])))+Bmag[1]*(0.7071067811865475*Bmag[1]*BmagInv[2]*geoY[2]+3.16227766016838*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+0.7071067811865475*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))+5.555838995037158*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2])*dfac_x2*wv)/q_-(1.414213562373095*(2.23606797749979*Gradpar[1]*Bmag[2]+Gradpar[0]*Bmag[1])*dfac_x)/m_))/dfac_m; 
  alphav[4] = 1.732050807568877*dfac_x2*(((Bmag[1]*(1.414213562373095*Bmag[2]*(3.51382110749967*BmagInv[2]*geoY[2]+4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])+0.7071067811865475*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[2]*(2.828427124746191*(2.874944542499729*BmagInv[1]*Bmag[2]+BmagInv[0]*Bmag[1])*geoY[2]+6.363961030678928*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2])+1.414213562373095*(0.4472135954999579*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]+(5.749889084999459*geoY[1]*Bmag[2]*Bmag[2]+Bmag[1]*(2.0*geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1]))*BmagInv[2]))*wm)/q_+((2.484646732989441*Bmag[1]*BmagInv[2]+1.414213562373095*(5.749889084999459*BmagInv[1]*Bmag[2]+BmagInv[0]*Bmag[1]))*Phi[2]+Phi[1]*(2.484646732989441*Bmag[2]*BmagInv[2]+1.414213562373095*(BmagInv[0]*Bmag[2]+0.4472135954999579*Bmag[1]*BmagInv[1])))*geoY[2]+(1.414213562373095*(5.749889084999459*geoY[1]*Bmag[2]+geoY[0]*Bmag[1])*BmagInv[2]+0.7071067811865475*(9.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])))*Phi[2]+Phi[1]*(1.414213562373095*(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2]+0.7071067811865475*((4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])*Bmag[2]+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))); 
  alphav[5] = (dfac_v*((1.732050807568877*(Bmag[1]*(1.414213562373095*Bmag[2]*(3.51382110749967*BmagInv[2]*geoY[2]+4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])+0.7071067811865475*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[2]*(2.828427124746191*(2.874944542499729*BmagInv[1]*Bmag[2]+BmagInv[0]*Bmag[1])*geoY[2]+6.363961030678928*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2])+1.414213562373095*(0.4472135954999579*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]+(5.749889084999459*geoY[1]*Bmag[2]*Bmag[2]+Bmag[1]*(2.0*geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1]))*BmagInv[2]))*dfac_x2*wv)/q_-(1.414213562373095*(Bmag[2]*(2.0*Gradpar[2]+2.23606797749979*Gradpar[0])+Bmag[1]*Gradpar[1])*dfac_x)/m_))/dfac_m; 
  alphav[6] = ((Bmag[2]*(1.414213562373095*(2.0*Bmag[1]*BmagInv[1]*geoY[2]+(2.23606797749979*geoY[0]*Bmag[2]+2.0*Bmag[1]*geoY[1])*BmagInv[2])+Bmag[2]*(3.16227766016838*BmagInv[0]*geoY[2]+0.7071067811865475*(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])))+0.7071067811865475*(7.857142857142857*Bmag[2]*Bmag[2]+Bmag[1]*Bmag[1])*BmagInv[2]*geoY[2]+Bmag[1]*(3.16227766016838*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+0.7071067811865475*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))*dfac_x2)/(dfac_m*q_); 
  alphav[7] = dfac_v*(1.414213562373095*(3.0*Phi[1]*(BmagInv[1]*(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])+BmagInv[0]*geoY[1]*Bmag[2])*dfac_x2*wv-(1.732050807568877*dfac_x*((Bmag[1]*Gradpar[2]+2.0*Gradpar[1]*Bmag[2])*wm+((2.0*Gradpar[1]*Phi[2]+Phi[1]*Gradpar[2])*q2)/q_))/m_)+dfac_x2*(((BmagInv[2]*(Bmag[2]*Bmag[2]*(27.10523708715755*geoY[2]+16.66751698511148*geoY[0])+Bmag[1]*(14.90788039793665*geoY[1]*Bmag[2]+2.121320343559642*geoY[0]*Bmag[1]))+(1.355261854357877*Bmag[1]*Bmag[1]*BmagInv[2]+16.66751698511148*BmagInv[0]*Bmag[2]*Bmag[2]+Bmag[1]*(14.90788039793665*BmagInv[1]*Bmag[2]+2.121320343559642*BmagInv[0]*Bmag[1]))*geoY[2]+1.414213562373095*(3.0*(2.23606797749979*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]+Bmag[1]*(BmagInv[1]*(2.0*geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])+2.0*BmagInv[0]*geoY[1]*Bmag[2]))+17.24966725499838*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]))*wm)/q_+Phi[2]*(0.7071067811865475*(23.57142857142857*BmagInv[0]*Bmag[2]*geoY[2]+10.54146332249901*Bmag[1]*geoY[1]*BmagInv[2])+1.414213562373095*(3.0*(2.23606797749979*BmagInv[0]*geoY[0]*Bmag[2]+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+17.24966725499838*BmagInv[1]*geoY[1]*Bmag[2]))+0.7071067811865475*(Phi[1]*(10.54146332249901*BmagInv[1]*Bmag[2]+3.0*BmagInv[0]*Bmag[1])*geoY[2]+BmagInv[2]*(23.57142857142857*geoY[0]*Bmag[2]*Phi[2]+Phi[1]*(10.54146332249901*geoY[1]*Bmag[2]+3.0*geoY[0]*Bmag[1])))+2.23606797749979*(Bmag[1]*(3.333503397022295*BmagInv[1]*Phi[2]+0.6060915267313265*Phi[1]*BmagInv[2])+12.12183053462653*Bmag[2]*BmagInv[2]*Phi[2])*geoY[2])*wv); 
  alphav[10] = ((1.414213562373095*((2.23606797749979*Bmag[2]*(1.571428571428571*Bmag[1]*BmagInv[2]+2.571428571428572*BmagInv[1]*Bmag[2])+Bmag[1]*(2.0*BmagInv[0]*Bmag[2]+0.4472135954999579*Bmag[1]*BmagInv[1]))*geoY[2]+(5.749889084999459*geoY[1]*Bmag[2]*Bmag[2]+Bmag[1]*(2.0*geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1]))*BmagInv[2])+6.363961030678928*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]*Bmag[2]+Bmag[1]*(1.414213562373095*(4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])*Bmag[2]+0.7071067811865475*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_x2)/(dfac_m*q_); 
  alphav[11] = 3.872983346207417*dfac_x2*(((BmagInv[2]*(Bmag[2]*Bmag[2]*(4.040610178208845*geoY[2]+2.484646732989441*geoY[0])+Bmag[1]*(2.222335598014864*geoY[1]*Bmag[2]+0.3162277660168379*geoY[0]*Bmag[1]))+(0.2020305089104422*Bmag[1]*Bmag[1]*BmagInv[2]+2.484646732989441*BmagInv[0]*Bmag[2]*Bmag[2]+Bmag[1]*(2.222335598014864*BmagInv[1]*Bmag[2]+0.3162277660168379*BmagInv[0]*Bmag[1]))*geoY[2]+1.414213562373095*((2.571428571428572*BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+Bmag[1]*(BmagInv[1]*(0.8944271909999159*geoY[0]*Bmag[2]+0.2*Bmag[1]*geoY[1])+0.8944271909999159*BmagInv[0]*geoY[1]*Bmag[2])))*wm)/q_+((4.040610178208845*Bmag[2]*BmagInv[2]+1.111167799007432*(2.23606797749979*BmagInv[0]*Bmag[2]+Bmag[1]*BmagInv[1]))*Phi[2]+Phi[1]*(0.2020305089104422*Bmag[1]*BmagInv[2]+0.7071067811865475*(1.571428571428571*BmagInv[1]*Bmag[2]+0.4472135954999579*BmagInv[0]*Bmag[1])))*geoY[2]+(1.111167799007432*(2.23606797749979*geoY[0]*Bmag[2]+Bmag[1]*geoY[1])*BmagInv[2]+1.414213562373095*((2.571428571428572*BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*Bmag[2]+0.4472135954999579*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*Phi[2]+Phi[1]*(0.7071067811865475*(1.571428571428571*geoY[1]*Bmag[2]+0.4472135954999579*geoY[0]*Bmag[1])*BmagInv[2]+1.414213562373095*(BmagInv[1]*(0.4472135954999579*geoY[0]*Bmag[2]+0.2*Bmag[1]*geoY[1])+0.4472135954999579*BmagInv[0]*geoY[1]*Bmag[2]))); 
  alphav[13] = (3.872983346207417*dfac_v*(1.414213562373095*(0.4472135954999579*((2.0*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x2*wv)/q_-(0.5773502691896258*(Bmag[1]*Gradpar[2]+2.0*Gradpar[1]*Bmag[2])*dfac_x)/m_)+(((2.571428571428572*BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+Bmag[1]*BmagInv[1]*(0.8944271909999159*geoY[0]*Bmag[2]+0.2*Bmag[1]*geoY[1]))*dfac_x2*wv)/q_)+((BmagInv[2]*(Bmag[2]*Bmag[2]*(4.040610178208845*geoY[2]+2.484646732989441*geoY[0])+Bmag[1]*(2.222335598014864*geoY[1]*Bmag[2]+0.3162277660168379*geoY[0]*Bmag[1]))+(0.2020305089104422*Bmag[1]*Bmag[1]*BmagInv[2]+2.484646732989441*BmagInv[0]*Bmag[2]*Bmag[2]+Bmag[1]*(2.222335598014864*BmagInv[1]*Bmag[2]+0.3162277660168379*BmagInv[0]*Bmag[1]))*geoY[2])*dfac_x2*wv)/q_))/dfac_m; 
  alphav[17] = (((0.4517539514526256*(20.0*Bmag[2]*Bmag[2]+Bmag[1]*Bmag[1])*BmagInv[2]+5.555838995037158*BmagInv[0]*Bmag[2]*Bmag[2]+Bmag[1]*(4.969293465978882*BmagInv[1]*Bmag[2]+0.7071067811865475*BmagInv[0]*Bmag[1]))*geoY[2]+(5.555838995037158*geoY[0]*Bmag[2]*Bmag[2]+Bmag[1]*(4.969293465978882*geoY[1]*Bmag[2]+0.7071067811865475*geoY[0]*Bmag[1]))*BmagInv[2]+1.414213562373095*(2.23606797749979*(2.571428571428572*BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+Bmag[1]*(BmagInv[1]*(2.0*geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])+2.0*BmagInv[0]*geoY[1]*Bmag[2])))*dfac_x2)/(dfac_m*q_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(2.449489742783178*alphav[2]-1.414213562373095*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(2.449489742783178*alphav[2]+1.414213562373095*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*(0.7348469228349533*alphav[17]-0.4242640687119285*alphav[13]-0.5477225575051661*alphav[11]-1.10227038425243*alphav[10]+0.3162277660168379*alphav[7]+0.8215838362577489*alphav[6]+0.6363961030678926*alphav[5]+0.8215838362577489*alphav[4]-0.4743416490252568*alphav[3]-0.6123724356957944*alphav[2]-0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.9185586535436913*alphav[17])+0.5303300858899104*alphav[13]+0.6846531968814574*alphav[11]-0.3952847075210473*alphav[7]+0.8215838362577489*alphav[6]-0.4743416490252568*alphav[3]-0.6123724356957944*alphav[2]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.7348469228349533*alphav[17]-0.4242640687119285*alphav[13]-0.5477225575051661*alphav[11]+1.10227038425243*alphav[10]+0.3162277660168379*alphav[7]+0.8215838362577489*alphav[6]-0.6363961030678926*alphav[5]-0.8215838362577489*alphav[4]-0.4743416490252568*alphav[3]-0.6123724356957944*alphav[2]+0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.5477225575051661*alphav[11])+0.3162277660168379*alphav[7]+0.8215838362577489*alphav[4]-0.6123724356957944*alphav[2]-0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.5477225575051661*alphav[11])+0.3162277660168379*alphav[7]-0.8215838362577489*alphav[4]-0.6123724356957944*alphav[2]+0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.7348469228349533*alphav[17])+0.4242640687119285*alphav[13]-0.5477225575051661*alphav[11]+1.10227038425243*alphav[10]+0.3162277660168379*alphav[7]-0.8215838362577489*alphav[6]-0.6363961030678926*alphav[5]+0.8215838362577489*alphav[4]+0.4743416490252568*alphav[3]-0.6123724356957944*alphav[2]-0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.9185586535436913*alphav[17]-0.5303300858899104*alphav[13]+0.6846531968814574*alphav[11]-0.3952847075210473*alphav[7]-0.8215838362577489*alphav[6]+0.4743416490252568*alphav[3]-0.6123724356957944*alphav[2]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.7348469228349533*alphav[17])+0.4242640687119285*alphav[13]-0.5477225575051661*alphav[11]-1.10227038425243*alphav[10]+0.3162277660168379*alphav[7]-0.8215838362577489*alphav[6]+0.6363961030678926*alphav[5]-0.8215838362577489*alphav[4]+0.4743416490252568*alphav[3]-0.6123724356957944*alphav[2]+0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*((-0.7348469228349533*alphav[17])-0.4242640687119285*alphav[13]+0.5477225575051661*alphav[11]+1.10227038425243*alphav[10]+0.3162277660168379*alphav[7]-0.8215838362577489*alphav[6]+0.6363961030678926*alphav[5]-0.8215838362577489*alphav[4]-0.4743416490252568*alphav[3]+0.6123724356957944*alphav[2]-0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.9185586535436913*alphav[17]+0.5303300858899104*alphav[13]-0.6846531968814574*alphav[11]-0.3952847075210473*alphav[7]-0.8215838362577489*alphav[6]-0.4743416490252568*alphav[3]+0.6123724356957944*alphav[2]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*((-0.7348469228349533*alphav[17])-0.4242640687119285*alphav[13]+0.5477225575051661*alphav[11]-1.10227038425243*alphav[10]+0.3162277660168379*alphav[7]-0.8215838362577489*alphav[6]-0.6363961030678926*alphav[5]+0.8215838362577489*alphav[4]-0.4743416490252568*alphav[3]+0.6123724356957944*alphav[2]+0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.5477225575051661*alphav[11]+0.3162277660168379*alphav[7]-0.8215838362577489*alphav[4]+0.6123724356957944*alphav[2]-0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.5477225575051661*alphav[11]+0.3162277660168379*alphav[7]+0.8215838362577489*alphav[4]+0.6123724356957944*alphav[2]+0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.7348469228349533*alphav[17]+0.4242640687119285*alphav[13]+0.5477225575051661*alphav[11]-1.10227038425243*alphav[10]+0.3162277660168379*alphav[7]+0.8215838362577489*alphav[6]-0.6363961030678926*alphav[5]-0.8215838362577489*alphav[4]+0.4743416490252568*alphav[3]+0.6123724356957944*alphav[2]-0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*((-0.9185586535436913*alphav[17])-0.5303300858899104*alphav[13]-0.6846531968814574*alphav[11]-0.3952847075210473*alphav[7]+0.8215838362577489*alphav[6]+0.4743416490252568*alphav[3]+0.6123724356957944*alphav[2]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.7348469228349533*alphav[17]+0.4242640687119285*alphav[13]+0.5477225575051661*alphav[11]+1.10227038425243*alphav[10]+0.3162277660168379*alphav[7]+0.8215838362577489*alphav[6]+0.6363961030678926*alphav[5]+0.8215838362577489*alphav[4]+0.4743416490252568*alphav[3]+0.6123724356957944*alphav[2]+0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  out[1] += 0.6123724356957944*(alphax[12]*f[12]+alphax[11]*f[11]+alphax[8]*f[8]+alphax[7]*f[7]+alphax[4]*f[4]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphav[17]*f[17]+alphav[13]*f[13]+alphav[11]*f[11]+alphav[10]*f[10]+alphav[7]*f[7]+alphav[6]*f[6]+alphav[5]*f[5]+alphav[4]*f[4]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[4] += 0.07071067811865474*(7.745966692414834*(alphav[10]*f[17]+f[10]*alphav[17]+alphav[5]*f[13]+f[5]*alphav[13]+alphax[4]*f[12]+f[4]*alphax[12])+(8.660254037844387*alphax[7]+7.745966692414834*alphav[4])*f[11]+8.660254037844387*f[7]*alphax[11]+7.745966692414834*f[4]*alphav[11]+8.660254037844386*(alphav[6]*f[10]+f[6]*alphav[10])+7.745966692414834*(alphax[2]*f[8]+f[2]*alphax[8]+alphav[1]*f[7]+f[1]*alphav[7])+8.660254037844386*(alphav[3]*f[5]+f[3]*alphav[5]+(alphav[2]+alphax[1])*f[4]+f[1]*alphax[4]+f[2]*(alphav[4]+alphax[0])+f[0]*alphax[2]+alphav[0]*f[1]+f[0]*alphav[1])); 
  out[5] += 0.07071067811865474*(8.660254037844387*(alphax[12]*f[18]+alphax[11]*f[17]+alphax[8]*f[14]+alphax[7]*f[13])+8.660254037844386*(alphax[4]*f[10]+alphax[2]*f[6]+alphax[1]*f[5]+alphax[0]*f[3])); 
  out[6] += 0.07071067811865474*(7.745966692414834*alphav[10]*f[19]+8.660254037844387*(alphav[11]*f[17]+f[11]*alphav[17])+7.745966692414834*(alphav[6]*f[16]+alphav[5]*f[15])+8.660254037844387*(alphav[7]*f[13]+f[7]*alphav[13])+8.660254037844386*(alphav[4]*f[10]+f[4]*alphav[10])+7.745966692414834*alphav[3]*f[9]+8.660254037844386*(alphav[2]*f[6]+f[2]*alphav[6]+alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[3]+f[0]*alphav[3])); 
  out[7] += 0.07071067811865474*(19.36491673103708*(alphax[8]*f[12]+f[8]*alphax[12])+17.32050807568877*(alphax[4]*f[11]+f[4]*alphax[11])+17.32050807568877*(alphax[1]*f[7]+f[1]*alphax[7])+19.36491673103709*(alphax[2]*f[4]+f[2]*alphax[4]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[8] += 0.07071067811865474*(17.32050807568877*alphav[10]*f[18]+19.36491673103708*(alphav[13]*f[17]+f[13]*alphav[17])+17.32050807568877*(alphav[6]*f[14]+alphav[4]*f[12])+19.36491673103708*(alphav[7]*f[11]+f[7]*alphav[11])+19.36491673103709*(alphav[5]*f[10]+f[5]*alphav[10])+17.32050807568877*alphav[2]*f[8]+19.36491673103709*(alphav[3]*f[6]+f[3]*alphav[6]+alphav[1]*f[4]+f[1]*alphav[4]+alphav[0]*f[2]+f[0]*alphav[2])); 
  out[10] += 0.07071067811865474*(6.928203230275509*alphav[17]*f[19]+7.745966692414834*(alphav[6]*f[19]+alphax[4]*f[18])+8.660254037844386*alphax[7]*f[17]+7.745966692414834*(alphav[4]*f[17]+f[4]*alphav[17]+alphav[10]*f[16])+6.928203230275509*alphav[13]*f[15]+7.745966692414834*(alphav[3]*f[15]+alphax[2]*f[14])+8.660254037844386*alphax[11]*f[13]+7.745966692414834*(alphav[1]*f[13]+f[1]*alphav[13]+f[10]*alphax[12]+alphav[10]*f[11]+f[10]*alphav[11])+8.660254037844386*((alphav[2]+alphax[1])*f[10]+f[2]*alphav[10])+7.745966692414834*(alphav[5]*f[9]+f[6]*alphax[8]+alphav[5]*f[7]+f[5]*alphav[7])+8.660254037844386*((alphav[4]+alphax[0])*f[6]+f[4]*alphav[6]+(alphax[4]+alphav[0])*f[5]+f[0]*alphav[5]+(alphax[2]+alphav[1])*f[3]+f[1]*alphav[3])); 
  out[11] += 0.01010152544552211*(38.72983346207417*alphav[17]*f[17]+60.62177826491071*(alphav[6]*f[17]+f[6]*alphav[17])+38.72983346207417*alphav[13]*f[13]+60.6217782649107*(alphav[3]*f[13]+f[3]*alphav[13])+(108.4435336938077*alphax[11]+121.2435565298214*alphax[2])*f[12]+(108.4435336938077*f[11]+121.2435565298214*f[2])*alphax[12]+(38.72983346207417*alphav[11]+60.6217782649107*alphav[2])*f[11]+121.2435565298214*(alphax[1]*f[11]+f[1]*alphax[11])+60.6217782649107*f[2]*alphav[11]+54.22176684690384*alphav[10]*f[10]+121.2435565298214*(alphax[4]*f[8]+f[4]*alphax[8])+(38.72983346207417*alphav[7]+121.2435565298214*alphax[4]+60.62177826491071*alphav[0])*f[7]+121.2435565298214*f[4]*alphax[7]+60.62177826491071*f[0]*alphav[7]+54.22176684690384*(alphav[5]*f[5]+alphav[4]*f[4])+135.5544171172596*(alphax[0]*f[4]+f[0]*alphax[4]+alphax[1]*f[2])+f[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphav[1])); 
  out[12] += 0.01010152544552211*(108.4435336938077*alphav[17]*f[18]+121.2435565298214*(alphav[6]*f[18]+alphav[5]*f[17]+f[5]*alphav[17])+121.2435565298214*(alphav[10]*(f[14]+f[13])+f[10]*alphav[13])+(38.72983346207417*alphax[12]+108.4435336938077*alphav[11]+121.2435565298214*alphav[2])*f[12]+60.6217782649107*(alphax[1]*f[12]+f[1]*alphax[12])+54.22176684690384*alphax[11]*f[11]+121.2435565298214*(alphav[1]*f[11]+f[1]*alphav[11])+135.5544171172596*(alphav[3]*f[10]+f[3]*alphav[10])+(38.72983346207417*alphax[8]+121.2435565298214*alphav[4])*f[8]+60.62177826491071*(alphax[0]*f[8]+f[0]*alphax[8])+121.2435565298214*(alphav[4]*f[7]+f[4]*alphav[7])+135.5544171172596*(alphav[5]*f[6]+f[5]*alphav[6])+54.22176684690384*alphax[4]*f[4]+135.5544171172596*(alphav[0]*f[4]+f[0]*alphav[4])+54.22176684690384*alphax[2]*f[2]+135.5544171172596*(alphav[1]*f[2]+f[1]*alphav[2])); 
  out[13] += 0.07071067811865474*(19.36491673103708*alphax[8]*f[18]+17.32050807568877*alphax[4]*f[17]+19.36491673103708*alphax[12]*f[14]+17.32050807568877*alphax[1]*f[13]+f[10]*(17.32050807568877*alphax[11]+19.36491673103708*alphax[2])+17.32050807568877*f[5]*alphax[7]+19.36491673103708*(alphax[4]*f[6]+alphax[0]*f[5]+alphax[1]*f[3])); 
  out[14] += 0.07071067811865474*(17.32050807568877*(alphav[5]*f[19]+alphav[4]*f[18])+19.36491673103708*(alphav[7]*f[17]+f[7]*alphav[17])+17.32050807568877*(alphav[3]*f[16]+alphav[10]*f[15]+alphav[2]*f[14])+19.36491673103708*(alphav[11]*f[13]+f[11]*alphav[13])+17.32050807568877*alphav[10]*f[12]+19.36491673103708*(alphav[1]*f[10]+f[1]*alphav[10])+17.32050807568877*alphav[6]*(f[9]+f[8])+19.36491673103708*(alphav[0]*f[6]+f[0]*alphav[6]+alphav[4]*f[5]+f[4]*alphav[5]+alphav[2]*f[3]+f[2]*alphav[3])); 
  out[15] += 0.07071067811865474*(8.660254037844387*alphax[4]*f[19]+8.660254037844386*(alphax[2]*f[16]+alphax[1]*f[15])+8.660254037844387*alphax[0]*f[9]); 
  out[16] += 0.07071067811865474*(8.660254037844387*alphav[4]*f[19]+7.745966692414834*alphav[17]*f[17]+8.660254037844386*(alphav[2]*f[16]+alphav[1]*f[15])+7.745966692414834*(alphav[13]*f[13]+alphav[10]*f[10])+8.660254037844387*alphav[0]*f[9]+7.745966692414834*(alphav[6]*f[6]+alphav[5]*f[5]+alphav[3]*f[3])); 
  out[17] += 0.002020305089104421*(242.4871130596428*alphav[10]*f[19]+(542.2176684690384*alphax[11]+606.217782649107*alphax[2])*f[18]+(542.2176684690384*alphax[12]+193.6491673103708*alphav[11]+303.1088913245535*alphav[2]+606.217782649107*alphax[1])*f[17]+(271.1088342345192*f[16]+193.6491673103708*f[11]+303.1088913245535*f[2])*alphav[17]+242.4871130596428*alphav[5]*f[15]+606.2177826491072*alphax[4]*f[14]+(193.6491673103708*alphav[7]+606.2177826491072*alphax[4]+303.1088913245536*alphav[0])*f[13]+(271.1088342345192*f[9]+193.6491673103708*f[7]+303.1088913245536*f[0])*alphav[13]+606.2177826491072*f[6]*alphax[12]+303.1088913245536*alphav[6]*f[11]+606.2177826491072*f[5]*alphax[11]+303.1088913245536*f[6]*alphav[11]+(606.217782649107*(alphax[8]+alphax[7])+271.1088342345192*alphav[4]+677.7720855862981*alphax[0])*f[10]+271.1088342345192*f[4]*alphav[10]+303.1088913245535*(alphav[3]*f[7]+f[3]*alphav[7])+677.7720855862981*(alphax[1]*f[6]+alphax[2]*f[5])+271.1088342345192*(alphav[1]*f[5]+f[1]*alphav[5])+677.7720855862981*f[3]*alphax[4]); 
  out[18] += 0.01010152544552211*((108.4435336938077*alphav[13]+121.2435565298214*alphav[3])*f[19]+(38.72983346207417*alphax[12]+108.4435336938077*alphav[11]+121.2435565298214*alphav[2]+60.6217782649107*alphax[1])*f[18]+(54.22176684690384*alphax[11]+121.2435565298214*alphav[1])*f[17]+(108.4435336938077*(f[15]+f[12])+121.2435565298214*f[1])*alphav[17]+121.2435565298214*(alphav[5]*f[16]+alphav[6]*f[15])+(38.72983346207417*alphax[8]+121.2435565298214*alphav[4]+60.62177826491071*alphax[0])*f[14]+121.2435565298214*(alphav[4]*f[13]+f[4]*alphav[13]+alphav[6]*f[12])+60.62177826491071*f[5]*alphax[12]+121.2435565298214*(alphav[5]*f[11]+f[5]*alphav[11])+(121.2435565298214*alphav[7]+54.22176684690384*alphax[4]+135.5544171172596*alphav[0])*f[10]+(121.2435565298214*(f[9]+f[8]+f[7])+135.5544171172596*f[0])*alphav[10]+60.6217782649107*f[3]*alphax[8]+54.22176684690384*alphax[2]*f[6]+135.5544171172596*(alphav[1]*f[6]+f[1]*alphav[6]+alphav[2]*f[5]+f[2]*alphav[5]+alphav[3]*f[4]+f[3]*alphav[4])); 
  out[19] += 0.01414213562373095*((38.72983346207417*(alphax[12]+alphav[11])+43.30127018922193*(alphav[2]+alphax[1]))*f[19]+34.64101615137754*(alphav[10]*f[17]+f[10]*alphav[17])+(38.72983346207417*alphax[8]+43.30127018922195*(alphav[4]+alphax[0]))*f[16]+(38.72983346207417*alphav[7]+43.30127018922195*(alphax[4]+alphav[0]))*f[15]+34.64101615137755*(alphav[5]*f[13]+f[5]*alphav[13])+38.72983346207418*(alphav[6]*f[10]+f[6]*alphav[10])+43.30127018922193*(alphax[2]+alphav[1])*f[9]+38.72983346207418*(alphav[3]*f[5]+f[3]*alphav[5])); 
  return cflFreq; 
} 
