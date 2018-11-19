#include <GyrokineticModDecl.h> 
double EmGyrokineticVol1x1vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[8]; 
  alphax[0] = (1.224744871391589*BmagInv[0]*geoY[0]*Apar[1]*dfac_x2+1.414213562373095*Gradpar[0]*dfac_x)*wv; 
  alphax[1] = 2.738612787525831*BmagInv[0]*geoY[0]*Apar[2]*dfac_x2*wv; 
  alphax[2] = (0.7071067811865475*BmagInv[0]*geoY[0]*Apar[1]*dfac_x2+0.8164965809277261*Gradpar[0]*dfac_x)/dfac_v; 
  alphax[3] = (1.58113883008419*BmagInv[0]*geoY[0]*Apar[2]*dfac_x2)/dfac_v; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.25*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(1.732050807568877*alphax[1]+alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(1.161895003862225*alphax[3]-0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.25*(0.5*alphax[0]-0.8660254037844386*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.25*((-1.161895003862225*alphax[3])+0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*((-1.161895003862225*alphax[3])-0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.25*(0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.25*(1.161895003862225*alphax[3]+0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphav[8]; 
  alphav[0] = (dfac_v*(Phi[1]*((-1.5*BmagInv[0]*geoY[0]*Apar[1]*dfac_x2)-1.732050807568877*Gradpar[0]*dfac_x)-7.5*BmagInv[0]*geoY[0]*Apar[2]*Phi[2]*dfac_x2)*q_)/m_; 
  alphav[1] = (2.23606797749979*dfac_v*((-1.5*BmagInv[0]*geoY[0]*(Apar[1]*Phi[2]+Phi[1]*Apar[2])*dfac_x2)-1.732050807568877*Gradpar[0]*Phi[2]*dfac_x)*q_)/m_; 
  alphav[4] = -(6.708203932499369*BmagInv[0]*geoY[0]*Apar[2]*Phi[2]*dfac_v*dfac_x2*q_)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.4472135954999579*alphav[4]-0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.25*(0.5*alphav[0]-0.5590169943749475*alphav[4]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.25*(0.4472135954999579*alphav[4]+0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*(0.4472135954999579*alphav[4]-0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.25*(0.5*alphav[0]-0.5590169943749475*alphav[4]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.25*(0.4472135954999579*alphav[4]+0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  out[1] += 0.8660254037844386*(alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphav[4]*f[4]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[3] += 0.1*(7.745966692414834*(alphax[3]*f[7]+alphax[2]*f[5]+alphav[1]*f[4]+f[1]*alphav[4])+8.660254037844386*(alphax[1]*f[3]+f[1]*alphax[3]+alphax[0]*f[2]+f[0]*alphax[2]+alphav[0]*f[1]+f[0]*alphav[1])); 
  out[4] += 0.1*(17.32050807568877*alphax[3]*f[6]+17.32050807568877*alphax[1]*f[4]+19.36491673103709*(alphax[2]*f[3]+f[2]*alphax[3]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[5] += 1.936491673103709*(alphav[4]*f[6]+alphav[1]*f[3]+alphav[0]*f[2]); 
  out[6] += 0.01428571428571429*(121.2435565298214*(alphax[2]*f[7]+alphax[1]*f[6])+121.2435565298214*alphax[3]*f[5]+(38.72983346207417*alphav[4]+121.2435565298214*alphax[3])*f[4]+60.62177826491071*(alphav[0]*f[4]+f[0]*alphav[4])+135.5544171172596*(alphax[0]*f[3]+f[0]*alphax[3]+alphax[1]*f[2])+f[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphav[1])); 
  out[7] += 0.1*(8.660254037844386*alphax[1]*f[7]+17.32050807568877*alphav[1]*f[6]+8.660254037844387*alphax[0]*f[5]+f[3]*(17.32050807568877*alphav[4]+7.745966692414834*alphax[3]+19.36491673103708*alphav[0])+(7.745966692414834*alphax[2]+19.36491673103708*alphav[1])*f[2]); 
  return cflFreq; 
} 
double EmGyrokineticStep2Vol1x1vSerP2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  double dvInv = 1.0/dxv[1]; 
  double dfac_v = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[2]*f[4]+dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  out[3] += -(0.2449489742783178*(4.47213595499958*dApardt[1]*f[4]+4.47213595499958*f[1]*dApardt[2]+5.0*dApardt[0]*f[1]+5.0*f[0]*dApardt[1])*dfac_v*q_)/m_; 
  out[5] += -(0.7071067811865475*(3.872983346207417*dApardt[2]*f[6]+3.872983346207417*dApardt[1]*f[3]+3.872983346207417*dApardt[0]*f[2])*dfac_v*q_)/m_; 
  out[6] += -(0.07824607964359516*(10.0*dApardt[2]*f[4]+15.65247584249853*dApardt[0]*f[4]+15.65247584249853*f[0]*dApardt[2]+14.0*dApardt[1]*f[1])*dfac_v*q_)/m_; 
  out[7] += -(0.1414213562373095*(17.32050807568877*dApardt[1]*f[6]+17.32050807568877*dApardt[2]*f[3]+19.36491673103708*dApardt[0]*f[3]+19.36491673103708*dApardt[1]*f[2])*dfac_v*q_)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -(0.3535533905932737*dApardt[0]*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = -(0.3535533905932737*dApardt[0]*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = -(0.25*(0.6324555320336759*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.25*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.25*(0.6324555320336759*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = -(0.25*(0.6324555320336759*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.25*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.25*(0.6324555320336759*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
#endif 
  return cflFreq; 
} 
double EmGyrokineticVol1x1vSerP2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[8]; 
  alphax[0] = -1.0*((1.732050807568877*(0.7071067811865475*(Bmag[1]*BmagInv[2]*geoY[2]+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))+1.414213562373095*Bmag[2]*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2]))*dfac_x2*m_*wv2)/q_+(1.732050807568877*(Bmag[2]*(Apar[0]*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2])+(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Apar[2]+1.118033988749895*(Apar[0]*BmagInv[0]*geoY[1]+geoY[0]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1])))+((-1.0*((-1.756910553749835*Apar[1]*Bmag[2])+Bmag[1]*((-0.31943828249997*Apar[2])-0.5*Apar[0])+0.7071067811865475*Apar[1])*BmagInv[2])+(1.756910553749835*BmagInv[1]*Apar[2]+BmagInv[0]*Apar[1])*Bmag[2]+Bmag[1]*(0.5*BmagInv[0]*Apar[2]+0.4472135954999579*Apar[1]*BmagInv[1])-1.414213562373095*BmagInv[1]*Apar[2])*geoY[2]+0.4472135954999579*geoY[1]*(Apar[1]*Bmag[1]*BmagInv[2]+BmagInv[1]*(4.5*Apar[1]*Bmag[2]+Bmag[1]*Apar[2]))+((1.756910553749835*geoY[1]*Apar[2]+geoY[0]*Apar[1])*Bmag[2]-1.0*(1.414213562373095*geoY[1]-0.5*geoY[0]*Bmag[1])*Apar[2])*BmagInv[2])*dfac_x2-1.0*(1.732050807568877*(1.58113883008419*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Apar[2]+(Apar[1]*(0.7071067811865475*BmagInv[1]-0.5*BmagInv[0]*Bmag[1])-0.5*Apar[0]*Bmag[1]*BmagInv[1])*geoY[1]+geoY[0]*(BmagInv[0]*(0.7071067811865475*Apar[1]-0.5*Apar[0]*Bmag[1])-0.5*Apar[1]*Bmag[1]*BmagInv[1]))*dfac_x2+1.414213562373095*Gradpar[0]*dfac_x))*wv+(0.5773502691896258*(0.7071067811865475*(Bmag[1]*BmagInv[2]*geoY[2]+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))+1.414213562373095*Bmag[2]*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2]))*dfac_x2*m_)/(dfac_v2*q_)); 
  alphax[1] = (1.732050807568877*((-1.0*(1.414213562373095*((BmagInv[0]*Bmag[2]+0.4472135954999579*Bmag[1]*BmagInv[1])*geoY[2]+(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2])+0.7071067811865475*((4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])*Bmag[2]+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))))-2.484646732989441*Bmag[2]*BmagInv[2]*geoY[2])*dfac_x2*m_*wv2)/q_+(1.414213562373095*(0.7745966692414833*Apar[1]*geoY[1]*BmagInv[2]*dfac_x2+Gradpar[1]*dfac_x)+1.732050807568877*(Apar[2]*((-0.7857142857142857*(Bmag[1]*BmagInv[1]*geoY[2]+2.23606797749979*geoY[0]*Bmag[2]*BmagInv[2]))+geoY[1]*((-2.571428571428571*BmagInv[1]*Bmag[2])-0.4472135954999579*BmagInv[0]*Bmag[1])+geoY[0]*(1.58113883008419*BmagInv[0]-0.4472135954999579*Bmag[1]*BmagInv[1]))+((2.23606797749979*(1.111167799007432*Apar[2]-0.7857142857142857*Apar[0]*Bmag[2])-2.857142857142857*Apar[2]*Bmag[2]-0.7857142857142857*Apar[1]*Bmag[1])*BmagInv[2]+BmagInv[0]*(1.414213562373095*Apar[2]-1.0*Apar[0]*Bmag[2])+((-1.756910553749835*BmagInv[0]*Apar[2])-2.571428571428571*Apar[1]*BmagInv[1])*Bmag[2]+0.4472135954999579*(Apar[1]*(1.414213562373095*BmagInv[1]-1.0*BmagInv[0]*Bmag[1])-1.0*Apar[0]*Bmag[1]*BmagInv[1]))*geoY[2]+geoY[1]*((-0.7857142857142857*Bmag[1]*Apar[2]*BmagInv[2])+9.0*BmagInv[1]*(0.4472135954999579*(0.7071067811865475*Apar[2]-0.5*Apar[0]*Bmag[2])-0.1*Apar[1]*Bmag[1])+BmagInv[0]*(0.7071067811865475*Apar[1]-0.5*Apar[0]*Bmag[1]))+(geoY[0]*(1.414213562373095*Apar[2]-1.0*Apar[0]*Bmag[2])-2.571428571428571*Apar[1]*geoY[1]*Bmag[2]-0.4472135954999579*Bmag[1]*(Apar[0]*geoY[1]+geoY[0]*Apar[1]))*BmagInv[2]+((-1.0*BmagInv[0]*geoY[0]*Apar[2])-0.5*(4.024922359499621*BmagInv[0]*Apar[1]*geoY[1]+geoY[0]*(4.024922359499621*Apar[1]*BmagInv[1]+2.23606797749979*Apar[0]*BmagInv[0])))*Bmag[2]+geoY[0]*(Apar[1]*(0.7071067811865475*BmagInv[1]-0.5*BmagInv[0]*Bmag[1])-0.5*Apar[0]*Bmag[1]*BmagInv[1]))*dfac_x2)*wv+(((-1.0*(0.8164965809277261*((BmagInv[0]*Bmag[2]+0.4472135954999579*Bmag[1]*BmagInv[1])*geoY[2]+(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2])+0.7071067811865475*(0.5773502691896258*(2.23606797749979*BmagInv[0]*geoY[0]*Bmag[2]+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+2.32379000772445*BmagInv[1]*geoY[1]*Bmag[2])))-1.434511460132578*Bmag[2]*BmagInv[2]*geoY[2])*dfac_x2*m_)/(dfac_v2*q_); 
  alphax[2] = -(1.0*(dfac_x2*(BmagInv[2]*(Bmag[1]*geoY[2]*((1.414213562373095*m_*wv)/q_+0.5*Apar[0])+(geoY[1]*(1.756910553749835*Apar[2]+Apar[0])+geoY[0]*Apar[1])*Bmag[2]+Bmag[1]*(0.5*geoY[0]*Apar[2]+0.4472135954999579*Apar[1]*geoY[1])-1.414213562373095*geoY[1]*Apar[2])+(1.414213562373095*(Bmag[2]*(2.0*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2])+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))*m_*wv)/q_+((-1.0*(2.23606797749979*((-0.7857142857142857*Apar[1]*Bmag[2])-0.1428571428571428*Bmag[1]*Apar[2])+0.7071067811865475*Apar[1])*BmagInv[2])+(BmagInv[1]*(1.756910553749835*Apar[2]+Apar[0])+BmagInv[0]*Apar[1])*Bmag[2]+Bmag[1]*(0.5*BmagInv[0]*Apar[2]+0.4472135954999579*Apar[1]*BmagInv[1])-1.414213562373095*BmagInv[1]*Apar[2])*geoY[2])-1.0*((geoY[0]*(BmagInv[1]*((-1.0*Apar[2]*Bmag[2])-0.5*Apar[1]*Bmag[1])+BmagInv[0]*(0.7071067811865475*Apar[1]-0.5*Apar[0]*Bmag[1]))+geoY[1]*(BmagInv[1]*(0.4472135954999579*((-4.5*Apar[1]*Bmag[2])-1.0*Bmag[1]*Apar[2])-0.5*Apar[0]*Bmag[1])+Apar[1]*(0.7071067811865475*BmagInv[1]-0.5*BmagInv[0]*Bmag[1]))+2.23606797749979*(0.7071067811865475*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Apar[2]-0.5*(Apar[0]*BmagInv[0]*geoY[1]+geoY[0]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))*Bmag[2])-1.0*BmagInv[0]*geoY[1]*Apar[2]*Bmag[2])*dfac_x2+0.8164965809277261*Gradpar[0]*dfac_x)))/dfac_v; 
  alphax[3] = (dfac_x2*((1.414213562373095*((-1.0*(2.0*(0.4472135954999579*Bmag[1]*BmagInv[1]*geoY[2]+(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2])+Bmag[2]*(2.0*BmagInv[0]*geoY[2]+4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))-3.51382110749967*Bmag[2]*BmagInv[2]*geoY[2])*m_*wv)/q_+Apar[2]*((-0.7857142857142857*Bmag[1]*BmagInv[1]*geoY[2])+geoY[1]*((-0.7857142857142857*Bmag[1]*BmagInv[2])-2.571428571428572*BmagInv[1]*Bmag[2]-0.4472135954999579*BmagInv[0]*Bmag[1])+geoY[0]*(1.58113883008419*BmagInv[0]-0.4472135954999579*Bmag[1]*BmagInv[1]))+((2.23606797749979*(1.111167799007432*Apar[2]-0.7857142857142857*Apar[0]*Bmag[2])-2.857142857142857*Apar[2]*Bmag[2]-0.7857142857142857*Apar[1]*Bmag[1])*BmagInv[2]+BmagInv[0]*(1.414213562373095*Apar[2]-1.0*Apar[0]*Bmag[2])+((-1.756910553749835*BmagInv[0]*Apar[2])-2.571428571428572*Apar[1]*BmagInv[1])*Bmag[2]+0.4472135954999579*(Apar[1]*(1.414213562373095*BmagInv[1]-1.0*BmagInv[0]*Bmag[1])-1.0*Apar[0]*Bmag[1]*BmagInv[1]))*geoY[2]+geoY[0]*((-1.756910553749835*Apar[2]*Bmag[2]*BmagInv[2])-0.5*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(0.7071067811865475*BmagInv[1]-0.5*BmagInv[0]*Bmag[1]))+Bmag[2]*((-2.571428571428572*Apar[1]*geoY[1]*BmagInv[2])+BmagInv[0]*((-1.0*geoY[0]*Apar[2])-2.012461179749811*Apar[1]*geoY[1])+geoY[0]*((-2.012461179749811*Apar[1]*BmagInv[1])-1.118033988749895*Apar[0]*BmagInv[0]))+geoY[1]*((-0.4472135954999579*Apar[0]*Bmag[1]*BmagInv[2])+BmagInv[1]*(0.4472135954999579*(6.363961030678928*Apar[2]-4.5*Apar[0]*Bmag[2])-0.9*Apar[1]*Bmag[1])+BmagInv[0]*(0.7071067811865475*Apar[1]-0.5*Apar[0]*Bmag[1]))+(geoY[0]*(1.414213562373095*Apar[2]-1.0*Apar[0]*Bmag[2])+0.4472135954999579*Apar[1]*(1.414213562373095*geoY[1]-1.0*geoY[0]*Bmag[1]))*BmagInv[2])+0.8164965809277261*Gradpar[1]*dfac_x)/dfac_v; 
  alphax[4] = (1.732050807568877*(2.23606797749979*((-1.111167799007432*Bmag[2]*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2]))-0.2020305089104422*Bmag[1]*BmagInv[2]*geoY[2])-1.0*(0.7071067811865475*Bmag[1]*(BmagInv[0]*geoY[2]+geoY[0]*BmagInv[2])+1.414213562373095*(BmagInv[1]*(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])+BmagInv[0]*geoY[1]*Bmag[2])))*dfac_x2*m_*wv2)/q_+(1.414213562373095*(0.7745966692414833*Apar[1]*BmagInv[1]*geoY[1]*dfac_x2+Gradpar[2]*dfac_x)+1.732050807568877*(Bmag[2]*((-1.756910553749835*(Apar[0]*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2])+(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Apar[2]))-1.0*(Apar[0]*BmagInv[0]*geoY[1]+geoY[0]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))-2.571428571428571*Apar[1]*BmagInv[1]*geoY[1])+Bmag[1]*(Apar[2]*((-1.071428571428571*BmagInv[2]*geoY[2])-0.5*BmagInv[0]*geoY[0])-0.4472135954999579*(Apar[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])+Apar[0]*BmagInv[1]*geoY[1]))+BmagInv[2]*((-0.31943828249997*Apar[0]*Bmag[1]*geoY[2])+((-2.857142857142857*geoY[1]*Apar[2])-1.756910553749835*geoY[0]*Apar[1])*Bmag[2]+Bmag[1]*((-0.31943828249997*geoY[0]*Apar[2])-0.7857142857142857*Apar[1]*geoY[1])+2.484646732989441*geoY[1]*Apar[2]+geoY[0]*(0.7071067811865475*Apar[1]-0.5*Apar[0]*Bmag[1]))+(Apar[1]*(0.4517539514526256-2.857142857142857*Bmag[2])*BmagInv[2]+((-2.857142857142857*BmagInv[1]*Apar[2])-1.756910553749835*BmagInv[0]*Apar[1])*Bmag[2]+Bmag[1]*((-0.31943828249997*BmagInv[0]*Apar[2])-0.7857142857142857*Apar[1]*BmagInv[1])+2.484646732989441*BmagInv[1]*Apar[2]+BmagInv[0]*(0.7071067811865475*Apar[1]-0.5*Apar[0]*Bmag[1]))*geoY[2]+(1.414213562373095*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])-0.7857142857142857*Bmag[1]*BmagInv[1]*geoY[1])*Apar[2])*dfac_x2)*wv+(0.5773502691896258*(2.23606797749979*((-1.111167799007432*Bmag[2]*(BmagInv[1]*geoY[2]+geoY[1]*BmagInv[2]))-0.2020305089104422*Bmag[1]*BmagInv[2]*geoY[2])-1.0*(0.7071067811865475*Bmag[1]*(BmagInv[0]*geoY[2]+geoY[0]*BmagInv[2])+1.414213562373095*(BmagInv[1]*(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])+BmagInv[0]*geoY[1]*Bmag[2])))*dfac_x2*m_)/(dfac_v2*q_); 
  alphax[5] = -(0.8164965809277261*(Bmag[2]*(0.8944271909999159*BmagInv[1]*geoY[2]+geoY[1]*(0.8944271909999159*BmagInv[2]+BmagInv[0])+geoY[0]*BmagInv[1])+0.4472135954999579*Bmag[1]*(BmagInv[2]*geoY[2]+BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))*dfac_x2*m_)/(dfac_v2*q_); 
  alphax[6] = (3.872983346207417*(dfac_x2*(Bmag[2]*(0.5773502691896258*((-(1.264911064067352*BmagInv[0]*geoY[1]*m_*wv)/q_)+Apar[2]*((-1.27775312999988*BmagInv[1]*geoY[2])-0.7857142857142857*BmagInv[0]*geoY[1])-0.7857142857142857*BmagInv[0]*Apar[1]*geoY[2]+((-1.27775312999988*geoY[1]*Apar[2])-0.7857142857142857*(Apar[0]*geoY[1]+geoY[0]*Apar[1]))*BmagInv[2])+BmagInv[1]*((-0.4536323543632774*(Apar[0]*geoY[2]+geoY[0]*Apar[2]))-0.6639400022069856*Apar[1]*geoY[1]))+0.5773502691896258*((1.414213562373095*(2.0*BmagInv[1]*((-0.4472135954999579*geoY[0]*Bmag[2])-0.2*Bmag[1]*geoY[1])-0.4472135954999579*geoY[0]*Bmag[1]*BmagInv[2])*m_*wv)/q_+0.4472135954999579*((-0.7857142857142857*Apar[1]*Bmag[1]*BmagInv[1]*geoY[2])+geoY[0]*(0.7071067811865475*Apar[1]-0.5*Apar[0]*Bmag[1])*BmagInv[2]-1.0*(Apar[0]*BmagInv[0]*geoY[1]+geoY[0]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))*Bmag[2]+((1.414213562373095*BmagInv[0]-0.7857142857142857*Bmag[1]*BmagInv[1])*geoY[1]+geoY[0]*(1.414213562373095*BmagInv[1]-0.5*BmagInv[0]*Bmag[1]))*Apar[2])+Apar[2]*(1.111167799007432*geoY[1]*BmagInv[2]-0.1428571428571428*BmagInv[0]*Bmag[1]*geoY[2])+Bmag[1]*((-0.223606797749979*Apar[0]*BmagInv[0]*geoY[2])-0.1428571428571428*geoY[0]*Apar[2]*BmagInv[2]-0.2*Apar[0]*BmagInv[1]*geoY[1])+Apar[1]*(0.3162277660168379*BmagInv[0]*geoY[2]+Bmag[1]*((-0.3513821107499669*geoY[1]*BmagInv[2])-0.2*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+0.2828427124746191*BmagInv[1]*geoY[1])+((0.2020305089104422*Apar[1]-0.1428571428571428*Apar[0]*Bmag[1])*BmagInv[2]+1.111167799007432*BmagInv[1]*Apar[2])*geoY[2])+BmagInv[2]*(2.23606797749979*((-0.3299144395369291*Apar[1]*Bmag[2])-0.1237179148263484*Bmag[1]*Apar[2])*geoY[2]-(1.28306605574357*geoY[1]*Bmag[2]*m_*wv)/q_)+(0.8164965809277261*((-0.2857142857142857*Bmag[1]*BmagInv[2])-1.571428571428571*BmagInv[1]*Bmag[2]-0.4472135954999579*BmagInv[0]*Bmag[1])*geoY[2]*m_*wv)/q_)+0.210818510677892*Gradpar[2]*dfac_x))/dfac_v; 
  alphax[7] = -(0.3651483716701107*(Bmag[2]*((3.51382110749967*BmagInv[2]+2.0*BmagInv[0])*geoY[2]+4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])+2.0*(0.4472135954999579*Bmag[1]*BmagInv[1]*geoY[2]+(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2])+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*dfac_x2*m_)/(dfac_v2*q_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*(2.23606797749979*alphax[4]-1.732050807568877*alphax[1]+alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(2.23606797749979*alphax[4]+1.732050807568877*alphax[1]+alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*((-0.7745966692414833*alphax[7])-1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]+1.161895003862225*alphax[3]-0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.25*(0.9682458365518543*alphax[7]-0.5590169943749475*alphax[5]+1.118033988749895*alphax[4]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.25*((-0.7745966692414833*alphax[7])+1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]-1.161895003862225*alphax[3]+0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*(0.7745966692414833*alphax[7]-1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]-1.161895003862225*alphax[3]-0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.25*((-0.9682458365518543*alphax[7])-0.5590169943749475*alphax[5]+1.118033988749895*alphax[4]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.25*(0.7745966692414833*alphax[7]+1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]+1.161895003862225*alphax[3]+0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphav[8]; 
  alphav[0] = dfac_v*((Phi[2]*(3.0*(2.23606797749979*BmagInv[0]*Bmag[2]*geoY[2]+(2.23606797749979*geoY[0]*Bmag[2]+Bmag[1]*geoY[1])*BmagInv[2])+(13.5*BmagInv[1]*geoY[1]+7.5*BmagInv[0]*geoY[0])*Bmag[2]+3.354101966249685*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+((11.78571428571429*Bmag[2]*BmagInv[2]+3.0*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*(1.5*Bmag[1]*BmagInv[2]+3.0*BmagInv[1]*Bmag[2]))*geoY[2]+Phi[1]*(3.0*geoY[1]*Bmag[2]*BmagInv[2]+1.5*(2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))))*dfac_x2*wv+(((Phi[2]*((-3.0*Apar[1]*BmagInv[1]*geoY[2])+2.23606797749979*(geoY[1]*(5.454823740581938*Apar[1]*Bmag[2]*BmagInv[2]+BmagInv[0]*(1.060660171779821*Apar[0]*Bmag[1]-1.5*Apar[1]))+geoY[0]*(1.060660171779821*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(1.060660171779821*BmagInv[0]*Bmag[1]-1.5*BmagInv[1])))+Apar[2]*(geoY[1]*(3.726970099484161*Bmag[1]*BmagInv[2]+BmagInv[1]*(12.19735668922089*Bmag[2]-13.5)+2.121320343559642*BmagInv[0]*Bmag[1])+geoY[0]*(2.121320343559642*Bmag[1]*BmagInv[1]-7.5*BmagInv[0]))+3.0*((2.23606797749979*geoY[0]*(0.7071067811865475*Apar[0]*Bmag[2]-1.0*Apar[2])+0.7071067811865475*Apar[0]*Bmag[1]*geoY[1]+Apar[1]*(0.7071067811865475*geoY[0]*Bmag[1]-1.0*geoY[1]))*BmagInv[2]+0.3535533905932737*(4.47213595499958*BmagInv[0]*geoY[0]*Apar[2]+9.0*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1])*geoY[1]+geoY[0]*(9.0*Apar[1]*BmagInv[1]+5.0*Apar[0]*BmagInv[0]))*Bmag[2])+0.3535533905932737*(23.57142857142857*geoY[0]*Apar[2]*Bmag[2]*BmagInv[2]+12.07476707849886*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]))+((((13.55261854357878*Apar[2]+8.333758492555738*Apar[0])*Bmag[2]-11.78571428571429*Apar[2]+3.726970099484161*Apar[1]*Bmag[1])*BmagInv[2]+3.0*(2.23606797749979*BmagInv[0]*(0.7071067811865475*Apar[0]*Bmag[2]-1.0*Apar[2])+0.7071067811865475*Bmag[1]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))+2.23606797749979*BmagInv[1]*(5.454823740581938*Apar[1]*Bmag[2]+1.666751698511147*Bmag[1]*Apar[2])+8.333758492555738*BmagInv[0]*Apar[2]*Bmag[2])*Phi[2]+Phi[1]*((0.3535533905932737*(10.54146332249901*Apar[1]*Bmag[2]+Bmag[1]*(1.91662969499982*Apar[2]+3.0*Apar[0]))-1.5*Apar[1])*BmagInv[2]+3.0*(0.3535533905932737*(2.0*BmagInv[0]*Apar[1]*Bmag[2]+Bmag[1]*(BmagInv[0]*Apar[2]+0.8944271909999159*Apar[1]*BmagInv[1]))+BmagInv[1]*(0.7071067811865475*Apar[0]*Bmag[2]-1.0*Apar[2]))+3.726970099484161*BmagInv[1]*Apar[2]*Bmag[2]))*geoY[2]+Phi[1]*(0.3535533905932737*(3.0*(2.0*(Apar[0]*geoY[1]+geoY[0]*Apar[1])*Bmag[2]+Bmag[1]*(geoY[0]*Apar[2]+0.8944271909999159*Apar[1]*geoY[1]))+10.54146332249901*geoY[1]*Apar[2]*Bmag[2])*BmagInv[2]+geoY[0]*(1.060660171779821*BmagInv[1]*(2.0*Apar[2]*Bmag[2]+Apar[1]*Bmag[1])+BmagInv[0]*(1.060660171779821*Apar[0]*Bmag[1]-1.5*Apar[1]))+geoY[1]*(1.060660171779821*BmagInv[1]*(4.024922359499621*Apar[1]*Bmag[2]+Bmag[1]*(0.8944271909999159*Apar[2]+Apar[0]))+Apar[1]*(1.060660171779821*BmagInv[0]*Bmag[1]-1.5*BmagInv[1]))+2.23606797749979*(1.060660171779821*(Apar[0]*BmagInv[0]*geoY[1]+geoY[0]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))*Bmag[2]-1.5*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Apar[2])+2.121320343559642*BmagInv[0]*geoY[1]*Apar[2]*Bmag[2]))*dfac_x2-1.0*(3.0*Phi[1]*geoY[1]*Apar[2]*BmagInv[2]*dfac_x2+1.732050807568877*(2.23606797749979*Gradpar[1]*Phi[2]+Gradpar[0]*Phi[1])*dfac_x))*q_)/m_); 
  alphav[1] = dfac_v*((Phi[2]*(17.24966725499838*BmagInv[1]*Bmag[2]*geoY[2]+(17.24966725499838*geoY[1]*Bmag[2]+3.0*geoY[0]*Bmag[1])*BmagInv[2]+13.5*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(6.037383539249432*BmagInv[1]*geoY[1]+3.354101966249685*BmagInv[0]*geoY[0]))+(Bmag[1]*(5.270731661249505*BmagInv[2]+3.0*BmagInv[0])*Phi[2]+Phi[1]*(5.270731661249505*Bmag[2]*BmagInv[2]+3.0*(BmagInv[0]*Bmag[2]+0.4472135954999579*Bmag[1]*BmagInv[1])))*geoY[2]+Phi[1]*(3.0*(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2]+1.5*(2.23606797749979*BmagInv[0]*geoY[0]*Bmag[2]+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+6.037383539249432*BmagInv[1]*geoY[1]*Bmag[2]))*dfac_x2*wv+(((Phi[2]*(2.121320343559642*Apar[0]*BmagInv[0]*Bmag[1]*geoY[2]+Apar[2]*(geoY[1]*(20.45558902718226*Bmag[2]*BmagInv[2]+BmagInv[0]*(12.19735668922089*Bmag[2]-13.5))+geoY[0]*(2.121320343559642*BmagInv[0]*Bmag[1]-13.5*BmagInv[1]))+geoY[1]*((-17.24966725499838*Apar[2]*BmagInv[2])+BmagInv[1]*(20.45558902718226*Apar[1]*Bmag[2]+5.454823740581938*Bmag[1]*Apar[2])+0.4472135954999579*(9.545941546018389*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(9.545941546018389*BmagInv[0]*Bmag[1]-13.5*BmagInv[1])))+(12.19735668922089*(Apar[0]*geoY[1]+geoY[0]*Apar[1])*Bmag[2]+Bmag[1]*(3.726970099484161*geoY[0]*Apar[2]+5.454823740581938*Apar[1]*geoY[1])+3.0*geoY[0]*(0.7071067811865475*Apar[0]*Bmag[1]-1.0*Apar[1]))*BmagInv[2]+(12.19735668922089*geoY[0]*BmagInv[1]*Apar[2]+9.545941546018389*(Apar[0]*BmagInv[0]*geoY[1]+geoY[0]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1])))*Bmag[2]+geoY[0]*(4.269074841227311*Apar[1]*Bmag[1]*BmagInv[1]+2.23606797749979*BmagInv[0]*(1.060660171779821*Apar[0]*Bmag[1]-1.5*Apar[1])))+((0.3535533905932737*(10.54146332249901*Apar[0]*Bmag[1]*BmagInv[2]+57.85714285714285*BmagInv[1]*Apar[2]*Bmag[2])+2.23606797749979*((-2.357142857142857*Apar[1]*BmagInv[2])+5.454823740581938*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1])*Bmag[2]+(1.666751698511147*BmagInv[0]*Bmag[1]-7.714285714285714*BmagInv[1])*Apar[2])+(20.45558902718226*Apar[1]*Bmag[2]+6.060915267313265*Bmag[1]*Apar[2])*BmagInv[2]+Apar[1]*(5.454823740581938*Bmag[1]*BmagInv[1]-3.0*BmagInv[0]))*Phi[2]+Phi[1]*((2.23606797749979*(1.666751698511147*Apar[0]*Bmag[2]-2.357142857142857*Apar[2])+6.060915267313265*Apar[2]*Bmag[2]+1.666751698511147*Apar[1]*Bmag[1])*BmagInv[2]+3.0*(BmagInv[0]*(0.7071067811865475*Apar[0]*Bmag[2]-1.0*Apar[2])+0.4472135954999579*(0.7071067811865475*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(0.7071067811865475*BmagInv[0]*Bmag[1]-1.0*BmagInv[1])))+BmagInv[1]*(5.454823740581938*Apar[1]*Bmag[2]+1.666751698511147*Bmag[1]*Apar[2])+3.726970099484161*BmagInv[0]*Apar[2]*Bmag[2]))*geoY[2]+Phi[1]*(Apar[2]*(3.726970099484161*geoY[0]*Bmag[2]*BmagInv[2]+geoY[1]*(5.454823740581938*BmagInv[1]*Bmag[2]+0.9486832980505137*BmagInv[0]*Bmag[1])+geoY[0]*(0.9486832980505137*Bmag[1]*BmagInv[1]-3.354101966249685*BmagInv[0]))+geoY[1]*(1.666751698511147*Bmag[1]*Apar[2]*BmagInv[2]+BmagInv[1]*(0.4472135954999579*(9.545941546018389*Apar[0]*Bmag[2]-13.5*Apar[2])+1.909188309203678*Apar[1]*Bmag[1])+BmagInv[0]*(1.060660171779821*Apar[0]*Bmag[1]-1.5*Apar[1]))+(3.0*(geoY[0]*(0.7071067811865475*Apar[0]*Bmag[2]-1.0*Apar[2])+0.3162277660168379*Bmag[1]*(Apar[0]*geoY[1]+geoY[0]*Apar[1]))+5.454823740581938*Apar[1]*geoY[1]*Bmag[2])*BmagInv[2]+1.060660171779821*(BmagInv[0]*(2.0*geoY[0]*Apar[2]+4.024922359499621*Apar[1]*geoY[1])+geoY[0]*(4.024922359499621*Apar[1]*BmagInv[1]+2.23606797749979*Apar[0]*BmagInv[0]))*Bmag[2]+geoY[0]*(1.060660171779821*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(1.060660171779821*BmagInv[0]*Bmag[1]-1.5*BmagInv[1]))))*dfac_x2-1.0*(1.341640786499874*Apar[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x2+1.732050807568877*((2.0*Gradpar[2]+2.23606797749979*Gradpar[0])*Phi[2]+Gradpar[1]*Phi[1])*dfac_x))*q_)/m_); 
  alphav[2] = 1.732050807568877*(((Bmag[2]*(3.928571428571428*BmagInv[2]+2.23606797749979*BmagInv[0])+Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*(0.5*Bmag[1]*BmagInv[2]+BmagInv[1]*Bmag[2]))*geoY[2]+((2.23606797749979*geoY[0]*Bmag[2]+Bmag[1]*geoY[1])*BmagInv[2]+(4.5*BmagInv[1]*geoY[1]+2.5*BmagInv[0]*geoY[0])*Bmag[2]+1.118033988749895*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(geoY[1]*Bmag[2]*BmagInv[2]+0.5*(2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))))*dfac_x2; 
  alphav[3] = 1.732050807568877*(((2.23606797749979*(0.7857142857142857*Bmag[1]*BmagInv[2]+2.571428571428571*BmagInv[1]*Bmag[2])+BmagInv[0]*Bmag[1])*Phi[2]+Phi[1]*(Bmag[2]*(1.756910553749835*BmagInv[2]+BmagInv[0])+0.4472135954999579*Bmag[1]*BmagInv[1]))*geoY[2]+((5.749889084999459*geoY[1]*Bmag[2]+geoY[0]*Bmag[1])*BmagInv[2]+0.5*(9.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])))*Phi[2]+Phi[1]*((geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])*BmagInv[2]+0.5*((4.024922359499621*BmagInv[1]*geoY[1]+2.23606797749979*BmagInv[0]*geoY[0])*Bmag[2]+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))))*dfac_x2; 
  alphav[4] = dfac_v*((Phi[2]*(11.78571428571429*BmagInv[0]*Bmag[2]*geoY[2]+(11.78571428571429*geoY[0]*Bmag[2]+5.270731661249505*Bmag[1]*geoY[1])*BmagInv[2]+3.0*(2.23606797749979*BmagInv[0]*geoY[0]*Bmag[2]+Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+17.24966725499838*BmagInv[1]*geoY[1]*Bmag[2])+(2.23606797749979*(8.571428571428571*Bmag[2]*BmagInv[2]+2.357142857142857*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*(2.23606797749979*(0.4285714285714285*Bmag[1]*BmagInv[2]+2.357142857142857*BmagInv[1]*Bmag[2])+1.5*BmagInv[0]*Bmag[1]))*geoY[2]+Phi[1]*((5.270731661249505*geoY[1]*Bmag[2]+1.5*geoY[0]*Bmag[1])*BmagInv[2]+3.0*(BmagInv[1]*(geoY[0]*Bmag[2]+0.4472135954999579*Bmag[1]*geoY[1])+BmagInv[0]*geoY[1]*Bmag[2])))*dfac_x2*wv+(((Phi[2]*((-5.270731661249505*Apar[1]*BmagInv[1]*geoY[2])+Bmag[2]*(0.3535533905932737*(57.85714285714285*Apar[1]*geoY[1]*BmagInv[2]+(57.85714285714285*BmagInv[1]*geoY[1]+23.57142857142857*BmagInv[0]*geoY[0])*Apar[2])+2.23606797749979*(5.454823740581938*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1])*geoY[1]+geoY[0]*(5.454823740581938*Apar[1]*BmagInv[1]+2.121320343559642*Apar[0]*BmagInv[0])))+(geoY[0]*(13.55261854357878*Apar[2]+8.333758492555738*Apar[0])*Bmag[2]+(6.060915267313265*Bmag[1]*geoY[1]-11.78571428571429*geoY[0])*Apar[2]+2.23606797749979*(1.666751698511147*Apar[0]*Bmag[1]*geoY[1]+Apar[1]*(1.666751698511147*geoY[0]*Bmag[1]-2.357142857142857*geoY[1])))*BmagInv[2]+Bmag[1]*(BmagInv[1]*(3.726970099484161*geoY[0]*Apar[2]+5.454823740581938*Apar[1]*geoY[1])+3.726970099484161*BmagInv[0]*geoY[1]*Apar[2])+3.0*(BmagInv[0]*((0.7071067811865475*Apar[0]*Bmag[1]-1.0*Apar[1])*geoY[1]-2.23606797749979*geoY[0]*Apar[2])+geoY[0]*(0.7071067811865475*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(0.7071067811865475*BmagInv[0]*Bmag[1]-1.0*BmagInv[1])))-17.24966725499838*BmagInv[1]*geoY[1]*Apar[2])+((Bmag[2]*(Apar[2]*(25.13902355192433*BmagInv[2]+13.55261854357878*BmagInv[0])+0.3535533905932737*(57.85714285714285*Apar[1]*BmagInv[1]+23.57142857142857*Apar[0]*BmagInv[0]))+(11.18033988749895*(1.212183053462653*Apar[0]*Bmag[2]-1.714285714285714*Apar[2])+6.060915267313265*Apar[1]*Bmag[1])*BmagInv[2]+(6.060915267313265*Bmag[1]*BmagInv[1]-11.78571428571429*BmagInv[0])*Apar[2]+3.726970099484161*Bmag[1]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))*Phi[2]+Phi[1]*((6.060915267313265*Apar[1]*Bmag[2]+2.272843225242474*Bmag[1]*Apar[2]+2.23606797749979*(0.3030457633656632*Apar[0]*Bmag[1]-0.4285714285714285*Apar[1]))*BmagInv[2]+BmagInv[1]*(2.23606797749979*(1.666751698511147*Apar[0]*Bmag[2]-2.357142857142857*Apar[2])+6.060915267313265*Apar[2]*Bmag[2])+0.3535533905932737*(10.54146332249901*BmagInv[0]*Apar[1]*Bmag[2]+Bmag[1]*(1.91662969499982*BmagInv[0]*Apar[2]+4.714285714285714*Apar[1]*BmagInv[1]))+BmagInv[0]*(1.060660171779821*Apar[0]*Bmag[1]-1.5*Apar[1])))*geoY[2]+Phi[1]*((geoY[1]*(2.23606797749979*(1.666751698511147*Apar[0]*Bmag[2]-2.357142857142857*Apar[2])+6.060915267313265*Apar[2]*Bmag[2])+0.3535533905932737*(10.54146332249901*geoY[0]*Apar[1]*Bmag[2]+Bmag[1]*(1.91662969499982*geoY[0]*Apar[2]+4.714285714285714*Apar[1]*geoY[1]))+geoY[0]*(1.060660171779821*Apar[0]*Bmag[1]-1.5*Apar[1]))*BmagInv[2]+3.0*(0.3535533905932737*(2.0*BmagInv[0]*geoY[0]*Apar[1]*Bmag[2]+Bmag[1]*(BmagInv[0]*geoY[0]*Apar[2]+0.8944271909999159*(Apar[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])+Apar[0]*BmagInv[1]*geoY[1])))-1.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Apar[2])+0.3535533905932737*(BmagInv[1]*(6.0*Apar[0]*geoY[0]*Bmag[2]+4.714285714285714*Bmag[1]*geoY[1]*Apar[2])+6.0*Apar[0]*BmagInv[0]*geoY[1]*Bmag[2])+(BmagInv[1]*(3.726970099484161*geoY[0]*Apar[2]+5.454823740581938*Apar[1]*geoY[1])+3.726970099484161*BmagInv[0]*geoY[1]*Apar[2])*Bmag[2]))*dfac_x2-1.0*(1.341640786499874*Apar[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x2+1.732050807568877*(2.0*Gradpar[1]*Phi[2]+Phi[1]*Gradpar[2])*dfac_x))*q_)/m_); 
  alphav[6] = 3.872983346207417*(((2.857142857142857*Bmag[2]*BmagInv[2]+0.7857142857142857*(2.23606797749979*BmagInv[0]*Bmag[2]+Bmag[1]*BmagInv[1]))*Phi[2]+Phi[1]*(0.1428571428571428*Bmag[1]*BmagInv[2]+0.7857142857142857*BmagInv[1]*Bmag[2]+0.223606797749979*BmagInv[0]*Bmag[1]))*geoY[2]+(0.7857142857142857*(2.23606797749979*geoY[0]*Bmag[2]+Bmag[1]*geoY[1])*BmagInv[2]+(2.571428571428572*BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*Bmag[2]+0.4472135954999579*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*((0.7857142857142857*geoY[1]*Bmag[2]+0.223606797749979*geoY[0]*Bmag[1])*BmagInv[2]+BmagInv[1]*(0.4472135954999579*geoY[0]*Bmag[2]+0.2*Bmag[1]*geoY[1])+0.4472135954999579*BmagInv[0]*geoY[1]*Bmag[2]))*dfac_x2; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.25*(1.732050807568877*alphav[2]-1.0*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(1.732050807568877*alphav[2]+alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*((-0.7745966692414833*alphav[6])+0.4472135954999579*alphav[4]+1.161895003862225*alphav[3]-0.8660254037844386*alphav[2]-0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.25*(0.9682458365518543*alphav[6]-0.5590169943749475*alphav[4]-0.8660254037844386*alphav[2]+0.5*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.25*((-0.7745966692414833*alphav[6])+0.4472135954999579*alphav[4]-1.161895003862225*alphav[3]-0.8660254037844386*alphav[2]+0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*(0.7745966692414833*alphav[6]+0.4472135954999579*alphav[4]-1.161895003862225*alphav[3]+0.8660254037844386*alphav[2]-0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.25*((-0.9682458365518543*alphav[6])-0.5590169943749475*alphav[4]+0.8660254037844386*alphav[2]+0.5*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.25*(0.7745966692414833*alphav[6]+0.4472135954999579*alphav[4]+1.161895003862225*alphav[3]+0.8660254037844386*alphav[2]+0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  out[1] += 0.8660254037844386*(alphax[7]*f[7]+alphax[6]*f[6]+alphax[5]*f[5]+alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphav[6]*f[6]+alphav[4]*f[4]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[3] += 0.1*(7.745966692414834*(alphax[3]*f[7]+f[3]*alphax[7])+(8.660254037844387*alphax[4]+7.745966692414834*alphav[3])*f[6]+8.660254037844387*f[4]*alphax[6]+7.745966692414834*(f[3]*alphav[6]+alphax[2]*f[5]+f[2]*alphax[5]+alphav[1]*f[4]+f[1]*alphav[4])+8.660254037844386*((alphav[2]+alphax[1])*f[3]+f[1]*alphax[3]+f[2]*(alphav[3]+alphax[0])+f[0]*alphax[2]+alphav[0]*f[1]+f[0]*alphav[1])); 
  out[4] += 0.1*(19.36491673103708*(alphax[5]*f[7]+f[5]*alphax[7])+17.32050807568877*(alphax[3]*f[6]+f[3]*alphax[6])+17.32050807568877*(alphax[1]*f[4]+f[1]*alphax[4])+19.36491673103709*(alphax[2]*f[3]+f[2]*alphax[3]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[5] += 0.1*(17.32050807568877*alphav[3]*f[7]+19.36491673103708*(alphav[4]*f[6]+f[4]*alphav[6])+17.32050807568877*alphav[2]*f[5]+19.36491673103709*(alphav[1]*f[3]+f[1]*alphav[3]+alphav[0]*f[2]+f[0]*alphav[2])); 
  out[6] += 0.01428571428571429*((108.4435336938077*alphax[6]+121.2435565298214*alphax[2])*f[7]+(108.4435336938077*f[6]+121.2435565298214*f[2])*alphax[7]+(38.72983346207417*alphav[6]+60.6217782649107*alphav[2])*f[6]+121.2435565298214*(alphax[1]*f[6]+f[1]*alphax[6])+60.6217782649107*f[2]*alphav[6]+121.2435565298214*(alphax[3]*f[5]+f[3]*alphax[5])+(38.72983346207417*alphav[4]+121.2435565298214*alphax[3]+60.62177826491071*alphav[0])*f[4]+121.2435565298214*f[3]*alphax[4]+60.62177826491071*f[0]*alphav[4]+54.22176684690384*alphav[3]*f[3]+135.5544171172596*(alphax[0]*f[3]+f[0]*alphax[3]+alphax[1]*f[2])+f[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphav[1])); 
  out[7] += 0.01428571428571429*((38.72983346207417*alphax[7]+108.4435336938077*alphav[6]+121.2435565298214*alphav[2])*f[7]+60.6217782649107*(alphax[1]*f[7]+f[1]*alphax[7])+54.22176684690384*alphax[6]*f[6]+121.2435565298214*(alphav[1]*f[6]+f[1]*alphav[6])+(38.72983346207417*alphax[5]+121.2435565298214*alphav[3])*f[5]+60.62177826491071*(alphax[0]*f[5]+f[0]*alphax[5])+121.2435565298214*alphav[3]*f[4]+f[3]*(121.2435565298214*alphav[4]+54.22176684690384*alphax[3])+135.5544171172596*(alphav[0]*f[3]+f[0]*alphav[3])+54.22176684690384*alphax[2]*f[2]+135.5544171172596*(alphav[1]*f[2]+f[1]*alphav[2])); 
  return cflFreq; 
} 
