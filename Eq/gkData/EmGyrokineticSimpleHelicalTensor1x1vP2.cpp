#include <GyrokineticModDecl.h> 
double EmGyrokineticSimpleHelicalVol1x1vTensorP2_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // f: Distribution function.
  // out: output increment.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wvpar = w[1];
  double rdvpar2 = 2.0/dxv[1];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;

  double hamil[9]; 
  hamil[0] = (0.3333333333333333*(3.0*rdvpar2Sq*(m_*wvparSq+1.414213562373095*phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 1.414213562373095*phi[1]*q_; 
  hamil[2] = (1.154700538379252*m_*wvpar)/rdvpar2; 
  hamil[4] = 1.414213562373095*phi[2]*q_; 
  hamil[5] = (0.2981423969999719*m_)/rdvpar2Sq; 

  double BstarZdBmag[9]; 
  BstarZdBmag[0] = 1.414213562373095*cmag[0]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[9]; 
  alphax[0] = (0.8660254037844386*BstarZdBmag[0]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[2] = (1.936491673103709*BstarZdBmag[0]*hamil[5]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.25*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.25*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*(0.5*alphax[0]-0.6708203932499369*alphax[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.6708203932499369*alphax[2]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*(0.5*alphax[0]-0.6708203932499369*alphax[2]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.6708203932499369*alphax[2]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[9]; 
  alphavpar[0] = -(0.8660254037844386*BstarZdBmag[0]*hamil[1]*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(1.936491673103709*BstarZdBmag[0]*hamil[4]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.25*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.25*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*(0.5*alphavpar[0]-0.6708203932499369*alphavpar[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*(0.5*alphavpar[0]-0.6708203932499369*alphavpar[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.8660254037844386*(alphax[2]*f[2]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[3] += 0.1732050807568877*(4.47213595499958*(alphax[2]*f[5]+alphavpar[1]*f[4])+5.0*(alphax[0]*f[2]+f[0]*alphax[2]+alphavpar[0]*f[1]+f[0]*alphavpar[1])); 
  out[4] += 1.936491673103709*(alphax[2]*f[3]+alphax[0]*f[1]); 
  out[5] += 1.936491673103709*(alphavpar[1]*f[3]+alphavpar[0]*f[2]); 
  out[6] += 0.1*(17.32050807568877*alphax[2]*f[7]+8.660254037844387*alphavpar[0]*f[4]+19.36491673103708*alphax[0]*f[3]+f[1]*(19.36491673103708*alphax[2]+7.745966692414834*alphavpar[1])); 
  out[7] += 0.1*(17.32050807568877*alphavpar[1]*f[6]+8.660254037844387*alphax[0]*f[5]+19.36491673103708*alphavpar[0]*f[3]+(7.745966692414834*alphax[2]+19.36491673103708*alphavpar[1])*f[2]); 
  out[8] += 0.5*(3.872983346207417*(alphax[0]*f[7]+alphavpar[0]*f[6])+3.464101615137754*(alphax[2]+alphavpar[1])*f[3]); 
  return cflFreq; 
} 
double EmGyrokineticSimpleHelicalVol1x1vTensorP2_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // f: Distribution function.
  // out: output increment.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wvpar = w[1];
  double rdvpar2 = 2.0/dxv[1];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;

  double hamil[9]; 
  hamil[0] = (0.3333333333333333*(3.0*rdvpar2Sq*(m_*wvparSq+1.414213562373095*phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 1.414213562373095*phi[1]*q_; 
  hamil[2] = (1.154700538379252*m_*wvpar)/rdvpar2; 
  hamil[4] = 1.414213562373095*phi[2]*q_; 
  hamil[5] = (0.2981423969999719*m_)/rdvpar2Sq; 

  double BstarZdBmag[9]; 
  BstarZdBmag[0] = 1.414213562373095*cmag[0]; 
  BstarZdBmag[1] = 1.414213562373095*cmag[1]; 
  BstarZdBmag[4] = 1.414213562373095*cmag[2]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[9]; 
  alphax[0] = (0.8660254037844386*BstarZdBmag[0]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.8660254037844386*BstarZdBmag[1]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[2] = (1.936491673103709*BstarZdBmag[0]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphax[3] = (1.936491673103709*BstarZdBmag[1]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphax[4] = (0.8660254037844386*hamil[2]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
  alphax[6] = (1.936491673103709*BstarZdBmag[4]*hamil[5]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.25*(2.23606797749979*alphax[4]-1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.25*(2.23606797749979*alphax[4]+1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*((-1.5*alphax[6])+1.118033988749895*alphax[4]+1.161895003862225*alphax[3]-0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(1.118033988749895*alphax[4]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(1.5*alphax[6]+1.118033988749895*alphax[4]-1.161895003862225*alphax[3]+0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*((-1.5*alphax[6])+1.118033988749895*alphax[4]-1.161895003862225*alphax[3]-0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(1.118033988749895*alphax[4]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(1.5*alphax[6]+1.118033988749895*alphax[4]+1.161895003862225*alphax[3]+0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[9]; 
  alphavpar[0] = -(0.8660254037844386*(2.23606797749979*BstarZdBmag[1]*hamil[4]+BstarZdBmag[0]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(1.732050807568877*(0.5*(2.23606797749979*BstarZdBmag[0]*hamil[4]+BstarZdBmag[1]*hamil[1])+BstarZdBmag[4]*hamil[4])*rdvpar2*rdx2)/m_; 
  alphavpar[4] = (1.732050807568877*((-1.0*BstarZdBmag[1]*hamil[4])-0.5*hamil[1]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.25*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.25*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*(0.4472135954999579*alphavpar[4]-0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.5*alphavpar[0]-0.5590169943749475*alphavpar[4]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.4472135954999579*alphavpar[4]+0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*(0.4472135954999579*alphavpar[4]-0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.5*alphavpar[0]-0.5590169943749475*alphavpar[4]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.4472135954999579*alphavpar[4]+0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.8660254037844386*(alphax[6]*f[6]+alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphavpar[4]*f[4]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[3] += 0.1*(7.745966692414834*(alphax[6]*f[8]+alphax[3]*f[7])+8.660254037844387*(alphax[4]*f[6]+f[4]*alphax[6])+7.745966692414834*(alphax[2]*f[5]+alphavpar[1]*f[4]+f[1]*alphavpar[4])+8.660254037844386*(alphax[1]*f[3]+f[1]*alphax[3]+alphax[0]*f[2]+f[0]*alphax[2]+alphavpar[0]*f[1]+f[0]*alphavpar[1])); 
  out[4] += 0.1*(17.32050807568877*(alphax[3]*f[6]+f[3]*alphax[6])+17.32050807568877*(alphax[1]*f[4]+f[1]*alphax[4])+19.36491673103709*(alphax[2]*f[3]+f[2]*alphax[3]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[5] += 1.936491673103709*(alphavpar[4]*f[6]+alphavpar[1]*f[3]+alphavpar[0]*f[2]); 
  out[6] += 0.01428571428571429*(108.4435336938077*(alphax[3]*f[8]+alphax[6]*f[7])+121.2435565298214*(alphax[2]*f[7]+alphax[1]*f[6]+f[1]*alphax[6])+121.2435565298214*alphax[3]*f[5]+(38.72983346207417*alphavpar[4]+121.2435565298214*alphax[3]+60.62177826491071*alphavpar[0])*f[4]+121.2435565298214*f[3]*alphax[4]+60.62177826491071*f[0]*alphavpar[4]+135.5544171172596*(alphax[0]*f[3]+f[0]*alphax[3]+alphax[1]*f[2])+f[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphavpar[1])); 
  out[7] += 0.1*(8.660254037844387*alphax[4]*f[8]+8.660254037844386*alphax[1]*f[7]+(7.745966692414834*alphax[6]+17.32050807568877*alphavpar[1])*f[6]+8.660254037844387*alphax[0]*f[5]+f[3]*(17.32050807568877*alphavpar[4]+7.745966692414834*alphax[3]+19.36491673103708*alphavpar[0])+(7.745966692414834*alphax[2]+19.36491673103708*alphavpar[1])*f[2]); 
  out[8] += 0.01428571428571429*(121.2435565298214*alphax[1]*f[8]+(121.2435565298214*alphax[4]+135.5544171172596*alphax[0])*f[7]+(86.60254037844389*alphavpar[4]+108.4435336938077*alphax[3]+135.5544171172596*alphavpar[0])*f[6]+108.4435336938077*f[3]*alphax[6]+135.5544171172596*(alphax[1]*f[5]+f[2]*alphavpar[4])+121.2435565298214*((alphax[2]+alphavpar[1])*f[3]+f[2]*alphax[3])); 
  return cflFreq; 
} 
double EmGyrokineticSimpleHelicalStep2Vol1x1vTensorP2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // dApardt: time derivative of Apar.
  // f: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[2]*f[4]+dApardt[1]*f[1]+dApardt[0]*f[0])*q_*rdvpar2)/m_; 
  out[3] += -(0.2449489742783178*(4.47213595499958*dApardt[1]*f[4]+4.47213595499958*f[1]*dApardt[2]+5.0*dApardt[0]*f[1]+5.0*f[0]*dApardt[1])*q_*rdvpar2)/m_; 
  out[5] += -(0.7071067811865475*(3.872983346207417*dApardt[2]*f[6]+3.872983346207417*dApardt[1]*f[3]+3.872983346207417*dApardt[0]*f[2])*q_*rdvpar2)/m_; 
  out[6] += -(0.07824607964359516*(10.0*dApardt[2]*f[4]+15.65247584249853*dApardt[0]*f[4]+15.65247584249853*f[0]*dApardt[2]+14.0*dApardt[1]*f[1])*q_*rdvpar2)/m_; 
  out[7] += -(0.1414213562373095*(17.32050807568877*dApardt[1]*f[6]+17.32050807568877*dApardt[2]*f[3]+19.36491673103708*dApardt[0]*f[3]+19.36491673103708*dApardt[1]*f[2])*q_*rdvpar2)/m_; 
  out[8] += -(0.1010152544552211*(17.32050807568877*dApardt[2]*f[6]+27.11088342345192*dApardt[0]*f[6]+24.24871130596428*dApardt[1]*f[3]+27.11088342345192*dApardt[2]*f[2])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -(0.3535533905932737*dApardt[0]*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = -(0.3535533905932737*dApardt[0]*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.25*(0.6324555320336759*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.25*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.25*(0.6324555320336759*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.25*(0.6324555320336759*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.25*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.25*(0.6324555320336759*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

return cflFreq; 
} 
