#include <GyrokineticModDecl.h> 
double GyrokineticGenGeoVol1x1vTensorP2_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f, double *out) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmag[0] = cmag[0]*jacobTotInv[0]; 

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
double GyrokineticGenGeoVol1x1vTensorP2_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f, double *out) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmag[0] = (1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmag[1] = (0.2*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2*wvpar+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmag[2] = ((2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[3] = ((b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (0.02857142857142857*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2*wvpar+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmag[6] = (1.0*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[9]; 
  alphax[0] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[2]*hamil[5]+BstarZdBmag[0]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[3]*hamil[5]+BstarZdBmag[1]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[2] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[0]*hamil[5]+BstarZdBmag[2]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[3] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[1]*hamil[5]+hamil[2]*BstarZdBmag[3])*rdvpar2*rdx2)/m_; 
  alphax[4] = (0.5*(3.872983346207417*hamil[5]*BstarZdBmag[6]+1.732050807568877*hamil[2]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
  alphax[5] = (1.732050807568877*BstarZdBmag[2]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphax[6] = (0.5*(1.732050807568877*hamil[2]*BstarZdBmag[6]+3.872983346207417*BstarZdBmag[4]*hamil[5])*rdvpar2*rdx2)/m_; 
  alphax[7] = (1.732050807568877*BstarZdBmag[3]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphax[8] = (1.732050807568877*hamil[5]*BstarZdBmag[6]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.25*(2.23606797749979*alphax[4]-1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.25*(2.23606797749979*alphax[4]+1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*(alphax[8]-0.7745966692414833*alphax[7]-1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]+1.161895003862225*alphax[3]-0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*((-1.25*alphax[8])+0.9682458365518543*alphax[7]-0.5590169943749475*alphax[5]+1.118033988749895*alphax[4]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(alphax[8]-0.7745966692414833*alphax[7]+1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]-1.161895003862225*alphax[3]+0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*(alphax[8]+0.7745966692414833*alphax[7]-1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]-1.161895003862225*alphax[3]-0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*((-1.25*alphax[8])-0.9682458365518543*alphax[7]-0.5590169943749475*alphax[5]+1.118033988749895*alphax[4]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(alphax[8]+0.7745966692414833*alphax[7]+1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]+1.161895003862225*alphax[3]+0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[9]; 
  alphavpar[0] = -(0.8660254037844386*(2.23606797749979*BstarZdBmag[1]*hamil[4]+BstarZdBmag[0]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(1.732050807568877*(0.5*(2.23606797749979*BstarZdBmag[0]*hamil[4]+BstarZdBmag[1]*hamil[1])+BstarZdBmag[4]*hamil[4])*rdvpar2*rdx2)/m_; 
  alphavpar[2] = -(0.8660254037844386*(2.23606797749979*BstarZdBmag[3]*hamil[4]+hamil[1]*BstarZdBmag[2])*rdvpar2*rdx2)/m_; 
  alphavpar[3] = -(1.0*(1.732050807568877*hamil[4]*BstarZdBmag[6]+0.8660254037844386*(2.23606797749979*BstarZdBmag[2]*hamil[4]+hamil[1]*BstarZdBmag[3]))*rdvpar2*rdx2)/m_; 
  alphavpar[4] = (1.732050807568877*((-1.0*BstarZdBmag[1]*hamil[4])-0.5*hamil[1]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
  alphavpar[6] = (((-0.8660254037844386*hamil[1]*BstarZdBmag[6])-1.732050807568877*BstarZdBmag[3]*hamil[4])*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.25*(1.732050807568877*alphavpar[2]-1.0*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.25*(1.732050807568877*alphavpar[2]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*((-0.7745966692414833*alphavpar[6])+0.4472135954999579*alphavpar[4]+1.161895003862225*alphavpar[3]-0.8660254037844386*alphavpar[2]-0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.9682458365518543*alphavpar[6]-0.5590169943749475*alphavpar[4]-0.8660254037844386*alphavpar[2]+0.5*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*((-0.7745966692414833*alphavpar[6])+0.4472135954999579*alphavpar[4]-1.161895003862225*alphavpar[3]-0.8660254037844386*alphavpar[2]+0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*(0.7745966692414833*alphavpar[6]+0.4472135954999579*alphavpar[4]-1.161895003862225*alphavpar[3]+0.8660254037844386*alphavpar[2]-0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*((-0.9682458365518543*alphavpar[6])-0.5590169943749475*alphavpar[4]+0.8660254037844386*alphavpar[2]+0.5*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.7745966692414833*alphavpar[6]+0.4472135954999579*alphavpar[4]+1.161895003862225*alphavpar[3]+0.8660254037844386*alphavpar[2]+0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.8660254037844386*(alphax[8]*f[8]+alphax[7]*f[7]+alphax[6]*f[6]+alphax[5]*f[5]+alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphavpar[6]*f[6]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[3] += 0.1*(7.745966692414834*(alphax[6]*f[8]+f[6]*alphax[8]+alphax[3]*f[7]+f[3]*alphax[7])+(8.660254037844387*alphax[4]+7.745966692414834*alphavpar[3])*f[6]+8.660254037844387*f[4]*alphax[6]+7.745966692414834*(f[3]*alphavpar[6]+alphax[2]*f[5]+f[2]*alphax[5]+alphavpar[1]*f[4]+f[1]*alphavpar[4])+8.660254037844386*((alphavpar[2]+alphax[1])*f[3]+f[1]*alphax[3]+f[2]*(alphavpar[3]+alphax[0])+f[0]*alphax[2]+alphavpar[0]*f[1]+f[0]*alphavpar[1])); 
  out[4] += 0.1*(17.32050807568877*(alphax[7]*f[8]+f[7]*alphax[8])+19.36491673103708*(alphax[5]*f[7]+f[5]*alphax[7])+17.32050807568877*(alphax[3]*f[6]+f[3]*alphax[6])+17.32050807568877*(alphax[1]*f[4]+f[1]*alphax[4])+19.36491673103709*(alphax[2]*f[3]+f[2]*alphax[3]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[5] += 0.1*(17.32050807568877*(alphavpar[6]*f[8]+alphavpar[3]*f[7])+19.36491673103708*(alphavpar[4]*f[6]+f[4]*alphavpar[6])+17.32050807568877*alphavpar[2]*f[5]+19.36491673103709*(alphavpar[1]*f[3]+f[1]*alphavpar[3]+alphavpar[0]*f[2]+f[0]*alphavpar[2])); 
  out[6] += 0.01428571428571429*(108.4435336938077*(alphax[3]*f[8]+f[3]*alphax[8])+(108.4435336938077*alphax[6]+121.2435565298214*alphax[2])*f[7]+(108.4435336938077*f[6]+121.2435565298214*f[2])*alphax[7]+(38.72983346207417*alphavpar[6]+60.6217782649107*alphavpar[2])*f[6]+121.2435565298214*(alphax[1]*f[6]+f[1]*alphax[6])+60.6217782649107*f[2]*alphavpar[6]+121.2435565298214*(alphax[3]*f[5]+f[3]*alphax[5])+(38.72983346207417*alphavpar[4]+121.2435565298214*alphax[3]+60.62177826491071*alphavpar[0])*f[4]+121.2435565298214*f[3]*alphax[4]+60.62177826491071*f[0]*alphavpar[4]+54.22176684690384*alphavpar[3]*f[3]+135.5544171172596*(alphax[0]*f[3]+f[0]*alphax[3]+alphax[1]*f[2])+f[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphavpar[1])); 
  out[7] += 0.01428571428571429*((38.72983346207417*alphax[8]+60.62177826491071*alphax[4]+108.4435336938077*alphavpar[3])*f[8]+60.62177826491071*f[4]*alphax[8]+(38.72983346207417*alphax[7]+108.4435336938077*alphavpar[6]+121.2435565298214*alphavpar[2])*f[7]+60.6217782649107*(alphax[1]*f[7]+f[1]*alphax[7])+54.22176684690384*alphax[6]*f[6]+121.2435565298214*(alphavpar[1]*f[6]+f[1]*alphavpar[6])+(38.72983346207417*alphax[5]+121.2435565298214*alphavpar[3])*f[5]+60.62177826491071*(alphax[0]*f[5]+f[0]*alphax[5])+121.2435565298214*alphavpar[3]*f[4]+f[3]*(121.2435565298214*alphavpar[4]+54.22176684690384*alphax[3])+135.5544171172596*(alphavpar[0]*f[3]+f[0]*alphavpar[3])+54.22176684690384*alphax[2]*f[2]+135.5544171172596*(alphavpar[1]*f[2]+f[1]*alphavpar[2])); 
  out[8] += 0.01428571428571429*((77.45966692414834*(alphax[7]+alphavpar[6])+121.2435565298214*(alphavpar[2]+alphax[1]))*f[8]+(77.45966692414834*f[7]+121.2435565298214*f[1])*alphax[8]+(86.60254037844389*alphax[5]+121.2435565298214*alphax[4]+108.4435336938077*alphavpar[3]+135.5544171172596*alphax[0])*f[7]+(86.60254037844389*f[5]+121.2435565298214*f[4]+135.5544171172596*f[0])*alphax[7]+(86.60254037844389*alphavpar[4]+108.4435336938077*alphax[3]+135.5544171172596*alphavpar[0])*f[6]+108.4435336938077*f[3]*alphax[6]+(121.2435565298214*f[5]+86.60254037844389*f[4]+135.5544171172596*f[0])*alphavpar[6]+135.5544171172596*(alphax[1]*f[5]+f[1]*alphax[5]+alphavpar[2]*f[4]+f[2]*alphavpar[4])+121.2435565298214*((alphax[2]+alphavpar[1])*f[3]+f[2]*alphax[3]+f[1]*alphavpar[3])); 
  return cflFreq; 
} 
double GyrokineticGenGeoVol1x1vTensorP2_Bvarsz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f, double *out) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmag[0] = cmag[0]*jacobTotInv[0]; 

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
double GyrokineticGenGeoVol1x1vTensorP2_Bvarsxz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f, double *out) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmag[0] = (1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmag[1] = (0.2*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2*wvpar+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmag[2] = ((2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[3] = ((b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (0.02857142857142857*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2*wvpar+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmag[6] = (1.0*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[9]; 
  alphax[0] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[2]*hamil[5]+BstarZdBmag[0]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[3]*hamil[5]+BstarZdBmag[1]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[2] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[0]*hamil[5]+BstarZdBmag[2]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[3] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[1]*hamil[5]+hamil[2]*BstarZdBmag[3])*rdvpar2*rdx2)/m_; 
  alphax[4] = (0.5*(3.872983346207417*hamil[5]*BstarZdBmag[6]+1.732050807568877*hamil[2]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
  alphax[5] = (1.732050807568877*BstarZdBmag[2]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphax[6] = (0.5*(1.732050807568877*hamil[2]*BstarZdBmag[6]+3.872983346207417*BstarZdBmag[4]*hamil[5])*rdvpar2*rdx2)/m_; 
  alphax[7] = (1.732050807568877*BstarZdBmag[3]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphax[8] = (1.732050807568877*hamil[5]*BstarZdBmag[6]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.25*(2.23606797749979*alphax[4]-1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.25*(2.23606797749979*alphax[4]+1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*(alphax[8]-0.7745966692414833*alphax[7]-1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]+1.161895003862225*alphax[3]-0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*((-1.25*alphax[8])+0.9682458365518543*alphax[7]-0.5590169943749475*alphax[5]+1.118033988749895*alphax[4]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(alphax[8]-0.7745966692414833*alphax[7]+1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]-1.161895003862225*alphax[3]+0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*(alphax[8]+0.7745966692414833*alphax[7]-1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]-1.161895003862225*alphax[3]-0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*((-1.25*alphax[8])-0.9682458365518543*alphax[7]-0.5590169943749475*alphax[5]+1.118033988749895*alphax[4]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(alphax[8]+0.7745966692414833*alphax[7]+1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]+1.161895003862225*alphax[3]+0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[9]; 
  alphavpar[0] = -(0.8660254037844386*(2.23606797749979*BstarZdBmag[1]*hamil[4]+BstarZdBmag[0]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(1.732050807568877*(0.5*(2.23606797749979*BstarZdBmag[0]*hamil[4]+BstarZdBmag[1]*hamil[1])+BstarZdBmag[4]*hamil[4])*rdvpar2*rdx2)/m_; 
  alphavpar[2] = -(0.8660254037844386*(2.23606797749979*BstarZdBmag[3]*hamil[4]+hamil[1]*BstarZdBmag[2])*rdvpar2*rdx2)/m_; 
  alphavpar[3] = -(1.0*(1.732050807568877*hamil[4]*BstarZdBmag[6]+0.8660254037844386*(2.23606797749979*BstarZdBmag[2]*hamil[4]+hamil[1]*BstarZdBmag[3]))*rdvpar2*rdx2)/m_; 
  alphavpar[4] = (1.732050807568877*((-1.0*BstarZdBmag[1]*hamil[4])-0.5*hamil[1]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
  alphavpar[6] = (((-0.8660254037844386*hamil[1]*BstarZdBmag[6])-1.732050807568877*BstarZdBmag[3]*hamil[4])*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.25*(1.732050807568877*alphavpar[2]-1.0*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.25*(1.732050807568877*alphavpar[2]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*((-0.7745966692414833*alphavpar[6])+0.4472135954999579*alphavpar[4]+1.161895003862225*alphavpar[3]-0.8660254037844386*alphavpar[2]-0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.9682458365518543*alphavpar[6]-0.5590169943749475*alphavpar[4]-0.8660254037844386*alphavpar[2]+0.5*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*((-0.7745966692414833*alphavpar[6])+0.4472135954999579*alphavpar[4]-1.161895003862225*alphavpar[3]-0.8660254037844386*alphavpar[2]+0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*(0.7745966692414833*alphavpar[6]+0.4472135954999579*alphavpar[4]-1.161895003862225*alphavpar[3]+0.8660254037844386*alphavpar[2]-0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*((-0.9682458365518543*alphavpar[6])-0.5590169943749475*alphavpar[4]+0.8660254037844386*alphavpar[2]+0.5*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.7745966692414833*alphavpar[6]+0.4472135954999579*alphavpar[4]+1.161895003862225*alphavpar[3]+0.8660254037844386*alphavpar[2]+0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.8660254037844386*(alphax[8]*f[8]+alphax[7]*f[7]+alphax[6]*f[6]+alphax[5]*f[5]+alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphavpar[6]*f[6]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[3] += 0.1*(7.745966692414834*(alphax[6]*f[8]+f[6]*alphax[8]+alphax[3]*f[7]+f[3]*alphax[7])+(8.660254037844387*alphax[4]+7.745966692414834*alphavpar[3])*f[6]+8.660254037844387*f[4]*alphax[6]+7.745966692414834*(f[3]*alphavpar[6]+alphax[2]*f[5]+f[2]*alphax[5]+alphavpar[1]*f[4]+f[1]*alphavpar[4])+8.660254037844386*((alphavpar[2]+alphax[1])*f[3]+f[1]*alphax[3]+f[2]*(alphavpar[3]+alphax[0])+f[0]*alphax[2]+alphavpar[0]*f[1]+f[0]*alphavpar[1])); 
  out[4] += 0.1*(17.32050807568877*(alphax[7]*f[8]+f[7]*alphax[8])+19.36491673103708*(alphax[5]*f[7]+f[5]*alphax[7])+17.32050807568877*(alphax[3]*f[6]+f[3]*alphax[6])+17.32050807568877*(alphax[1]*f[4]+f[1]*alphax[4])+19.36491673103709*(alphax[2]*f[3]+f[2]*alphax[3]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[5] += 0.1*(17.32050807568877*(alphavpar[6]*f[8]+alphavpar[3]*f[7])+19.36491673103708*(alphavpar[4]*f[6]+f[4]*alphavpar[6])+17.32050807568877*alphavpar[2]*f[5]+19.36491673103709*(alphavpar[1]*f[3]+f[1]*alphavpar[3]+alphavpar[0]*f[2]+f[0]*alphavpar[2])); 
  out[6] += 0.01428571428571429*(108.4435336938077*(alphax[3]*f[8]+f[3]*alphax[8])+(108.4435336938077*alphax[6]+121.2435565298214*alphax[2])*f[7]+(108.4435336938077*f[6]+121.2435565298214*f[2])*alphax[7]+(38.72983346207417*alphavpar[6]+60.6217782649107*alphavpar[2])*f[6]+121.2435565298214*(alphax[1]*f[6]+f[1]*alphax[6])+60.6217782649107*f[2]*alphavpar[6]+121.2435565298214*(alphax[3]*f[5]+f[3]*alphax[5])+(38.72983346207417*alphavpar[4]+121.2435565298214*alphax[3]+60.62177826491071*alphavpar[0])*f[4]+121.2435565298214*f[3]*alphax[4]+60.62177826491071*f[0]*alphavpar[4]+54.22176684690384*alphavpar[3]*f[3]+135.5544171172596*(alphax[0]*f[3]+f[0]*alphax[3]+alphax[1]*f[2])+f[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphavpar[1])); 
  out[7] += 0.01428571428571429*((38.72983346207417*alphax[8]+60.62177826491071*alphax[4]+108.4435336938077*alphavpar[3])*f[8]+60.62177826491071*f[4]*alphax[8]+(38.72983346207417*alphax[7]+108.4435336938077*alphavpar[6]+121.2435565298214*alphavpar[2])*f[7]+60.6217782649107*(alphax[1]*f[7]+f[1]*alphax[7])+54.22176684690384*alphax[6]*f[6]+121.2435565298214*(alphavpar[1]*f[6]+f[1]*alphavpar[6])+(38.72983346207417*alphax[5]+121.2435565298214*alphavpar[3])*f[5]+60.62177826491071*(alphax[0]*f[5]+f[0]*alphax[5])+121.2435565298214*alphavpar[3]*f[4]+f[3]*(121.2435565298214*alphavpar[4]+54.22176684690384*alphax[3])+135.5544171172596*(alphavpar[0]*f[3]+f[0]*alphavpar[3])+54.22176684690384*alphax[2]*f[2]+135.5544171172596*(alphavpar[1]*f[2]+f[1]*alphavpar[2])); 
  out[8] += 0.01428571428571429*((77.45966692414834*(alphax[7]+alphavpar[6])+121.2435565298214*(alphavpar[2]+alphax[1]))*f[8]+(77.45966692414834*f[7]+121.2435565298214*f[1])*alphax[8]+(86.60254037844389*alphax[5]+121.2435565298214*alphax[4]+108.4435336938077*alphavpar[3]+135.5544171172596*alphax[0])*f[7]+(86.60254037844389*f[5]+121.2435565298214*f[4]+135.5544171172596*f[0])*alphax[7]+(86.60254037844389*alphavpar[4]+108.4435336938077*alphax[3]+135.5544171172596*alphavpar[0])*f[6]+108.4435336938077*f[3]*alphax[6]+(121.2435565298214*f[5]+86.60254037844389*f[4]+135.5544171172596*f[0])*alphavpar[6]+135.5544171172596*(alphax[1]*f[5]+f[1]*alphax[5]+alphavpar[2]*f[4]+f[2]*alphavpar[4])+121.2435565298214*((alphax[2]+alphavpar[1])*f[3]+f[2]*alphax[3]+f[1]*alphavpar[3])); 
  return cflFreq; 
} 
