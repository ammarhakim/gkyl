#include <DeltaFGyrokineticModDecl.h> 
double LinearDeltaFGyrokineticGenGeoVol1x1vSerP2_Bvarsxz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0, const double *f1, double *out) 
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

  double hamil0[8]; 
  hamil0[0] = (0.3333333333333333*m_*(3.0*rdvpar2Sq*wvparSq+1.0))/rdvpar2Sq; 
  hamil0[2] = (1.154700538379252*m_*wvpar)/rdvpar2; 
  hamil0[5] = (0.2981423969999719*m_)/rdvpar2Sq; 

  double hamil1[8]; 
  hamil1[0] = 1.414213562373095*phi[0]*q_; 
  hamil1[1] = 1.414213562373095*phi[1]*q_; 
  hamil1[4] = 1.414213562373095*phi[2]*q_; 

  double BstarZdBmag[8]; 
  BstarZdBmag[0] = (1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmag[1] = (0.2*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2*wvpar+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmag[2] = ((2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[3] = ((b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (0.02857142857142857*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2*wvpar+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmag[6] = (1.0*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alpha0x[8]; 
  alpha0x[0] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[2]*hamil0[5]+BstarZdBmag[0]*hamil0[2])*rdvpar2*rdx2)/m_; 
  alpha0x[1] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[3]*hamil0[5]+BstarZdBmag[1]*hamil0[2])*rdvpar2*rdx2)/m_; 
  alpha0x[2] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[0]*hamil0[5]+BstarZdBmag[2]*hamil0[2])*rdvpar2*rdx2)/m_; 
  alpha0x[3] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[1]*hamil0[5]+hamil0[2]*BstarZdBmag[3])*rdvpar2*rdx2)/m_; 
  alpha0x[4] = (0.5*(3.872983346207417*hamil0[5]*BstarZdBmag[6]+1.732050807568877*hamil0[2]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
  alpha0x[5] = (1.732050807568877*BstarZdBmag[2]*hamil0[5]*rdvpar2*rdx2)/m_; 
  alpha0x[6] = (0.5*(1.732050807568877*hamil0[2]*BstarZdBmag[6]+3.872983346207417*BstarZdBmag[4]*hamil0[5])*rdvpar2*rdx2)/m_; 
  alpha0x[7] = (1.732050807568877*BstarZdBmag[3]*hamil0[5]*rdvpar2*rdx2)/m_; 
  double alpha1x[8]; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*((-0.7745966692414833*alpha0x[7])-1.5*alpha0x[6]+0.4472135954999579*alpha0x[5]+1.118033988749895*alpha0x[4]+1.161895003862225*alpha0x[3]-0.6708203932499369*alpha0x[2]-0.8660254037844386*alpha0x[1]+0.5*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.9682458365518543*alpha0x[7]-0.5590169943749475*alpha0x[5]+1.118033988749895*alpha0x[4]-0.8660254037844386*alpha0x[1]+0.5*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*((-0.7745966692414833*alpha0x[7])+1.5*alpha0x[6]+0.4472135954999579*alpha0x[5]+1.118033988749895*alpha0x[4]-1.161895003862225*alpha0x[3]+0.6708203932499369*alpha0x[2]-0.8660254037844386*alpha0x[1]+0.5*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*(0.7745966692414833*alpha0x[7]-1.5*alpha0x[6]+0.4472135954999579*alpha0x[5]+1.118033988749895*alpha0x[4]-1.161895003862225*alpha0x[3]-0.6708203932499369*alpha0x[2]+0.8660254037844386*alpha0x[1]+0.5*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*((-0.9682458365518543*alpha0x[7])-0.5590169943749475*alpha0x[5]+1.118033988749895*alpha0x[4]+0.8660254037844386*alpha0x[1]+0.5*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.7745966692414833*alpha0x[7]+1.5*alpha0x[6]+0.4472135954999579*alpha0x[5]+1.118033988749895*alpha0x[4]+1.161895003862225*alpha0x[3]+0.6708203932499369*alpha0x[2]+0.8660254037844386*alpha0x[1]+0.5*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0vpar[8]; 
  double alpha1vpar[8]; 
  alpha1vpar[0] = -(0.8660254037844386*(2.23606797749979*BstarZdBmag[1]*hamil1[4]+BstarZdBmag[0]*hamil1[1])*rdvpar2*rdx2)/m_; 
  alpha1vpar[1] = -(1.732050807568877*(0.5*(2.23606797749979*BstarZdBmag[0]*hamil1[4]+BstarZdBmag[1]*hamil1[1])+BstarZdBmag[4]*hamil1[4])*rdvpar2*rdx2)/m_; 
  alpha1vpar[2] = -(0.8660254037844386*(2.23606797749979*BstarZdBmag[3]*hamil1[4]+hamil1[1]*BstarZdBmag[2])*rdvpar2*rdx2)/m_; 
  alpha1vpar[3] = -(1.0*(1.732050807568877*hamil1[4]*BstarZdBmag[6]+0.8660254037844386*(2.23606797749979*BstarZdBmag[2]*hamil1[4]+hamil1[1]*BstarZdBmag[3]))*rdvpar2*rdx2)/m_; 
  alpha1vpar[4] = (1.732050807568877*((-1.0*BstarZdBmag[1]*hamil1[4])-0.5*hamil1[1]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
  alpha1vpar[6] = (((-0.8660254037844386*hamil1[1]*BstarZdBmag[6])-1.732050807568877*BstarZdBmag[3]*hamil1[4])*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*((-0.7745966692414833*alpha1vpar[6])+0.4472135954999579*alpha1vpar[4]+1.161895003862225*alpha1vpar[3]-0.8660254037844386*alpha1vpar[2]-0.6708203932499369*alpha1vpar[1]+0.5*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.9682458365518543*alpha1vpar[6]-0.5590169943749475*alpha1vpar[4]-0.8660254037844386*alpha1vpar[2]+0.5*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*((-0.7745966692414833*alpha1vpar[6])+0.4472135954999579*alpha1vpar[4]-1.161895003862225*alpha1vpar[3]-0.8660254037844386*alpha1vpar[2]+0.6708203932499369*alpha1vpar[1]+0.5*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*(0.7745966692414833*alpha1vpar[6]+0.4472135954999579*alpha1vpar[4]-1.161895003862225*alpha1vpar[3]+0.8660254037844386*alpha1vpar[2]-0.6708203932499369*alpha1vpar[1]+0.5*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*((-0.9682458365518543*alpha1vpar[6])-0.5590169943749475*alpha1vpar[4]+0.8660254037844386*alpha1vpar[2]+0.5*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.7745966692414833*alpha1vpar[6]+0.4472135954999579*alpha1vpar[4]+1.161895003862225*alpha1vpar[3]+0.8660254037844386*alpha1vpar[2]+0.6708203932499369*alpha1vpar[1]+0.5*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[1] += 0.8660254037844386*(alpha0x[7]*f1[7]+alpha0x[6]*f1[6]+alpha0x[5]*f1[5]+alpha0x[4]*f1[4]+alpha0x[3]*f1[3]+alpha0x[2]*f1[2]+alpha0x[1]*f1[1]+alpha0x[0]*f1[0]); 
  out[2] += 0.8660254037844386*(alpha1vpar[6]*f0[6]+alpha1vpar[4]*f0[4]+alpha1vpar[3]*f0[3]+alpha1vpar[2]*f0[2]+alpha1vpar[1]*f0[1]+alpha1vpar[0]*f0[0]); 
  out[3] += 0.1*(7.745966692414834*(alpha0x[3]*f1[7]+f1[3]*alpha0x[7])+8.660254037844387*alpha0x[4]*f1[6]+7.745966692414834*(alpha1vpar[3]*f0[6]+f0[3]*alpha1vpar[6])+8.660254037844387*f1[4]*alpha0x[6]+7.745966692414834*(alpha0x[2]*f1[5]+f1[2]*alpha0x[5]+alpha1vpar[1]*f0[4]+f0[1]*alpha1vpar[4])+8.660254037844386*(alpha0x[1]*f1[3]+alpha1vpar[2]*f0[3]+f0[2]*alpha1vpar[3]+f1[1]*alpha0x[3]+alpha0x[0]*f1[2]+f1[0]*alpha0x[2]+alpha1vpar[0]*f0[1]+f0[0]*alpha1vpar[1])); 
  out[4] += 0.1*(19.36491673103708*(alpha0x[5]*f1[7]+f1[5]*alpha0x[7])+17.32050807568877*(alpha0x[3]*f1[6]+f1[3]*alpha0x[6])+17.32050807568877*(alpha0x[1]*f1[4]+f1[1]*alpha0x[4])+19.36491673103709*(alpha0x[2]*f1[3]+f1[2]*alpha0x[3]+alpha0x[0]*f1[1]+f1[0]*alpha0x[1])); 
  out[5] += 0.1*(17.32050807568877*alpha1vpar[3]*f0[7]+19.36491673103708*(alpha1vpar[4]*f0[6]+f0[4]*alpha1vpar[6])+17.32050807568877*alpha1vpar[2]*f0[5]+19.36491673103709*(alpha1vpar[1]*f0[3]+f0[1]*alpha1vpar[3]+alpha1vpar[0]*f0[2]+f0[0]*alpha1vpar[2])); 
  out[6] += 0.01428571428571429*((108.4435336938077*alpha0x[6]+121.2435565298214*alpha0x[2])*f1[7]+108.4435336938077*f1[6]*alpha0x[7]+121.2435565298214*(f1[2]*alpha0x[7]+alpha0x[1]*f1[6])+38.72983346207417*alpha1vpar[6]*f0[6]+60.6217782649107*(alpha1vpar[2]*f0[6]+f0[2]*alpha1vpar[6])+121.2435565298214*f1[1]*alpha0x[6]+121.2435565298214*(alpha0x[3]*f1[5]+f1[3]*alpha0x[5]+alpha0x[3]*f1[4])+38.72983346207417*alpha1vpar[4]*f0[4]+60.62177826491071*(alpha1vpar[0]*f0[4]+f0[0]*alpha1vpar[4])+f1[3]*(121.2435565298214*alpha0x[4]+135.5544171172596*alpha0x[0])+54.22176684690384*alpha1vpar[3]*f0[3]+135.5544171172596*(f1[0]*alpha0x[3]+alpha0x[1]*f1[2]+f1[1]*alpha0x[2])+54.22176684690384*alpha1vpar[1]*f0[1]); 
  out[7] += 0.01428571428571429*((38.72983346207417*alpha0x[7]+60.6217782649107*alpha0x[1])*f1[7]+(108.4435336938077*alpha1vpar[6]+121.2435565298214*alpha1vpar[2])*f0[7]+60.6217782649107*f1[1]*alpha0x[7]+54.22176684690384*alpha0x[6]*f1[6]+121.2435565298214*(alpha1vpar[1]*f0[6]+f0[1]*alpha1vpar[6])+(38.72983346207417*alpha0x[5]+60.62177826491071*alpha0x[0])*f1[5]+121.2435565298214*alpha1vpar[3]*f0[5]+60.62177826491071*f1[0]*alpha0x[5]+121.2435565298214*(alpha1vpar[3]*f0[4]+f0[3]*alpha1vpar[4])+54.22176684690384*alpha0x[3]*f1[3]+135.5544171172596*(alpha1vpar[0]*f0[3]+f0[0]*alpha1vpar[3])+54.22176684690384*alpha0x[2]*f1[2]+135.5544171172596*(alpha1vpar[1]*f0[2]+f0[1]*alpha1vpar[2])); 
  return cflFreq; 
} 
