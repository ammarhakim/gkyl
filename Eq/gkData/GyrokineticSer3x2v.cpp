#include <GyrokineticModDecl.h> 
double GyrokineticVol3x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dyInv = 1.0/dxv[1]; 
  double dzInv = 1.0/dxv[2]; 
  double dvInv = 1.0/dxv[3]; 
  double dmInv = 1.0/dxv[4]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_z = 2.0/dxv[2]; 
  double dfac_v = 2.0/dxv[3]; 
  double dfac_m = 2.0/dxv[4]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wz = w[2]; 
  double wv = w[3]; 
  double wm = w[4]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  out[1] += -0.375*(BmagInv[0]*Phi[7]*f[7]+BmagInv[1]*Phi[6]*f[7]+BmagInv[1]*f[3]*Phi[7]+BmagInv[0]*f[3]*Phi[6]+BmagInv[0]*f[1]*Phi[4]+f[0]*BmagInv[1]*Phi[4]+BmagInv[1]*f[1]*Phi[2]+BmagInv[0]*f[0]*Phi[2])*dfac_x*dfac_y; 
  out[2] += (0.0883883476483184*(6.928203230275509*(BcurvY[1]*f[1]+BcurvY[0]*f[0])*dfac_m*dfac_v*dfac_y*m_*wv2+4.0*(BcurvY[1]*f[9]+BcurvY[0]*f[4])*dfac_m*dfac_y*m_*wv+4.242640687119286*dfac_m*dfac_v*dfac_x*dfac_y*(Bmag[1]*BmagInv[1]*f[1]*wm+BmagInv[0]*f[0]*Bmag[1]*wm+BmagInv[1]*Phi[7]*f[16]*q_+BmagInv[0]*Phi[7]*f[8]*q_+BmagInv[1]*Phi[5]*f[7]*q_+BmagInv[1]*Phi[4]*f[6]*q_+BmagInv[0]*f[3]*Phi[5]*q_+BmagInv[0]*f[2]*Phi[4]*q_+BmagInv[1]*Phi[1]*f[1]*q_+BmagInv[0]*f[0]*Phi[1]*q_)+2.449489742783178*Bmag[1]*(BmagInv[1]*f[12]+BmagInv[0]*f[5])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[3] += 1.732050807568877*f[0]*dfac_z*wv; 
  out[4] += -(0.0883883476483184*(4.242640687119286*(BcurvY[0]*Phi[7]*f[7]+BcurvY[1]*Phi[6]*f[7]+BcurvY[1]*f[3]*Phi[7]+BcurvY[0]*f[3]*Phi[6]+BcurvY[0]*f[1]*Phi[4]+f[0]*BcurvY[1]*Phi[4]+BcurvY[1]*f[1]*Phi[2]+BcurvY[0]*f[0]*Phi[2])*dfac_v*dfac_y*m_*wv+12.0*(f[6]*Phi[7]+f[2]*Phi[6]+f[1]*Phi[5]+f[0]*Phi[3])*dfac_v*dfac_z*q_+2.449489742783178*(BcurvY[0]*Phi[7]*f[18]+BcurvY[1]*Phi[6]*f[18]+BcurvY[1]*Phi[7]*f[11]+BcurvY[0]*Phi[6]*f[11]+BcurvY[0]*Phi[4]*f[9]+BcurvY[1]*Phi[2]*f[9]+BcurvY[1]*Phi[4]*f[4]+BcurvY[0]*Phi[2]*f[4])*dfac_y*m_))/m_; 
  out[6] += (0.125*(4.898979485566357*(BcurvY[0]*f[1]+f[0]*BcurvY[1])*dfac_m*dfac_v*dfac_y*m_*wv2+2.828427124746191*(BcurvY[0]*f[9]+BcurvY[1]*f[4])*dfac_m*dfac_y*m_*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*(BmagInv[0]*Bmag[1]*f[1]*wm+f[0]*Bmag[1]*BmagInv[1]*wm-1.0*BmagInv[1]*Phi[6]*f[16]*q_-1.0*BmagInv[0]*Phi[6]*f[8]*q_+BmagInv[0]*Phi[5]*f[7]*q_-1.0*BmagInv[1]*Phi[2]*f[6]*q_+BmagInv[1]*f[3]*Phi[5]*q_-1.0*BmagInv[0]*Phi[2]*f[2]*q_+BmagInv[0]*Phi[1]*f[1]*q_+f[0]*BmagInv[1]*Phi[1]*q_)+1.732050807568877*Bmag[1]*(BmagInv[0]*f[12]+BmagInv[1]*f[5])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[7] += 0.125*(13.85640646055102*f[1]*dfac_z*wv-3.0*(BmagInv[0]*Phi[4]*f[7]+BmagInv[1]*Phi[2]*f[7]+BmagInv[0]*f[1]*Phi[7]+f[0]*BmagInv[1]*Phi[7]+BmagInv[1]*f[1]*Phi[6]+BmagInv[0]*f[0]*Phi[6]+BmagInv[1]*f[3]*Phi[4]+BmagInv[0]*Phi[2]*f[3])*dfac_x*dfac_y); 
  out[8] += (0.125*(4.898979485566357*(BcurvY[1]*f[7]+BcurvY[0]*f[3])*dfac_m*dfac_v*dfac_y*m_*wv2+13.85640646055102*f[2]*dfac_m*dfac_v*dfac_z*q_*wv+2.828427124746191*(BcurvY[1]*f[18]+BcurvY[0]*f[11])*dfac_m*dfac_y*m_*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*(Bmag[1]*BmagInv[1]*f[7]*wm+BmagInv[0]*Bmag[1]*f[3]*wm+BmagInv[1]*Phi[4]*f[16]*q_+BmagInv[0]*Phi[4]*f[8]*q_+BmagInv[1]*Phi[1]*f[7]*q_+BmagInv[1]*f[6]*Phi[7]*q_+BmagInv[0]*f[2]*Phi[7]*q_+BmagInv[1]*f[1]*Phi[5]*q_+BmagInv[0]*f[0]*Phi[5]*q_+BmagInv[0]*Phi[1]*f[3]*q_)+1.732050807568877*Bmag[1]*(BmagInv[1]*f[21]+BmagInv[0]*f[14])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[9] += -(0.01767766952966368*(4.242640687119286*(9.0*BcurvY[1]*Phi[7]*f[7]+5.0*BcurvY[0]*Phi[6]*f[7]+5.0*BcurvY[0]*f[3]*Phi[7]+5.0*BcurvY[1]*f[3]*Phi[6]+9.0*BcurvY[1]*f[1]*Phi[4]+5.0*BcurvY[0]*f[0]*Phi[4]+5.0*BcurvY[0]*f[1]*Phi[2]+5.0*f[0]*BcurvY[1]*Phi[2])*dfac_v*dfac_y*m_*wv+60.0*(f[2]*Phi[7]+Phi[6]*f[6]+f[0]*Phi[5]+f[1]*Phi[3])*dfac_v*dfac_z*q_+21.21320343559643*(BmagInv[0]*Phi[7]*f[18]+BmagInv[1]*Phi[6]*f[18]+BmagInv[1]*Phi[7]*f[11]+BmagInv[0]*Phi[6]*f[11]+BmagInv[0]*Phi[4]*f[9]+BmagInv[1]*Phi[2]*f[9]+BmagInv[1]*Phi[4]*f[4]+BmagInv[0]*Phi[2]*f[4])*dfac_x*dfac_y*m_+2.449489742783178*(9.0*BcurvY[1]*Phi[7]*f[18]+5.0*BcurvY[0]*Phi[6]*f[18]+5.0*BcurvY[0]*Phi[7]*f[11]+5.0*BcurvY[1]*Phi[6]*f[11]+9.0*BcurvY[1]*Phi[4]*f[9]+5.0*BcurvY[0]*Phi[2]*f[9]+5.0*BcurvY[0]*Phi[4]*f[4]+5.0*BcurvY[1]*Phi[2]*f[4])*dfac_y*m_))/m_; 
  out[10] += (0.125*(1.732050807568877*dfac_m*dfac_v*dfac_y*m_*(2.828427124746191*BcurvY[1]*f[9]*m_*wv2+2.828427124746191*BcurvY[0]*f[4]*m_*wv2-1.0*BcurvY[0]*Phi[7]*f[26]*q_-1.0*BcurvY[1]*Phi[6]*f[26]*q_-1.0*BcurvY[1]*Phi[7]*f[19]*q_-1.0*BcurvY[0]*Phi[6]*f[19]*q_-1.0*BcurvY[0]*Phi[4]*f[17]*q_-1.0*BcurvY[1]*Phi[2]*f[17]*q_-1.0*BcurvY[1]*Phi[4]*f[10]*q_-1.0*BcurvY[0]*Phi[2]*f[10]*q_)-3.0*(BcurvY[0]*Phi[7]*f[16]+BcurvY[1]*Phi[6]*f[16]+BcurvY[1]*Phi[7]*f[8]+BcurvY[0]*Phi[6]*f[8]+BcurvY[0]*Phi[4]*f[6]+BcurvY[1]*Phi[2]*f[6]+BcurvY[1]*f[2]*Phi[4]+BcurvY[0]*Phi[2]*f[2])*dfac_m*dfac_v2*dfac_y*m_*q_*wv+2.828427124746191*(BcurvY[1]*f[1]+BcurvY[0]*f[0])*dfac_m*dfac_y*m2*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*m_*(Bmag[1]*BmagInv[1]*f[9]*wm+BmagInv[0]*Bmag[1]*f[4]*wm+BmagInv[1]*Phi[7]*f[26]*q_+BmagInv[0]*Phi[7]*f[19]*q_+BmagInv[1]*Phi[5]*f[18]*q_+BmagInv[1]*Phi[4]*f[17]*q_+BmagInv[0]*Phi[5]*f[11]*q_+BmagInv[0]*Phi[4]*f[10]*q_+BmagInv[1]*Phi[1]*f[9]*q_+BmagInv[0]*Phi[1]*f[4]*q_)-8.485281374238571*(f[1]*Phi[7]+Phi[5]*f[6]+f[0]*Phi[6]+f[2]*Phi[3])*dfac_m*dfac_v2*dfac_z*q2+1.732050807568877*Bmag[1]*(BmagInv[1]*f[23]+BmagInv[0]*f[15])*dfac_v*dfac_x*dfac_y*m_))/(dfac_m*dfac_v*m_*q_); 
  out[11] += (0.125*(13.85640646055102*f[4]*dfac_z*m_*wv-3.0*(BcurvY[0]*Phi[4]*f[7]+BcurvY[1]*Phi[2]*f[7]+BcurvY[0]*f[1]*Phi[7]+f[0]*BcurvY[1]*Phi[7]+BcurvY[1]*f[1]*Phi[6]+BcurvY[0]*f[0]*Phi[6]+BcurvY[1]*f[3]*Phi[4]+BcurvY[0]*Phi[2]*f[3])*dfac_v*dfac_y*m_*wv-8.485281374238571*(Phi[7]*f[16]+Phi[6]*f[8]+Phi[5]*f[7]+Phi[3]*f[3])*dfac_v*dfac_z*q_-1.732050807568877*(BcurvY[0]*Phi[4]*f[18]+BcurvY[1]*Phi[2]*f[18]+BcurvY[1]*Phi[4]*f[11]+BcurvY[0]*Phi[2]*f[11]+BcurvY[0]*Phi[7]*f[9]+BcurvY[1]*Phi[6]*f[9]+BcurvY[1]*f[4]*Phi[7]+BcurvY[0]*f[4]*Phi[6])*dfac_y*m_))/m_; 
  out[12] += -0.375*(BmagInv[0]*Phi[7]*f[21]+BmagInv[1]*Phi[6]*f[21]+BmagInv[1]*Phi[7]*f[14]+BmagInv[0]*Phi[6]*f[14]+BmagInv[0]*Phi[4]*f[12]+BmagInv[1]*Phi[2]*f[12]+BmagInv[1]*Phi[4]*f[5]+BmagInv[0]*Phi[2]*f[5])*dfac_x*dfac_y; 
  out[13] += (0.0883883476483184*(6.928203230275509*(BcurvY[1]*f[12]+BcurvY[0]*f[5])*dfac_m*dfac_v*dfac_y*m_*wv2+4.0*(BcurvY[1]*f[23]+BcurvY[0]*f[15])*dfac_m*dfac_y*m_*wv+4.242640687119286*dfac_m*dfac_v*dfac_x*dfac_y*(Bmag[1]*BmagInv[1]*f[12]*wm+BmagInv[0]*Bmag[1]*f[5]*wm+BmagInv[1]*Phi[7]*f[27]*q_+BmagInv[0]*Phi[7]*f[22]*q_+BmagInv[1]*Phi[5]*f[21]*q_+BmagInv[1]*Phi[4]*f[20]*q_+BmagInv[0]*Phi[5]*f[14]*q_+BmagInv[0]*Phi[4]*f[13]*q_+BmagInv[1]*Phi[1]*f[12]*q_+BmagInv[0]*Phi[1]*f[5]*q_)+2.449489742783178*Bmag[1]*(BmagInv[1]*f[1]+BmagInv[0]*f[0])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[14] += 1.732050807568877*f[5]*dfac_z*wv; 
  out[15] += -(0.0883883476483184*(4.242640687119286*(BcurvY[0]*Phi[7]*f[21]+BcurvY[1]*Phi[6]*f[21]+BcurvY[1]*Phi[7]*f[14]+BcurvY[0]*Phi[6]*f[14]+BcurvY[0]*Phi[4]*f[12]+BcurvY[1]*Phi[2]*f[12]+BcurvY[1]*Phi[4]*f[5]+BcurvY[0]*Phi[2]*f[5])*dfac_v*dfac_y*m_*wv+12.0*(Phi[7]*f[20]+Phi[6]*f[13]+Phi[5]*f[12]+Phi[3]*f[5])*dfac_v*dfac_z*q_+2.449489742783178*(BcurvY[0]*Phi[7]*f[29]+BcurvY[1]*Phi[6]*f[29]+BcurvY[1]*Phi[7]*f[25]+BcurvY[0]*Phi[6]*f[25]+BcurvY[0]*Phi[4]*f[23]+BcurvY[1]*Phi[2]*f[23]+BcurvY[1]*Phi[4]*f[15]+BcurvY[0]*Phi[2]*f[15])*dfac_y*m_))/m_; 
  out[16] += (0.0883883476483184*(6.928203230275509*(BcurvY[0]*f[7]+BcurvY[1]*f[3])*dfac_m*dfac_v*dfac_y*m_*wv2+19.59591794226543*f[6]*dfac_m*dfac_v*dfac_z*q_*wv+4.0*(BcurvY[0]*f[18]+BcurvY[1]*f[11])*dfac_m*dfac_y*m_*wv+4.242640687119286*dfac_m*dfac_v*dfac_x*dfac_y*(BmagInv[0]*Bmag[1]*f[7]*wm+Bmag[1]*BmagInv[1]*f[3]*wm-1.0*BmagInv[1]*Phi[2]*f[16]*q_-1.0*BmagInv[0]*Phi[2]*f[8]*q_+BmagInv[0]*Phi[1]*f[7]*q_-1.0*BmagInv[1]*Phi[6]*f[6]*q_-1.0*BmagInv[0]*f[2]*Phi[6]*q_+BmagInv[0]*f[1]*Phi[5]*q_+f[0]*BmagInv[1]*Phi[5]*q_+BmagInv[1]*Phi[1]*f[3]*q_)+2.449489742783178*Bmag[1]*(BmagInv[0]*f[21]+BmagInv[1]*f[14])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[17] += (0.01767766952966368*(1.732050807568877*dfac_m*dfac_v*dfac_y*m_*(20.0*BcurvY[0]*f[9]*m_*wv2+20.0*BcurvY[1]*f[4]*m_*wv2-12.72792206135786*BcurvY[1]*Phi[7]*f[26]*q_-7.071067811865476*BcurvY[0]*Phi[6]*f[26]*q_-7.071067811865476*BcurvY[0]*Phi[7]*f[19]*q_-7.071067811865476*BcurvY[1]*Phi[6]*f[19]*q_-12.72792206135786*BcurvY[1]*Phi[4]*f[17]*q_-7.071067811865476*BcurvY[0]*Phi[2]*f[17]*q_-7.071067811865476*BcurvY[0]*Phi[4]*f[10]*q_-7.071067811865476*BcurvY[1]*Phi[2]*f[10]*q_)-4.242640687119286*(9.0*BcurvY[1]*Phi[7]*f[16]+5.0*BcurvY[0]*Phi[6]*f[16]+5.0*BcurvY[0]*Phi[7]*f[8]+5.0*BcurvY[1]*Phi[6]*f[8]+9.0*BcurvY[1]*Phi[4]*f[6]+5.0*BcurvY[0]*Phi[2]*f[6]+5.0*BcurvY[0]*f[2]*Phi[4]+5.0*BcurvY[1]*Phi[2]*f[2])*dfac_m*dfac_v2*dfac_y*m_*q_*wv+20.0*(BcurvY[0]*f[1]+f[0]*BcurvY[1])*dfac_m*dfac_y*m2*wv+21.21320343559643*dfac_m*dfac_v*dfac_x*dfac_y*m_*(BmagInv[0]*Bmag[1]*f[9]*wm+Bmag[1]*BmagInv[1]*f[4]*wm-1.0*BmagInv[1]*Phi[6]*f[26]*q_-1.0*BmagInv[0]*Phi[6]*f[19]*q_+BmagInv[0]*Phi[5]*f[18]*q_-1.0*BmagInv[1]*Phi[2]*f[17]*q_+BmagInv[1]*Phi[5]*f[11]*q_-1.0*BmagInv[0]*Phi[2]*f[10]*q_+BmagInv[0]*Phi[1]*f[9]*q_+BmagInv[1]*Phi[1]*f[4]*q_)-60.0*(f[0]*Phi[7]+Phi[3]*f[6]+f[1]*Phi[6]+f[2]*Phi[5])*dfac_m*dfac_v2*dfac_z*q2+12.24744871391589*Bmag[1]*(BmagInv[0]*f[23]+BmagInv[1]*f[15])*dfac_v*dfac_x*dfac_y*m_))/(dfac_m*dfac_v*m_*q_); 
  out[18] += (0.025*(69.28203230275508*f[9]*dfac_z*m_*wv-3.0*(9.0*BcurvY[1]*Phi[4]*f[7]+5.0*BcurvY[0]*Phi[2]*f[7]+9.0*BcurvY[1]*f[1]*Phi[7]+5.0*BcurvY[0]*f[0]*Phi[7]+5.0*BcurvY[0]*f[1]*Phi[6]+5.0*f[0]*BcurvY[1]*Phi[6]+5.0*BcurvY[0]*f[3]*Phi[4]+5.0*BcurvY[1]*Phi[2]*f[3])*dfac_v*dfac_y*m_*wv-42.42640687119286*(Phi[6]*f[16]+Phi[7]*f[8]+Phi[3]*f[7]+f[3]*Phi[5])*dfac_v*dfac_z*q_-15.0*(BmagInv[0]*Phi[4]*f[18]+BmagInv[1]*Phi[2]*f[18]+BmagInv[1]*Phi[4]*f[11]+BmagInv[0]*Phi[2]*f[11]+BmagInv[0]*Phi[7]*f[9]+BmagInv[1]*Phi[6]*f[9]+BmagInv[1]*f[4]*Phi[7]+BmagInv[0]*f[4]*Phi[6])*dfac_x*dfac_y*m_-1.732050807568877*(9.0*BcurvY[1]*Phi[4]*f[18]+5.0*BcurvY[0]*Phi[2]*f[18]+5.0*BcurvY[0]*Phi[4]*f[11]+5.0*BcurvY[1]*Phi[2]*f[11]+9.0*BcurvY[1]*Phi[7]*f[9]+5.0*BcurvY[0]*Phi[6]*f[9]+5.0*BcurvY[0]*f[4]*Phi[7]+5.0*BcurvY[1]*f[4]*Phi[6])*dfac_y*m_))/m_; 
  out[19] += (0.0883883476483184*(1.732050807568877*dfac_m*dfac_v*dfac_y*m_*(4.0*BcurvY[1]*f[18]*m_*wv2+4.0*BcurvY[0]*f[11]*m_*wv2-1.414213562373095*BcurvY[0]*Phi[4]*f[26]*q_-1.414213562373095*BcurvY[1]*Phi[2]*f[26]*q_-1.414213562373095*BcurvY[1]*Phi[4]*f[19]*q_-1.414213562373095*BcurvY[0]*Phi[2]*f[19]*q_-1.414213562373095*BcurvY[0]*Phi[7]*f[17]*q_-1.414213562373095*BcurvY[1]*Phi[6]*f[17]*q_-1.414213562373095*BcurvY[1]*Phi[7]*f[10]*q_-1.414213562373095*BcurvY[0]*Phi[6]*f[10]*q_)+19.59591794226543*f[10]*dfac_m*dfac_v*dfac_z*m_*q_*wv-4.242640687119286*(BcurvY[0]*Phi[4]*f[16]+BcurvY[1]*Phi[2]*f[16]+BcurvY[1]*Phi[4]*f[8]+BcurvY[0]*Phi[2]*f[8]+BcurvY[0]*f[6]*Phi[7]+BcurvY[1]*f[2]*Phi[7]+BcurvY[1]*Phi[6]*f[6]+BcurvY[0]*f[2]*Phi[6])*dfac_m*dfac_v2*dfac_y*m_*q_*wv+4.0*(BcurvY[1]*f[7]+BcurvY[0]*f[3])*dfac_m*dfac_y*m2*wv+4.242640687119286*dfac_m*dfac_v*dfac_x*dfac_y*m_*(Bmag[1]*BmagInv[1]*f[18]*wm+BmagInv[0]*Bmag[1]*f[11]*wm+BmagInv[1]*Phi[4]*f[26]*q_+BmagInv[0]*Phi[4]*f[19]*q_+BmagInv[1]*Phi[1]*f[18]*q_+BmagInv[1]*Phi[7]*f[17]*q_+BmagInv[0]*Phi[1]*f[11]*q_+BmagInv[0]*Phi[7]*f[10]*q_+BmagInv[1]*Phi[5]*f[9]*q_+BmagInv[0]*f[4]*Phi[5]*q_)-12.0*(Phi[5]*f[16]+Phi[3]*f[8]+Phi[7]*f[7]+f[3]*Phi[6])*dfac_m*dfac_v2*dfac_z*q2+2.449489742783178*Bmag[1]*(BmagInv[1]*f[29]+BmagInv[0]*f[25])*dfac_v*dfac_x*dfac_y*m_))/(dfac_m*dfac_v*m_*q_); 
  out[20] += (0.125*(4.898979485566357*(BcurvY[0]*f[12]+BcurvY[1]*f[5])*dfac_m*dfac_v*dfac_y*m_*wv2+2.828427124746191*(BcurvY[0]*f[23]+BcurvY[1]*f[15])*dfac_m*dfac_y*m_*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*(BmagInv[0]*Bmag[1]*f[12]*wm+Bmag[1]*BmagInv[1]*f[5]*wm-1.0*BmagInv[1]*Phi[6]*f[27]*q_-1.0*BmagInv[0]*Phi[6]*f[22]*q_+BmagInv[0]*Phi[5]*f[21]*q_-1.0*BmagInv[1]*Phi[2]*f[20]*q_+BmagInv[1]*Phi[5]*f[14]*q_-1.0*BmagInv[0]*Phi[2]*f[13]*q_+BmagInv[0]*Phi[1]*f[12]*q_+BmagInv[1]*Phi[1]*f[5]*q_)+1.732050807568877*Bmag[1]*(BmagInv[0]*f[1]+f[0]*BmagInv[1])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[21] += 0.125*(13.85640646055102*f[12]*dfac_z*wv-3.0*(BmagInv[0]*Phi[4]*f[21]+BmagInv[1]*Phi[2]*f[21]+BmagInv[1]*Phi[4]*f[14]+BmagInv[0]*Phi[2]*f[14]+BmagInv[0]*Phi[7]*f[12]+BmagInv[1]*Phi[6]*f[12]+BmagInv[1]*f[5]*Phi[7]+BmagInv[0]*f[5]*Phi[6])*dfac_x*dfac_y); 
  out[22] += (0.125*(4.898979485566357*(BcurvY[1]*f[21]+BcurvY[0]*f[14])*dfac_m*dfac_v*dfac_y*m_*wv2+13.85640646055102*f[13]*dfac_m*dfac_v*dfac_z*q_*wv+2.828427124746191*(BcurvY[1]*f[29]+BcurvY[0]*f[25])*dfac_m*dfac_y*m_*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*(Bmag[1]*BmagInv[1]*f[21]*wm+BmagInv[0]*Bmag[1]*f[14]*wm+BmagInv[1]*Phi[4]*f[27]*q_+BmagInv[0]*Phi[4]*f[22]*q_+BmagInv[1]*Phi[1]*f[21]*q_+BmagInv[1]*Phi[7]*f[20]*q_+BmagInv[0]*Phi[1]*f[14]*q_+BmagInv[0]*Phi[7]*f[13]*q_+BmagInv[1]*Phi[5]*f[12]*q_+BmagInv[0]*Phi[5]*f[5]*q_)+1.732050807568877*Bmag[1]*(BmagInv[1]*f[7]+BmagInv[0]*f[3])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[23] += -(0.01767766952966368*(4.242640687119286*(9.0*BcurvY[1]*Phi[7]*f[21]+5.0*BcurvY[0]*Phi[6]*f[21]+5.0*BcurvY[0]*Phi[7]*f[14]+5.0*BcurvY[1]*Phi[6]*f[14]+9.0*BcurvY[1]*Phi[4]*f[12]+5.0*BcurvY[0]*Phi[2]*f[12]+5.0*BcurvY[0]*Phi[4]*f[5]+5.0*BcurvY[1]*Phi[2]*f[5])*dfac_v*dfac_y*m_*wv+60.0*(Phi[6]*f[20]+Phi[7]*f[13]+Phi[3]*f[12]+Phi[5]*f[5])*dfac_v*dfac_z*q_+21.21320343559643*(BmagInv[0]*Phi[7]*f[29]+BmagInv[1]*Phi[6]*f[29]+BmagInv[1]*Phi[7]*f[25]+BmagInv[0]*Phi[6]*f[25]+BmagInv[0]*Phi[4]*f[23]+BmagInv[1]*Phi[2]*f[23]+BmagInv[1]*Phi[4]*f[15]+BmagInv[0]*Phi[2]*f[15])*dfac_x*dfac_y*m_+2.449489742783178*(9.0*BcurvY[1]*Phi[7]*f[29]+5.0*BcurvY[0]*Phi[6]*f[29]+5.0*BcurvY[0]*Phi[7]*f[25]+5.0*BcurvY[1]*Phi[6]*f[25]+9.0*BcurvY[1]*Phi[4]*f[23]+5.0*BcurvY[0]*Phi[2]*f[23]+5.0*BcurvY[0]*Phi[4]*f[15]+5.0*BcurvY[1]*Phi[2]*f[15])*dfac_y*m_))/m_; 
  out[24] += (0.125*(1.732050807568877*dfac_m*dfac_v*dfac_y*m_*(2.828427124746191*BcurvY[1]*f[23]*m_*wv2+2.828427124746191*BcurvY[0]*f[15]*m_*wv2-1.0*BcurvY[0]*Phi[7]*f[31]*q_-1.0*BcurvY[1]*Phi[6]*f[31]*q_-1.0*BcurvY[1]*Phi[7]*f[30]*q_-1.0*BcurvY[0]*Phi[6]*f[30]*q_-1.0*BcurvY[0]*Phi[4]*f[28]*q_-1.0*BcurvY[1]*Phi[2]*f[28]*q_-1.0*BcurvY[1]*Phi[4]*f[24]*q_-1.0*BcurvY[0]*Phi[2]*f[24]*q_)-3.0*(BcurvY[0]*Phi[7]*f[27]+BcurvY[1]*Phi[6]*f[27]+BcurvY[1]*Phi[7]*f[22]+BcurvY[0]*Phi[6]*f[22]+BcurvY[0]*Phi[4]*f[20]+BcurvY[1]*Phi[2]*f[20]+BcurvY[1]*Phi[4]*f[13]+BcurvY[0]*Phi[2]*f[13])*dfac_m*dfac_v2*dfac_y*m_*q_*wv+2.828427124746191*(BcurvY[1]*f[12]+BcurvY[0]*f[5])*dfac_m*dfac_y*m2*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*m_*(Bmag[1]*BmagInv[1]*f[23]*wm+BmagInv[0]*Bmag[1]*f[15]*wm+BmagInv[1]*Phi[7]*f[31]*q_+BmagInv[0]*Phi[7]*f[30]*q_+BmagInv[1]*Phi[5]*f[29]*q_+BmagInv[1]*Phi[4]*f[28]*q_+BmagInv[0]*Phi[5]*f[25]*q_+BmagInv[0]*Phi[4]*f[24]*q_+BmagInv[1]*Phi[1]*f[23]*q_+BmagInv[0]*Phi[1]*f[15]*q_)-8.485281374238571*(Phi[5]*f[20]+Phi[3]*f[13]+Phi[7]*f[12]+f[5]*Phi[6])*dfac_m*dfac_v2*dfac_z*q2+1.732050807568877*Bmag[1]*(BmagInv[1]*f[9]+BmagInv[0]*f[4])*dfac_v*dfac_x*dfac_y*m_))/(dfac_m*dfac_v*m_*q_); 
  out[25] += (0.125*(13.85640646055102*f[15]*dfac_z*m_*wv-3.0*(BcurvY[0]*Phi[4]*f[21]+BcurvY[1]*Phi[2]*f[21]+BcurvY[1]*Phi[4]*f[14]+BcurvY[0]*Phi[2]*f[14]+BcurvY[0]*Phi[7]*f[12]+BcurvY[1]*Phi[6]*f[12]+BcurvY[1]*f[5]*Phi[7]+BcurvY[0]*f[5]*Phi[6])*dfac_v*dfac_y*m_*wv-8.485281374238571*(Phi[7]*f[27]+Phi[6]*f[22]+Phi[5]*f[21]+Phi[3]*f[14])*dfac_v*dfac_z*q_-1.732050807568877*(BcurvY[0]*Phi[4]*f[29]+BcurvY[1]*Phi[2]*f[29]+BcurvY[1]*Phi[4]*f[25]+BcurvY[0]*Phi[2]*f[25]+BcurvY[0]*Phi[7]*f[23]+BcurvY[1]*Phi[6]*f[23]+BcurvY[1]*Phi[7]*f[15]+BcurvY[0]*Phi[6]*f[15])*dfac_y*m_))/m_; 
  out[26] += (0.025*(1.732050807568877*dfac_m*dfac_v*dfac_y*m_*(14.14213562373095*BcurvY[0]*f[18]*m_*wv2+14.14213562373095*BcurvY[1]*f[11]*m_*wv2-9.0*BcurvY[1]*Phi[4]*f[26]*q_-5.0*BcurvY[0]*Phi[2]*f[26]*q_-5.0*BcurvY[0]*Phi[4]*f[19]*q_-5.0*BcurvY[1]*Phi[2]*f[19]*q_-9.0*BcurvY[1]*Phi[7]*f[17]*q_-5.0*BcurvY[0]*Phi[6]*f[17]*q_-5.0*BcurvY[0]*Phi[7]*f[10]*q_-5.0*BcurvY[1]*Phi[6]*f[10]*q_)+69.28203230275508*f[17]*dfac_m*dfac_v*dfac_z*m_*q_*wv-3.0*(9.0*BcurvY[1]*Phi[4]*f[16]+5.0*BcurvY[0]*Phi[2]*f[16]+5.0*BcurvY[0]*Phi[4]*f[8]+5.0*BcurvY[1]*Phi[2]*f[8]+9.0*BcurvY[1]*f[6]*Phi[7]+5.0*BcurvY[0]*f[2]*Phi[7]+5.0*BcurvY[0]*Phi[6]*f[6]+5.0*BcurvY[1]*f[2]*Phi[6])*dfac_m*dfac_v2*dfac_y*m_*q_*wv+14.14213562373095*(BcurvY[0]*f[7]+BcurvY[1]*f[3])*dfac_m*dfac_y*m2*wv+15.0*dfac_m*dfac_v*dfac_x*dfac_y*m_*(BmagInv[0]*Bmag[1]*f[18]*wm+Bmag[1]*BmagInv[1]*f[11]*wm-1.0*BmagInv[1]*Phi[2]*f[26]*q_-1.0*BmagInv[0]*Phi[2]*f[19]*q_+BmagInv[0]*Phi[1]*f[18]*q_-1.0*BmagInv[1]*Phi[6]*f[17]*q_+BmagInv[1]*Phi[1]*f[11]*q_-1.0*BmagInv[0]*Phi[6]*f[10]*q_+BmagInv[0]*Phi[5]*f[9]*q_+BmagInv[1]*f[4]*Phi[5]*q_)-42.42640687119286*(Phi[3]*f[16]+Phi[5]*f[8]+Phi[6]*f[7]+f[3]*Phi[7])*dfac_m*dfac_v2*dfac_z*q2+8.660254037844386*Bmag[1]*(BmagInv[0]*f[29]+BmagInv[1]*f[25])*dfac_v*dfac_x*dfac_y*m_))/(dfac_m*dfac_v*m_*q_); 
  out[27] += (0.0883883476483184*(6.928203230275509*(BcurvY[0]*f[21]+BcurvY[1]*f[14])*dfac_m*dfac_v*dfac_y*m_*wv2+19.59591794226543*f[20]*dfac_m*dfac_v*dfac_z*q_*wv+4.0*(BcurvY[0]*f[29]+BcurvY[1]*f[25])*dfac_m*dfac_y*m_*wv+4.242640687119286*dfac_m*dfac_v*dfac_x*dfac_y*(BmagInv[0]*Bmag[1]*f[21]*wm+Bmag[1]*BmagInv[1]*f[14]*wm-1.0*BmagInv[1]*Phi[2]*f[27]*q_-1.0*BmagInv[0]*Phi[2]*f[22]*q_+BmagInv[0]*Phi[1]*f[21]*q_-1.0*BmagInv[1]*Phi[6]*f[20]*q_+BmagInv[1]*Phi[1]*f[14]*q_-1.0*BmagInv[0]*Phi[6]*f[13]*q_+BmagInv[0]*Phi[5]*f[12]*q_+BmagInv[1]*Phi[5]*f[5]*q_)+2.449489742783178*Bmag[1]*(BmagInv[0]*f[7]+BmagInv[1]*f[3])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[28] += (0.01767766952966368*(1.732050807568877*dfac_m*dfac_v*dfac_y*m_*(20.0*BcurvY[0]*f[23]*m_*wv2+20.0*BcurvY[1]*f[15]*m_*wv2-12.72792206135786*BcurvY[1]*Phi[7]*f[31]*q_-7.071067811865476*BcurvY[0]*Phi[6]*f[31]*q_-7.071067811865476*BcurvY[0]*Phi[7]*f[30]*q_-7.071067811865476*BcurvY[1]*Phi[6]*f[30]*q_-12.72792206135786*BcurvY[1]*Phi[4]*f[28]*q_-7.071067811865476*BcurvY[0]*Phi[2]*f[28]*q_-7.071067811865476*BcurvY[0]*Phi[4]*f[24]*q_-7.071067811865476*BcurvY[1]*Phi[2]*f[24]*q_)-4.242640687119286*(9.0*BcurvY[1]*Phi[7]*f[27]+5.0*BcurvY[0]*Phi[6]*f[27]+5.0*BcurvY[0]*Phi[7]*f[22]+5.0*BcurvY[1]*Phi[6]*f[22]+9.0*BcurvY[1]*Phi[4]*f[20]+5.0*BcurvY[0]*Phi[2]*f[20]+5.0*BcurvY[0]*Phi[4]*f[13]+5.0*BcurvY[1]*Phi[2]*f[13])*dfac_m*dfac_v2*dfac_y*m_*q_*wv+20.0*(BcurvY[0]*f[12]+BcurvY[1]*f[5])*dfac_m*dfac_y*m2*wv+21.21320343559643*dfac_m*dfac_v*dfac_x*dfac_y*m_*(BmagInv[0]*Bmag[1]*f[23]*wm+Bmag[1]*BmagInv[1]*f[15]*wm-1.0*BmagInv[1]*Phi[6]*f[31]*q_-1.0*BmagInv[0]*Phi[6]*f[30]*q_+BmagInv[0]*Phi[5]*f[29]*q_-1.0*BmagInv[1]*Phi[2]*f[28]*q_+BmagInv[1]*Phi[5]*f[25]*q_-1.0*BmagInv[0]*Phi[2]*f[24]*q_+BmagInv[0]*Phi[1]*f[23]*q_+BmagInv[1]*Phi[1]*f[15]*q_)-60.0*(Phi[3]*f[20]+Phi[5]*f[13]+Phi[6]*f[12]+f[5]*Phi[7])*dfac_m*dfac_v2*dfac_z*q2+12.24744871391589*Bmag[1]*(BmagInv[0]*f[9]+BmagInv[1]*f[4])*dfac_v*dfac_x*dfac_y*m_))/(dfac_m*dfac_v*m_*q_); 
  out[29] += (0.025*(69.28203230275508*f[23]*dfac_z*m_*wv-3.0*(9.0*BcurvY[1]*Phi[4]*f[21]+5.0*BcurvY[0]*Phi[2]*f[21]+5.0*BcurvY[0]*Phi[4]*f[14]+5.0*BcurvY[1]*Phi[2]*f[14]+9.0*BcurvY[1]*Phi[7]*f[12]+5.0*BcurvY[0]*Phi[6]*f[12]+5.0*BcurvY[0]*f[5]*Phi[7]+5.0*BcurvY[1]*f[5]*Phi[6])*dfac_v*dfac_y*m_*wv-42.42640687119286*(Phi[6]*f[27]+Phi[7]*f[22]+Phi[3]*f[21]+Phi[5]*f[14])*dfac_v*dfac_z*q_-15.0*(BmagInv[0]*Phi[4]*f[29]+BmagInv[1]*Phi[2]*f[29]+BmagInv[1]*Phi[4]*f[25]+BmagInv[0]*Phi[2]*f[25]+BmagInv[0]*Phi[7]*f[23]+BmagInv[1]*Phi[6]*f[23]+BmagInv[1]*Phi[7]*f[15]+BmagInv[0]*Phi[6]*f[15])*dfac_x*dfac_y*m_-1.732050807568877*(9.0*BcurvY[1]*Phi[4]*f[29]+5.0*BcurvY[0]*Phi[2]*f[29]+5.0*BcurvY[0]*Phi[4]*f[25]+5.0*BcurvY[1]*Phi[2]*f[25]+9.0*BcurvY[1]*Phi[7]*f[23]+5.0*BcurvY[0]*Phi[6]*f[23]+5.0*BcurvY[0]*Phi[7]*f[15]+5.0*BcurvY[1]*Phi[6]*f[15])*dfac_y*m_))/m_; 
  out[30] += (0.0883883476483184*(1.732050807568877*dfac_m*dfac_v*dfac_y*m_*(4.0*BcurvY[1]*f[29]*m_*wv2+4.0*BcurvY[0]*f[25]*m_*wv2-1.414213562373095*BcurvY[0]*Phi[4]*f[31]*q_-1.414213562373095*BcurvY[1]*Phi[2]*f[31]*q_-1.414213562373095*BcurvY[1]*Phi[4]*f[30]*q_-1.414213562373095*BcurvY[0]*Phi[2]*f[30]*q_-1.414213562373095*BcurvY[0]*Phi[7]*f[28]*q_-1.414213562373095*BcurvY[1]*Phi[6]*f[28]*q_-1.414213562373095*BcurvY[1]*Phi[7]*f[24]*q_-1.414213562373095*BcurvY[0]*Phi[6]*f[24]*q_)+19.59591794226543*f[24]*dfac_m*dfac_v*dfac_z*m_*q_*wv-4.242640687119286*(BcurvY[0]*Phi[4]*f[27]+BcurvY[1]*Phi[2]*f[27]+BcurvY[1]*Phi[4]*f[22]+BcurvY[0]*Phi[2]*f[22]+BcurvY[0]*Phi[7]*f[20]+BcurvY[1]*Phi[6]*f[20]+BcurvY[1]*Phi[7]*f[13]+BcurvY[0]*Phi[6]*f[13])*dfac_m*dfac_v2*dfac_y*m_*q_*wv+4.0*(BcurvY[1]*f[21]+BcurvY[0]*f[14])*dfac_m*dfac_y*m2*wv+4.242640687119286*dfac_m*dfac_v*dfac_x*dfac_y*m_*(Bmag[1]*BmagInv[1]*f[29]*wm+BmagInv[0]*Bmag[1]*f[25]*wm+BmagInv[1]*Phi[4]*f[31]*q_+BmagInv[0]*Phi[4]*f[30]*q_+BmagInv[1]*Phi[1]*f[29]*q_+BmagInv[1]*Phi[7]*f[28]*q_+BmagInv[0]*Phi[1]*f[25]*q_+BmagInv[0]*Phi[7]*f[24]*q_+BmagInv[1]*Phi[5]*f[23]*q_+BmagInv[0]*Phi[5]*f[15]*q_)-12.0*(Phi[5]*f[27]+Phi[3]*f[22]+Phi[7]*f[21]+Phi[6]*f[14])*dfac_m*dfac_v2*dfac_z*q2+2.449489742783178*Bmag[1]*(BmagInv[1]*f[18]+BmagInv[0]*f[11])*dfac_v*dfac_x*dfac_y*m_))/(dfac_m*dfac_v*m_*q_); 
  out[31] += (0.025*(1.732050807568877*dfac_m*dfac_v*dfac_y*m_*(14.14213562373095*BcurvY[0]*f[29]*m_*wv2+14.14213562373095*BcurvY[1]*f[25]*m_*wv2-9.0*BcurvY[1]*Phi[4]*f[31]*q_-5.0*BcurvY[0]*Phi[2]*f[31]*q_-5.0*BcurvY[0]*Phi[4]*f[30]*q_-5.0*BcurvY[1]*Phi[2]*f[30]*q_-9.0*BcurvY[1]*Phi[7]*f[28]*q_-5.0*BcurvY[0]*Phi[6]*f[28]*q_-5.0*BcurvY[0]*Phi[7]*f[24]*q_-5.0*BcurvY[1]*Phi[6]*f[24]*q_)+69.28203230275508*f[28]*dfac_m*dfac_v*dfac_z*m_*q_*wv-3.0*(9.0*BcurvY[1]*Phi[4]*f[27]+5.0*BcurvY[0]*Phi[2]*f[27]+5.0*BcurvY[0]*Phi[4]*f[22]+5.0*BcurvY[1]*Phi[2]*f[22]+9.0*BcurvY[1]*Phi[7]*f[20]+5.0*BcurvY[0]*Phi[6]*f[20]+5.0*BcurvY[0]*Phi[7]*f[13]+5.0*BcurvY[1]*Phi[6]*f[13])*dfac_m*dfac_v2*dfac_y*m_*q_*wv+14.14213562373095*(BcurvY[0]*f[21]+BcurvY[1]*f[14])*dfac_m*dfac_y*m2*wv+15.0*dfac_m*dfac_v*dfac_x*dfac_y*m_*(BmagInv[0]*Bmag[1]*f[29]*wm+Bmag[1]*BmagInv[1]*f[25]*wm-1.0*BmagInv[1]*Phi[2]*f[31]*q_-1.0*BmagInv[0]*Phi[2]*f[30]*q_+BmagInv[0]*Phi[1]*f[29]*q_-1.0*BmagInv[1]*Phi[6]*f[28]*q_+BmagInv[1]*Phi[1]*f[25]*q_-1.0*BmagInv[0]*Phi[6]*f[24]*q_+BmagInv[0]*Phi[5]*f[23]*q_+BmagInv[1]*Phi[5]*f[15]*q_)-42.42640687119286*(Phi[3]*f[27]+Phi[5]*f[22]+Phi[6]*f[21]+Phi[7]*f[14])*dfac_m*dfac_v2*dfac_z*q2+8.660254037844386*Bmag[1]*(BmagInv[0]*f[18]+BmagInv[1]*f[11])*dfac_v*dfac_x*dfac_y*m_))/(dfac_m*dfac_v*m_*q_); 
  double cflFreq = 0.0; 
  cflFreq += fabs(-0.2165063509461096*BmagInv[0]*Phi[2]*dfac_y)*dxInv; 
  cflFreq += fabs((0.0441941738241592*(8.0*BcurvY[0]*m_*wv2+4.898979485566357*BmagInv[0]*Bmag[1]*dfac_x*wm+4.898979485566357*BmagInv[0]*Phi[1]*dfac_x*q_))/q_)*dyInv; 
  cflFreq += fabs(wv)*dzInv; 
  cflFreq += fabs((0.0441941738241592*((-4.898979485566357*BcurvY[0]*Phi[2]*dfac_y*m_*wv)-13.85640646055102*Phi[3]*dfac_z*q_))/m_)*dvInv; 
  return cflFreq; 
} 
