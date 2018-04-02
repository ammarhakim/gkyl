#include <GyrokineticModDecl.h> 
double GyrokineticVol2x2vSerP1(const double mcByq, const double *w, const double *dxv, const double *Bmag, const double *BcurvY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dyInv = 1.0/dxv[1]; 
  double dvInv = 1.0/dxv[2]; 
  double dmInv = 1.0/dxv[3]; 
  double dfac_x = 2.0*dxInv; 
  double dfac_y = 2.0*dyInv; 
  double dfac_v = 2.0*dvInv; 
  double dfac_m = 2.0*dmInv; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = 4.0*dvInv*dvInv; 
  out[1] += -1.5*(f[1]*Phi[3]+f[0]*Phi[2])*dfac_x*dfac_y*mcByq; 
  out[2] += (0.5*(1.732050807568877*(BcurvY[1]*f[1]+BcurvY[0]*f[0])*dfac_m*dfac_v*dfac_y*mcByq*wv2+(BcurvY[1]*f[6]+BcurvY[0]*f[3])*dfac_m*dfac_y*mcByq*wv+3.0*f[0]*Bmag[1]*dfac_m*dfac_v*dfac_x*dfac_y*mcByq*wm+3.0*(f[2]*Phi[3]+f[0]*Phi[1])*dfac_m*dfac_v*dfac_x*dfac_y*mcByq+1.732050807568877*Bmag[1]*f[4]*dfac_v*dfac_x*dfac_y*mcByq))/(dfac_m*dfac_v); 
  out[3] += -0.25*(3.0*(BcurvY[0]*f[1]*Phi[3]+f[0]*BcurvY[1]*Phi[3]+BcurvY[1]*f[1]*Phi[2]+BcurvY[0]*f[0]*Phi[2])*dfac_v*dfac_y*mcByq*wv+1.732050807568877*(BcurvY[0]*Phi[3]*f[6]+BcurvY[1]*Phi[2]*f[6]+BcurvY[1]*Phi[3]*f[3]+BcurvY[0]*Phi[2]*f[3])*dfac_y*mcByq); 
  out[5] += (0.5*(1.732050807568877*(BcurvY[0]*f[1]+f[0]*BcurvY[1])*dfac_m*dfac_v*dfac_y*mcByq*wv2+(BcurvY[0]*f[6]+BcurvY[1]*f[3])*dfac_m*dfac_y*mcByq*wv+3.0*Bmag[1]*f[1]*dfac_m*dfac_v*dfac_x*dfac_y*mcByq*wm-3.0*(Phi[2]*f[2]-1.0*Phi[1]*f[1])*dfac_m*dfac_v*dfac_x*dfac_y*mcByq+1.732050807568877*Bmag[1]*f[8]*dfac_v*dfac_x*dfac_y*mcByq))/(dfac_m*dfac_v); 
  out[6] += -0.05*(3.0*(9.0*BcurvY[1]*f[1]*Phi[3]+5.0*BcurvY[0]*f[0]*Phi[3]+5.0*BcurvY[0]*f[1]*Phi[2]+5.0*f[0]*BcurvY[1]*Phi[2])*dfac_v*dfac_y*mcByq*wv+30.0*(Phi[3]*f[6]+Phi[2]*f[3])*dfac_x*dfac_y*mcByq+1.732050807568877*(9.0*BcurvY[1]*Phi[3]*f[6]+5.0*BcurvY[0]*Phi[2]*f[6]+5.0*BcurvY[0]*Phi[3]*f[3]+5.0*BcurvY[1]*Phi[2]*f[3])*dfac_y*mcByq); 
  out[7] += (0.25*(3.464101615137754*(BcurvY[1]*f[6]+BcurvY[0]*f[3])*dfac_m*dfac_v*dfac_y*mcByq*wv2-3.0*(BcurvY[0]*Phi[3]*f[5]+BcurvY[1]*Phi[2]*f[5]+BcurvY[1]*f[2]*Phi[3]+BcurvY[0]*Phi[2]*f[2])*dfac_m*dfac_v2*dfac_y*mcByq*wv+2.0*(BcurvY[1]*f[1]+BcurvY[0]*f[0])*dfac_m*dfac_y*mcByq*wv+6.0*Bmag[1]*f[3]*dfac_m*dfac_v*dfac_x*dfac_y*mcByq*wm+6.0*(Phi[3]*f[7]+Phi[1]*f[3])*dfac_m*dfac_v*dfac_x*dfac_y*mcByq+3.464101615137754*Bmag[1]*f[10]*dfac_v*dfac_x*dfac_y*mcByq-1.732050807568877*(BcurvY[0]*Phi[3]*f[11]+BcurvY[1]*Phi[2]*f[11]+BcurvY[1]*Phi[3]*f[7]+BcurvY[0]*Phi[2]*f[7])*dfac_m*dfac_v*dfac_y*mcByq))/(dfac_m*dfac_v); 
  out[8] += -1.5*(Phi[3]*f[8]+Phi[2]*f[4])*dfac_x*dfac_y*mcByq; 
  out[9] += (0.5*(1.732050807568877*(BcurvY[1]*f[8]+BcurvY[0]*f[4])*dfac_m*dfac_v*dfac_y*mcByq*wv2+(BcurvY[1]*f[13]+BcurvY[0]*f[10])*dfac_m*dfac_y*mcByq*wv+3.0*Bmag[1]*f[4]*dfac_m*dfac_v*dfac_x*dfac_y*mcByq*wm+3.0*(Phi[3]*f[9]+Phi[1]*f[4])*dfac_m*dfac_v*dfac_x*dfac_y*mcByq+1.732050807568877*f[0]*Bmag[1]*dfac_v*dfac_x*dfac_y*mcByq))/(dfac_m*dfac_v); 
  out[10] += -0.25*(3.0*(BcurvY[0]*Phi[3]*f[8]+BcurvY[1]*Phi[2]*f[8]+BcurvY[1]*Phi[3]*f[4]+BcurvY[0]*Phi[2]*f[4])*dfac_v*dfac_y*mcByq*wv+1.732050807568877*(BcurvY[0]*Phi[3]*f[13]+BcurvY[1]*Phi[2]*f[13]+BcurvY[1]*Phi[3]*f[10]+BcurvY[0]*Phi[2]*f[10])*dfac_y*mcByq); 
  out[11] += (0.05*(17.32050807568877*(BcurvY[0]*f[6]+BcurvY[1]*f[3])*dfac_m*dfac_v*dfac_y*mcByq*wv2-3.0*(9.0*BcurvY[1]*Phi[3]*f[5]+5.0*BcurvY[0]*Phi[2]*f[5]+5.0*BcurvY[0]*f[2]*Phi[3]+5.0*BcurvY[1]*Phi[2]*f[2])*dfac_m*dfac_v2*dfac_y*mcByq*wv+10.0*(BcurvY[0]*f[1]+f[0]*BcurvY[1])*dfac_m*dfac_y*mcByq*wv+30.0*Bmag[1]*f[6]*dfac_m*dfac_v*dfac_x*dfac_y*mcByq*wm-30.0*(Phi[2]*f[7]-1.0*Phi[1]*f[6])*dfac_m*dfac_v*dfac_x*dfac_y*mcByq+17.32050807568877*Bmag[1]*f[13]*dfac_v*dfac_x*dfac_y*mcByq-1.732050807568877*(9.0*BcurvY[1]*Phi[3]*f[11]+5.0*BcurvY[0]*Phi[2]*f[11]+5.0*BcurvY[0]*Phi[3]*f[7]+5.0*BcurvY[1]*Phi[2]*f[7])*dfac_m*dfac_v*dfac_y*mcByq))/(dfac_m*dfac_v); 
  out[12] += (0.5*(1.732050807568877*(BcurvY[0]*f[8]+BcurvY[1]*f[4])*dfac_m*dfac_v*dfac_y*mcByq*wv2+(BcurvY[0]*f[13]+BcurvY[1]*f[10])*dfac_m*dfac_y*mcByq*wv+3.0*Bmag[1]*f[8]*dfac_m*dfac_v*dfac_x*dfac_y*mcByq*wm-3.0*(Phi[2]*f[9]-1.0*Phi[1]*f[8])*dfac_m*dfac_v*dfac_x*dfac_y*mcByq+1.732050807568877*Bmag[1]*f[1]*dfac_v*dfac_x*dfac_y*mcByq))/(dfac_m*dfac_v); 
  out[13] += -0.05*(3.0*(9.0*BcurvY[1]*Phi[3]*f[8]+5.0*BcurvY[0]*Phi[2]*f[8]+5.0*BcurvY[0]*Phi[3]*f[4]+5.0*BcurvY[1]*Phi[2]*f[4])*dfac_v*dfac_y*mcByq*wv+30.0*(Phi[3]*f[13]+Phi[2]*f[10])*dfac_x*dfac_y*mcByq+1.732050807568877*(9.0*BcurvY[1]*Phi[3]*f[13]+5.0*BcurvY[0]*Phi[2]*f[13]+5.0*BcurvY[0]*Phi[3]*f[10]+5.0*BcurvY[1]*Phi[2]*f[10])*dfac_y*mcByq); 
  out[14] += (0.25*(3.464101615137754*(BcurvY[1]*f[13]+BcurvY[0]*f[10])*dfac_m*dfac_v*dfac_y*mcByq*wv2-3.0*(BcurvY[0]*Phi[3]*f[12]+BcurvY[1]*Phi[2]*f[12]+BcurvY[1]*Phi[3]*f[9]+BcurvY[0]*Phi[2]*f[9])*dfac_m*dfac_v2*dfac_y*mcByq*wv+2.0*(BcurvY[1]*f[8]+BcurvY[0]*f[4])*dfac_m*dfac_y*mcByq*wv+6.0*Bmag[1]*f[10]*dfac_m*dfac_v*dfac_x*dfac_y*mcByq*wm+6.0*(Phi[3]*f[14]+Phi[1]*f[10])*dfac_m*dfac_v*dfac_x*dfac_y*mcByq+3.464101615137754*Bmag[1]*f[3]*dfac_v*dfac_x*dfac_y*mcByq-1.732050807568877*(BcurvY[0]*Phi[3]*f[15]+BcurvY[1]*Phi[2]*f[15]+BcurvY[1]*Phi[3]*f[14]+BcurvY[0]*Phi[2]*f[14])*dfac_m*dfac_v*dfac_y*mcByq))/(dfac_m*dfac_v); 
  out[15] += (0.05*(17.32050807568877*(BcurvY[0]*f[13]+BcurvY[1]*f[10])*dfac_m*dfac_v*dfac_y*mcByq*wv2-3.0*(9.0*BcurvY[1]*Phi[3]*f[12]+5.0*BcurvY[0]*Phi[2]*f[12]+5.0*BcurvY[0]*Phi[3]*f[9]+5.0*BcurvY[1]*Phi[2]*f[9])*dfac_m*dfac_v2*dfac_y*mcByq*wv+10.0*(BcurvY[0]*f[8]+BcurvY[1]*f[4])*dfac_m*dfac_y*mcByq*wv+30.0*Bmag[1]*f[13]*dfac_m*dfac_v*dfac_x*dfac_y*mcByq*wm-30.0*(Phi[2]*f[14]-1.0*Phi[1]*f[13])*dfac_m*dfac_v*dfac_x*dfac_y*mcByq+17.32050807568877*Bmag[1]*f[6]*dfac_v*dfac_x*dfac_y*mcByq-1.732050807568877*(9.0*BcurvY[1]*Phi[3]*f[15]+5.0*BcurvY[0]*Phi[2]*f[15]+5.0*BcurvY[0]*Phi[3]*f[14]+5.0*BcurvY[1]*Phi[2]*f[14])*dfac_m*dfac_v*dfac_y*mcByq))/(dfac_m*dfac_v); 
  double cflFreq = 0.0; 
  cflFreq += fabs(-0.8660254037844386*Phi[2]*dfac_y*mcByq)*dxInv; 
  cflFreq += fabs(0.5*BcurvY[0]*mcByq*wv2+0.4330127018922193*dfac_x*mcByq*(2.0*Bmag[1]*wm+2.0*Phi[1]))*dyInv; 
  cflFreq += fabs(-0.4330127018922193*BcurvY[0]*Phi[2]*dfac_y*mcByq*wv)*dvInv; 
  return cflFreq; 
} 
