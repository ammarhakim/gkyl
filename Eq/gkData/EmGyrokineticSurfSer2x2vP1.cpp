#include <GyrokineticModDecl.h> 
double EmGyrokineticSurf2x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*(2.0*BdriftX[0]*m_*wv2+q_*(1.732050807568877*((BmagInv[0]*Apar[2]*dfac_y-1.0*BdriftX[0]*Apar[1])*wv-1.0*BmagInv[0]*Phi[2]*dfac_y)+(Apar[0]*BdriftX[0]-3.0*BmagInv[0]*Apar[3]*dfac_y)*wv+3.0*BmagInv[0]*Phi[3]*dfac_y)))/q_; 

  double alpha[8]; 
  alpha[0] = (0.7071067811865475*(2.0*BdriftX[0]*m_*wv2+q_*(1.732050807568877*((BmagInv[0]*Apar[2]*dfac_y-1.0*BdriftX[0]*Apar[1])*wv-1.0*BmagInv[0]*Phi[2]*dfac_y)+(Apar[0]*BdriftX[0]-3.0*BmagInv[0]*Apar[3]*dfac_y)*wv+3.0*BmagInv[0]*Phi[3]*dfac_y)))/q_; 
  alpha[1] = 0.408248290463863*BdriftX[0]*(1.732050807568877*Apar[2]-3.0*Apar[3])*wv; 
  alpha[2] = (0.8164965809277261*BdriftX[0]*m_*wv)/(dfac_v*q_); 
  double alphaUp[8]; 
  alphaUp[0] = (0.7071067811865475*(2.0*BdriftX[0]*m_*wv2+q_*(1.732050807568877*((BmagInv[0]*Apar[2]*dfac_y-1.0*BdriftX[0]*Apar[1])*wv-1.0*BmagInv[0]*Phi[2]*dfac_y)+(Apar[0]*BdriftX[0]-3.0*BmagInv[0]*Apar[3]*dfac_y)*wv+3.0*BmagInv[0]*Phi[3]*dfac_y)))/q_; 
  alphaUp[1] = 0.408248290463863*BdriftX[0]*(1.732050807568877*Apar[2]-3.0*Apar[3])*wv; 
  alphaUp[2] = (0.8164965809277261*BdriftX[0]*m_*wv)/(dfac_v*q_); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[6]+alpha[1]*fl[5])+alpha[2]*fl[3]+alpha[1]*fl[2]+alpha[0]*(1.732050807568877*fl[1]+fl[0]))*dfac_x; 
  incr[1] = -0.1767766952966368*(3.0*(alpha[2]*fl[6]+alpha[1]*fl[5])+1.732050807568877*(alpha[2]*fl[3]+alpha[1]*fl[2])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[2] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[11]+fl[7])+alpha[0]*(1.732050807568877*fl[5]+fl[2])+alpha[1]*(1.732050807568877*fl[1]+fl[0]))*dfac_x; 
  incr[3] = 0.1767766952966368*(alpha[1]*(1.732050807568877*fl[11]+fl[7])+alpha[0]*(1.732050807568877*fl[6]+fl[3])+(1.732050807568877*fl[1]+fl[0])*alpha[2])*dfac_x; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[13]+alpha[1]*fl[12])+alpha[2]*fl[10]+alpha[1]*fl[9]+alpha[0]*(1.732050807568877*fl[8]+fl[4]))*dfac_x; 
  incr[5] = -0.1767766952966368*(alpha[2]*(3.0*fl[11]+1.732050807568877*fl[7])+alpha[0]*(3.0*fl[5]+1.732050807568877*fl[2])+alpha[1]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[6] = -0.1767766952966368*(alpha[1]*(3.0*fl[11]+1.732050807568877*fl[7])+alpha[0]*(3.0*fl[6]+1.732050807568877*fl[3])+(3.0*fl[1]+1.732050807568877*fl[0])*alpha[2])*dfac_x; 
  incr[7] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[11]+fl[7])+1.732050807568877*(alpha[1]*fl[6]+alpha[2]*fl[5])+alpha[1]*fl[3]+alpha[2]*fl[2])*dfac_x; 
  incr[8] = -0.1767766952966368*(3.0*(alpha[2]*fl[13]+alpha[1]*fl[12])+1.732050807568877*(alpha[2]*fl[10]+alpha[1]*fl[9])+alpha[0]*(3.0*fl[8]+1.732050807568877*fl[4]))*dfac_x; 
  incr[9] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[15]+fl[14])+alpha[0]*(1.732050807568877*fl[12]+fl[9])+alpha[1]*(1.732050807568877*fl[8]+fl[4]))*dfac_x; 
  incr[10] = 0.1767766952966368*(alpha[1]*(1.732050807568877*fl[15]+fl[14])+alpha[0]*(1.732050807568877*fl[13]+fl[10])+alpha[2]*(1.732050807568877*fl[8]+fl[4]))*dfac_x; 
  incr[11] = -0.1767766952966368*(alpha[0]*(3.0*fl[11]+1.732050807568877*fl[7])+3.0*(alpha[1]*fl[6]+alpha[2]*fl[5])+1.732050807568877*(alpha[1]*fl[3]+alpha[2]*fl[2]))*dfac_x; 
  incr[12] = -0.1767766952966368*(alpha[2]*(3.0*fl[15]+1.732050807568877*fl[14])+alpha[0]*(3.0*fl[12]+1.732050807568877*fl[9])+alpha[1]*(3.0*fl[8]+1.732050807568877*fl[4]))*dfac_x; 
  incr[13] = -0.1767766952966368*(alpha[1]*(3.0*fl[15]+1.732050807568877*fl[14])+alpha[0]*(3.0*fl[13]+1.732050807568877*fl[10])+alpha[2]*(3.0*fl[8]+1.732050807568877*fl[4]))*dfac_x; 
  incr[14] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[15]+fl[14])+1.732050807568877*(alpha[1]*fl[13]+alpha[2]*fl[12])+alpha[1]*fl[10]+alpha[2]*fl[9])*dfac_x; 
  incr[15] = -0.1767766952966368*(alpha[0]*(3.0*fl[15]+1.732050807568877*fl[14])+3.0*(alpha[1]*fl[13]+alpha[2]*fl[12])+1.732050807568877*(alpha[1]*fl[10]+alpha[2]*fl[9]))*dfac_x; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[6]+alpha[1]*fr[5])-1.0*(alpha[2]*fr[3]+alpha[1]*fr[2])+alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac_x; 
  incr[1] = 0.1767766952966368*(3.0*(alpha[2]*fr[6]+alpha[1]*fr[5])-1.732050807568877*(alpha[2]*fr[3]+alpha[1]*fr[2])+alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac_x; 
  incr[2] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[11]-1.0*fr[7])+alpha[0]*(1.732050807568877*fr[5]-1.0*fr[2])+alpha[1]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac_x; 
  incr[3] = -0.1767766952966368*(alpha[1]*(1.732050807568877*fr[11]-1.0*fr[7])+alpha[0]*(1.732050807568877*fr[6]-1.0*fr[3])+(1.732050807568877*fr[1]-1.0*fr[0])*alpha[2])*dfac_x; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[13]+alpha[1]*fr[12])-1.0*(alpha[2]*fr[10]+alpha[1]*fr[9])+alpha[0]*(1.732050807568877*fr[8]-1.0*fr[4]))*dfac_x; 
  incr[5] = 0.1767766952966368*(alpha[2]*(3.0*fr[11]-1.732050807568877*fr[7])+alpha[0]*(3.0*fr[5]-1.732050807568877*fr[2])+alpha[1]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac_x; 
  incr[6] = 0.1767766952966368*(alpha[1]*(3.0*fr[11]-1.732050807568877*fr[7])+alpha[0]*(3.0*fr[6]-1.732050807568877*fr[3])+(3.0*fr[1]-1.732050807568877*fr[0])*alpha[2])*dfac_x; 
  incr[7] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[11]-1.0*fr[7])+1.732050807568877*(alpha[1]*fr[6]+alpha[2]*fr[5])-1.0*(alpha[1]*fr[3]+alpha[2]*fr[2]))*dfac_x; 
  incr[8] = 0.1767766952966368*(3.0*(alpha[2]*fr[13]+alpha[1]*fr[12])-1.732050807568877*(alpha[2]*fr[10]+alpha[1]*fr[9])+alpha[0]*(3.0*fr[8]-1.732050807568877*fr[4]))*dfac_x; 
  incr[9] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[15]-1.0*fr[14])+alpha[0]*(1.732050807568877*fr[12]-1.0*fr[9])+alpha[1]*(1.732050807568877*fr[8]-1.0*fr[4]))*dfac_x; 
  incr[10] = -0.1767766952966368*(alpha[1]*(1.732050807568877*fr[15]-1.0*fr[14])+alpha[0]*(1.732050807568877*fr[13]-1.0*fr[10])+alpha[2]*(1.732050807568877*fr[8]-1.0*fr[4]))*dfac_x; 
  incr[11] = 0.1767766952966368*(alpha[0]*(3.0*fr[11]-1.732050807568877*fr[7])+3.0*(alpha[1]*fr[6]+alpha[2]*fr[5])-1.732050807568877*(alpha[1]*fr[3]+alpha[2]*fr[2]))*dfac_x; 
  incr[12] = 0.1767766952966368*(alpha[2]*(3.0*fr[15]-1.732050807568877*fr[14])+alpha[0]*(3.0*fr[12]-1.732050807568877*fr[9])+alpha[1]*(3.0*fr[8]-1.732050807568877*fr[4]))*dfac_x; 
  incr[13] = 0.1767766952966368*(alpha[1]*(3.0*fr[15]-1.732050807568877*fr[14])+alpha[0]*(3.0*fr[13]-1.732050807568877*fr[10])+alpha[2]*(3.0*fr[8]-1.732050807568877*fr[4]))*dfac_x; 
  incr[14] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[15]-1.0*fr[14])+1.732050807568877*(alpha[1]*fr[13]+alpha[2]*fr[12])-1.0*(alpha[1]*fr[10]+alpha[2]*fr[9]))*dfac_x; 
  incr[15] = 0.1767766952966368*(alpha[0]*(3.0*fr[15]-1.732050807568877*fr[14])+3.0*(alpha[1]*fr[13]+alpha[2]*fr[12])-1.732050807568877*(alpha[1]*fr[10]+alpha[2]*fr[9]))*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[16];
double fupwindQuad[16];
double limQuad[16];
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[0] = (-0.25*(fl[15]+fl[14]))+0.25*(fl[13]+fl[12]+fl[11]+fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5]+fl[4]+fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[0] = 0.25*fr[15]-0.25*(fr[14]+fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.25*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[1] = 0.25*fl[15]-0.25*(fl[14]+fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])-0.25*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[1] = (-0.25*(fr[15]+fr[14]))+0.25*(fr[13]+fr[12]+fr[11]+fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5]+fr[4]+fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[2] = 0.25*(fl[15]+fl[14]+fl[13])-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*fl[5]-0.25*(fl[4]+fl[3])+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[2] = (-0.25*fr[15])+0.25*fr[14]-0.25*fr[13]+0.25*(fr[12]+fr[11]+fr[10])-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*(fr[5]+fr[4]+fr[3])+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[3] = (-0.25*fl[15])+0.25*fl[14]-0.25*fl[13]+0.25*(fl[12]+fl[11]+fl[10])-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*(fl[5]+fl[4]+fl[3])+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[3] = 0.25*(fr[15]+fr[14]+fr[13])-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*fr[5]-0.25*(fr[4]+fr[3])+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[4] = 0.25*(fl[15]+fl[14])-0.25*fl[13]+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*(fl[5]+fl[4])+0.25*fl[3]-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[4] = (-0.25*fr[15])+0.25*(fr[14]+fr[13])-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*fr[5]-0.25*fr[4]+0.25*fr[3]-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[5] = (-0.25*fl[15])+0.25*(fl[14]+fl[13])-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*fl[5]-0.25*fl[4]+0.25*fl[3]-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[5] = 0.25*(fr[15]+fr[14])-0.25*fr[13]+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*(fr[5]+fr[4])+0.25*fr[3]-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[6] = (-0.25*(fl[15]+fl[14]+fl[13]+fl[12]))+0.25*fl[11]-0.25*(fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5])-0.25*fl[4]+0.25*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[6] = 0.25*fr[15]-0.25*fr[14]+0.25*(fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5]+fr[4])+0.25*(fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[7] = 0.25*fl[15]-0.25*fl[14]+0.25*(fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5]+fl[4])+0.25*(fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[7] = (-0.25*(fr[15]+fr[14]+fr[13]+fr[12]))+0.25*fr[11]-0.25*(fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5])-0.25*fr[4]+0.25*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[8] = 0.25*(fl[15]+fl[14])-0.25*(fl[13]+fl[12])+0.25*fl[11]-0.25*(fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5])+0.25*fl[4]-0.25*(fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[8] = (-0.25*fr[15])+0.25*(fr[14]+fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5]+fr[4])-0.25*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[9] = (-0.25*fl[15])+0.25*(fl[14]+fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5]+fl[4])-0.25*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[9] = 0.25*(fr[15]+fr[14])-0.25*(fr[13]+fr[12])+0.25*fr[11]-0.25*(fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5])+0.25*fr[4]-0.25*(fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[10] = (-0.25*(fl[15]+fl[14]+fl[13]))+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*(fl[5]+fl[4])-0.25*fl[3]+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[10] = 0.25*fr[15]-0.25*fr[14]+0.25*fr[13]-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*fr[5]+0.25*fr[4]-0.25*fr[3]+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[11] = 0.25*fl[15]-0.25*fl[14]+0.25*fl[13]-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*fl[5]+0.25*fl[4]-0.25*fl[3]+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[11] = (-0.25*(fr[15]+fr[14]+fr[13]))+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*(fr[5]+fr[4])-0.25*fr[3]+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[12] = (-0.25*(fl[15]+fl[14]))+0.25*fl[13]-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*fl[5]+0.25*(fl[4]+fl[3])-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[12] = 0.25*fr[15]-0.25*(fr[14]+fr[13])+0.25*(fr[12]+fr[11]+fr[10])-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*(fr[5]+fr[4]+fr[3])-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[13] = 0.25*fl[15]-0.25*(fl[14]+fl[13])+0.25*(fl[12]+fl[11]+fl[10])-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*(fl[5]+fl[4]+fl[3])-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[13] = (-0.25*(fr[15]+fr[14]))+0.25*fr[13]-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*fr[5]+0.25*(fr[4]+fr[3])-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[14] = 0.25*(fl[15]+fl[14]+fl[13]+fl[12]+fl[11]+fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5]+fl[4]+fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[14] = (-0.25*fr[15])+0.25*fr[14]-0.25*(fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5])+0.25*(fr[4]+fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[15] = (-0.25*fl[15])+0.25*fl[14]-0.25*(fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5])+0.25*(fl[4]+fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[15] = 0.25*(fr[15]+fr[14]+fr[13]+fr[12]+fr[11]+fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5]+fr[4]+fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[4] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[5] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  fupwind[6] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[7] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[8] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[9] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[10] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[11] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[12] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[13] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[14] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[15] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[6]+alpha[1]*fupwind[5])-1.0*(alpha[2]*fupwind[3]+alpha[1]*fupwind[2])+alpha[0]*(1.732050807568877*fupwind[1]-1.0*fupwind[0]))*dfac_x; 
  incr[1] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[6]+alpha[1]*fupwind[5])-1.732050807568877*(alpha[2]*fupwind[3]+alpha[1]*fupwind[2])+alpha[0]*(3.0*fupwind[1]-1.732050807568877*fupwind[0]))*dfac_x; 
  incr[2] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fupwind[11]-1.0*fupwind[7])+alpha[0]*(1.732050807568877*fupwind[5]-1.0*fupwind[2])+alpha[1]*(1.732050807568877*fupwind[1]-1.0*fupwind[0]))*dfac_x; 
  incr[3] = -0.1767766952966368*(alpha[1]*(1.732050807568877*fupwind[11]-1.0*fupwind[7])+alpha[0]*(1.732050807568877*fupwind[6]-1.0*fupwind[3])+(1.732050807568877*fupwind[1]-1.0*fupwind[0])*alpha[2])*dfac_x; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[13]+alpha[1]*fupwind[12])-1.0*(alpha[2]*fupwind[10]+alpha[1]*fupwind[9])+alpha[0]*(1.732050807568877*fupwind[8]-1.0*fupwind[4]))*dfac_x; 
  incr[5] = 0.1767766952966368*(alpha[2]*(3.0*fupwind[11]-1.732050807568877*fupwind[7])+alpha[0]*(3.0*fupwind[5]-1.732050807568877*fupwind[2])+alpha[1]*(3.0*fupwind[1]-1.732050807568877*fupwind[0]))*dfac_x; 
  incr[6] = 0.1767766952966368*(alpha[1]*(3.0*fupwind[11]-1.732050807568877*fupwind[7])+alpha[0]*(3.0*fupwind[6]-1.732050807568877*fupwind[3])+(3.0*fupwind[1]-1.732050807568877*fupwind[0])*alpha[2])*dfac_x; 
  incr[7] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fupwind[11]-1.0*fupwind[7])+1.732050807568877*(alpha[1]*fupwind[6]+alpha[2]*fupwind[5])-1.0*(alpha[1]*fupwind[3]+alpha[2]*fupwind[2]))*dfac_x; 
  incr[8] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[13]+alpha[1]*fupwind[12])-1.732050807568877*(alpha[2]*fupwind[10]+alpha[1]*fupwind[9])+alpha[0]*(3.0*fupwind[8]-1.732050807568877*fupwind[4]))*dfac_x; 
  incr[9] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fupwind[15]-1.0*fupwind[14])+alpha[0]*(1.732050807568877*fupwind[12]-1.0*fupwind[9])+alpha[1]*(1.732050807568877*fupwind[8]-1.0*fupwind[4]))*dfac_x; 
  incr[10] = -0.1767766952966368*(alpha[1]*(1.732050807568877*fupwind[15]-1.0*fupwind[14])+alpha[0]*(1.732050807568877*fupwind[13]-1.0*fupwind[10])+alpha[2]*(1.732050807568877*fupwind[8]-1.0*fupwind[4]))*dfac_x; 
  incr[11] = 0.1767766952966368*(alpha[0]*(3.0*fupwind[11]-1.732050807568877*fupwind[7])+3.0*(alpha[1]*fupwind[6]+alpha[2]*fupwind[5])-1.732050807568877*(alpha[1]*fupwind[3]+alpha[2]*fupwind[2]))*dfac_x; 
  incr[12] = 0.1767766952966368*(alpha[2]*(3.0*fupwind[15]-1.732050807568877*fupwind[14])+alpha[0]*(3.0*fupwind[12]-1.732050807568877*fupwind[9])+alpha[1]*(3.0*fupwind[8]-1.732050807568877*fupwind[4]))*dfac_x; 
  incr[13] = 0.1767766952966368*(alpha[1]*(3.0*fupwind[15]-1.732050807568877*fupwind[14])+alpha[0]*(3.0*fupwind[13]-1.732050807568877*fupwind[10])+alpha[2]*(3.0*fupwind[8]-1.732050807568877*fupwind[4]))*dfac_x; 
  incr[14] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fupwind[15]-1.0*fupwind[14])+1.732050807568877*(alpha[1]*fupwind[13]+alpha[2]*fupwind[12])-1.0*(alpha[1]*fupwind[10]+alpha[2]*fupwind[9]))*dfac_x; 
  incr[15] = 0.1767766952966368*(alpha[0]*(3.0*fupwind[15]-1.732050807568877*fupwind[14])+3.0*(alpha[1]*fupwind[13]+alpha[2]*fupwind[12])-1.732050807568877*(alpha[1]*fupwind[10]+alpha[2]*fupwind[9]))*dfac_x; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSer_Y_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*(2.0*BdriftY[0]*m_*wv2+q_*((-1.732050807568877*((BmagInv[0]*Apar[1]*dfac_x+BdriftY[0]*Apar[2])*wv-1.0*BmagInv[0]*Phi[1]*dfac_x))+(3.0*BmagInv[0]*Apar[3]*dfac_x+Apar[0]*BdriftY[0])*wv-3.0*BmagInv[0]*Phi[3]*dfac_x)))/q_; 

  double alpha[8]; 
  alpha[0] = -(0.7071067811865475*(q_*(1.732050807568877*((BmagInv[0]*Apar[1]*dfac_x+BdriftY[0]*Apar[2])*wv-1.0*BmagInv[0]*Phi[1]*dfac_x)+((-3.0*BmagInv[0]*Apar[3]*dfac_x)-1.0*Apar[0]*BdriftY[0])*wv+3.0*BmagInv[0]*Phi[3]*dfac_x)-2.0*BdriftY[0]*m_*wv2))/q_; 
  alpha[1] = 0.408248290463863*BdriftY[0]*(1.732050807568877*Apar[1]-3.0*Apar[3])*wv; 
  alpha[2] = (0.8164965809277261*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  double alphaUp[8]; 
  alphaUp[0] = -(0.7071067811865475*(q_*(1.732050807568877*((BmagInv[0]*Apar[1]*dfac_x+BdriftY[0]*Apar[2])*wv-1.0*BmagInv[0]*Phi[1]*dfac_x)+((-3.0*BmagInv[0]*Apar[3]*dfac_x)-1.0*Apar[0]*BdriftY[0])*wv+3.0*BmagInv[0]*Phi[3]*dfac_x)-2.0*BdriftY[0]*m_*wv2))/q_; 
  alphaUp[1] = 0.408248290463863*BdriftY[0]*(1.732050807568877*Apar[1]-3.0*Apar[3])*wv; 
  alphaUp[2] = (0.8164965809277261*BdriftY[0]*m_*wv)/(dfac_v*q_); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[7]+alpha[1]*fl[5])+alpha[2]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_y; 
  incr[1] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[11]+fl[6])+1.732050807568877*(alpha[0]*fl[5]+alpha[1]*fl[2])+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_y; 
  incr[2] = -0.1767766952966368*(3.0*(alpha[2]*fl[7]+alpha[1]*fl[5])+1.732050807568877*alpha[2]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_y; 
  incr[3] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[11]+alpha[0]*fl[7])+alpha[1]*fl[6]+alpha[0]*fl[3]+alpha[2]*(1.732050807568877*fl[2]+fl[0]))*dfac_y; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[14]+alpha[1]*fl[12])+alpha[2]*fl[10]+1.732050807568877*alpha[0]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4])*dfac_y; 
  incr[5] = -0.1767766952966368*(alpha[2]*(3.0*fl[11]+1.732050807568877*fl[6])+3.0*(alpha[0]*fl[5]+alpha[1]*fl[2])+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_y; 
  incr[6] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[11]+alpha[1]*fl[7])+alpha[0]*fl[6]+1.732050807568877*alpha[2]*fl[5]+alpha[1]*fl[3]+fl[1]*alpha[2])*dfac_y; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[1]*fl[11]+alpha[0]*fl[7])+1.732050807568877*(alpha[1]*fl[6]+alpha[0]*fl[3])+alpha[2]*(3.0*fl[2]+1.732050807568877*fl[0]))*dfac_y; 
  incr[8] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[15]+fl[13])+1.732050807568877*(alpha[0]*fl[12]+alpha[1]*fl[9])+alpha[0]*fl[8]+alpha[1]*fl[4])*dfac_y; 
  incr[9] = -0.1767766952966368*(3.0*(alpha[2]*fl[14]+alpha[1]*fl[12])+1.732050807568877*alpha[2]*fl[10]+3.0*alpha[0]*fl[9]+1.732050807568877*(alpha[1]*fl[8]+alpha[0]*fl[4]))*dfac_y; 
  incr[10] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14])+alpha[1]*fl[13]+alpha[0]*fl[10]+alpha[2]*(1.732050807568877*fl[9]+fl[4]))*dfac_y; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[0]*fl[11]+alpha[1]*fl[7])+1.732050807568877*alpha[0]*fl[6]+3.0*alpha[2]*fl[5]+1.732050807568877*(alpha[1]*fl[3]+fl[1]*alpha[2]))*dfac_y; 
  incr[12] = -0.1767766952966368*(alpha[2]*(3.0*fl[15]+1.732050807568877*fl[13])+3.0*(alpha[0]*fl[12]+alpha[1]*fl[9])+1.732050807568877*(alpha[0]*fl[8]+alpha[1]*fl[4]))*dfac_y; 
  incr[13] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14])+alpha[0]*fl[13]+1.732050807568877*alpha[2]*fl[12]+alpha[1]*fl[10]+alpha[2]*fl[8])*dfac_y; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14])+1.732050807568877*(alpha[1]*fl[13]+alpha[0]*fl[10])+alpha[2]*(3.0*fl[9]+1.732050807568877*fl[4]))*dfac_y; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14])+1.732050807568877*alpha[0]*fl[13]+3.0*alpha[2]*fl[12]+1.732050807568877*(alpha[1]*fl[10]+alpha[2]*fl[8]))*dfac_y; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[7]+alpha[1]*fr[5])-1.0*alpha[2]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_y; 
  incr[1] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[11]-1.0*fr[6])+1.732050807568877*(alpha[0]*fr[5]+alpha[1]*fr[2])-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_y; 
  incr[2] = 0.1767766952966368*(3.0*(alpha[2]*fr[7]+alpha[1]*fr[5])-1.732050807568877*alpha[2]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_y; 
  incr[3] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[11]+alpha[0]*fr[7])-1.0*(alpha[1]*fr[6]+alpha[0]*fr[3])+alpha[2]*(1.732050807568877*fr[2]-1.0*fr[0]))*dfac_y; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[14]+alpha[1]*fr[12])-1.0*alpha[2]*fr[10]+1.732050807568877*alpha[0]*fr[9]-1.0*(alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_y; 
  incr[5] = 0.1767766952966368*(alpha[2]*(3.0*fr[11]-1.732050807568877*fr[6])+3.0*(alpha[0]*fr[5]+alpha[1]*fr[2])-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_y; 
  incr[6] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[11]+alpha[1]*fr[7])-1.0*alpha[0]*fr[6]+1.732050807568877*alpha[2]*fr[5]-1.0*(alpha[1]*fr[3]+fr[1]*alpha[2]))*dfac_y; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[1]*fr[11]+alpha[0]*fr[7])-1.732050807568877*(alpha[1]*fr[6]+alpha[0]*fr[3])+alpha[2]*(3.0*fr[2]-1.732050807568877*fr[0]))*dfac_y; 
  incr[8] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[15]-1.0*fr[13])+1.732050807568877*(alpha[0]*fr[12]+alpha[1]*fr[9])-1.0*(alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_y; 
  incr[9] = 0.1767766952966368*(3.0*(alpha[2]*fr[14]+alpha[1]*fr[12])-1.732050807568877*alpha[2]*fr[10]+3.0*alpha[0]*fr[9]-1.732050807568877*(alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_y; 
  incr[10] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14])-1.0*(alpha[1]*fr[13]+alpha[0]*fr[10])+alpha[2]*(1.732050807568877*fr[9]-1.0*fr[4]))*dfac_y; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[0]*fr[11]+alpha[1]*fr[7])-1.732050807568877*alpha[0]*fr[6]+3.0*alpha[2]*fr[5]-1.732050807568877*(alpha[1]*fr[3]+fr[1]*alpha[2]))*dfac_y; 
  incr[12] = 0.1767766952966368*(alpha[2]*(3.0*fr[15]-1.732050807568877*fr[13])+3.0*(alpha[0]*fr[12]+alpha[1]*fr[9])-1.732050807568877*(alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_y; 
  incr[13] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14])-1.0*alpha[0]*fr[13]+1.732050807568877*alpha[2]*fr[12]-1.0*(alpha[1]*fr[10]+alpha[2]*fr[8]))*dfac_y; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14])-1.732050807568877*(alpha[1]*fr[13]+alpha[0]*fr[10])+alpha[2]*(3.0*fr[9]-1.732050807568877*fr[4]))*dfac_y; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14])-1.732050807568877*alpha[0]*fr[13]+3.0*alpha[2]*fr[12]-1.732050807568877*(alpha[1]*fr[10]+alpha[2]*fr[8]))*dfac_y; 
  }
#elif upwindType == QUAD 
double fupwind[16];
double fupwindQuad[16];
double limQuad[16];
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[0] = (-0.25*fl[15])+0.25*fl[14]-0.25*fl[13]+0.25*(fl[12]+fl[11]+fl[10])-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*(fl[5]+fl[4]+fl[3])+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[0] = 0.25*fr[15]-0.25*(fr[14]+fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.25*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[1] = 0.25*(fl[15]+fl[14]+fl[13])-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*fl[5]-0.25*(fl[4]+fl[3])+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[1] = (-0.25*(fr[15]+fr[14]))+0.25*(fr[13]+fr[12]+fr[11]+fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5]+fr[4]+fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[2] = 0.25*fl[15]-0.25*(fl[14]+fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])-0.25*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[2] = (-0.25*fr[15])+0.25*fr[14]-0.25*fr[13]+0.25*(fr[12]+fr[11]+fr[10])-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*(fr[5]+fr[4]+fr[3])+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[3] = (-0.25*(fl[15]+fl[14]))+0.25*(fl[13]+fl[12]+fl[11]+fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5]+fl[4]+fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[3] = 0.25*(fr[15]+fr[14]+fr[13])-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*fr[5]-0.25*(fr[4]+fr[3])+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[4] = 0.25*fl[15]-0.25*fl[14]+0.25*(fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5]+fl[4])+0.25*(fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[4] = (-0.25*fr[15])+0.25*(fr[14]+fr[13])-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*fr[5]-0.25*fr[4]+0.25*fr[3]-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[5] = (-0.25*(fl[15]+fl[14]+fl[13]+fl[12]))+0.25*fl[11]-0.25*(fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5])-0.25*fl[4]+0.25*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[5] = 0.25*(fr[15]+fr[14])-0.25*fr[13]+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*(fr[5]+fr[4])+0.25*fr[3]-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[6] = (-0.25*fl[15])+0.25*(fl[14]+fl[13])-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*fl[5]-0.25*fl[4]+0.25*fl[3]-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[6] = 0.25*fr[15]-0.25*fr[14]+0.25*(fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5]+fr[4])+0.25*(fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[7] = 0.25*(fl[15]+fl[14])-0.25*fl[13]+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*(fl[5]+fl[4])+0.25*fl[3]-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[7] = (-0.25*(fr[15]+fr[14]+fr[13]+fr[12]))+0.25*fr[11]-0.25*(fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5])-0.25*fr[4]+0.25*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[8] = 0.25*fl[15]-0.25*fl[14]+0.25*fl[13]-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*fl[5]+0.25*fl[4]-0.25*fl[3]+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[8] = (-0.25*fr[15])+0.25*(fr[14]+fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5]+fr[4])-0.25*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[9] = (-0.25*(fl[15]+fl[14]+fl[13]))+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*(fl[5]+fl[4])-0.25*fl[3]+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[9] = 0.25*(fr[15]+fr[14])-0.25*(fr[13]+fr[12])+0.25*fr[11]-0.25*(fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5])+0.25*fr[4]-0.25*(fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[10] = (-0.25*fl[15])+0.25*(fl[14]+fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5]+fl[4])-0.25*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[10] = 0.25*fr[15]-0.25*fr[14]+0.25*fr[13]-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*fr[5]+0.25*fr[4]-0.25*fr[3]+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[11] = 0.25*(fl[15]+fl[14])-0.25*(fl[13]+fl[12])+0.25*fl[11]-0.25*(fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5])+0.25*fl[4]-0.25*(fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[11] = (-0.25*(fr[15]+fr[14]+fr[13]))+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*(fr[5]+fr[4])-0.25*fr[3]+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[12] = (-0.25*fl[15])+0.25*fl[14]-0.25*(fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5])+0.25*(fl[4]+fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[12] = 0.25*fr[15]-0.25*(fr[14]+fr[13])+0.25*(fr[12]+fr[11]+fr[10])-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*(fr[5]+fr[4]+fr[3])-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[13] = 0.25*(fl[15]+fl[14]+fl[13]+fl[12]+fl[11]+fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5]+fl[4]+fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[13] = (-0.25*(fr[15]+fr[14]))+0.25*fr[13]-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*fr[5]+0.25*(fr[4]+fr[3])-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[14] = 0.25*fl[15]-0.25*(fl[14]+fl[13])+0.25*(fl[12]+fl[11]+fl[10])-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*(fl[5]+fl[4]+fl[3])-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[14] = (-0.25*fr[15])+0.25*fr[14]-0.25*(fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5])+0.25*(fr[4]+fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[15] = (-0.25*(fl[15]+fl[14]))+0.25*fl[13]-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*fl[5]+0.25*(fl[4]+fl[3])-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[15] = 0.25*(fr[15]+fr[14]+fr[13]+fr[12]+fr[11]+fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5]+fr[4]+fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[4] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[5] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  fupwind[6] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[7] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[8] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[9] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[10] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[11] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[12] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[13] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[14] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[15] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[7]+alpha[1]*fupwind[5])-1.0*alpha[2]*fupwind[3]+1.732050807568877*alpha[0]*fupwind[2]-1.0*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_y; 
  incr[1] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fupwind[11]-1.0*fupwind[6])+1.732050807568877*(alpha[0]*fupwind[5]+alpha[1]*fupwind[2])-1.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_y; 
  incr[2] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[7]+alpha[1]*fupwind[5])-1.732050807568877*alpha[2]*fupwind[3]+3.0*alpha[0]*fupwind[2]-1.732050807568877*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_y; 
  incr[3] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fupwind[11]+alpha[0]*fupwind[7])-1.0*(alpha[1]*fupwind[6]+alpha[0]*fupwind[3])+alpha[2]*(1.732050807568877*fupwind[2]-1.0*fupwind[0]))*dfac_y; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[14]+alpha[1]*fupwind[12])-1.0*alpha[2]*fupwind[10]+1.732050807568877*alpha[0]*fupwind[9]-1.0*(alpha[1]*fupwind[8]+alpha[0]*fupwind[4]))*dfac_y; 
  incr[5] = 0.1767766952966368*(alpha[2]*(3.0*fupwind[11]-1.732050807568877*fupwind[6])+3.0*(alpha[0]*fupwind[5]+alpha[1]*fupwind[2])-1.732050807568877*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_y; 
  incr[6] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fupwind[11]+alpha[1]*fupwind[7])-1.0*alpha[0]*fupwind[6]+1.732050807568877*alpha[2]*fupwind[5]-1.0*(alpha[1]*fupwind[3]+fupwind[1]*alpha[2]))*dfac_y; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[1]*fupwind[11]+alpha[0]*fupwind[7])-1.732050807568877*(alpha[1]*fupwind[6]+alpha[0]*fupwind[3])+alpha[2]*(3.0*fupwind[2]-1.732050807568877*fupwind[0]))*dfac_y; 
  incr[8] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fupwind[15]-1.0*fupwind[13])+1.732050807568877*(alpha[0]*fupwind[12]+alpha[1]*fupwind[9])-1.0*(alpha[0]*fupwind[8]+alpha[1]*fupwind[4]))*dfac_y; 
  incr[9] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[14]+alpha[1]*fupwind[12])-1.732050807568877*alpha[2]*fupwind[10]+3.0*alpha[0]*fupwind[9]-1.732050807568877*(alpha[1]*fupwind[8]+alpha[0]*fupwind[4]))*dfac_y; 
  incr[10] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14])-1.0*(alpha[1]*fupwind[13]+alpha[0]*fupwind[10])+alpha[2]*(1.732050807568877*fupwind[9]-1.0*fupwind[4]))*dfac_y; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[0]*fupwind[11]+alpha[1]*fupwind[7])-1.732050807568877*alpha[0]*fupwind[6]+3.0*alpha[2]*fupwind[5]-1.732050807568877*(alpha[1]*fupwind[3]+fupwind[1]*alpha[2]))*dfac_y; 
  incr[12] = 0.1767766952966368*(alpha[2]*(3.0*fupwind[15]-1.732050807568877*fupwind[13])+3.0*(alpha[0]*fupwind[12]+alpha[1]*fupwind[9])-1.732050807568877*(alpha[0]*fupwind[8]+alpha[1]*fupwind[4]))*dfac_y; 
  incr[13] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14])-1.0*alpha[0]*fupwind[13]+1.732050807568877*alpha[2]*fupwind[12]-1.0*(alpha[1]*fupwind[10]+alpha[2]*fupwind[8]))*dfac_y; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14])-1.732050807568877*(alpha[1]*fupwind[13]+alpha[0]*fupwind[10])+alpha[2]*(3.0*fupwind[9]-1.732050807568877*fupwind[4]))*dfac_y; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14])-1.732050807568877*alpha[0]*fupwind[13]+3.0*alpha[2]*fupwind[12]-1.732050807568877*(alpha[1]*fupwind[10]+alpha[2]*fupwind[8]))*dfac_y; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  double incrEmMod[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.0625*(1.732050807568877*(dfac_v*(2.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*wv+(BdriftY[0]*(Apar[1]*Phi[3]+Apar[0]*Phi[2])*dfac_y+BdriftX[0]*(Apar[2]*Phi[3]+Apar[0]*Phi[1])*dfac_x)*q_)-2.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_)+dfac_v*(3.0*BmagInv[0]*(Phi[1]*Apar[2]-1.0*Apar[1]*Phi[2])*dfac_x*dfac_y+4.0*dApardtPrev[0])*q_))/(dfac_v*m_); 

  double alpha[8]; 
  alpha[0] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(2.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*wv+(BdriftY[0]*(Apar[1]*Phi[3]+Apar[0]*Phi[2])*dfac_y+BdriftX[0]*(Apar[2]*Phi[3]+Apar[0]*Phi[1])*dfac_x)*q_)-2.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_)+BmagInv[0]*(3.0*Phi[1]*Apar[2]-3.0*Apar[1]*Phi[2])*dfac_v*dfac_x*dfac_y*q_))/(dfac_v*m_); 
  alpha[1] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(2.0*BdriftY[0]*Phi[3]*dfac_y*m_*wv+(BdriftY[0]*(Apar[0]*Phi[3]+Apar[1]*Phi[2])*dfac_y+BdriftX[0]*(Apar[3]*Phi[3]+Apar[1]*Phi[1])*dfac_x)*q_)-2.0*BdriftY[0]*Phi[3]*dfac_y*m_)+BmagInv[0]*(3.0*Phi[1]*Apar[3]-3.0*Apar[1]*Phi[3])*dfac_v*dfac_x*dfac_y*q_))/(dfac_v*m_); 
  alpha[2] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(2.0*BdriftX[0]*Phi[3]*dfac_x*m_*wv+(BdriftY[0]*(Apar[3]*Phi[3]+Apar[2]*Phi[2])*dfac_y+BdriftX[0]*(Apar[0]*Phi[3]+Phi[1]*Apar[2])*dfac_x)*q_)-2.0*BdriftX[0]*Phi[3]*dfac_x*m_)+BmagInv[0]*(3.0*Apar[2]*Phi[3]-3.0*Phi[2]*Apar[3])*dfac_v*dfac_x*dfac_y*q_))/(dfac_v*m_); 
  alpha[4] = -(0.6123724356957944*(BdriftY[0]*(Apar[2]*Phi[3]+Phi[2]*Apar[3])*dfac_y+BdriftX[0]*(Apar[1]*Phi[3]+Phi[1]*Apar[3])*dfac_x)*q_)/m_; 
  double alphaUp[8]; 
  alphaUp[0] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(2.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*wv+(BdriftY[0]*(Apar[1]*Phi[3]+Apar[0]*Phi[2])*dfac_y+BdriftX[0]*(Apar[2]*Phi[3]+Apar[0]*Phi[1])*dfac_x)*q_)-2.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_)+dfac_v*(BmagInv[0]*(3.0*Phi[1]*Apar[2]-3.0*Apar[1]*Phi[2])*dfac_x*dfac_y+4.0*dApardtPrev[0])*q_))/(dfac_v*m_); 
  alphaUp[1] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(2.0*BdriftY[0]*Phi[3]*dfac_y*m_*wv+(BdriftY[0]*(Apar[0]*Phi[3]+Apar[1]*Phi[2])*dfac_y+BdriftX[0]*(Apar[3]*Phi[3]+Apar[1]*Phi[1])*dfac_x)*q_)-2.0*BdriftY[0]*Phi[3]*dfac_y*m_)+dfac_v*(BmagInv[0]*(3.0*Phi[1]*Apar[3]-3.0*Apar[1]*Phi[3])*dfac_x*dfac_y+4.0*dApardtPrev[1])*q_))/(dfac_v*m_); 
  alphaUp[2] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(2.0*BdriftX[0]*Phi[3]*dfac_x*m_*wv+(BdriftY[0]*(Apar[3]*Phi[3]+Apar[2]*Phi[2])*dfac_y+BdriftX[0]*(Apar[0]*Phi[3]+Phi[1]*Apar[2])*dfac_x)*q_)-2.0*BdriftX[0]*Phi[3]*dfac_x*m_)+dfac_v*(BmagInv[0]*(3.0*Apar[2]*Phi[3]-3.0*Phi[2]*Apar[3])*dfac_x*dfac_y+4.0*dApardtPrev[2])*q_))/(dfac_v*m_); 
  alphaUp[4] = -(0.3535533905932737*(1.732050807568877*(BdriftY[0]*(Apar[2]*Phi[3]+Phi[2]*Apar[3])*dfac_y+BdriftX[0]*(Apar[1]*Phi[3]+Phi[1]*Apar[3])*dfac_x)+4.0*dApardtPrev[3])*q_)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  emModR[0] += 0.5*(1.732050807568877*fl[3]+fl[0])*dfac_v; 
  emModR[1] += 0.5*(1.732050807568877*fl[6]+fl[1])*dfac_v; 
  emModR[2] += 0.5*(1.732050807568877*fl[7]+fl[2])*dfac_v; 
  emModR[3] += -0.5*(3.0*fl[3]+1.732050807568877*fl[0])*dfac_v; 
  emModR[4] += 0.5*(1.732050807568877*fl[10]+fl[4])*dfac_v; 
  emModR[5] += 0.5*(1.732050807568877*fl[11]+fl[5])*dfac_v; 
  emModR[6] += -0.5*(3.0*fl[6]+1.732050807568877*fl[1])*dfac_v; 
  emModR[7] += -0.5*(3.0*fl[7]+1.732050807568877*fl[2])*dfac_v; 
  emModR[8] += 0.5*(1.732050807568877*fl[13]+fl[8])*dfac_v; 
  emModR[9] += 0.5*(1.732050807568877*fl[14]+fl[9])*dfac_v; 
  emModR[10] += -0.5*(3.0*fl[10]+1.732050807568877*fl[4])*dfac_v; 
  emModR[11] += -0.5*(3.0*fl[11]+1.732050807568877*fl[5])*dfac_v; 
  emModR[12] += 0.5*(1.732050807568877*fl[15]+fl[12])*dfac_v; 
  emModR[13] += -0.5*(3.0*fl[13]+1.732050807568877*fl[8])*dfac_v; 
  emModR[14] += -0.5*(3.0*fl[14]+1.732050807568877*fl[9])*dfac_v; 
  emModR[15] += -0.5*(3.0*fl[15]+1.732050807568877*fl[12])*dfac_v; 
  emModL[0] += -0.5*(1.732050807568877*fl[3]+fl[0])*dfac_v; 
  emModL[1] += -0.5*(1.732050807568877*fl[6]+fl[1])*dfac_v; 
  emModL[2] += -0.5*(1.732050807568877*fl[7]+fl[2])*dfac_v; 
  emModL[3] += -0.5*(3.0*fl[3]+1.732050807568877*fl[0])*dfac_v; 
  emModL[4] += -0.5*(1.732050807568877*fl[10]+fl[4])*dfac_v; 
  emModL[5] += -0.5*(1.732050807568877*fl[11]+fl[5])*dfac_v; 
  emModL[6] += -0.5*(3.0*fl[6]+1.732050807568877*fl[1])*dfac_v; 
  emModL[7] += -0.5*(3.0*fl[7]+1.732050807568877*fl[2])*dfac_v; 
  emModL[8] += -0.5*(1.732050807568877*fl[13]+fl[8])*dfac_v; 
  emModL[9] += -0.5*(1.732050807568877*fl[14]+fl[9])*dfac_v; 
  emModL[10] += -0.5*(3.0*fl[10]+1.732050807568877*fl[4])*dfac_v; 
  emModL[11] += -0.5*(3.0*fl[11]+1.732050807568877*fl[5])*dfac_v; 
  emModL[12] += -0.5*(1.732050807568877*fl[15]+fl[12])*dfac_v; 
  emModL[13] += -0.5*(3.0*fl[13]+1.732050807568877*fl[8])*dfac_v; 
  emModL[14] += -0.5*(3.0*fl[14]+1.732050807568877*fl[9])*dfac_v; 
  emModL[15] += -0.5*(3.0*fl[15]+1.732050807568877*fl[12])*dfac_v; 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[11]+alpha[2]*fl[7]+alpha[1]*fl[6])+alpha[4]*fl[5]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[11]+alpha[4]*fl[7]+alpha[0]*fl[6])+alpha[2]*fl[5]+fl[2]*alpha[4]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[11]+alpha[0]*fl[7]+alpha[4]*fl[6])+alpha[1]*fl[5]+fl[1]*alpha[4]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.1767766952966368*(3.0*(alpha[4]*fl[11]+alpha[2]*fl[7]+alpha[1]*fl[6])+1.732050807568877*alpha[4]*fl[5]+3.0*alpha[0]*fl[3]+1.732050807568877*(alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+alpha[4]*fl[12]+1.732050807568877*alpha[0]*fl[10]+alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[11]+alpha[1]*fl[7]+alpha[2]*fl[6])+alpha[0]*fl[5]+(1.732050807568877*fl[3]+fl[0])*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*(alpha[2]*fl[11]+alpha[4]*fl[7]+alpha[0]*fl[6])+1.732050807568877*(alpha[2]*fl[5]+fl[2]*alpha[4])+3.0*alpha[1]*fl[3]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[1]*fl[11]+alpha[0]*fl[7]+alpha[4]*fl[6])+1.732050807568877*(alpha[1]*fl[5]+fl[1]*alpha[4])+3.0*alpha[2]*fl[3]+1.732050807568877*(alpha[0]*fl[2]+fl[0]*alpha[2]))*dfac_v; 
  incr[8] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+alpha[2]*fl[12]+1.732050807568877*alpha[1]*fl[10]+alpha[4]*fl[9]+alpha[0]*fl[8]+alpha[1]*fl[4])*dfac_v; 
  incr[9] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+alpha[1]*fl[12]+1.732050807568877*alpha[2]*fl[10]+alpha[0]*fl[9]+alpha[4]*fl[8]+alpha[2]*fl[4])*dfac_v; 
  incr[10] = -0.1767766952966368*(3.0*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+1.732050807568877*alpha[4]*fl[12]+3.0*alpha[0]*fl[10]+1.732050807568877*(alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4]))*dfac_v; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[0]*fl[11]+alpha[1]*fl[7]+alpha[2]*fl[6])+1.732050807568877*alpha[0]*fl[5]+3.0*fl[3]*alpha[4]+1.732050807568877*(fl[0]*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2]))*dfac_v; 
  incr[12] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+alpha[0]*fl[12]+1.732050807568877*alpha[4]*fl[10]+alpha[1]*fl[9]+alpha[2]*fl[8]+alpha[4]*fl[4])*dfac_v; 
  incr[13] = -0.1767766952966368*(3.0*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+1.732050807568877*alpha[2]*fl[12]+3.0*alpha[1]*fl[10]+1.732050807568877*(alpha[4]*fl[9]+alpha[0]*fl[8]+alpha[1]*fl[4]))*dfac_v; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+1.732050807568877*alpha[1]*fl[12]+3.0*alpha[2]*fl[10]+1.732050807568877*(alpha[0]*fl[9]+alpha[4]*fl[8]+alpha[2]*fl[4]))*dfac_v; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+1.732050807568877*alpha[0]*fl[12]+3.0*alpha[4]*fl[10]+1.732050807568877*(alpha[1]*fl[9]+alpha[2]*fl[8]+alpha[4]*fl[4]))*dfac_v; 
  } else { 
  emModR[0] += -0.5*(1.732050807568877*fr[3]-1.0*fr[0])*dfac_v; 
  emModR[1] += -0.5*(1.732050807568877*fr[6]-1.0*fr[1])*dfac_v; 
  emModR[2] += -0.5*(1.732050807568877*fr[7]-1.0*fr[2])*dfac_v; 
  emModR[3] += 0.5*(3.0*fr[3]-1.732050807568877*fr[0])*dfac_v; 
  emModR[4] += -0.5*(1.732050807568877*fr[10]-1.0*fr[4])*dfac_v; 
  emModR[5] += -0.5*(1.732050807568877*fr[11]-1.0*fr[5])*dfac_v; 
  emModR[6] += 0.5*(3.0*fr[6]-1.732050807568877*fr[1])*dfac_v; 
  emModR[7] += 0.5*(3.0*fr[7]-1.732050807568877*fr[2])*dfac_v; 
  emModR[8] += -0.5*(1.732050807568877*fr[13]-1.0*fr[8])*dfac_v; 
  emModR[9] += -0.5*(1.732050807568877*fr[14]-1.0*fr[9])*dfac_v; 
  emModR[10] += 0.5*(3.0*fr[10]-1.732050807568877*fr[4])*dfac_v; 
  emModR[11] += 0.5*(3.0*fr[11]-1.732050807568877*fr[5])*dfac_v; 
  emModR[12] += -0.5*(1.732050807568877*fr[15]-1.0*fr[12])*dfac_v; 
  emModR[13] += 0.5*(3.0*fr[13]-1.732050807568877*fr[8])*dfac_v; 
  emModR[14] += 0.5*(3.0*fr[14]-1.732050807568877*fr[9])*dfac_v; 
  emModR[15] += 0.5*(3.0*fr[15]-1.732050807568877*fr[12])*dfac_v; 
  emModL[0] += 0.5*(1.732050807568877*fr[3]-1.0*fr[0])*dfac_v; 
  emModL[1] += 0.5*(1.732050807568877*fr[6]-1.0*fr[1])*dfac_v; 
  emModL[2] += 0.5*(1.732050807568877*fr[7]-1.0*fr[2])*dfac_v; 
  emModL[3] += 0.5*(3.0*fr[3]-1.732050807568877*fr[0])*dfac_v; 
  emModL[4] += 0.5*(1.732050807568877*fr[10]-1.0*fr[4])*dfac_v; 
  emModL[5] += 0.5*(1.732050807568877*fr[11]-1.0*fr[5])*dfac_v; 
  emModL[6] += 0.5*(3.0*fr[6]-1.732050807568877*fr[1])*dfac_v; 
  emModL[7] += 0.5*(3.0*fr[7]-1.732050807568877*fr[2])*dfac_v; 
  emModL[8] += 0.5*(1.732050807568877*fr[13]-1.0*fr[8])*dfac_v; 
  emModL[9] += 0.5*(1.732050807568877*fr[14]-1.0*fr[9])*dfac_v; 
  emModL[10] += 0.5*(3.0*fr[10]-1.732050807568877*fr[4])*dfac_v; 
  emModL[11] += 0.5*(3.0*fr[11]-1.732050807568877*fr[5])*dfac_v; 
  emModL[12] += 0.5*(1.732050807568877*fr[15]-1.0*fr[12])*dfac_v; 
  emModL[13] += 0.5*(3.0*fr[13]-1.732050807568877*fr[8])*dfac_v; 
  emModL[14] += 0.5*(3.0*fr[14]-1.732050807568877*fr[9])*dfac_v; 
  emModL[15] += 0.5*(3.0*fr[15]-1.732050807568877*fr[12])*dfac_v; 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[11]+alpha[2]*fr[7]+alpha[1]*fr[6])-1.0*alpha[4]*fr[5]+1.732050807568877*alpha[0]*fr[3]-1.0*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[11]+alpha[4]*fr[7]+alpha[0]*fr[6])-1.0*(alpha[2]*fr[5]+fr[2]*alpha[4])+1.732050807568877*alpha[1]*fr[3]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[11]+alpha[0]*fr[7]+alpha[4]*fr[6])-1.0*(alpha[1]*fr[5]+fr[1]*alpha[4])+1.732050807568877*alpha[2]*fr[3]-1.0*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[4]*fr[11]+alpha[2]*fr[7]+alpha[1]*fr[6])-1.732050807568877*alpha[4]*fr[5]+3.0*alpha[0]*fr[3]-1.732050807568877*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.0*alpha[4]*fr[12]+1.732050807568877*alpha[0]*fr[10]-1.0*(alpha[2]*fr[9]+alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[11]+alpha[1]*fr[7]+alpha[2]*fr[6])-1.0*alpha[0]*fr[5]+1.732050807568877*fr[3]*alpha[4]-1.0*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[2]*fr[11]+alpha[4]*fr[7]+alpha[0]*fr[6])-1.732050807568877*(alpha[2]*fr[5]+fr[2]*alpha[4])+3.0*alpha[1]*fr[3]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[1]*fr[11]+alpha[0]*fr[7]+alpha[4]*fr[6])-1.732050807568877*(alpha[1]*fr[5]+fr[1]*alpha[4])+3.0*alpha[2]*fr[3]-1.732050807568877*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.0*alpha[2]*fr[12]+1.732050807568877*alpha[1]*fr[10]-1.0*(alpha[4]*fr[9]+alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.0*alpha[1]*fr[12]+1.732050807568877*alpha[2]*fr[10]-1.0*(alpha[0]*fr[9]+alpha[4]*fr[8]+alpha[2]*fr[4]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.732050807568877*alpha[4]*fr[12]+3.0*alpha[0]*fr[10]-1.732050807568877*(alpha[2]*fr[9]+alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[0]*fr[11]+alpha[1]*fr[7]+alpha[2]*fr[6])-1.732050807568877*alpha[0]*fr[5]+3.0*fr[3]*alpha[4]-1.732050807568877*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.0*alpha[0]*fr[12]+1.732050807568877*alpha[4]*fr[10]-1.0*(alpha[1]*fr[9]+alpha[2]*fr[8]+alpha[4]*fr[4]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.732050807568877*alpha[2]*fr[12]+3.0*alpha[1]*fr[10]-1.732050807568877*(alpha[4]*fr[9]+alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.732050807568877*alpha[1]*fr[12]+3.0*alpha[2]*fr[10]-1.732050807568877*(alpha[0]*fr[9]+alpha[4]*fr[8]+alpha[2]*fr[4]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.732050807568877*alpha[0]*fr[12]+3.0*alpha[4]*fr[10]-1.732050807568877*(alpha[1]*fr[9]+alpha[2]*fr[8]+alpha[4]*fr[4]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[16];
double fupwindQuad[16];
double limQuad[16];
  if(0.3535533905932737*alphaUp[4]-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[0] = (-0.25*fl[15])+0.25*(fl[14]+fl[13])-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*fl[5]-0.25*fl[4]+0.25*fl[3]-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[0] = 0.25*fr[15]-0.25*(fr[14]+fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.25*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*(alphaUp[4]+alphaUp[2]) > 0) {
  fupwindQuad[1] = 0.25*(fl[15]+fl[14])-0.25*fl[13]+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*(fl[5]+fl[4])+0.25*fl[3]-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[1] = (-0.25*(fr[15]+fr[14]))+0.25*(fr[13]+fr[12]+fr[11]+fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5]+fr[4]+fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[4])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[2] = 0.25*fl[15]-0.25*fl[14]+0.25*(fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5]+fl[4])+0.25*(fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[2] = (-0.25*fr[15])+0.25*fr[14]-0.25*fr[13]+0.25*(fr[12]+fr[11]+fr[10])-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*(fr[5]+fr[4]+fr[3])+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[4]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[3] = (-0.25*(fl[15]+fl[14]+fl[13]+fl[12]))+0.25*fl[11]-0.25*(fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5])-0.25*fl[4]+0.25*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[3] = 0.25*(fr[15]+fr[14]+fr[13])-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*fr[5]-0.25*(fr[4]+fr[3])+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[4]-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[4] = 0.25*fl[15]-0.25*(fl[14]+fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])-0.25*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[4] = (-0.25*fr[15])+0.25*(fr[14]+fr[13])-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*fr[5]-0.25*fr[4]+0.25*fr[3]-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*(alphaUp[4]+alphaUp[2]) > 0) {
  fupwindQuad[5] = (-0.25*(fl[15]+fl[14]))+0.25*(fl[13]+fl[12]+fl[11]+fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5]+fl[4]+fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[5] = 0.25*(fr[15]+fr[14])-0.25*fr[13]+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*(fr[5]+fr[4])+0.25*fr[3]-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[4])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[6] = (-0.25*fl[15])+0.25*fl[14]-0.25*fl[13]+0.25*(fl[12]+fl[11]+fl[10])-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*(fl[5]+fl[4]+fl[3])+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[6] = 0.25*fr[15]-0.25*fr[14]+0.25*(fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5]+fr[4])+0.25*(fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[4]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[7] = 0.25*(fl[15]+fl[14]+fl[13])-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*fl[5]-0.25*(fl[4]+fl[3])+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[7] = (-0.25*(fr[15]+fr[14]+fr[13]+fr[12]))+0.25*fr[11]-0.25*(fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5])-0.25*fr[4]+0.25*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[4]-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[8] = 0.25*fl[15]-0.25*(fl[14]+fl[13])+0.25*(fl[12]+fl[11]+fl[10])-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*(fl[5]+fl[4]+fl[3])-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[8] = (-0.25*fr[15])+0.25*(fr[14]+fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5]+fr[4])-0.25*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*(alphaUp[4]+alphaUp[2]) > 0) {
  fupwindQuad[9] = (-0.25*(fl[15]+fl[14]))+0.25*fl[13]-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*fl[5]+0.25*(fl[4]+fl[3])-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[9] = 0.25*(fr[15]+fr[14])-0.25*(fr[13]+fr[12])+0.25*fr[11]-0.25*(fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5])+0.25*fr[4]-0.25*(fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[4])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[10] = (-0.25*fl[15])+0.25*fl[14]-0.25*(fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5])+0.25*(fl[4]+fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[10] = 0.25*fr[15]-0.25*fr[14]+0.25*fr[13]-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*fr[5]+0.25*fr[4]-0.25*fr[3]+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[4]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[11] = 0.25*(fl[15]+fl[14]+fl[13]+fl[12]+fl[11]+fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5]+fl[4]+fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[11] = (-0.25*(fr[15]+fr[14]+fr[13]))+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*(fr[5]+fr[4])-0.25*fr[3]+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[4]-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[12] = (-0.25*fl[15])+0.25*(fl[14]+fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5]+fl[4])-0.25*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[12] = 0.25*fr[15]-0.25*(fr[14]+fr[13])+0.25*(fr[12]+fr[11]+fr[10])-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*(fr[5]+fr[4]+fr[3])-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*(alphaUp[4]+alphaUp[2]) > 0) {
  fupwindQuad[13] = 0.25*(fl[15]+fl[14])-0.25*(fl[13]+fl[12])+0.25*fl[11]-0.25*(fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5])+0.25*fl[4]-0.25*(fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[13] = (-0.25*(fr[15]+fr[14]))+0.25*fr[13]-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*fr[5]+0.25*(fr[4]+fr[3])-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[4])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[14] = 0.25*fl[15]-0.25*fl[14]+0.25*fl[13]-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*fl[5]+0.25*fl[4]-0.25*fl[3]+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[14] = (-0.25*fr[15])+0.25*fr[14]-0.25*(fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5])+0.25*(fr[4]+fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[4]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[15] = (-0.25*(fl[15]+fl[14]+fl[13]))+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*(fl[5]+fl[4])-0.25*fl[3]+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[15] = 0.25*(fr[15]+fr[14]+fr[13]+fr[12]+fr[11]+fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5]+fr[4]+fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[4] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[5] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  fupwind[6] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[7] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[8] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[9] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[10] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[11] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[12] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[13] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[14] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[15] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  emModR[0] += -0.5*(1.732050807568877*fupwind[3]-1.0*fupwind[0])*dfac_v; 
  emModR[1] += -0.5*(1.732050807568877*fupwind[6]-1.0*fupwind[1])*dfac_v; 
  emModR[2] += -0.5*(1.732050807568877*fupwind[7]-1.0*fupwind[2])*dfac_v; 
  emModR[3] += 0.5*(3.0*fupwind[3]-1.732050807568877*fupwind[0])*dfac_v; 
  emModR[4] += -0.5*(1.732050807568877*fupwind[10]-1.0*fupwind[4])*dfac_v; 
  emModR[5] += -0.5*(1.732050807568877*fupwind[11]-1.0*fupwind[5])*dfac_v; 
  emModR[6] += 0.5*(3.0*fupwind[6]-1.732050807568877*fupwind[1])*dfac_v; 
  emModR[7] += 0.5*(3.0*fupwind[7]-1.732050807568877*fupwind[2])*dfac_v; 
  emModR[8] += -0.5*(1.732050807568877*fupwind[13]-1.0*fupwind[8])*dfac_v; 
  emModR[9] += -0.5*(1.732050807568877*fupwind[14]-1.0*fupwind[9])*dfac_v; 
  emModR[10] += 0.5*(3.0*fupwind[10]-1.732050807568877*fupwind[4])*dfac_v; 
  emModR[11] += 0.5*(3.0*fupwind[11]-1.732050807568877*fupwind[5])*dfac_v; 
  emModR[12] += -0.5*(1.732050807568877*fupwind[15]-1.0*fupwind[12])*dfac_v; 
  emModR[13] += 0.5*(3.0*fupwind[13]-1.732050807568877*fupwind[8])*dfac_v; 
  emModR[14] += 0.5*(3.0*fupwind[14]-1.732050807568877*fupwind[9])*dfac_v; 
  emModR[15] += 0.5*(3.0*fupwind[15]-1.732050807568877*fupwind[12])*dfac_v; 
  emModL[0] += 0.5*(1.732050807568877*fupwind[3]-1.0*fupwind[0])*dfac_v; 
  emModL[1] += 0.5*(1.732050807568877*fupwind[6]-1.0*fupwind[1])*dfac_v; 
  emModL[2] += 0.5*(1.732050807568877*fupwind[7]-1.0*fupwind[2])*dfac_v; 
  emModL[3] += 0.5*(3.0*fupwind[3]-1.732050807568877*fupwind[0])*dfac_v; 
  emModL[4] += 0.5*(1.732050807568877*fupwind[10]-1.0*fupwind[4])*dfac_v; 
  emModL[5] += 0.5*(1.732050807568877*fupwind[11]-1.0*fupwind[5])*dfac_v; 
  emModL[6] += 0.5*(3.0*fupwind[6]-1.732050807568877*fupwind[1])*dfac_v; 
  emModL[7] += 0.5*(3.0*fupwind[7]-1.732050807568877*fupwind[2])*dfac_v; 
  emModL[8] += 0.5*(1.732050807568877*fupwind[13]-1.0*fupwind[8])*dfac_v; 
  emModL[9] += 0.5*(1.732050807568877*fupwind[14]-1.0*fupwind[9])*dfac_v; 
  emModL[10] += 0.5*(3.0*fupwind[10]-1.732050807568877*fupwind[4])*dfac_v; 
  emModL[11] += 0.5*(3.0*fupwind[11]-1.732050807568877*fupwind[5])*dfac_v; 
  emModL[12] += 0.5*(1.732050807568877*fupwind[15]-1.0*fupwind[12])*dfac_v; 
  emModL[13] += 0.5*(3.0*fupwind[13]-1.732050807568877*fupwind[8])*dfac_v; 
  emModL[14] += 0.5*(3.0*fupwind[14]-1.732050807568877*fupwind[9])*dfac_v; 
  emModL[15] += 0.5*(3.0*fupwind[15]-1.732050807568877*fupwind[12])*dfac_v; 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fupwind[11]+alpha[2]*fupwind[7]+alpha[1]*fupwind[6])-1.0*alpha[4]*fupwind[5]+1.732050807568877*alpha[0]*fupwind[3]-1.0*(alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[11]+alpha[4]*fupwind[7]+alpha[0]*fupwind[6])-1.0*(alpha[2]*fupwind[5]+fupwind[2]*alpha[4])+1.732050807568877*alpha[1]*fupwind[3]-1.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fupwind[11]+alpha[0]*fupwind[7]+alpha[4]*fupwind[6])-1.0*(alpha[1]*fupwind[5]+fupwind[1]*alpha[4])+1.732050807568877*alpha[2]*fupwind[3]-1.0*(alpha[0]*fupwind[2]+fupwind[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[4]*fupwind[11]+alpha[2]*fupwind[7]+alpha[1]*fupwind[6])-1.732050807568877*alpha[4]*fupwind[5]+3.0*alpha[0]*fupwind[3]-1.732050807568877*(alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fupwind[15]+alpha[2]*fupwind[14]+alpha[1]*fupwind[13])-1.0*alpha[4]*fupwind[12]+1.732050807568877*alpha[0]*fupwind[10]-1.0*(alpha[2]*fupwind[9]+alpha[1]*fupwind[8]+alpha[0]*fupwind[4]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fupwind[11]+alpha[1]*fupwind[7]+alpha[2]*fupwind[6])-1.0*alpha[0]*fupwind[5]+1.732050807568877*fupwind[3]*alpha[4]-1.0*(fupwind[0]*alpha[4]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[11]+alpha[4]*fupwind[7]+alpha[0]*fupwind[6])-1.732050807568877*(alpha[2]*fupwind[5]+fupwind[2]*alpha[4])+3.0*alpha[1]*fupwind[3]-1.732050807568877*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[1]*fupwind[11]+alpha[0]*fupwind[7]+alpha[4]*fupwind[6])-1.732050807568877*(alpha[1]*fupwind[5]+fupwind[1]*alpha[4])+3.0*alpha[2]*fupwind[3]-1.732050807568877*(alpha[0]*fupwind[2]+fupwind[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[15]+alpha[4]*fupwind[14]+alpha[0]*fupwind[13])-1.0*alpha[2]*fupwind[12]+1.732050807568877*alpha[1]*fupwind[10]-1.0*(alpha[4]*fupwind[9]+alpha[0]*fupwind[8]+alpha[1]*fupwind[4]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14]+alpha[4]*fupwind[13])-1.0*alpha[1]*fupwind[12]+1.732050807568877*alpha[2]*fupwind[10]-1.0*(alpha[0]*fupwind[9]+alpha[4]*fupwind[8]+alpha[2]*fupwind[4]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[4]*fupwind[15]+alpha[2]*fupwind[14]+alpha[1]*fupwind[13])-1.732050807568877*alpha[4]*fupwind[12]+3.0*alpha[0]*fupwind[10]-1.732050807568877*(alpha[2]*fupwind[9]+alpha[1]*fupwind[8]+alpha[0]*fupwind[4]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[0]*fupwind[11]+alpha[1]*fupwind[7]+alpha[2]*fupwind[6])-1.732050807568877*alpha[0]*fupwind[5]+3.0*fupwind[3]*alpha[4]-1.732050807568877*(fupwind[0]*alpha[4]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14]+alpha[2]*fupwind[13])-1.0*alpha[0]*fupwind[12]+1.732050807568877*alpha[4]*fupwind[10]-1.0*(alpha[1]*fupwind[9]+alpha[2]*fupwind[8]+alpha[4]*fupwind[4]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[15]+alpha[4]*fupwind[14]+alpha[0]*fupwind[13])-1.732050807568877*alpha[2]*fupwind[12]+3.0*alpha[1]*fupwind[10]-1.732050807568877*(alpha[4]*fupwind[9]+alpha[0]*fupwind[8]+alpha[1]*fupwind[4]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14]+alpha[4]*fupwind[13])-1.732050807568877*alpha[1]*fupwind[12]+3.0*alpha[2]*fupwind[10]-1.732050807568877*(alpha[0]*fupwind[9]+alpha[4]*fupwind[8]+alpha[2]*fupwind[4]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14]+alpha[2]*fupwind[13])-1.732050807568877*alpha[0]*fupwind[12]+3.0*alpha[4]*fupwind[10]-1.732050807568877*(alpha[1]*fupwind[9]+alpha[2]*fupwind[8]+alpha[4]*fupwind[4]))*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSerStep2_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.0625*(1.732050807568877*(dfac_v*(2.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*wv+(BdriftY[0]*(Apar[1]*Phi[3]+Apar[0]*Phi[2])*dfac_y+BdriftX[0]*(Apar[2]*Phi[3]+Apar[0]*Phi[1])*dfac_x)*q_)-2.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_)+dfac_v*(3.0*BmagInv[0]*(Phi[1]*Apar[2]-1.0*Apar[1]*Phi[2])*dfac_x*dfac_y+4.0*dApardtPrev[0])*q_))/(dfac_v*m_); 

  double alpha[8]; 
  alpha[0] = -(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 
  alpha[2] = -(1.414213562373095*dApardt[2]*q_)/m_; 
  alpha[4] = -(1.414213562373095*dApardt[3]*q_)/m_; 
  double alphaUp[8]; 
  alphaUp[0] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(2.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*wv+(BdriftY[0]*(Apar[1]*Phi[3]+Apar[0]*Phi[2])*dfac_y+BdriftX[0]*(Apar[2]*Phi[3]+Apar[0]*Phi[1])*dfac_x)*q_)-2.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_)+dfac_v*(BmagInv[0]*(3.0*Phi[1]*Apar[2]-3.0*Apar[1]*Phi[2])*dfac_x*dfac_y+4.0*dApardtPrev[0])*q_))/(dfac_v*m_); 
  alphaUp[1] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(2.0*BdriftY[0]*Phi[3]*dfac_y*m_*wv+(BdriftY[0]*(Apar[0]*Phi[3]+Apar[1]*Phi[2])*dfac_y+BdriftX[0]*(Apar[3]*Phi[3]+Apar[1]*Phi[1])*dfac_x)*q_)-2.0*BdriftY[0]*Phi[3]*dfac_y*m_)+dfac_v*(BmagInv[0]*(3.0*Phi[1]*Apar[3]-3.0*Apar[1]*Phi[3])*dfac_x*dfac_y+4.0*dApardtPrev[1])*q_))/(dfac_v*m_); 
  alphaUp[2] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(2.0*BdriftX[0]*Phi[3]*dfac_x*m_*wv+(BdriftY[0]*(Apar[3]*Phi[3]+Apar[2]*Phi[2])*dfac_y+BdriftX[0]*(Apar[0]*Phi[3]+Phi[1]*Apar[2])*dfac_x)*q_)-2.0*BdriftX[0]*Phi[3]*dfac_x*m_)+dfac_v*(BmagInv[0]*(3.0*Apar[2]*Phi[3]-3.0*Phi[2]*Apar[3])*dfac_x*dfac_y+4.0*dApardtPrev[2])*q_))/(dfac_v*m_); 
  alphaUp[4] = -(0.3535533905932737*(1.732050807568877*(BdriftY[0]*(Apar[2]*Phi[3]+Phi[2]*Apar[3])*dfac_y+BdriftX[0]*(Apar[1]*Phi[3]+Phi[1]*Apar[3])*dfac_x)+4.0*dApardtPrev[3])*q_)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[11]+alpha[2]*fl[7]+alpha[1]*fl[6])+alpha[4]*fl[5]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[11]+alpha[4]*fl[7]+alpha[0]*fl[6])+alpha[2]*fl[5]+fl[2]*alpha[4]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[11]+alpha[0]*fl[7]+alpha[4]*fl[6])+alpha[1]*fl[5]+fl[1]*alpha[4]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.1767766952966368*(3.0*(alpha[4]*fl[11]+alpha[2]*fl[7]+alpha[1]*fl[6])+1.732050807568877*alpha[4]*fl[5]+3.0*alpha[0]*fl[3]+1.732050807568877*(alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+alpha[4]*fl[12]+1.732050807568877*alpha[0]*fl[10]+alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[11]+alpha[1]*fl[7]+alpha[2]*fl[6])+alpha[0]*fl[5]+(1.732050807568877*fl[3]+fl[0])*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*(alpha[2]*fl[11]+alpha[4]*fl[7]+alpha[0]*fl[6])+1.732050807568877*(alpha[2]*fl[5]+fl[2]*alpha[4])+3.0*alpha[1]*fl[3]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[1]*fl[11]+alpha[0]*fl[7]+alpha[4]*fl[6])+1.732050807568877*(alpha[1]*fl[5]+fl[1]*alpha[4])+3.0*alpha[2]*fl[3]+1.732050807568877*(alpha[0]*fl[2]+fl[0]*alpha[2]))*dfac_v; 
  incr[8] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+alpha[2]*fl[12]+1.732050807568877*alpha[1]*fl[10]+alpha[4]*fl[9]+alpha[0]*fl[8]+alpha[1]*fl[4])*dfac_v; 
  incr[9] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+alpha[1]*fl[12]+1.732050807568877*alpha[2]*fl[10]+alpha[0]*fl[9]+alpha[4]*fl[8]+alpha[2]*fl[4])*dfac_v; 
  incr[10] = -0.1767766952966368*(3.0*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+1.732050807568877*alpha[4]*fl[12]+3.0*alpha[0]*fl[10]+1.732050807568877*(alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4]))*dfac_v; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[0]*fl[11]+alpha[1]*fl[7]+alpha[2]*fl[6])+1.732050807568877*alpha[0]*fl[5]+3.0*fl[3]*alpha[4]+1.732050807568877*(fl[0]*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2]))*dfac_v; 
  incr[12] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+alpha[0]*fl[12]+1.732050807568877*alpha[4]*fl[10]+alpha[1]*fl[9]+alpha[2]*fl[8]+alpha[4]*fl[4])*dfac_v; 
  incr[13] = -0.1767766952966368*(3.0*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+1.732050807568877*alpha[2]*fl[12]+3.0*alpha[1]*fl[10]+1.732050807568877*(alpha[4]*fl[9]+alpha[0]*fl[8]+alpha[1]*fl[4]))*dfac_v; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+1.732050807568877*alpha[1]*fl[12]+3.0*alpha[2]*fl[10]+1.732050807568877*(alpha[0]*fl[9]+alpha[4]*fl[8]+alpha[2]*fl[4]))*dfac_v; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+1.732050807568877*alpha[0]*fl[12]+3.0*alpha[4]*fl[10]+1.732050807568877*(alpha[1]*fl[9]+alpha[2]*fl[8]+alpha[4]*fl[4]))*dfac_v; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[11]+alpha[2]*fr[7]+alpha[1]*fr[6])-1.0*alpha[4]*fr[5]+1.732050807568877*alpha[0]*fr[3]-1.0*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[11]+alpha[4]*fr[7]+alpha[0]*fr[6])-1.0*(alpha[2]*fr[5]+fr[2]*alpha[4])+1.732050807568877*alpha[1]*fr[3]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[11]+alpha[0]*fr[7]+alpha[4]*fr[6])-1.0*(alpha[1]*fr[5]+fr[1]*alpha[4])+1.732050807568877*alpha[2]*fr[3]-1.0*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[4]*fr[11]+alpha[2]*fr[7]+alpha[1]*fr[6])-1.732050807568877*alpha[4]*fr[5]+3.0*alpha[0]*fr[3]-1.732050807568877*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.0*alpha[4]*fr[12]+1.732050807568877*alpha[0]*fr[10]-1.0*(alpha[2]*fr[9]+alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[11]+alpha[1]*fr[7]+alpha[2]*fr[6])-1.0*alpha[0]*fr[5]+1.732050807568877*fr[3]*alpha[4]-1.0*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[2]*fr[11]+alpha[4]*fr[7]+alpha[0]*fr[6])-1.732050807568877*(alpha[2]*fr[5]+fr[2]*alpha[4])+3.0*alpha[1]*fr[3]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[1]*fr[11]+alpha[0]*fr[7]+alpha[4]*fr[6])-1.732050807568877*(alpha[1]*fr[5]+fr[1]*alpha[4])+3.0*alpha[2]*fr[3]-1.732050807568877*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.0*alpha[2]*fr[12]+1.732050807568877*alpha[1]*fr[10]-1.0*(alpha[4]*fr[9]+alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.0*alpha[1]*fr[12]+1.732050807568877*alpha[2]*fr[10]-1.0*(alpha[0]*fr[9]+alpha[4]*fr[8]+alpha[2]*fr[4]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.732050807568877*alpha[4]*fr[12]+3.0*alpha[0]*fr[10]-1.732050807568877*(alpha[2]*fr[9]+alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[0]*fr[11]+alpha[1]*fr[7]+alpha[2]*fr[6])-1.732050807568877*alpha[0]*fr[5]+3.0*fr[3]*alpha[4]-1.732050807568877*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.0*alpha[0]*fr[12]+1.732050807568877*alpha[4]*fr[10]-1.0*(alpha[1]*fr[9]+alpha[2]*fr[8]+alpha[4]*fr[4]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.732050807568877*alpha[2]*fr[12]+3.0*alpha[1]*fr[10]-1.732050807568877*(alpha[4]*fr[9]+alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.732050807568877*alpha[1]*fr[12]+3.0*alpha[2]*fr[10]-1.732050807568877*(alpha[0]*fr[9]+alpha[4]*fr[8]+alpha[2]*fr[4]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.732050807568877*alpha[0]*fr[12]+3.0*alpha[4]*fr[10]-1.732050807568877*(alpha[1]*fr[9]+alpha[2]*fr[8]+alpha[4]*fr[4]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[16];
double fupwindQuad[16];
double limQuad[16];
  if(0.3535533905932737*alphaUp[4]-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[0] = (-0.25*fl[15])+0.25*(fl[14]+fl[13])-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*fl[5]-0.25*fl[4]+0.25*fl[3]-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[0] = 0.25*fr[15]-0.25*(fr[14]+fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.25*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*(alphaUp[4]+alphaUp[2]) > 0) {
  fupwindQuad[1] = 0.25*(fl[15]+fl[14])-0.25*fl[13]+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*(fl[5]+fl[4])+0.25*fl[3]-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[1] = (-0.25*(fr[15]+fr[14]))+0.25*(fr[13]+fr[12]+fr[11]+fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5]+fr[4]+fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[4])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[2] = 0.25*fl[15]-0.25*fl[14]+0.25*(fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5]+fl[4])+0.25*(fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[2] = (-0.25*fr[15])+0.25*fr[14]-0.25*fr[13]+0.25*(fr[12]+fr[11]+fr[10])-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*(fr[5]+fr[4]+fr[3])+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[4]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[3] = (-0.25*(fl[15]+fl[14]+fl[13]+fl[12]))+0.25*fl[11]-0.25*(fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5])-0.25*fl[4]+0.25*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[3] = 0.25*(fr[15]+fr[14]+fr[13])-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*fr[5]-0.25*(fr[4]+fr[3])+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[4]-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[4] = 0.25*fl[15]-0.25*(fl[14]+fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])-0.25*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[4] = (-0.25*fr[15])+0.25*(fr[14]+fr[13])-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*fr[5]-0.25*fr[4]+0.25*fr[3]-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*(alphaUp[4]+alphaUp[2]) > 0) {
  fupwindQuad[5] = (-0.25*(fl[15]+fl[14]))+0.25*(fl[13]+fl[12]+fl[11]+fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5]+fl[4]+fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[5] = 0.25*(fr[15]+fr[14])-0.25*fr[13]+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*(fr[5]+fr[4])+0.25*fr[3]-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[4])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[6] = (-0.25*fl[15])+0.25*fl[14]-0.25*fl[13]+0.25*(fl[12]+fl[11]+fl[10])-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*(fl[5]+fl[4]+fl[3])+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[6] = 0.25*fr[15]-0.25*fr[14]+0.25*(fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5]+fr[4])+0.25*(fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[4]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[7] = 0.25*(fl[15]+fl[14]+fl[13])-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*fl[5]-0.25*(fl[4]+fl[3])+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[7] = (-0.25*(fr[15]+fr[14]+fr[13]+fr[12]))+0.25*fr[11]-0.25*(fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5])-0.25*fr[4]+0.25*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[4]-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[8] = 0.25*fl[15]-0.25*(fl[14]+fl[13])+0.25*(fl[12]+fl[11]+fl[10])-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*(fl[5]+fl[4]+fl[3])-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[8] = (-0.25*fr[15])+0.25*(fr[14]+fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5]+fr[4])-0.25*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*(alphaUp[4]+alphaUp[2]) > 0) {
  fupwindQuad[9] = (-0.25*(fl[15]+fl[14]))+0.25*fl[13]-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*fl[5]+0.25*(fl[4]+fl[3])-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[9] = 0.25*(fr[15]+fr[14])-0.25*(fr[13]+fr[12])+0.25*fr[11]-0.25*(fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5])+0.25*fr[4]-0.25*(fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[4])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[10] = (-0.25*fl[15])+0.25*fl[14]-0.25*(fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5])+0.25*(fl[4]+fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[10] = 0.25*fr[15]-0.25*fr[14]+0.25*fr[13]-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*fr[5]+0.25*fr[4]-0.25*fr[3]+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[4]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[11] = 0.25*(fl[15]+fl[14]+fl[13]+fl[12]+fl[11]+fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5]+fl[4]+fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[11] = (-0.25*(fr[15]+fr[14]+fr[13]))+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*(fr[5]+fr[4])-0.25*fr[3]+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[4]-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[12] = (-0.25*fl[15])+0.25*(fl[14]+fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5]+fl[4])-0.25*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[12] = 0.25*fr[15]-0.25*(fr[14]+fr[13])+0.25*(fr[12]+fr[11]+fr[10])-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*(fr[5]+fr[4]+fr[3])-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*(alphaUp[4]+alphaUp[2]) > 0) {
  fupwindQuad[13] = 0.25*(fl[15]+fl[14])-0.25*(fl[13]+fl[12])+0.25*fl[11]-0.25*(fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5])+0.25*fl[4]-0.25*(fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[13] = (-0.25*(fr[15]+fr[14]))+0.25*fr[13]-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*fr[5]+0.25*(fr[4]+fr[3])-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[4])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[14] = 0.25*fl[15]-0.25*fl[14]+0.25*fl[13]-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*fl[5]+0.25*fl[4]-0.25*fl[3]+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[14] = (-0.25*fr[15])+0.25*fr[14]-0.25*(fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5])+0.25*(fr[4]+fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[4]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[15] = (-0.25*(fl[15]+fl[14]+fl[13]))+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*(fl[5]+fl[4])-0.25*fl[3]+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[15] = 0.25*(fr[15]+fr[14]+fr[13]+fr[12]+fr[11]+fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5]+fr[4]+fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[4] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[5] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  fupwind[6] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[7] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[8] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[9] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[10] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[11] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[12] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[13] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[14] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[15] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fupwind[11]+alpha[2]*fupwind[7]+alpha[1]*fupwind[6])-1.0*alpha[4]*fupwind[5]+1.732050807568877*alpha[0]*fupwind[3]-1.0*(alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[11]+alpha[4]*fupwind[7]+alpha[0]*fupwind[6])-1.0*(alpha[2]*fupwind[5]+fupwind[2]*alpha[4])+1.732050807568877*alpha[1]*fupwind[3]-1.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fupwind[11]+alpha[0]*fupwind[7]+alpha[4]*fupwind[6])-1.0*(alpha[1]*fupwind[5]+fupwind[1]*alpha[4])+1.732050807568877*alpha[2]*fupwind[3]-1.0*(alpha[0]*fupwind[2]+fupwind[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[4]*fupwind[11]+alpha[2]*fupwind[7]+alpha[1]*fupwind[6])-1.732050807568877*alpha[4]*fupwind[5]+3.0*alpha[0]*fupwind[3]-1.732050807568877*(alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fupwind[15]+alpha[2]*fupwind[14]+alpha[1]*fupwind[13])-1.0*alpha[4]*fupwind[12]+1.732050807568877*alpha[0]*fupwind[10]-1.0*(alpha[2]*fupwind[9]+alpha[1]*fupwind[8]+alpha[0]*fupwind[4]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fupwind[11]+alpha[1]*fupwind[7]+alpha[2]*fupwind[6])-1.0*alpha[0]*fupwind[5]+1.732050807568877*fupwind[3]*alpha[4]-1.0*(fupwind[0]*alpha[4]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[11]+alpha[4]*fupwind[7]+alpha[0]*fupwind[6])-1.732050807568877*(alpha[2]*fupwind[5]+fupwind[2]*alpha[4])+3.0*alpha[1]*fupwind[3]-1.732050807568877*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[1]*fupwind[11]+alpha[0]*fupwind[7]+alpha[4]*fupwind[6])-1.732050807568877*(alpha[1]*fupwind[5]+fupwind[1]*alpha[4])+3.0*alpha[2]*fupwind[3]-1.732050807568877*(alpha[0]*fupwind[2]+fupwind[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[15]+alpha[4]*fupwind[14]+alpha[0]*fupwind[13])-1.0*alpha[2]*fupwind[12]+1.732050807568877*alpha[1]*fupwind[10]-1.0*(alpha[4]*fupwind[9]+alpha[0]*fupwind[8]+alpha[1]*fupwind[4]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14]+alpha[4]*fupwind[13])-1.0*alpha[1]*fupwind[12]+1.732050807568877*alpha[2]*fupwind[10]-1.0*(alpha[0]*fupwind[9]+alpha[4]*fupwind[8]+alpha[2]*fupwind[4]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[4]*fupwind[15]+alpha[2]*fupwind[14]+alpha[1]*fupwind[13])-1.732050807568877*alpha[4]*fupwind[12]+3.0*alpha[0]*fupwind[10]-1.732050807568877*(alpha[2]*fupwind[9]+alpha[1]*fupwind[8]+alpha[0]*fupwind[4]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[0]*fupwind[11]+alpha[1]*fupwind[7]+alpha[2]*fupwind[6])-1.732050807568877*alpha[0]*fupwind[5]+3.0*fupwind[3]*alpha[4]-1.732050807568877*(fupwind[0]*alpha[4]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14]+alpha[2]*fupwind[13])-1.0*alpha[0]*fupwind[12]+1.732050807568877*alpha[4]*fupwind[10]-1.0*(alpha[1]*fupwind[9]+alpha[2]*fupwind[8]+alpha[4]*fupwind[4]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[15]+alpha[4]*fupwind[14]+alpha[0]*fupwind[13])-1.732050807568877*alpha[2]*fupwind[12]+3.0*alpha[1]*fupwind[10]-1.732050807568877*(alpha[4]*fupwind[9]+alpha[0]*fupwind[8]+alpha[1]*fupwind[4]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14]+alpha[4]*fupwind[13])-1.732050807568877*alpha[1]*fupwind[12]+3.0*alpha[2]*fupwind[10]-1.732050807568877*(alpha[0]*fupwind[9]+alpha[4]*fupwind[8]+alpha[2]*fupwind[4]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14]+alpha[2]*fupwind[13])-1.732050807568877*alpha[0]*fupwind[12]+3.0*alpha[4]*fupwind[10]-1.732050807568877*(alpha[1]*fupwind[9]+alpha[2]*fupwind[8]+alpha[4]*fupwind[4]))*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.125*(1.732050807568877*(2.0*BdriftX[1]*m_*wv2+q_*((3.0*BmagInv[1]*Phi[3]+BmagInv[0]*Phi[2])*dfac_y-1.0*((3.0*BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_y-1.0*(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1]))*wv))-2.0*BdriftX[0]*m_*wv2+q_*((3.0*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_y-1.0*(3.0*Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0]))*wv-3.0*(BmagInv[0]*Phi[3]+BmagInv[1]*Phi[2])*dfac_y)))/q_; 

  double alpha[8]; 
  alpha[0] = -(0.7071067811865475*(1.732050807568877*(2.0*BdriftX[1]*m_*wv2+q_*((((-3.0*BmagInv[1]*Apar[3])-1.0*BmagInv[0]*Apar[2])*dfac_y+Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*wv+(3.0*BmagInv[1]*Phi[3]+BmagInv[0]*Phi[2])*dfac_y))-2.0*BdriftX[0]*m_*wv2+q_*((3.0*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_y-3.0*Apar[1]*BdriftX[1]-1.0*Apar[0]*BdriftX[0])*wv-3.0*(BmagInv[0]*Phi[3]+BmagInv[1]*Phi[2])*dfac_y)))/q_; 
  alpha[1] = 0.408248290463863*(1.732050807568877*(3.0*BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])-3.0*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2]))*wv; 
  alpha[2] = (0.2357022603955158*(3.464101615137754*BdriftX[0]-6.0*BdriftX[1])*m_*wv)/(dfac_v*q_); 
  double alphaUp[8]; 
  alphaUp[0] = -(0.7071067811865475*(1.732050807568877*(2.0*BdriftX[1]*m_*wv2+q_*((((-3.0*BmagInv[1]*Apar[3])-1.0*BmagInv[0]*Apar[2])*dfac_y+Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*wv+(3.0*BmagInv[1]*Phi[3]+BmagInv[0]*Phi[2])*dfac_y))-2.0*BdriftX[0]*m_*wv2+q_*((3.0*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_y-3.0*Apar[1]*BdriftX[1]-1.0*Apar[0]*BdriftX[0])*wv-3.0*(BmagInv[0]*Phi[3]+BmagInv[1]*Phi[2])*dfac_y)))/q_; 
  alphaUp[1] = 0.408248290463863*(1.732050807568877*(3.0*BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])-3.0*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2]))*wv; 
  alphaUp[2] = (0.2357022603955158*(3.464101615137754*BdriftX[0]-6.0*BdriftX[1])*m_*wv)/(dfac_v*q_); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[6]+alpha[1]*fl[5])+alpha[2]*fl[3]+alpha[1]*fl[2]+alpha[0]*(1.732050807568877*fl[1]+fl[0]))*dfac_x; 
  incr[1] = -0.1767766952966368*(3.0*(alpha[2]*fl[6]+alpha[1]*fl[5])+1.732050807568877*(alpha[2]*fl[3]+alpha[1]*fl[2])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[2] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[11]+fl[7])+alpha[0]*(1.732050807568877*fl[5]+fl[2])+alpha[1]*(1.732050807568877*fl[1]+fl[0]))*dfac_x; 
  incr[3] = 0.1767766952966368*(alpha[1]*(1.732050807568877*fl[11]+fl[7])+alpha[0]*(1.732050807568877*fl[6]+fl[3])+(1.732050807568877*fl[1]+fl[0])*alpha[2])*dfac_x; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[13]+alpha[1]*fl[12])+alpha[2]*fl[10]+alpha[1]*fl[9]+alpha[0]*(1.732050807568877*fl[8]+fl[4]))*dfac_x; 
  incr[5] = -0.1767766952966368*(alpha[2]*(3.0*fl[11]+1.732050807568877*fl[7])+alpha[0]*(3.0*fl[5]+1.732050807568877*fl[2])+alpha[1]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[6] = -0.1767766952966368*(alpha[1]*(3.0*fl[11]+1.732050807568877*fl[7])+alpha[0]*(3.0*fl[6]+1.732050807568877*fl[3])+(3.0*fl[1]+1.732050807568877*fl[0])*alpha[2])*dfac_x; 
  incr[7] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[11]+fl[7])+1.732050807568877*(alpha[1]*fl[6]+alpha[2]*fl[5])+alpha[1]*fl[3]+alpha[2]*fl[2])*dfac_x; 
  incr[8] = -0.1767766952966368*(3.0*(alpha[2]*fl[13]+alpha[1]*fl[12])+1.732050807568877*(alpha[2]*fl[10]+alpha[1]*fl[9])+alpha[0]*(3.0*fl[8]+1.732050807568877*fl[4]))*dfac_x; 
  incr[9] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[15]+fl[14])+alpha[0]*(1.732050807568877*fl[12]+fl[9])+alpha[1]*(1.732050807568877*fl[8]+fl[4]))*dfac_x; 
  incr[10] = 0.1767766952966368*(alpha[1]*(1.732050807568877*fl[15]+fl[14])+alpha[0]*(1.732050807568877*fl[13]+fl[10])+alpha[2]*(1.732050807568877*fl[8]+fl[4]))*dfac_x; 
  incr[11] = -0.1767766952966368*(alpha[0]*(3.0*fl[11]+1.732050807568877*fl[7])+3.0*(alpha[1]*fl[6]+alpha[2]*fl[5])+1.732050807568877*(alpha[1]*fl[3]+alpha[2]*fl[2]))*dfac_x; 
  incr[12] = -0.1767766952966368*(alpha[2]*(3.0*fl[15]+1.732050807568877*fl[14])+alpha[0]*(3.0*fl[12]+1.732050807568877*fl[9])+alpha[1]*(3.0*fl[8]+1.732050807568877*fl[4]))*dfac_x; 
  incr[13] = -0.1767766952966368*(alpha[1]*(3.0*fl[15]+1.732050807568877*fl[14])+alpha[0]*(3.0*fl[13]+1.732050807568877*fl[10])+alpha[2]*(3.0*fl[8]+1.732050807568877*fl[4]))*dfac_x; 
  incr[14] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[15]+fl[14])+1.732050807568877*(alpha[1]*fl[13]+alpha[2]*fl[12])+alpha[1]*fl[10]+alpha[2]*fl[9])*dfac_x; 
  incr[15] = -0.1767766952966368*(alpha[0]*(3.0*fl[15]+1.732050807568877*fl[14])+3.0*(alpha[1]*fl[13]+alpha[2]*fl[12])+1.732050807568877*(alpha[1]*fl[10]+alpha[2]*fl[9]))*dfac_x; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[6]+alpha[1]*fr[5])-1.0*(alpha[2]*fr[3]+alpha[1]*fr[2])+alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac_x; 
  incr[1] = 0.1767766952966368*(3.0*(alpha[2]*fr[6]+alpha[1]*fr[5])-1.732050807568877*(alpha[2]*fr[3]+alpha[1]*fr[2])+alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac_x; 
  incr[2] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[11]-1.0*fr[7])+alpha[0]*(1.732050807568877*fr[5]-1.0*fr[2])+alpha[1]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac_x; 
  incr[3] = -0.1767766952966368*(alpha[1]*(1.732050807568877*fr[11]-1.0*fr[7])+alpha[0]*(1.732050807568877*fr[6]-1.0*fr[3])+(1.732050807568877*fr[1]-1.0*fr[0])*alpha[2])*dfac_x; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[13]+alpha[1]*fr[12])-1.0*(alpha[2]*fr[10]+alpha[1]*fr[9])+alpha[0]*(1.732050807568877*fr[8]-1.0*fr[4]))*dfac_x; 
  incr[5] = 0.1767766952966368*(alpha[2]*(3.0*fr[11]-1.732050807568877*fr[7])+alpha[0]*(3.0*fr[5]-1.732050807568877*fr[2])+alpha[1]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac_x; 
  incr[6] = 0.1767766952966368*(alpha[1]*(3.0*fr[11]-1.732050807568877*fr[7])+alpha[0]*(3.0*fr[6]-1.732050807568877*fr[3])+(3.0*fr[1]-1.732050807568877*fr[0])*alpha[2])*dfac_x; 
  incr[7] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[11]-1.0*fr[7])+1.732050807568877*(alpha[1]*fr[6]+alpha[2]*fr[5])-1.0*(alpha[1]*fr[3]+alpha[2]*fr[2]))*dfac_x; 
  incr[8] = 0.1767766952966368*(3.0*(alpha[2]*fr[13]+alpha[1]*fr[12])-1.732050807568877*(alpha[2]*fr[10]+alpha[1]*fr[9])+alpha[0]*(3.0*fr[8]-1.732050807568877*fr[4]))*dfac_x; 
  incr[9] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[15]-1.0*fr[14])+alpha[0]*(1.732050807568877*fr[12]-1.0*fr[9])+alpha[1]*(1.732050807568877*fr[8]-1.0*fr[4]))*dfac_x; 
  incr[10] = -0.1767766952966368*(alpha[1]*(1.732050807568877*fr[15]-1.0*fr[14])+alpha[0]*(1.732050807568877*fr[13]-1.0*fr[10])+alpha[2]*(1.732050807568877*fr[8]-1.0*fr[4]))*dfac_x; 
  incr[11] = 0.1767766952966368*(alpha[0]*(3.0*fr[11]-1.732050807568877*fr[7])+3.0*(alpha[1]*fr[6]+alpha[2]*fr[5])-1.732050807568877*(alpha[1]*fr[3]+alpha[2]*fr[2]))*dfac_x; 
  incr[12] = 0.1767766952966368*(alpha[2]*(3.0*fr[15]-1.732050807568877*fr[14])+alpha[0]*(3.0*fr[12]-1.732050807568877*fr[9])+alpha[1]*(3.0*fr[8]-1.732050807568877*fr[4]))*dfac_x; 
  incr[13] = 0.1767766952966368*(alpha[1]*(3.0*fr[15]-1.732050807568877*fr[14])+alpha[0]*(3.0*fr[13]-1.732050807568877*fr[10])+alpha[2]*(3.0*fr[8]-1.732050807568877*fr[4]))*dfac_x; 
  incr[14] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[15]-1.0*fr[14])+1.732050807568877*(alpha[1]*fr[13]+alpha[2]*fr[12])-1.0*(alpha[1]*fr[10]+alpha[2]*fr[9]))*dfac_x; 
  incr[15] = 0.1767766952966368*(alpha[0]*(3.0*fr[15]-1.732050807568877*fr[14])+3.0*(alpha[1]*fr[13]+alpha[2]*fr[12])-1.732050807568877*(alpha[1]*fr[10]+alpha[2]*fr[9]))*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[16];
double fupwindQuad[16];
double limQuad[16];
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[0] = (-0.25*(fl[15]+fl[14]))+0.25*(fl[13]+fl[12]+fl[11]+fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5]+fl[4]+fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[0] = 0.25*fr[15]-0.25*(fr[14]+fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.25*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[1] = 0.25*fl[15]-0.25*(fl[14]+fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])-0.25*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[1] = (-0.25*(fr[15]+fr[14]))+0.25*(fr[13]+fr[12]+fr[11]+fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5]+fr[4]+fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[2] = 0.25*(fl[15]+fl[14]+fl[13])-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*fl[5]-0.25*(fl[4]+fl[3])+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[2] = (-0.25*fr[15])+0.25*fr[14]-0.25*fr[13]+0.25*(fr[12]+fr[11]+fr[10])-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*(fr[5]+fr[4]+fr[3])+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[3] = (-0.25*fl[15])+0.25*fl[14]-0.25*fl[13]+0.25*(fl[12]+fl[11]+fl[10])-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*(fl[5]+fl[4]+fl[3])+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[3] = 0.25*(fr[15]+fr[14]+fr[13])-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*fr[5]-0.25*(fr[4]+fr[3])+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[4] = 0.25*(fl[15]+fl[14])-0.25*fl[13]+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*(fl[5]+fl[4])+0.25*fl[3]-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[4] = (-0.25*fr[15])+0.25*(fr[14]+fr[13])-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*fr[5]-0.25*fr[4]+0.25*fr[3]-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[5] = (-0.25*fl[15])+0.25*(fl[14]+fl[13])-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*fl[5]-0.25*fl[4]+0.25*fl[3]-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[5] = 0.25*(fr[15]+fr[14])-0.25*fr[13]+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*(fr[5]+fr[4])+0.25*fr[3]-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[6] = (-0.25*(fl[15]+fl[14]+fl[13]+fl[12]))+0.25*fl[11]-0.25*(fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5])-0.25*fl[4]+0.25*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[6] = 0.25*fr[15]-0.25*fr[14]+0.25*(fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5]+fr[4])+0.25*(fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[7] = 0.25*fl[15]-0.25*fl[14]+0.25*(fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5]+fl[4])+0.25*(fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[7] = (-0.25*(fr[15]+fr[14]+fr[13]+fr[12]))+0.25*fr[11]-0.25*(fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5])-0.25*fr[4]+0.25*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[8] = 0.25*(fl[15]+fl[14])-0.25*(fl[13]+fl[12])+0.25*fl[11]-0.25*(fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5])+0.25*fl[4]-0.25*(fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[8] = (-0.25*fr[15])+0.25*(fr[14]+fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5]+fr[4])-0.25*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*alphaUp[0]-0.3535533905932737*(alphaUp[2]+alphaUp[1]) > 0) {
  fupwindQuad[9] = (-0.25*fl[15])+0.25*(fl[14]+fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5]+fl[4])-0.25*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[9] = 0.25*(fr[15]+fr[14])-0.25*(fr[13]+fr[12])+0.25*fr[11]-0.25*(fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5])+0.25*fr[4]-0.25*(fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[10] = (-0.25*(fl[15]+fl[14]+fl[13]))+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*(fl[5]+fl[4])-0.25*fl[3]+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[10] = 0.25*fr[15]-0.25*fr[14]+0.25*fr[13]-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*fr[5]+0.25*fr[4]-0.25*fr[3]+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*alphaUp[2] > 0) {
  fupwindQuad[11] = 0.25*fl[15]-0.25*fl[14]+0.25*fl[13]-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*fl[5]+0.25*fl[4]-0.25*fl[3]+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[11] = (-0.25*(fr[15]+fr[14]+fr[13]))+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*(fr[5]+fr[4])-0.25*fr[3]+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[12] = (-0.25*(fl[15]+fl[14]))+0.25*fl[13]-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*fl[5]+0.25*(fl[4]+fl[3])-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[12] = 0.25*fr[15]-0.25*(fr[14]+fr[13])+0.25*(fr[12]+fr[11]+fr[10])-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*(fr[5]+fr[4]+fr[3])-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[13] = 0.25*fl[15]-0.25*(fl[14]+fl[13])+0.25*(fl[12]+fl[11]+fl[10])-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*(fl[5]+fl[4]+fl[3])-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[13] = (-0.25*(fr[15]+fr[14]))+0.25*fr[13]-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*fr[5]+0.25*(fr[4]+fr[3])-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[14] = 0.25*(fl[15]+fl[14]+fl[13]+fl[12]+fl[11]+fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5]+fl[4]+fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[14] = (-0.25*fr[15])+0.25*fr[14]-0.25*(fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5])+0.25*(fr[4]+fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[15] = (-0.25*fl[15])+0.25*fl[14]-0.25*(fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5])+0.25*(fl[4]+fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[15] = 0.25*(fr[15]+fr[14]+fr[13]+fr[12]+fr[11]+fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5]+fr[4]+fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[4] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[5] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  fupwind[6] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[7] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[8] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[9] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[10] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[11] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[12] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[13] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[14] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[15] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[6]+alpha[1]*fupwind[5])-1.0*(alpha[2]*fupwind[3]+alpha[1]*fupwind[2])+alpha[0]*(1.732050807568877*fupwind[1]-1.0*fupwind[0]))*dfac_x; 
  incr[1] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[6]+alpha[1]*fupwind[5])-1.732050807568877*(alpha[2]*fupwind[3]+alpha[1]*fupwind[2])+alpha[0]*(3.0*fupwind[1]-1.732050807568877*fupwind[0]))*dfac_x; 
  incr[2] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fupwind[11]-1.0*fupwind[7])+alpha[0]*(1.732050807568877*fupwind[5]-1.0*fupwind[2])+alpha[1]*(1.732050807568877*fupwind[1]-1.0*fupwind[0]))*dfac_x; 
  incr[3] = -0.1767766952966368*(alpha[1]*(1.732050807568877*fupwind[11]-1.0*fupwind[7])+alpha[0]*(1.732050807568877*fupwind[6]-1.0*fupwind[3])+(1.732050807568877*fupwind[1]-1.0*fupwind[0])*alpha[2])*dfac_x; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[13]+alpha[1]*fupwind[12])-1.0*(alpha[2]*fupwind[10]+alpha[1]*fupwind[9])+alpha[0]*(1.732050807568877*fupwind[8]-1.0*fupwind[4]))*dfac_x; 
  incr[5] = 0.1767766952966368*(alpha[2]*(3.0*fupwind[11]-1.732050807568877*fupwind[7])+alpha[0]*(3.0*fupwind[5]-1.732050807568877*fupwind[2])+alpha[1]*(3.0*fupwind[1]-1.732050807568877*fupwind[0]))*dfac_x; 
  incr[6] = 0.1767766952966368*(alpha[1]*(3.0*fupwind[11]-1.732050807568877*fupwind[7])+alpha[0]*(3.0*fupwind[6]-1.732050807568877*fupwind[3])+(3.0*fupwind[1]-1.732050807568877*fupwind[0])*alpha[2])*dfac_x; 
  incr[7] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fupwind[11]-1.0*fupwind[7])+1.732050807568877*(alpha[1]*fupwind[6]+alpha[2]*fupwind[5])-1.0*(alpha[1]*fupwind[3]+alpha[2]*fupwind[2]))*dfac_x; 
  incr[8] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[13]+alpha[1]*fupwind[12])-1.732050807568877*(alpha[2]*fupwind[10]+alpha[1]*fupwind[9])+alpha[0]*(3.0*fupwind[8]-1.732050807568877*fupwind[4]))*dfac_x; 
  incr[9] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fupwind[15]-1.0*fupwind[14])+alpha[0]*(1.732050807568877*fupwind[12]-1.0*fupwind[9])+alpha[1]*(1.732050807568877*fupwind[8]-1.0*fupwind[4]))*dfac_x; 
  incr[10] = -0.1767766952966368*(alpha[1]*(1.732050807568877*fupwind[15]-1.0*fupwind[14])+alpha[0]*(1.732050807568877*fupwind[13]-1.0*fupwind[10])+alpha[2]*(1.732050807568877*fupwind[8]-1.0*fupwind[4]))*dfac_x; 
  incr[11] = 0.1767766952966368*(alpha[0]*(3.0*fupwind[11]-1.732050807568877*fupwind[7])+3.0*(alpha[1]*fupwind[6]+alpha[2]*fupwind[5])-1.732050807568877*(alpha[1]*fupwind[3]+alpha[2]*fupwind[2]))*dfac_x; 
  incr[12] = 0.1767766952966368*(alpha[2]*(3.0*fupwind[15]-1.732050807568877*fupwind[14])+alpha[0]*(3.0*fupwind[12]-1.732050807568877*fupwind[9])+alpha[1]*(3.0*fupwind[8]-1.732050807568877*fupwind[4]))*dfac_x; 
  incr[13] = 0.1767766952966368*(alpha[1]*(3.0*fupwind[15]-1.732050807568877*fupwind[14])+alpha[0]*(3.0*fupwind[13]-1.732050807568877*fupwind[10])+alpha[2]*(3.0*fupwind[8]-1.732050807568877*fupwind[4]))*dfac_x; 
  incr[14] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fupwind[15]-1.0*fupwind[14])+1.732050807568877*(alpha[1]*fupwind[13]+alpha[2]*fupwind[12])-1.0*(alpha[1]*fupwind[10]+alpha[2]*fupwind[9]))*dfac_x; 
  incr[15] = 0.1767766952966368*(alpha[0]*(3.0*fupwind[15]-1.732050807568877*fupwind[14])+3.0*(alpha[1]*fupwind[13]+alpha[2]*fupwind[12])-1.732050807568877*(alpha[1]*fupwind[10]+alpha[2]*fupwind[9]))*dfac_x; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSer_Y_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*(2.0*BdriftY[0]*m_*wv2-1.732050807568877*((BmagInv[0]*Apar[1]*dfac_x+BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2])*q_*wv-1.0*BmagInv[0]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))+q_*((3.0*BmagInv[0]*Apar[3]*dfac_x+Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*wv-3.0*BmagInv[0]*Phi[3]*dfac_x)))/q_; 

  double alpha[8]; 
  alpha[0] = -(0.7071067811865475*((-2.0*BdriftY[0]*m_*wv2)+1.732050807568877*((BmagInv[0]*Apar[1]*dfac_x+BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2])*q_*wv-1.0*BmagInv[0]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))+q_*(((-3.0*BmagInv[0]*Apar[3]*dfac_x)-1.0*(Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0]))*wv+3.0*BmagInv[0]*Phi[3]*dfac_x)))/q_; 
  alpha[1] = -(0.7071067811865475*((-2.0*BdriftY[1]*m_*wv2)+1.732050807568877*((Apar[1]*BmagInv[1]*dfac_x+BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2])*q_*wv-1.0*BmagInv[1]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))+q_*(((-3.0*BmagInv[1]*Apar[3]*dfac_x)-1.0*(Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1]))*wv+3.0*BmagInv[1]*Phi[3]*dfac_x)))/q_; 
  alpha[2] = (0.8164965809277261*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  alpha[3] = (0.7071067811865475*BmagInv[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alpha[4] = (0.8164965809277261*BdriftY[1]*m_*wv)/(dfac_v*q_); 
  alpha[5] = (0.7071067811865475*Bmag[1]*BmagInv[1]*dfac_x)/(dfac_m*q_); 
  double alphaUp[8]; 
  alphaUp[0] = -(0.7071067811865475*((-2.0*BdriftY[0]*m_*wv2)+1.732050807568877*((BmagInv[0]*Apar[1]*dfac_x+BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2])*q_*wv-1.0*BmagInv[0]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))+q_*(((-3.0*BmagInv[0]*Apar[3]*dfac_x)-1.0*(Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0]))*wv+3.0*BmagInv[0]*Phi[3]*dfac_x)))/q_; 
  alphaUp[1] = -(0.7071067811865475*((-2.0*BdriftY[1]*m_*wv2)+1.732050807568877*((Apar[1]*BmagInv[1]*dfac_x+BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2])*q_*wv-1.0*BmagInv[1]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))+q_*(((-3.0*BmagInv[1]*Apar[3]*dfac_x)-1.0*(Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1]))*wv+3.0*BmagInv[1]*Phi[3]*dfac_x)))/q_; 
  alphaUp[2] = (0.8164965809277261*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  alphaUp[3] = (0.7071067811865475*BmagInv[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alphaUp[4] = (0.8164965809277261*BdriftY[1]*m_*wv)/(dfac_v*q_); 
  alphaUp[5] = (0.7071067811865475*Bmag[1]*BmagInv[1]*dfac_x)/(dfac_m*q_); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[5]*fl[12]+alpha[4]*fl[11]+alpha[3]*fl[9])+alpha[5]*fl[8]+1.732050807568877*alpha[2]*fl[7]+alpha[4]*fl[6]+1.732050807568877*alpha[1]*fl[5]+alpha[3]*fl[4]+alpha[2]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_y; 
  incr[1] = 0.1767766952966368*(1.732050807568877*(alpha[3]*fl[12]+alpha[2]*fl[11]+alpha[5]*fl[9])+alpha[3]*fl[8]+1.732050807568877*alpha[4]*fl[7]+alpha[2]*fl[6]+1.732050807568877*alpha[0]*fl[5]+fl[4]*alpha[5]+fl[3]*alpha[4]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_y; 
  incr[2] = -0.1767766952966368*(3.0*(alpha[5]*fl[12]+alpha[4]*fl[11]+alpha[3]*fl[9])+1.732050807568877*alpha[5]*fl[8]+3.0*alpha[2]*fl[7]+1.732050807568877*alpha[4]*fl[6]+3.0*alpha[1]*fl[5]+1.732050807568877*(alpha[3]*fl[4]+alpha[2]*fl[3])+3.0*alpha[0]*fl[2]+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_y; 
  incr[3] = 0.1767766952966368*(1.732050807568877*(alpha[5]*fl[15]+alpha[3]*fl[14])+alpha[5]*fl[13]+1.732050807568877*alpha[1]*fl[11]+alpha[3]*fl[10]+1.732050807568877*alpha[0]*fl[7]+alpha[1]*fl[6]+alpha[4]*(1.732050807568877*fl[5]+fl[1])+alpha[0]*fl[3]+alpha[2]*(1.732050807568877*fl[2]+fl[0]))*dfac_y; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[15]+alpha[2]*fl[14])+alpha[4]*fl[13]+1.732050807568877*alpha[1]*fl[12]+alpha[2]*fl[10]+1.732050807568877*alpha[0]*fl[9]+alpha[1]*fl[8]+alpha[5]*(1.732050807568877*fl[5]+fl[1])+alpha[0]*fl[4]+(1.732050807568877*fl[2]+fl[0])*alpha[3])*dfac_y; 
  incr[5] = -0.1767766952966368*(3.0*(alpha[3]*fl[12]+alpha[2]*fl[11]+alpha[5]*fl[9])+1.732050807568877*alpha[3]*fl[8]+3.0*alpha[4]*fl[7]+1.732050807568877*alpha[2]*fl[6]+3.0*alpha[0]*fl[5]+1.732050807568877*(fl[4]*alpha[5]+fl[3]*alpha[4])+3.0*alpha[1]*fl[2]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_y; 
  incr[6] = 0.1767766952966368*(1.732050807568877*(alpha[3]*fl[15]+alpha[5]*fl[14])+alpha[3]*fl[13]+1.732050807568877*alpha[0]*fl[11]+alpha[5]*fl[10]+1.732050807568877*alpha[1]*fl[7]+alpha[0]*fl[6]+1.732050807568877*alpha[2]*fl[5]+(1.732050807568877*fl[2]+fl[0])*alpha[4]+alpha[1]*fl[3]+fl[1]*alpha[2])*dfac_y; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[5]*fl[15]+alpha[3]*fl[14])+1.732050807568877*alpha[5]*fl[13]+3.0*alpha[1]*fl[11]+1.732050807568877*alpha[3]*fl[10]+3.0*alpha[0]*fl[7]+1.732050807568877*alpha[1]*fl[6]+3.0*alpha[4]*fl[5]+1.732050807568877*(fl[1]*alpha[4]+alpha[0]*fl[3])+alpha[2]*(3.0*fl[2]+1.732050807568877*fl[0]))*dfac_y; 
  incr[8] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[15]+alpha[4]*fl[14])+alpha[2]*fl[13]+1.732050807568877*alpha[0]*fl[12]+alpha[4]*fl[10]+1.732050807568877*alpha[1]*fl[9]+alpha[0]*fl[8]+1.732050807568877*alpha[3]*fl[5]+(1.732050807568877*fl[2]+fl[0])*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3])*dfac_y; 
  incr[9] = -0.1767766952966368*(3.0*(alpha[4]*fl[15]+alpha[2]*fl[14])+1.732050807568877*alpha[4]*fl[13]+3.0*alpha[1]*fl[12]+1.732050807568877*alpha[2]*fl[10]+3.0*alpha[0]*fl[9]+1.732050807568877*alpha[1]*fl[8]+3.0*alpha[5]*fl[5]+1.732050807568877*(fl[1]*alpha[5]+alpha[0]*fl[4])+(3.0*fl[2]+1.732050807568877*fl[0])*alpha[3])*dfac_y; 
  incr[10] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14])+alpha[1]*fl[13]+1.732050807568877*(alpha[4]*fl[12]+alpha[5]*fl[11])+alpha[0]*fl[10]+1.732050807568877*alpha[2]*fl[9]+alpha[4]*fl[8]+1.732050807568877*alpha[3]*fl[7]+alpha[5]*fl[6]+alpha[2]*fl[4]+alpha[3]*fl[3])*dfac_y; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[3]*fl[15]+alpha[5]*fl[14])+1.732050807568877*alpha[3]*fl[13]+3.0*alpha[0]*fl[11]+1.732050807568877*alpha[5]*fl[10]+3.0*alpha[1]*fl[7]+1.732050807568877*alpha[0]*fl[6]+3.0*(alpha[2]*fl[5]+fl[2]*alpha[4])+1.732050807568877*(fl[0]*alpha[4]+alpha[1]*fl[3]+fl[1]*alpha[2]))*dfac_y; 
  incr[12] = -0.1767766952966368*(3.0*(alpha[2]*fl[15]+alpha[4]*fl[14])+1.732050807568877*alpha[2]*fl[13]+3.0*alpha[0]*fl[12]+1.732050807568877*alpha[4]*fl[10]+3.0*alpha[1]*fl[9]+1.732050807568877*alpha[0]*fl[8]+3.0*(alpha[3]*fl[5]+fl[2]*alpha[5])+1.732050807568877*(fl[0]*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3]))*dfac_y; 
  incr[13] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14])+alpha[0]*fl[13]+1.732050807568877*(alpha[2]*fl[12]+alpha[3]*fl[11])+alpha[1]*fl[10]+1.732050807568877*alpha[4]*fl[9]+alpha[2]*fl[8]+1.732050807568877*alpha[5]*fl[7]+alpha[3]*fl[6]+fl[3]*alpha[5]+alpha[4]*fl[4])*dfac_y; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14])+1.732050807568877*alpha[1]*fl[13]+3.0*(alpha[4]*fl[12]+alpha[5]*fl[11])+1.732050807568877*alpha[0]*fl[10]+3.0*alpha[2]*fl[9]+1.732050807568877*alpha[4]*fl[8]+3.0*alpha[3]*fl[7]+1.732050807568877*(alpha[5]*fl[6]+alpha[2]*fl[4]+alpha[3]*fl[3]))*dfac_y; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14])+1.732050807568877*alpha[0]*fl[13]+3.0*(alpha[2]*fl[12]+alpha[3]*fl[11])+1.732050807568877*alpha[1]*fl[10]+3.0*alpha[4]*fl[9]+1.732050807568877*alpha[2]*fl[8]+3.0*alpha[5]*fl[7]+1.732050807568877*(alpha[3]*fl[6]+fl[3]*alpha[5]+alpha[4]*fl[4]))*dfac_y; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fr[12]+alpha[4]*fr[11]+alpha[3]*fr[9])-1.0*alpha[5]*fr[8]+1.732050807568877*alpha[2]*fr[7]-1.0*alpha[4]*fr[6]+1.732050807568877*alpha[1]*fr[5]-1.0*(alpha[3]*fr[4]+alpha[2]*fr[3])+1.732050807568877*alpha[0]*fr[2]-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_y; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fr[12]+alpha[2]*fr[11]+alpha[5]*fr[9])-1.0*alpha[3]*fr[8]+1.732050807568877*alpha[4]*fr[7]-1.0*alpha[2]*fr[6]+1.732050807568877*alpha[0]*fr[5]-1.0*(fr[4]*alpha[5]+fr[3]*alpha[4])+1.732050807568877*alpha[1]*fr[2]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_y; 
  incr[2] = 0.1767766952966368*(3.0*(alpha[5]*fr[12]+alpha[4]*fr[11]+alpha[3]*fr[9])-1.732050807568877*alpha[5]*fr[8]+3.0*alpha[2]*fr[7]-1.732050807568877*alpha[4]*fr[6]+3.0*alpha[1]*fr[5]-1.732050807568877*(alpha[3]*fr[4]+alpha[2]*fr[3])+3.0*alpha[0]*fr[2]-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_y; 
  incr[3] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fr[15]+alpha[3]*fr[14])-1.0*alpha[5]*fr[13]+1.732050807568877*alpha[1]*fr[11]-1.0*alpha[3]*fr[10]+1.732050807568877*alpha[0]*fr[7]-1.0*alpha[1]*fr[6]+1.732050807568877*alpha[4]*fr[5]-1.0*(fr[1]*alpha[4]+alpha[0]*fr[3])+alpha[2]*(1.732050807568877*fr[2]-1.0*fr[0]))*dfac_y; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[15]+alpha[2]*fr[14])-1.0*alpha[4]*fr[13]+1.732050807568877*alpha[1]*fr[12]-1.0*alpha[2]*fr[10]+1.732050807568877*alpha[0]*fr[9]-1.0*alpha[1]*fr[8]+1.732050807568877*alpha[5]*fr[5]-1.0*(fr[1]*alpha[5]+alpha[0]*fr[4])+(1.732050807568877*fr[2]-1.0*fr[0])*alpha[3])*dfac_y; 
  incr[5] = 0.1767766952966368*(3.0*(alpha[3]*fr[12]+alpha[2]*fr[11]+alpha[5]*fr[9])-1.732050807568877*alpha[3]*fr[8]+3.0*alpha[4]*fr[7]-1.732050807568877*alpha[2]*fr[6]+3.0*alpha[0]*fr[5]-1.732050807568877*(fr[4]*alpha[5]+fr[3]*alpha[4])+3.0*alpha[1]*fr[2]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_y; 
  incr[6] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fr[15]+alpha[5]*fr[14])-1.0*alpha[3]*fr[13]+1.732050807568877*alpha[0]*fr[11]-1.0*alpha[5]*fr[10]+1.732050807568877*alpha[1]*fr[7]-1.0*alpha[0]*fr[6]+1.732050807568877*(alpha[2]*fr[5]+fr[2]*alpha[4])-1.0*(fr[0]*alpha[4]+alpha[1]*fr[3]+fr[1]*alpha[2]))*dfac_y; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[5]*fr[15]+alpha[3]*fr[14])-1.732050807568877*alpha[5]*fr[13]+3.0*alpha[1]*fr[11]-1.732050807568877*alpha[3]*fr[10]+3.0*alpha[0]*fr[7]-1.732050807568877*alpha[1]*fr[6]+3.0*alpha[4]*fr[5]-1.732050807568877*(fr[1]*alpha[4]+alpha[0]*fr[3])+alpha[2]*(3.0*fr[2]-1.732050807568877*fr[0]))*dfac_y; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[15]+alpha[4]*fr[14])-1.0*alpha[2]*fr[13]+1.732050807568877*alpha[0]*fr[12]-1.0*alpha[4]*fr[10]+1.732050807568877*alpha[1]*fr[9]-1.0*alpha[0]*fr[8]+1.732050807568877*(alpha[3]*fr[5]+fr[2]*alpha[5])-1.0*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_y; 
  incr[9] = 0.1767766952966368*(3.0*(alpha[4]*fr[15]+alpha[2]*fr[14])-1.732050807568877*alpha[4]*fr[13]+3.0*alpha[1]*fr[12]-1.732050807568877*alpha[2]*fr[10]+3.0*alpha[0]*fr[9]-1.732050807568877*alpha[1]*fr[8]+3.0*alpha[5]*fr[5]-1.732050807568877*(fr[1]*alpha[5]+alpha[0]*fr[4])+(3.0*fr[2]-1.732050807568877*fr[0])*alpha[3])*dfac_y; 
  incr[10] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14])-1.0*alpha[1]*fr[13]+1.732050807568877*(alpha[4]*fr[12]+alpha[5]*fr[11])-1.0*alpha[0]*fr[10]+1.732050807568877*alpha[2]*fr[9]-1.0*alpha[4]*fr[8]+1.732050807568877*alpha[3]*fr[7]-1.0*(alpha[5]*fr[6]+alpha[2]*fr[4]+alpha[3]*fr[3]))*dfac_y; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[3]*fr[15]+alpha[5]*fr[14])-1.732050807568877*alpha[3]*fr[13]+3.0*alpha[0]*fr[11]-1.732050807568877*alpha[5]*fr[10]+3.0*alpha[1]*fr[7]-1.732050807568877*alpha[0]*fr[6]+3.0*(alpha[2]*fr[5]+fr[2]*alpha[4])-1.732050807568877*(fr[0]*alpha[4]+alpha[1]*fr[3]+fr[1]*alpha[2]))*dfac_y; 
  incr[12] = 0.1767766952966368*(3.0*(alpha[2]*fr[15]+alpha[4]*fr[14])-1.732050807568877*alpha[2]*fr[13]+3.0*alpha[0]*fr[12]-1.732050807568877*alpha[4]*fr[10]+3.0*alpha[1]*fr[9]-1.732050807568877*alpha[0]*fr[8]+3.0*(alpha[3]*fr[5]+fr[2]*alpha[5])-1.732050807568877*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_y; 
  incr[13] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14])-1.0*alpha[0]*fr[13]+1.732050807568877*(alpha[2]*fr[12]+alpha[3]*fr[11])-1.0*alpha[1]*fr[10]+1.732050807568877*alpha[4]*fr[9]-1.0*alpha[2]*fr[8]+1.732050807568877*alpha[5]*fr[7]-1.0*(alpha[3]*fr[6]+fr[3]*alpha[5]+alpha[4]*fr[4]))*dfac_y; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14])-1.732050807568877*alpha[1]*fr[13]+3.0*(alpha[4]*fr[12]+alpha[5]*fr[11])-1.732050807568877*alpha[0]*fr[10]+3.0*alpha[2]*fr[9]-1.732050807568877*alpha[4]*fr[8]+3.0*alpha[3]*fr[7]-1.732050807568877*(alpha[5]*fr[6]+alpha[2]*fr[4]+alpha[3]*fr[3]))*dfac_y; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14])-1.732050807568877*alpha[0]*fr[13]+3.0*(alpha[2]*fr[12]+alpha[3]*fr[11])-1.732050807568877*alpha[1]*fr[10]+3.0*alpha[4]*fr[9]-1.732050807568877*alpha[2]*fr[8]+3.0*alpha[5]*fr[7]-1.732050807568877*(alpha[3]*fr[6]+fr[3]*alpha[5]+alpha[4]*fr[4]))*dfac_y; 
  }
#elif upwindType == QUAD 
double fupwind[16];
double fupwindQuad[16];
double limQuad[16];
  if(0.3535533905932737*(alphaUp[5]+alphaUp[4])-0.3535533905932737*(alphaUp[3]+alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[0] = (-0.25*fl[15])+0.25*fl[14]-0.25*fl[13]+0.25*(fl[12]+fl[11]+fl[10])-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*(fl[5]+fl[4]+fl[3])+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[0] = 0.25*fr[15]-0.25*(fr[14]+fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.25*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*(alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2]) > 0) {
  fupwindQuad[1] = 0.25*(fl[15]+fl[14]+fl[13])-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*fl[5]-0.25*(fl[4]+fl[3])+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[1] = (-0.25*(fr[15]+fr[14]))+0.25*(fr[13]+fr[12]+fr[11]+fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5]+fr[4]+fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*(alphaUp[5]+alphaUp[4])-0.3535533905932737*(alphaUp[3]+alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[2] = 0.25*fl[15]-0.25*(fl[14]+fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])-0.25*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[2] = (-0.25*fr[15])+0.25*fr[14]-0.25*fr[13]+0.25*(fr[12]+fr[11]+fr[10])-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*(fr[5]+fr[4]+fr[3])+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[1]+alphaUp[0])-0.3535533905932737*(alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2]) > 0) {
  fupwindQuad[3] = (-0.25*(fl[15]+fl[14]))+0.25*(fl[13]+fl[12]+fl[11]+fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5]+fl[4]+fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[3] = 0.25*(fr[15]+fr[14]+fr[13])-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*fr[5]-0.25*(fr[4]+fr[3])+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[5]-0.3535533905932737*(alphaUp[4]+alphaUp[3])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[4] = 0.25*fl[15]-0.25*fl[14]+0.25*(fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5]+fl[4])+0.25*(fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[4] = (-0.25*fr[15])+0.25*(fr[14]+fr[13])-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*fr[5]-0.25*fr[4]+0.25*fr[3]-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if((-0.3535533905932737*alphaUp[5])+0.3535533905932737*alphaUp[4]-0.3535533905932737*alphaUp[3]+0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[5] = (-0.25*(fl[15]+fl[14]+fl[13]+fl[12]))+0.25*fl[11]-0.25*(fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5])-0.25*fl[4]+0.25*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[5] = 0.25*(fr[15]+fr[14])-0.25*fr[13]+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*(fr[5]+fr[4])+0.25*fr[3]-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[5]-0.3535533905932737*(alphaUp[4]+alphaUp[3])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[6] = (-0.25*fl[15])+0.25*(fl[14]+fl[13])-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*fl[5]-0.25*fl[4]+0.25*fl[3]-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[6] = 0.25*fr[15]-0.25*fr[14]+0.25*(fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5]+fr[4])+0.25*(fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if((-0.3535533905932737*alphaUp[5])+0.3535533905932737*alphaUp[4]-0.3535533905932737*alphaUp[3]+0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[7] = 0.25*(fl[15]+fl[14])-0.25*fl[13]+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*(fl[5]+fl[4])+0.25*fl[3]-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[7] = (-0.25*(fr[15]+fr[14]+fr[13]+fr[12]))+0.25*fr[11]-0.25*(fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5])-0.25*fr[4]+0.25*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[5])+0.3535533905932737*(alphaUp[4]+alphaUp[3])-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[8] = 0.25*fl[15]-0.25*fl[14]+0.25*fl[13]-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*fl[5]+0.25*fl[4]-0.25*fl[3]+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[8] = (-0.25*fr[15])+0.25*(fr[14]+fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5]+fr[4])-0.25*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*alphaUp[5]-0.3535533905932737*alphaUp[4]+0.3535533905932737*alphaUp[3]-0.3535533905932737*alphaUp[2]+0.3535533905932737*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[9] = (-0.25*(fl[15]+fl[14]+fl[13]))+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*(fl[5]+fl[4])-0.25*fl[3]+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[9] = 0.25*(fr[15]+fr[14])-0.25*(fr[13]+fr[12])+0.25*fr[11]-0.25*(fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5])+0.25*fr[4]-0.25*(fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[5])+0.3535533905932737*(alphaUp[4]+alphaUp[3])-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[10] = (-0.25*fl[15])+0.25*(fl[14]+fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5]+fl[4])-0.25*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[10] = 0.25*fr[15]-0.25*fr[14]+0.25*fr[13]-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*fr[5]+0.25*fr[4]-0.25*fr[3]+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*alphaUp[5]-0.3535533905932737*alphaUp[4]+0.3535533905932737*alphaUp[3]-0.3535533905932737*alphaUp[2]+0.3535533905932737*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[11] = 0.25*(fl[15]+fl[14])-0.25*(fl[13]+fl[12])+0.25*fl[11]-0.25*(fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5])+0.25*fl[4]-0.25*(fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[11] = (-0.25*(fr[15]+fr[14]+fr[13]))+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*(fr[5]+fr[4])-0.25*fr[3]+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*(alphaUp[5]+alphaUp[4]))+0.3535533905932737*(alphaUp[3]+alphaUp[2])-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[12] = (-0.25*fl[15])+0.25*fl[14]-0.25*(fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5])+0.25*(fl[4]+fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[12] = 0.25*fr[15]-0.25*(fr[14]+fr[13])+0.25*(fr[12]+fr[11]+fr[10])-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*(fr[5]+fr[4]+fr[3])-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[13] = 0.25*(fl[15]+fl[14]+fl[13]+fl[12]+fl[11]+fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5]+fl[4]+fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[13] = (-0.25*(fr[15]+fr[14]))+0.25*fr[13]-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*fr[5]+0.25*(fr[4]+fr[3])-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*(alphaUp[5]+alphaUp[4]))+0.3535533905932737*(alphaUp[3]+alphaUp[2])-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[14] = 0.25*fl[15]-0.25*(fl[14]+fl[13])+0.25*(fl[12]+fl[11]+fl[10])-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*(fl[5]+fl[4]+fl[3])-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[14] = (-0.25*fr[15])+0.25*fr[14]-0.25*(fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5])+0.25*(fr[4]+fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[15] = (-0.25*(fl[15]+fl[14]))+0.25*fl[13]-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*fl[5]+0.25*(fl[4]+fl[3])-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[15] = 0.25*(fr[15]+fr[14]+fr[13]+fr[12]+fr[11]+fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5]+fr[4]+fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[4] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[5] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  fupwind[6] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[7] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[8] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[9] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[10] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[11] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[12] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[13] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[14] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[15] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fupwind[12]+alpha[4]*fupwind[11]+alpha[3]*fupwind[9])-1.0*alpha[5]*fupwind[8]+1.732050807568877*alpha[2]*fupwind[7]-1.0*alpha[4]*fupwind[6]+1.732050807568877*alpha[1]*fupwind[5]-1.0*(alpha[3]*fupwind[4]+alpha[2]*fupwind[3])+1.732050807568877*alpha[0]*fupwind[2]-1.0*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_y; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fupwind[12]+alpha[2]*fupwind[11]+alpha[5]*fupwind[9])-1.0*alpha[3]*fupwind[8]+1.732050807568877*alpha[4]*fupwind[7]-1.0*alpha[2]*fupwind[6]+1.732050807568877*alpha[0]*fupwind[5]-1.0*(fupwind[4]*alpha[5]+fupwind[3]*alpha[4])+1.732050807568877*alpha[1]*fupwind[2]-1.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_y; 
  incr[2] = 0.1767766952966368*(3.0*(alpha[5]*fupwind[12]+alpha[4]*fupwind[11]+alpha[3]*fupwind[9])-1.732050807568877*alpha[5]*fupwind[8]+3.0*alpha[2]*fupwind[7]-1.732050807568877*alpha[4]*fupwind[6]+3.0*alpha[1]*fupwind[5]-1.732050807568877*(alpha[3]*fupwind[4]+alpha[2]*fupwind[3])+3.0*alpha[0]*fupwind[2]-1.732050807568877*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_y; 
  incr[3] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fupwind[15]+alpha[3]*fupwind[14])-1.0*alpha[5]*fupwind[13]+1.732050807568877*alpha[1]*fupwind[11]-1.0*alpha[3]*fupwind[10]+1.732050807568877*alpha[0]*fupwind[7]-1.0*alpha[1]*fupwind[6]+1.732050807568877*alpha[4]*fupwind[5]-1.0*(fupwind[1]*alpha[4]+alpha[0]*fupwind[3])+alpha[2]*(1.732050807568877*fupwind[2]-1.0*fupwind[0]))*dfac_y; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fupwind[15]+alpha[2]*fupwind[14])-1.0*alpha[4]*fupwind[13]+1.732050807568877*alpha[1]*fupwind[12]-1.0*alpha[2]*fupwind[10]+1.732050807568877*alpha[0]*fupwind[9]-1.0*alpha[1]*fupwind[8]+1.732050807568877*alpha[5]*fupwind[5]-1.0*(fupwind[1]*alpha[5]+alpha[0]*fupwind[4])+(1.732050807568877*fupwind[2]-1.0*fupwind[0])*alpha[3])*dfac_y; 
  incr[5] = 0.1767766952966368*(3.0*(alpha[3]*fupwind[12]+alpha[2]*fupwind[11]+alpha[5]*fupwind[9])-1.732050807568877*alpha[3]*fupwind[8]+3.0*alpha[4]*fupwind[7]-1.732050807568877*alpha[2]*fupwind[6]+3.0*alpha[0]*fupwind[5]-1.732050807568877*(fupwind[4]*alpha[5]+fupwind[3]*alpha[4])+3.0*alpha[1]*fupwind[2]-1.732050807568877*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_y; 
  incr[6] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fupwind[15]+alpha[5]*fupwind[14])-1.0*alpha[3]*fupwind[13]+1.732050807568877*alpha[0]*fupwind[11]-1.0*alpha[5]*fupwind[10]+1.732050807568877*alpha[1]*fupwind[7]-1.0*alpha[0]*fupwind[6]+1.732050807568877*(alpha[2]*fupwind[5]+fupwind[2]*alpha[4])-1.0*(fupwind[0]*alpha[4]+alpha[1]*fupwind[3]+fupwind[1]*alpha[2]))*dfac_y; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[5]*fupwind[15]+alpha[3]*fupwind[14])-1.732050807568877*alpha[5]*fupwind[13]+3.0*alpha[1]*fupwind[11]-1.732050807568877*alpha[3]*fupwind[10]+3.0*alpha[0]*fupwind[7]-1.732050807568877*alpha[1]*fupwind[6]+3.0*alpha[4]*fupwind[5]-1.732050807568877*(fupwind[1]*alpha[4]+alpha[0]*fupwind[3])+alpha[2]*(3.0*fupwind[2]-1.732050807568877*fupwind[0]))*dfac_y; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[15]+alpha[4]*fupwind[14])-1.0*alpha[2]*fupwind[13]+1.732050807568877*alpha[0]*fupwind[12]-1.0*alpha[4]*fupwind[10]+1.732050807568877*alpha[1]*fupwind[9]-1.0*alpha[0]*fupwind[8]+1.732050807568877*(alpha[3]*fupwind[5]+fupwind[2]*alpha[5])-1.0*(fupwind[0]*alpha[5]+alpha[1]*fupwind[4]+fupwind[1]*alpha[3]))*dfac_y; 
  incr[9] = 0.1767766952966368*(3.0*(alpha[4]*fupwind[15]+alpha[2]*fupwind[14])-1.732050807568877*alpha[4]*fupwind[13]+3.0*alpha[1]*fupwind[12]-1.732050807568877*alpha[2]*fupwind[10]+3.0*alpha[0]*fupwind[9]-1.732050807568877*alpha[1]*fupwind[8]+3.0*alpha[5]*fupwind[5]-1.732050807568877*(fupwind[1]*alpha[5]+alpha[0]*fupwind[4])+(3.0*fupwind[2]-1.732050807568877*fupwind[0])*alpha[3])*dfac_y; 
  incr[10] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14])-1.0*alpha[1]*fupwind[13]+1.732050807568877*(alpha[4]*fupwind[12]+alpha[5]*fupwind[11])-1.0*alpha[0]*fupwind[10]+1.732050807568877*alpha[2]*fupwind[9]-1.0*alpha[4]*fupwind[8]+1.732050807568877*alpha[3]*fupwind[7]-1.0*(alpha[5]*fupwind[6]+alpha[2]*fupwind[4]+alpha[3]*fupwind[3]))*dfac_y; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[3]*fupwind[15]+alpha[5]*fupwind[14])-1.732050807568877*alpha[3]*fupwind[13]+3.0*alpha[0]*fupwind[11]-1.732050807568877*alpha[5]*fupwind[10]+3.0*alpha[1]*fupwind[7]-1.732050807568877*alpha[0]*fupwind[6]+3.0*(alpha[2]*fupwind[5]+fupwind[2]*alpha[4])-1.732050807568877*(fupwind[0]*alpha[4]+alpha[1]*fupwind[3]+fupwind[1]*alpha[2]))*dfac_y; 
  incr[12] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[15]+alpha[4]*fupwind[14])-1.732050807568877*alpha[2]*fupwind[13]+3.0*alpha[0]*fupwind[12]-1.732050807568877*alpha[4]*fupwind[10]+3.0*alpha[1]*fupwind[9]-1.732050807568877*alpha[0]*fupwind[8]+3.0*(alpha[3]*fupwind[5]+fupwind[2]*alpha[5])-1.732050807568877*(fupwind[0]*alpha[5]+alpha[1]*fupwind[4]+fupwind[1]*alpha[3]))*dfac_y; 
  incr[13] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14])-1.0*alpha[0]*fupwind[13]+1.732050807568877*(alpha[2]*fupwind[12]+alpha[3]*fupwind[11])-1.0*alpha[1]*fupwind[10]+1.732050807568877*alpha[4]*fupwind[9]-1.0*alpha[2]*fupwind[8]+1.732050807568877*alpha[5]*fupwind[7]-1.0*(alpha[3]*fupwind[6]+fupwind[3]*alpha[5]+alpha[4]*fupwind[4]))*dfac_y; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14])-1.732050807568877*alpha[1]*fupwind[13]+3.0*(alpha[4]*fupwind[12]+alpha[5]*fupwind[11])-1.732050807568877*alpha[0]*fupwind[10]+3.0*alpha[2]*fupwind[9]-1.732050807568877*alpha[4]*fupwind[8]+3.0*alpha[3]*fupwind[7]-1.732050807568877*(alpha[5]*fupwind[6]+alpha[2]*fupwind[4]+alpha[3]*fupwind[3]))*dfac_y; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14])-1.732050807568877*alpha[0]*fupwind[13]+3.0*(alpha[2]*fupwind[12]+alpha[3]*fupwind[11])-1.732050807568877*alpha[1]*fupwind[10]+3.0*alpha[4]*fupwind[9]-1.732050807568877*alpha[2]*fupwind[8]+3.0*alpha[5]*fupwind[7]-1.732050807568877*(alpha[3]*fupwind[6]+fupwind[3]*alpha[5]+alpha[4]*fupwind[4]))*dfac_y; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  double incrEmMod[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.0625*(1.732050807568877*(2.0*dfac_v*m_*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*wv+Bmag[1]*dfac_x*((Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*dfac_v*q_-2.0*BdriftX[0]*m_)*wm+q_*(dfac_v*(((Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[3]+(Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*Phi[2])*dfac_y+((BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*Phi[3]+(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[1])*dfac_x)*q_-2.0*((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_))+dfac_v*q_*(3.0*Bmag[1]*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_x*dfac_y*wm+(4.0*dApardtPrev[0]-3.0*(BmagInv[1]*(Apar[1]*Phi[3]-1.0*Phi[1]*Apar[3])+BmagInv[0]*(Apar[1]*Phi[2]-1.0*Phi[1]*Apar[2]))*dfac_x*dfac_y)*q_)))/(dfac_v*m_*q_); 

  double alpha[8]; 
  alpha[0] = -(0.3535533905932737*(1.732050807568877*(2.0*dfac_v*m_*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*wv+Bmag[1]*dfac_x*((Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*dfac_v*q_-2.0*BdriftX[0]*m_)*wm-2.0*((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*q_+dfac_v*(((Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[3]+(Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*Phi[2])*dfac_y+((BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*Phi[3]+(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[1])*dfac_x)*q2)+dfac_v*dfac_x*dfac_y*(3.0*Bmag[1]*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*q_*wm+(BmagInv[1]*(3.0*Phi[1]*Apar[3]-3.0*Apar[1]*Phi[3])+BmagInv[0]*(3.0*Phi[1]*Apar[2]-3.0*Apar[1]*Phi[2]))*q2)))/(dfac_v*m_*q_); 
  alpha[1] = -(0.07071067811865474*(1.732050807568877*(10.0*dfac_v*m_*(BdriftX[1]*Bmag[1]*dfac_x*wm+((BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*q_)*wv+Bmag[1]*dfac_x*(5.0*(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*dfac_v*q_-10.0*BdriftX[1]*m_)*wm-10.0*((BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*m_*q_+dfac_v*(((9.0*Apar[1]*BdriftY[1]+5.0*Apar[0]*BdriftY[0])*Phi[3]+5.0*(Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[2])*dfac_y+5.0*((BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*Phi[3]+(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*Phi[1])*dfac_x)*q2)+dfac_v*dfac_x*dfac_y*(15.0*Bmag[1]*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*q_*wm+(BmagInv[0]*(15.0*Phi[1]*Apar[3]-15.0*Apar[1]*Phi[3])+BmagInv[1]*(15.0*Phi[1]*Apar[2]-15.0*Apar[1]*Phi[2]))*q2)))/(dfac_v*m_*q_); 
  alpha[2] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(dfac_x*(2.0*BdriftX[0]*Phi[3]*m_*wv+Bmag[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*wm)+(((BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2])*Phi[3]+Phi[2]*(BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2]))*dfac_y+((Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[3]+Phi[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2]))*dfac_x)*q_)-2.0*BdriftX[0]*Phi[3]*dfac_x*m_)+BmagInv[0]*(3.0*Apar[2]*Phi[3]-3.0*Phi[2]*Apar[3])*dfac_v*dfac_x*dfac_y*q_))/(dfac_v*m_); 
  alpha[3] = -(0.3535533905932737*Bmag[1]*dfac_x*(dfac_v*(2.0*BdriftX[0]*m_*wv+(1.732050807568877*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_y+Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*q_)-2.0*BdriftX[0]*m_))/(dfac_m*dfac_v*m_*q_); 
  alpha[4] = -(0.07071067811865474*(1.732050807568877*(dfac_v*(dfac_x*(10.0*BdriftX[1]*Phi[3]*m_*wv+5.0*Bmag[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*wm)+(((9.0*BdriftY[1]*Apar[3]+5.0*BdriftY[0]*Apar[2])*Phi[3]+5.0*Phi[2]*(BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2]))*dfac_y+5.0*((Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*Phi[3]+Phi[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2]))*dfac_x)*q_)-10.0*BdriftX[1]*Phi[3]*dfac_x*m_)+BmagInv[1]*(15.0*Apar[2]*Phi[3]-15.0*Phi[2]*Apar[3])*dfac_v*dfac_x*dfac_y*q_))/(dfac_v*m_); 
  alpha[5] = -(0.3535533905932737*Bmag[1]*dfac_x*(dfac_v*(2.0*BdriftX[1]*m_*wv+(1.732050807568877*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_y+Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*q_)-2.0*BdriftX[1]*m_))/(dfac_m*dfac_v*m_*q_); 
  alpha[6] = -(0.3535533905932737*Bmag[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*dfac_x)/(dfac_m*m_); 
  alpha[7] = -(0.3535533905932737*Bmag[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*dfac_x)/(dfac_m*m_); 
  double alphaUp[8]; 
  alphaUp[0] = -(0.3535533905932737*(1.732050807568877*(2.0*dfac_v*m_*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*wv+Bmag[1]*dfac_x*((Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*dfac_v*q_-2.0*BdriftX[0]*m_)*wm-2.0*((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*q_+dfac_v*(((Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[3]+(Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*Phi[2])*dfac_y+((BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*Phi[3]+(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[1])*dfac_x)*q2)+dfac_v*(3.0*Bmag[1]*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_x*dfac_y*q_*wm+((BmagInv[1]*(3.0*Phi[1]*Apar[3]-3.0*Apar[1]*Phi[3])+BmagInv[0]*(3.0*Phi[1]*Apar[2]-3.0*Apar[1]*Phi[2]))*dfac_x*dfac_y+4.0*dApardtPrev[0])*q2)))/(dfac_v*m_*q_); 
  alphaUp[1] = -(0.07071067811865474*(1.732050807568877*(10.0*dfac_v*m_*(BdriftX[1]*Bmag[1]*dfac_x*wm+((BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*q_)*wv+Bmag[1]*dfac_x*(5.0*(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*dfac_v*q_-10.0*BdriftX[1]*m_)*wm-10.0*((BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*m_*q_+dfac_v*(((9.0*Apar[1]*BdriftY[1]+5.0*Apar[0]*BdriftY[0])*Phi[3]+5.0*(Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[2])*dfac_y+5.0*((BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*Phi[3]+(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*Phi[1])*dfac_x)*q2)+dfac_v*(15.0*Bmag[1]*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_x*dfac_y*q_*wm+((BmagInv[0]*(15.0*Phi[1]*Apar[3]-15.0*Apar[1]*Phi[3])+BmagInv[1]*(15.0*Phi[1]*Apar[2]-15.0*Apar[1]*Phi[2]))*dfac_x*dfac_y+20.0*dApardtPrev[1])*q2)))/(dfac_v*m_*q_); 
  alphaUp[2] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(dfac_x*(2.0*BdriftX[0]*Phi[3]*m_*wv+Bmag[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*wm)+(((BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2])*Phi[3]+Phi[2]*(BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2]))*dfac_y+((Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[3]+Phi[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2]))*dfac_x)*q_)-2.0*BdriftX[0]*Phi[3]*dfac_x*m_)+dfac_v*(BmagInv[0]*(3.0*Apar[2]*Phi[3]-3.0*Phi[2]*Apar[3])*dfac_x*dfac_y+4.0*dApardtPrev[2])*q_))/(dfac_v*m_); 
  alphaUp[3] = -(0.3535533905932737*Bmag[1]*dfac_x*(dfac_v*(2.0*BdriftX[0]*m_*wv+(1.732050807568877*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_y+Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*q_)-2.0*BdriftX[0]*m_))/(dfac_m*dfac_v*m_*q_); 
  alphaUp[4] = -(0.07071067811865474*(1.732050807568877*(dfac_v*(dfac_x*(10.0*BdriftX[1]*Phi[3]*m_*wv+5.0*Bmag[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*wm)+(((9.0*BdriftY[1]*Apar[3]+5.0*BdriftY[0]*Apar[2])*Phi[3]+5.0*Phi[2]*(BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2]))*dfac_y+5.0*((Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*Phi[3]+Phi[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2]))*dfac_x)*q_)-10.0*BdriftX[1]*Phi[3]*dfac_x*m_)+dfac_v*(BmagInv[1]*(15.0*Apar[2]*Phi[3]-15.0*Phi[2]*Apar[3])*dfac_x*dfac_y+20.0*dApardtPrev[3])*q_))/(dfac_v*m_); 
  alphaUp[5] = -(0.3535533905932737*Bmag[1]*dfac_x*(dfac_v*(2.0*BdriftX[1]*m_*wv+(1.732050807568877*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_y+Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*q_)-2.0*BdriftX[1]*m_))/(dfac_m*dfac_v*m_*q_); 
  alphaUp[6] = -(0.3535533905932737*Bmag[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*dfac_x)/(dfac_m*m_); 
  alphaUp[7] = -(0.3535533905932737*Bmag[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*dfac_x)/(dfac_m*m_); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  emModR[0] += 0.5*(1.732050807568877*fl[3]+fl[0])*dfac_v; 
  emModR[1] += 0.5*(1.732050807568877*fl[6]+fl[1])*dfac_v; 
  emModR[2] += 0.5*(1.732050807568877*fl[7]+fl[2])*dfac_v; 
  emModR[3] += -0.5*(3.0*fl[3]+1.732050807568877*fl[0])*dfac_v; 
  emModR[4] += 0.5*(1.732050807568877*fl[10]+fl[4])*dfac_v; 
  emModR[5] += 0.5*(1.732050807568877*fl[11]+fl[5])*dfac_v; 
  emModR[6] += -0.5*(3.0*fl[6]+1.732050807568877*fl[1])*dfac_v; 
  emModR[7] += -0.5*(3.0*fl[7]+1.732050807568877*fl[2])*dfac_v; 
  emModR[8] += 0.5*(1.732050807568877*fl[13]+fl[8])*dfac_v; 
  emModR[9] += 0.5*(1.732050807568877*fl[14]+fl[9])*dfac_v; 
  emModR[10] += -0.5*(3.0*fl[10]+1.732050807568877*fl[4])*dfac_v; 
  emModR[11] += -0.5*(3.0*fl[11]+1.732050807568877*fl[5])*dfac_v; 
  emModR[12] += 0.5*(1.732050807568877*fl[15]+fl[12])*dfac_v; 
  emModR[13] += -0.5*(3.0*fl[13]+1.732050807568877*fl[8])*dfac_v; 
  emModR[14] += -0.5*(3.0*fl[14]+1.732050807568877*fl[9])*dfac_v; 
  emModR[15] += -0.5*(3.0*fl[15]+1.732050807568877*fl[12])*dfac_v; 
  emModL[0] += -0.5*(1.732050807568877*fl[3]+fl[0])*dfac_v; 
  emModL[1] += -0.5*(1.732050807568877*fl[6]+fl[1])*dfac_v; 
  emModL[2] += -0.5*(1.732050807568877*fl[7]+fl[2])*dfac_v; 
  emModL[3] += -0.5*(3.0*fl[3]+1.732050807568877*fl[0])*dfac_v; 
  emModL[4] += -0.5*(1.732050807568877*fl[10]+fl[4])*dfac_v; 
  emModL[5] += -0.5*(1.732050807568877*fl[11]+fl[5])*dfac_v; 
  emModL[6] += -0.5*(3.0*fl[6]+1.732050807568877*fl[1])*dfac_v; 
  emModL[7] += -0.5*(3.0*fl[7]+1.732050807568877*fl[2])*dfac_v; 
  emModL[8] += -0.5*(1.732050807568877*fl[13]+fl[8])*dfac_v; 
  emModL[9] += -0.5*(1.732050807568877*fl[14]+fl[9])*dfac_v; 
  emModL[10] += -0.5*(3.0*fl[10]+1.732050807568877*fl[4])*dfac_v; 
  emModL[11] += -0.5*(3.0*fl[11]+1.732050807568877*fl[5])*dfac_v; 
  emModL[12] += -0.5*(1.732050807568877*fl[15]+fl[12])*dfac_v; 
  emModL[13] += -0.5*(3.0*fl[13]+1.732050807568877*fl[8])*dfac_v; 
  emModL[14] += -0.5*(3.0*fl[14]+1.732050807568877*fl[9])*dfac_v; 
  emModL[15] += -0.5*(3.0*fl[15]+1.732050807568877*fl[12])*dfac_v; 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[7]*fl[15]+alpha[6]*fl[14]+alpha[5]*fl[13])+alpha[7]*fl[12]+1.732050807568877*(alpha[4]*fl[11]+alpha[3]*fl[10])+alpha[6]*fl[9]+alpha[5]*fl[8]+1.732050807568877*(alpha[2]*fl[7]+alpha[1]*fl[6])+alpha[4]*fl[5]+alpha[3]*fl[4]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*(alpha[6]*fl[15]+alpha[7]*fl[14]+alpha[3]*fl[13])+alpha[6]*fl[12]+1.732050807568877*(alpha[2]*fl[11]+alpha[5]*fl[10])+alpha[7]*fl[9]+alpha[3]*fl[8]+1.732050807568877*(alpha[4]*fl[7]+alpha[0]*fl[6])+alpha[2]*fl[5]+fl[4]*alpha[5]+fl[2]*alpha[4]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(1.732050807568877*(alpha[5]*fl[15]+alpha[3]*fl[14]+alpha[7]*fl[13])+alpha[5]*fl[12]+1.732050807568877*(alpha[1]*fl[11]+alpha[6]*fl[10])+alpha[3]*fl[9]+alpha[7]*fl[8]+1.732050807568877*(alpha[0]*fl[7]+alpha[4]*fl[6])+fl[4]*alpha[6]+alpha[1]*fl[5]+fl[1]*alpha[4]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.1767766952966368*(3.0*(alpha[7]*fl[15]+alpha[6]*fl[14]+alpha[5]*fl[13])+1.732050807568877*alpha[7]*fl[12]+3.0*(alpha[4]*fl[11]+alpha[3]*fl[10])+1.732050807568877*(alpha[6]*fl[9]+alpha[5]*fl[8])+3.0*(alpha[2]*fl[7]+alpha[1]*fl[6])+1.732050807568877*(alpha[4]*fl[5]+alpha[3]*fl[4])+3.0*alpha[0]*fl[3]+1.732050807568877*(alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+alpha[4]*fl[12]+1.732050807568877*(alpha[7]*fl[11]+alpha[0]*fl[10])+alpha[2]*fl[9]+alpha[1]*fl[8]+1.732050807568877*alpha[6]*fl[7]+fl[5]*alpha[7]+1.732050807568877*alpha[5]*fl[6]+fl[2]*alpha[6]+fl[1]*alpha[5]+alpha[0]*fl[4]+alpha[3]*(1.732050807568877*fl[3]+fl[0]))*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*(alpha[3]*fl[15]+alpha[5]*fl[14]+alpha[6]*fl[13])+alpha[3]*fl[12]+1.732050807568877*(alpha[0]*fl[11]+alpha[7]*fl[10])+alpha[5]*fl[9]+alpha[6]*fl[8]+1.732050807568877*alpha[1]*fl[7]+fl[4]*alpha[7]+1.732050807568877*alpha[2]*fl[6]+alpha[0]*fl[5]+(1.732050807568877*fl[3]+fl[0])*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*(alpha[6]*fl[15]+alpha[7]*fl[14]+alpha[3]*fl[13])+1.732050807568877*alpha[6]*fl[12]+3.0*(alpha[2]*fl[11]+alpha[5]*fl[10])+1.732050807568877*(alpha[7]*fl[9]+alpha[3]*fl[8])+3.0*(alpha[4]*fl[7]+alpha[0]*fl[6])+1.732050807568877*(alpha[2]*fl[5]+fl[4]*alpha[5]+fl[2]*alpha[4])+3.0*alpha[1]*fl[3]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[5]*fl[15]+alpha[3]*fl[14]+alpha[7]*fl[13])+1.732050807568877*alpha[5]*fl[12]+3.0*(alpha[1]*fl[11]+alpha[6]*fl[10])+1.732050807568877*(alpha[3]*fl[9]+alpha[7]*fl[8])+3.0*(alpha[0]*fl[7]+alpha[4]*fl[6])+1.732050807568877*(fl[4]*alpha[6]+alpha[1]*fl[5]+fl[1]*alpha[4])+3.0*alpha[2]*fl[3]+1.732050807568877*(alpha[0]*fl[2]+fl[0]*alpha[2]))*dfac_v; 
  incr[8] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+alpha[2]*fl[12]+1.732050807568877*(alpha[6]*fl[11]+alpha[1]*fl[10])+alpha[4]*fl[9]+alpha[0]*fl[8]+alpha[7]*(1.732050807568877*fl[7]+fl[2])+1.732050807568877*alpha[3]*fl[6]+fl[5]*alpha[6]+(1.732050807568877*fl[3]+fl[0])*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3])*dfac_v; 
  incr[9] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+alpha[1]*fl[12]+1.732050807568877*(alpha[5]*fl[11]+alpha[2]*fl[10])+alpha[0]*fl[9]+alpha[4]*fl[8]+1.732050807568877*alpha[3]*fl[7]+(1.732050807568877*fl[6]+fl[1])*alpha[7]+(1.732050807568877*fl[3]+fl[0])*alpha[6]+alpha[5]*fl[5]+alpha[2]*fl[4]+fl[2]*alpha[3])*dfac_v; 
  incr[10] = -0.1767766952966368*(3.0*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+1.732050807568877*alpha[4]*fl[12]+3.0*(alpha[7]*fl[11]+alpha[0]*fl[10])+1.732050807568877*(alpha[2]*fl[9]+alpha[1]*fl[8])+3.0*alpha[6]*fl[7]+1.732050807568877*fl[5]*alpha[7]+3.0*alpha[5]*fl[6]+1.732050807568877*(fl[2]*alpha[6]+fl[1]*alpha[5]+alpha[0]*fl[4])+alpha[3]*(3.0*fl[3]+1.732050807568877*fl[0]))*dfac_v; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[3]*fl[15]+alpha[5]*fl[14]+alpha[6]*fl[13])+1.732050807568877*alpha[3]*fl[12]+3.0*(alpha[0]*fl[11]+alpha[7]*fl[10])+1.732050807568877*(alpha[5]*fl[9]+alpha[6]*fl[8])+3.0*alpha[1]*fl[7]+1.732050807568877*fl[4]*alpha[7]+3.0*alpha[2]*fl[6]+1.732050807568877*alpha[0]*fl[5]+3.0*fl[3]*alpha[4]+1.732050807568877*(fl[0]*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2]))*dfac_v; 
  incr[12] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+alpha[0]*fl[12]+1.732050807568877*(alpha[3]*fl[11]+alpha[4]*fl[10])+alpha[1]*fl[9]+alpha[2]*fl[8]+1.732050807568877*alpha[5]*fl[7]+(1.732050807568877*fl[3]+fl[0])*alpha[7]+alpha[6]*(1.732050807568877*fl[6]+fl[1])+alpha[3]*fl[5]+fl[2]*alpha[5]+alpha[4]*fl[4])*dfac_v; 
  incr[13] = -0.1767766952966368*(3.0*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+1.732050807568877*alpha[2]*fl[12]+3.0*(alpha[6]*fl[11]+alpha[1]*fl[10])+1.732050807568877*(alpha[4]*fl[9]+alpha[0]*fl[8])+alpha[7]*(3.0*fl[7]+1.732050807568877*fl[2])+3.0*alpha[3]*fl[6]+1.732050807568877*fl[5]*alpha[6]+3.0*fl[3]*alpha[5]+1.732050807568877*(fl[0]*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3]))*dfac_v; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+1.732050807568877*alpha[1]*fl[12]+3.0*(alpha[5]*fl[11]+alpha[2]*fl[10])+1.732050807568877*(alpha[0]*fl[9]+alpha[4]*fl[8])+3.0*alpha[3]*fl[7]+(3.0*fl[6]+1.732050807568877*fl[1])*alpha[7]+3.0*fl[3]*alpha[6]+1.732050807568877*(fl[0]*alpha[6]+alpha[5]*fl[5]+alpha[2]*fl[4]+fl[2]*alpha[3]))*dfac_v; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+1.732050807568877*alpha[0]*fl[12]+3.0*(alpha[3]*fl[11]+alpha[4]*fl[10])+1.732050807568877*(alpha[1]*fl[9]+alpha[2]*fl[8])+3.0*alpha[5]*fl[7]+(3.0*fl[3]+1.732050807568877*fl[0])*alpha[7]+3.0*alpha[6]*fl[6]+1.732050807568877*(fl[1]*alpha[6]+alpha[3]*fl[5]+fl[2]*alpha[5]+alpha[4]*fl[4]))*dfac_v; 
  } else { 
  emModR[0] += -0.5*(1.732050807568877*fr[3]-1.0*fr[0])*dfac_v; 
  emModR[1] += -0.5*(1.732050807568877*fr[6]-1.0*fr[1])*dfac_v; 
  emModR[2] += -0.5*(1.732050807568877*fr[7]-1.0*fr[2])*dfac_v; 
  emModR[3] += 0.5*(3.0*fr[3]-1.732050807568877*fr[0])*dfac_v; 
  emModR[4] += -0.5*(1.732050807568877*fr[10]-1.0*fr[4])*dfac_v; 
  emModR[5] += -0.5*(1.732050807568877*fr[11]-1.0*fr[5])*dfac_v; 
  emModR[6] += 0.5*(3.0*fr[6]-1.732050807568877*fr[1])*dfac_v; 
  emModR[7] += 0.5*(3.0*fr[7]-1.732050807568877*fr[2])*dfac_v; 
  emModR[8] += -0.5*(1.732050807568877*fr[13]-1.0*fr[8])*dfac_v; 
  emModR[9] += -0.5*(1.732050807568877*fr[14]-1.0*fr[9])*dfac_v; 
  emModR[10] += 0.5*(3.0*fr[10]-1.732050807568877*fr[4])*dfac_v; 
  emModR[11] += 0.5*(3.0*fr[11]-1.732050807568877*fr[5])*dfac_v; 
  emModR[12] += -0.5*(1.732050807568877*fr[15]-1.0*fr[12])*dfac_v; 
  emModR[13] += 0.5*(3.0*fr[13]-1.732050807568877*fr[8])*dfac_v; 
  emModR[14] += 0.5*(3.0*fr[14]-1.732050807568877*fr[9])*dfac_v; 
  emModR[15] += 0.5*(3.0*fr[15]-1.732050807568877*fr[12])*dfac_v; 
  emModL[0] += 0.5*(1.732050807568877*fr[3]-1.0*fr[0])*dfac_v; 
  emModL[1] += 0.5*(1.732050807568877*fr[6]-1.0*fr[1])*dfac_v; 
  emModL[2] += 0.5*(1.732050807568877*fr[7]-1.0*fr[2])*dfac_v; 
  emModL[3] += 0.5*(3.0*fr[3]-1.732050807568877*fr[0])*dfac_v; 
  emModL[4] += 0.5*(1.732050807568877*fr[10]-1.0*fr[4])*dfac_v; 
  emModL[5] += 0.5*(1.732050807568877*fr[11]-1.0*fr[5])*dfac_v; 
  emModL[6] += 0.5*(3.0*fr[6]-1.732050807568877*fr[1])*dfac_v; 
  emModL[7] += 0.5*(3.0*fr[7]-1.732050807568877*fr[2])*dfac_v; 
  emModL[8] += 0.5*(1.732050807568877*fr[13]-1.0*fr[8])*dfac_v; 
  emModL[9] += 0.5*(1.732050807568877*fr[14]-1.0*fr[9])*dfac_v; 
  emModL[10] += 0.5*(3.0*fr[10]-1.732050807568877*fr[4])*dfac_v; 
  emModL[11] += 0.5*(3.0*fr[11]-1.732050807568877*fr[5])*dfac_v; 
  emModL[12] += 0.5*(1.732050807568877*fr[15]-1.0*fr[12])*dfac_v; 
  emModL[13] += 0.5*(3.0*fr[13]-1.732050807568877*fr[8])*dfac_v; 
  emModL[14] += 0.5*(3.0*fr[14]-1.732050807568877*fr[9])*dfac_v; 
  emModL[15] += 0.5*(3.0*fr[15]-1.732050807568877*fr[12])*dfac_v; 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[7]*fr[15]+alpha[6]*fr[14]+alpha[5]*fr[13])-1.0*alpha[7]*fr[12]+1.732050807568877*(alpha[4]*fr[11]+alpha[3]*fr[10])-1.0*(alpha[6]*fr[9]+alpha[5]*fr[8])+1.732050807568877*(alpha[2]*fr[7]+alpha[1]*fr[6])-1.0*(alpha[4]*fr[5]+alpha[3]*fr[4])+1.732050807568877*alpha[0]*fr[3]-1.0*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[6]*fr[15]+alpha[7]*fr[14]+alpha[3]*fr[13])-1.0*alpha[6]*fr[12]+1.732050807568877*(alpha[2]*fr[11]+alpha[5]*fr[10])-1.0*(alpha[7]*fr[9]+alpha[3]*fr[8])+1.732050807568877*(alpha[4]*fr[7]+alpha[0]*fr[6])-1.0*(alpha[2]*fr[5]+fr[4]*alpha[5]+fr[2]*alpha[4])+1.732050807568877*alpha[1]*fr[3]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fr[15]+alpha[3]*fr[14]+alpha[7]*fr[13])-1.0*alpha[5]*fr[12]+1.732050807568877*(alpha[1]*fr[11]+alpha[6]*fr[10])-1.0*(alpha[3]*fr[9]+alpha[7]*fr[8])+1.732050807568877*(alpha[0]*fr[7]+alpha[4]*fr[6])-1.0*(fr[4]*alpha[6]+alpha[1]*fr[5]+fr[1]*alpha[4])+1.732050807568877*alpha[2]*fr[3]-1.0*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[7]*fr[15]+alpha[6]*fr[14]+alpha[5]*fr[13])-1.732050807568877*alpha[7]*fr[12]+3.0*(alpha[4]*fr[11]+alpha[3]*fr[10])-1.732050807568877*(alpha[6]*fr[9]+alpha[5]*fr[8])+3.0*(alpha[2]*fr[7]+alpha[1]*fr[6])-1.732050807568877*(alpha[4]*fr[5]+alpha[3]*fr[4])+3.0*alpha[0]*fr[3]-1.732050807568877*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.0*alpha[4]*fr[12]+1.732050807568877*(alpha[7]*fr[11]+alpha[0]*fr[10])-1.0*(alpha[2]*fr[9]+alpha[1]*fr[8])+1.732050807568877*alpha[6]*fr[7]-1.0*fr[5]*alpha[7]+1.732050807568877*alpha[5]*fr[6]-1.0*(fr[2]*alpha[6]+fr[1]*alpha[5]+alpha[0]*fr[4])+alpha[3]*(1.732050807568877*fr[3]-1.0*fr[0]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fr[15]+alpha[5]*fr[14]+alpha[6]*fr[13])-1.0*alpha[3]*fr[12]+1.732050807568877*(alpha[0]*fr[11]+alpha[7]*fr[10])-1.0*(alpha[5]*fr[9]+alpha[6]*fr[8])+1.732050807568877*alpha[1]*fr[7]-1.0*fr[4]*alpha[7]+1.732050807568877*alpha[2]*fr[6]-1.0*alpha[0]*fr[5]+1.732050807568877*fr[3]*alpha[4]-1.0*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[6]*fr[15]+alpha[7]*fr[14]+alpha[3]*fr[13])-1.732050807568877*alpha[6]*fr[12]+3.0*(alpha[2]*fr[11]+alpha[5]*fr[10])-1.732050807568877*(alpha[7]*fr[9]+alpha[3]*fr[8])+3.0*(alpha[4]*fr[7]+alpha[0]*fr[6])-1.732050807568877*(alpha[2]*fr[5]+fr[4]*alpha[5]+fr[2]*alpha[4])+3.0*alpha[1]*fr[3]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[5]*fr[15]+alpha[3]*fr[14]+alpha[7]*fr[13])-1.732050807568877*alpha[5]*fr[12]+3.0*(alpha[1]*fr[11]+alpha[6]*fr[10])-1.732050807568877*(alpha[3]*fr[9]+alpha[7]*fr[8])+3.0*(alpha[0]*fr[7]+alpha[4]*fr[6])-1.732050807568877*(fr[4]*alpha[6]+alpha[1]*fr[5]+fr[1]*alpha[4])+3.0*alpha[2]*fr[3]-1.732050807568877*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.0*alpha[2]*fr[12]+1.732050807568877*(alpha[6]*fr[11]+alpha[1]*fr[10])-1.0*(alpha[4]*fr[9]+alpha[0]*fr[8])+alpha[7]*(1.732050807568877*fr[7]-1.0*fr[2])+1.732050807568877*alpha[3]*fr[6]-1.0*fr[5]*alpha[6]+1.732050807568877*fr[3]*alpha[5]-1.0*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.0*alpha[1]*fr[12]+1.732050807568877*(alpha[5]*fr[11]+alpha[2]*fr[10])-1.0*(alpha[0]*fr[9]+alpha[4]*fr[8])+1.732050807568877*alpha[3]*fr[7]+(1.732050807568877*fr[6]-1.0*fr[1])*alpha[7]+1.732050807568877*fr[3]*alpha[6]-1.0*(fr[0]*alpha[6]+alpha[5]*fr[5]+alpha[2]*fr[4]+fr[2]*alpha[3]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.732050807568877*alpha[4]*fr[12]+3.0*(alpha[7]*fr[11]+alpha[0]*fr[10])-1.732050807568877*(alpha[2]*fr[9]+alpha[1]*fr[8])+3.0*alpha[6]*fr[7]-1.732050807568877*fr[5]*alpha[7]+3.0*alpha[5]*fr[6]-1.732050807568877*(fr[2]*alpha[6]+fr[1]*alpha[5]+alpha[0]*fr[4])+alpha[3]*(3.0*fr[3]-1.732050807568877*fr[0]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[3]*fr[15]+alpha[5]*fr[14]+alpha[6]*fr[13])-1.732050807568877*alpha[3]*fr[12]+3.0*(alpha[0]*fr[11]+alpha[7]*fr[10])-1.732050807568877*(alpha[5]*fr[9]+alpha[6]*fr[8])+3.0*alpha[1]*fr[7]-1.732050807568877*fr[4]*alpha[7]+3.0*alpha[2]*fr[6]-1.732050807568877*alpha[0]*fr[5]+3.0*fr[3]*alpha[4]-1.732050807568877*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.0*alpha[0]*fr[12]+1.732050807568877*(alpha[3]*fr[11]+alpha[4]*fr[10])-1.0*(alpha[1]*fr[9]+alpha[2]*fr[8])+1.732050807568877*alpha[5]*fr[7]+(1.732050807568877*fr[3]-1.0*fr[0])*alpha[7]+1.732050807568877*alpha[6]*fr[6]-1.0*(fr[1]*alpha[6]+alpha[3]*fr[5]+fr[2]*alpha[5]+alpha[4]*fr[4]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.732050807568877*alpha[2]*fr[12]+3.0*(alpha[6]*fr[11]+alpha[1]*fr[10])-1.732050807568877*(alpha[4]*fr[9]+alpha[0]*fr[8])+alpha[7]*(3.0*fr[7]-1.732050807568877*fr[2])+3.0*alpha[3]*fr[6]-1.732050807568877*fr[5]*alpha[6]+3.0*fr[3]*alpha[5]-1.732050807568877*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.732050807568877*alpha[1]*fr[12]+3.0*(alpha[5]*fr[11]+alpha[2]*fr[10])-1.732050807568877*(alpha[0]*fr[9]+alpha[4]*fr[8])+3.0*alpha[3]*fr[7]+(3.0*fr[6]-1.732050807568877*fr[1])*alpha[7]+3.0*fr[3]*alpha[6]-1.732050807568877*(fr[0]*alpha[6]+alpha[5]*fr[5]+alpha[2]*fr[4]+fr[2]*alpha[3]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.732050807568877*alpha[0]*fr[12]+3.0*(alpha[3]*fr[11]+alpha[4]*fr[10])-1.732050807568877*(alpha[1]*fr[9]+alpha[2]*fr[8])+3.0*alpha[5]*fr[7]+(3.0*fr[3]-1.732050807568877*fr[0])*alpha[7]+3.0*alpha[6]*fr[6]-1.732050807568877*(fr[1]*alpha[6]+alpha[3]*fr[5]+fr[2]*alpha[5]+alpha[4]*fr[4]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[16];
double fupwindQuad[16];
double limQuad[16];
  if((-0.3535533905932737*alphaUp[7])+0.3535533905932737*(alphaUp[6]+alphaUp[5]+alphaUp[4])-0.3535533905932737*(alphaUp[3]+alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[0] = (-0.25*fl[15])+0.25*(fl[14]+fl[13])-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*fl[5]-0.25*fl[4]+0.25*fl[3]-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[0] = 0.25*fr[15]-0.25*(fr[14]+fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.25*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[7]+alphaUp[6])-0.3535533905932737*(alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2])+0.3535533905932737*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[1] = 0.25*(fl[15]+fl[14])-0.25*fl[13]+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*(fl[5]+fl[4])+0.25*fl[3]-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[1] = (-0.25*(fr[15]+fr[14]))+0.25*(fr[13]+fr[12]+fr[11]+fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5]+fr[4]+fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[7]-0.3535533905932737*alphaUp[6]+0.3535533905932737*alphaUp[5]-0.3535533905932737*(alphaUp[4]+alphaUp[3])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[2] = 0.25*fl[15]-0.25*fl[14]+0.25*(fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5]+fl[4])+0.25*(fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[2] = (-0.25*fr[15])+0.25*fr[14]-0.25*fr[13]+0.25*(fr[12]+fr[11]+fr[10])-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*(fr[5]+fr[4]+fr[3])+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if((-0.3535533905932737*(alphaUp[7]+alphaUp[6]+alphaUp[5]))+0.3535533905932737*alphaUp[4]-0.3535533905932737*alphaUp[3]+0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[3] = (-0.25*(fl[15]+fl[14]+fl[13]+fl[12]))+0.25*fl[11]-0.25*(fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5])-0.25*fl[4]+0.25*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[3] = 0.25*(fr[15]+fr[14]+fr[13])-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*fr[5]-0.25*(fr[4]+fr[3])+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[7])+0.3535533905932737*(alphaUp[6]+alphaUp[5]+alphaUp[4])-0.3535533905932737*(alphaUp[3]+alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[4] = 0.25*fl[15]-0.25*(fl[14]+fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])-0.25*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[4] = (-0.25*fr[15])+0.25*(fr[14]+fr[13])-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*fr[5]-0.25*fr[4]+0.25*fr[3]-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[7]+alphaUp[6])-0.3535533905932737*(alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2])+0.3535533905932737*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[5] = (-0.25*(fl[15]+fl[14]))+0.25*(fl[13]+fl[12]+fl[11]+fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5]+fl[4]+fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[5] = 0.25*(fr[15]+fr[14])-0.25*fr[13]+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*(fr[5]+fr[4])+0.25*fr[3]-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[7]-0.3535533905932737*alphaUp[6]+0.3535533905932737*alphaUp[5]-0.3535533905932737*(alphaUp[4]+alphaUp[3])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[6] = (-0.25*fl[15])+0.25*fl[14]-0.25*fl[13]+0.25*(fl[12]+fl[11]+fl[10])-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*(fl[5]+fl[4]+fl[3])+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[6] = 0.25*fr[15]-0.25*fr[14]+0.25*(fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5]+fr[4])+0.25*(fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if((-0.3535533905932737*(alphaUp[7]+alphaUp[6]+alphaUp[5]))+0.3535533905932737*alphaUp[4]-0.3535533905932737*alphaUp[3]+0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[7] = 0.25*(fl[15]+fl[14]+fl[13])-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*fl[5]-0.25*(fl[4]+fl[3])+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[7] = (-0.25*(fr[15]+fr[14]+fr[13]+fr[12]))+0.25*fr[11]-0.25*(fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5])-0.25*fr[4]+0.25*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[7]-0.3535533905932737*(alphaUp[6]+alphaUp[5])+0.3535533905932737*(alphaUp[4]+alphaUp[3])-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[8] = 0.25*fl[15]-0.25*(fl[14]+fl[13])+0.25*(fl[12]+fl[11]+fl[10])-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*(fl[5]+fl[4]+fl[3])-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[8] = (-0.25*fr[15])+0.25*(fr[14]+fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5]+fr[4])-0.25*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if((-0.3535533905932737*(alphaUp[7]+alphaUp[6]))+0.3535533905932737*alphaUp[5]-0.3535533905932737*alphaUp[4]+0.3535533905932737*alphaUp[3]-0.3535533905932737*alphaUp[2]+0.3535533905932737*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[9] = (-0.25*(fl[15]+fl[14]))+0.25*fl[13]-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*fl[5]+0.25*(fl[4]+fl[3])-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[9] = 0.25*(fr[15]+fr[14])-0.25*(fr[13]+fr[12])+0.25*fr[11]-0.25*(fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5])+0.25*fr[4]-0.25*(fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[7])+0.3535533905932737*alphaUp[6]-0.3535533905932737*(alphaUp[5]+alphaUp[4])+0.3535533905932737*(alphaUp[3]+alphaUp[2])-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[10] = (-0.25*fl[15])+0.25*fl[14]-0.25*(fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5])+0.25*(fl[4]+fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[10] = 0.25*fr[15]-0.25*fr[14]+0.25*fr[13]-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*fr[5]+0.25*fr[4]-0.25*fr[3]+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[7]+alphaUp[6]+alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[11] = 0.25*(fl[15]+fl[14]+fl[13]+fl[12]+fl[11]+fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5]+fl[4]+fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[11] = (-0.25*(fr[15]+fr[14]+fr[13]))+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*(fr[5]+fr[4])-0.25*fr[3]+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[7]-0.3535533905932737*(alphaUp[6]+alphaUp[5])+0.3535533905932737*(alphaUp[4]+alphaUp[3])-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[12] = (-0.25*fl[15])+0.25*(fl[14]+fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5]+fl[4])-0.25*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[12] = 0.25*fr[15]-0.25*(fr[14]+fr[13])+0.25*(fr[12]+fr[11]+fr[10])-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*(fr[5]+fr[4]+fr[3])-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if((-0.3535533905932737*(alphaUp[7]+alphaUp[6]))+0.3535533905932737*alphaUp[5]-0.3535533905932737*alphaUp[4]+0.3535533905932737*alphaUp[3]-0.3535533905932737*alphaUp[2]+0.3535533905932737*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[13] = 0.25*(fl[15]+fl[14])-0.25*(fl[13]+fl[12])+0.25*fl[11]-0.25*(fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5])+0.25*fl[4]-0.25*(fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[13] = (-0.25*(fr[15]+fr[14]))+0.25*fr[13]-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*fr[5]+0.25*(fr[4]+fr[3])-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[7])+0.3535533905932737*alphaUp[6]-0.3535533905932737*(alphaUp[5]+alphaUp[4])+0.3535533905932737*(alphaUp[3]+alphaUp[2])-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[14] = 0.25*fl[15]-0.25*fl[14]+0.25*fl[13]-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*fl[5]+0.25*fl[4]-0.25*fl[3]+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[14] = (-0.25*fr[15])+0.25*fr[14]-0.25*(fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5])+0.25*(fr[4]+fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[7]+alphaUp[6]+alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[15] = (-0.25*(fl[15]+fl[14]+fl[13]))+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*(fl[5]+fl[4])-0.25*fl[3]+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[15] = 0.25*(fr[15]+fr[14]+fr[13]+fr[12]+fr[11]+fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5]+fr[4]+fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[4] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[5] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  fupwind[6] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[7] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[8] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[9] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[10] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[11] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[12] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[13] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[14] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[15] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  emModR[0] += -0.5*(1.732050807568877*fupwind[3]-1.0*fupwind[0])*dfac_v; 
  emModR[1] += -0.5*(1.732050807568877*fupwind[6]-1.0*fupwind[1])*dfac_v; 
  emModR[2] += -0.5*(1.732050807568877*fupwind[7]-1.0*fupwind[2])*dfac_v; 
  emModR[3] += 0.5*(3.0*fupwind[3]-1.732050807568877*fupwind[0])*dfac_v; 
  emModR[4] += -0.5*(1.732050807568877*fupwind[10]-1.0*fupwind[4])*dfac_v; 
  emModR[5] += -0.5*(1.732050807568877*fupwind[11]-1.0*fupwind[5])*dfac_v; 
  emModR[6] += 0.5*(3.0*fupwind[6]-1.732050807568877*fupwind[1])*dfac_v; 
  emModR[7] += 0.5*(3.0*fupwind[7]-1.732050807568877*fupwind[2])*dfac_v; 
  emModR[8] += -0.5*(1.732050807568877*fupwind[13]-1.0*fupwind[8])*dfac_v; 
  emModR[9] += -0.5*(1.732050807568877*fupwind[14]-1.0*fupwind[9])*dfac_v; 
  emModR[10] += 0.5*(3.0*fupwind[10]-1.732050807568877*fupwind[4])*dfac_v; 
  emModR[11] += 0.5*(3.0*fupwind[11]-1.732050807568877*fupwind[5])*dfac_v; 
  emModR[12] += -0.5*(1.732050807568877*fupwind[15]-1.0*fupwind[12])*dfac_v; 
  emModR[13] += 0.5*(3.0*fupwind[13]-1.732050807568877*fupwind[8])*dfac_v; 
  emModR[14] += 0.5*(3.0*fupwind[14]-1.732050807568877*fupwind[9])*dfac_v; 
  emModR[15] += 0.5*(3.0*fupwind[15]-1.732050807568877*fupwind[12])*dfac_v; 
  emModL[0] += 0.5*(1.732050807568877*fupwind[3]-1.0*fupwind[0])*dfac_v; 
  emModL[1] += 0.5*(1.732050807568877*fupwind[6]-1.0*fupwind[1])*dfac_v; 
  emModL[2] += 0.5*(1.732050807568877*fupwind[7]-1.0*fupwind[2])*dfac_v; 
  emModL[3] += 0.5*(3.0*fupwind[3]-1.732050807568877*fupwind[0])*dfac_v; 
  emModL[4] += 0.5*(1.732050807568877*fupwind[10]-1.0*fupwind[4])*dfac_v; 
  emModL[5] += 0.5*(1.732050807568877*fupwind[11]-1.0*fupwind[5])*dfac_v; 
  emModL[6] += 0.5*(3.0*fupwind[6]-1.732050807568877*fupwind[1])*dfac_v; 
  emModL[7] += 0.5*(3.0*fupwind[7]-1.732050807568877*fupwind[2])*dfac_v; 
  emModL[8] += 0.5*(1.732050807568877*fupwind[13]-1.0*fupwind[8])*dfac_v; 
  emModL[9] += 0.5*(1.732050807568877*fupwind[14]-1.0*fupwind[9])*dfac_v; 
  emModL[10] += 0.5*(3.0*fupwind[10]-1.732050807568877*fupwind[4])*dfac_v; 
  emModL[11] += 0.5*(3.0*fupwind[11]-1.732050807568877*fupwind[5])*dfac_v; 
  emModL[12] += 0.5*(1.732050807568877*fupwind[15]-1.0*fupwind[12])*dfac_v; 
  emModL[13] += 0.5*(3.0*fupwind[13]-1.732050807568877*fupwind[8])*dfac_v; 
  emModL[14] += 0.5*(3.0*fupwind[14]-1.732050807568877*fupwind[9])*dfac_v; 
  emModL[15] += 0.5*(3.0*fupwind[15]-1.732050807568877*fupwind[12])*dfac_v; 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[7]*fupwind[15]+alpha[6]*fupwind[14]+alpha[5]*fupwind[13])-1.0*alpha[7]*fupwind[12]+1.732050807568877*(alpha[4]*fupwind[11]+alpha[3]*fupwind[10])-1.0*(alpha[6]*fupwind[9]+alpha[5]*fupwind[8])+1.732050807568877*(alpha[2]*fupwind[7]+alpha[1]*fupwind[6])-1.0*(alpha[4]*fupwind[5]+alpha[3]*fupwind[4])+1.732050807568877*alpha[0]*fupwind[3]-1.0*(alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[6]*fupwind[15]+alpha[7]*fupwind[14]+alpha[3]*fupwind[13])-1.0*alpha[6]*fupwind[12]+1.732050807568877*(alpha[2]*fupwind[11]+alpha[5]*fupwind[10])-1.0*(alpha[7]*fupwind[9]+alpha[3]*fupwind[8])+1.732050807568877*(alpha[4]*fupwind[7]+alpha[0]*fupwind[6])-1.0*(alpha[2]*fupwind[5]+fupwind[4]*alpha[5]+fupwind[2]*alpha[4])+1.732050807568877*alpha[1]*fupwind[3]-1.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fupwind[15]+alpha[3]*fupwind[14]+alpha[7]*fupwind[13])-1.0*alpha[5]*fupwind[12]+1.732050807568877*(alpha[1]*fupwind[11]+alpha[6]*fupwind[10])-1.0*(alpha[3]*fupwind[9]+alpha[7]*fupwind[8])+1.732050807568877*(alpha[0]*fupwind[7]+alpha[4]*fupwind[6])-1.0*(fupwind[4]*alpha[6]+alpha[1]*fupwind[5]+fupwind[1]*alpha[4])+1.732050807568877*alpha[2]*fupwind[3]-1.0*(alpha[0]*fupwind[2]+fupwind[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[7]*fupwind[15]+alpha[6]*fupwind[14]+alpha[5]*fupwind[13])-1.732050807568877*alpha[7]*fupwind[12]+3.0*(alpha[4]*fupwind[11]+alpha[3]*fupwind[10])-1.732050807568877*(alpha[6]*fupwind[9]+alpha[5]*fupwind[8])+3.0*(alpha[2]*fupwind[7]+alpha[1]*fupwind[6])-1.732050807568877*(alpha[4]*fupwind[5]+alpha[3]*fupwind[4])+3.0*alpha[0]*fupwind[3]-1.732050807568877*(alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fupwind[15]+alpha[2]*fupwind[14]+alpha[1]*fupwind[13])-1.0*alpha[4]*fupwind[12]+1.732050807568877*(alpha[7]*fupwind[11]+alpha[0]*fupwind[10])-1.0*(alpha[2]*fupwind[9]+alpha[1]*fupwind[8])+1.732050807568877*alpha[6]*fupwind[7]-1.0*fupwind[5]*alpha[7]+1.732050807568877*alpha[5]*fupwind[6]-1.0*(fupwind[2]*alpha[6]+fupwind[1]*alpha[5]+alpha[0]*fupwind[4])+alpha[3]*(1.732050807568877*fupwind[3]-1.0*fupwind[0]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fupwind[15]+alpha[5]*fupwind[14]+alpha[6]*fupwind[13])-1.0*alpha[3]*fupwind[12]+1.732050807568877*(alpha[0]*fupwind[11]+alpha[7]*fupwind[10])-1.0*(alpha[5]*fupwind[9]+alpha[6]*fupwind[8])+1.732050807568877*alpha[1]*fupwind[7]-1.0*fupwind[4]*alpha[7]+1.732050807568877*alpha[2]*fupwind[6]-1.0*alpha[0]*fupwind[5]+1.732050807568877*fupwind[3]*alpha[4]-1.0*(fupwind[0]*alpha[4]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[6]*fupwind[15]+alpha[7]*fupwind[14]+alpha[3]*fupwind[13])-1.732050807568877*alpha[6]*fupwind[12]+3.0*(alpha[2]*fupwind[11]+alpha[5]*fupwind[10])-1.732050807568877*(alpha[7]*fupwind[9]+alpha[3]*fupwind[8])+3.0*(alpha[4]*fupwind[7]+alpha[0]*fupwind[6])-1.732050807568877*(alpha[2]*fupwind[5]+fupwind[4]*alpha[5]+fupwind[2]*alpha[4])+3.0*alpha[1]*fupwind[3]-1.732050807568877*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[5]*fupwind[15]+alpha[3]*fupwind[14]+alpha[7]*fupwind[13])-1.732050807568877*alpha[5]*fupwind[12]+3.0*(alpha[1]*fupwind[11]+alpha[6]*fupwind[10])-1.732050807568877*(alpha[3]*fupwind[9]+alpha[7]*fupwind[8])+3.0*(alpha[0]*fupwind[7]+alpha[4]*fupwind[6])-1.732050807568877*(fupwind[4]*alpha[6]+alpha[1]*fupwind[5]+fupwind[1]*alpha[4])+3.0*alpha[2]*fupwind[3]-1.732050807568877*(alpha[0]*fupwind[2]+fupwind[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[15]+alpha[4]*fupwind[14]+alpha[0]*fupwind[13])-1.0*alpha[2]*fupwind[12]+1.732050807568877*(alpha[6]*fupwind[11]+alpha[1]*fupwind[10])-1.0*(alpha[4]*fupwind[9]+alpha[0]*fupwind[8])+alpha[7]*(1.732050807568877*fupwind[7]-1.0*fupwind[2])+1.732050807568877*alpha[3]*fupwind[6]-1.0*fupwind[5]*alpha[6]+1.732050807568877*fupwind[3]*alpha[5]-1.0*(fupwind[0]*alpha[5]+alpha[1]*fupwind[4]+fupwind[1]*alpha[3]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14]+alpha[4]*fupwind[13])-1.0*alpha[1]*fupwind[12]+1.732050807568877*(alpha[5]*fupwind[11]+alpha[2]*fupwind[10])-1.0*(alpha[0]*fupwind[9]+alpha[4]*fupwind[8])+1.732050807568877*alpha[3]*fupwind[7]+(1.732050807568877*fupwind[6]-1.0*fupwind[1])*alpha[7]+1.732050807568877*fupwind[3]*alpha[6]-1.0*(fupwind[0]*alpha[6]+alpha[5]*fupwind[5]+alpha[2]*fupwind[4]+fupwind[2]*alpha[3]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[4]*fupwind[15]+alpha[2]*fupwind[14]+alpha[1]*fupwind[13])-1.732050807568877*alpha[4]*fupwind[12]+3.0*(alpha[7]*fupwind[11]+alpha[0]*fupwind[10])-1.732050807568877*(alpha[2]*fupwind[9]+alpha[1]*fupwind[8])+3.0*alpha[6]*fupwind[7]-1.732050807568877*fupwind[5]*alpha[7]+3.0*alpha[5]*fupwind[6]-1.732050807568877*(fupwind[2]*alpha[6]+fupwind[1]*alpha[5]+alpha[0]*fupwind[4])+alpha[3]*(3.0*fupwind[3]-1.732050807568877*fupwind[0]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[3]*fupwind[15]+alpha[5]*fupwind[14]+alpha[6]*fupwind[13])-1.732050807568877*alpha[3]*fupwind[12]+3.0*(alpha[0]*fupwind[11]+alpha[7]*fupwind[10])-1.732050807568877*(alpha[5]*fupwind[9]+alpha[6]*fupwind[8])+3.0*alpha[1]*fupwind[7]-1.732050807568877*fupwind[4]*alpha[7]+3.0*alpha[2]*fupwind[6]-1.732050807568877*alpha[0]*fupwind[5]+3.0*fupwind[3]*alpha[4]-1.732050807568877*(fupwind[0]*alpha[4]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14]+alpha[2]*fupwind[13])-1.0*alpha[0]*fupwind[12]+1.732050807568877*(alpha[3]*fupwind[11]+alpha[4]*fupwind[10])-1.0*(alpha[1]*fupwind[9]+alpha[2]*fupwind[8])+1.732050807568877*alpha[5]*fupwind[7]+(1.732050807568877*fupwind[3]-1.0*fupwind[0])*alpha[7]+1.732050807568877*alpha[6]*fupwind[6]-1.0*(fupwind[1]*alpha[6]+alpha[3]*fupwind[5]+fupwind[2]*alpha[5]+alpha[4]*fupwind[4]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[15]+alpha[4]*fupwind[14]+alpha[0]*fupwind[13])-1.732050807568877*alpha[2]*fupwind[12]+3.0*(alpha[6]*fupwind[11]+alpha[1]*fupwind[10])-1.732050807568877*(alpha[4]*fupwind[9]+alpha[0]*fupwind[8])+alpha[7]*(3.0*fupwind[7]-1.732050807568877*fupwind[2])+3.0*alpha[3]*fupwind[6]-1.732050807568877*fupwind[5]*alpha[6]+3.0*fupwind[3]*alpha[5]-1.732050807568877*(fupwind[0]*alpha[5]+alpha[1]*fupwind[4]+fupwind[1]*alpha[3]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14]+alpha[4]*fupwind[13])-1.732050807568877*alpha[1]*fupwind[12]+3.0*(alpha[5]*fupwind[11]+alpha[2]*fupwind[10])-1.732050807568877*(alpha[0]*fupwind[9]+alpha[4]*fupwind[8])+3.0*alpha[3]*fupwind[7]+(3.0*fupwind[6]-1.732050807568877*fupwind[1])*alpha[7]+3.0*fupwind[3]*alpha[6]-1.732050807568877*(fupwind[0]*alpha[6]+alpha[5]*fupwind[5]+alpha[2]*fupwind[4]+fupwind[2]*alpha[3]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14]+alpha[2]*fupwind[13])-1.732050807568877*alpha[0]*fupwind[12]+3.0*(alpha[3]*fupwind[11]+alpha[4]*fupwind[10])-1.732050807568877*(alpha[1]*fupwind[9]+alpha[2]*fupwind[8])+3.0*alpha[5]*fupwind[7]+(3.0*fupwind[3]-1.732050807568877*fupwind[0])*alpha[7]+3.0*alpha[6]*fupwind[6]-1.732050807568877*(fupwind[1]*alpha[6]+alpha[3]*fupwind[5]+fupwind[2]*alpha[5]+alpha[4]*fupwind[4]))*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSerStep2_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.0625*(1.732050807568877*(2.0*dfac_v*m_*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*wv+Bmag[1]*dfac_x*((Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*dfac_v*q_-2.0*BdriftX[0]*m_)*wm+q_*(dfac_v*(((Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[3]+(Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*Phi[2])*dfac_y+((BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*Phi[3]+(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[1])*dfac_x)*q_-2.0*((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_))+dfac_v*q_*(3.0*Bmag[1]*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_x*dfac_y*wm+(4.0*dApardtPrev[0]-3.0*(BmagInv[1]*(Apar[1]*Phi[3]-1.0*Phi[1]*Apar[3])+BmagInv[0]*(Apar[1]*Phi[2]-1.0*Phi[1]*Apar[2]))*dfac_x*dfac_y)*q_)))/(dfac_v*m_*q_); 

  double alpha[8]; 
  alpha[0] = -(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 
  alpha[2] = -(1.414213562373095*dApardt[2]*q_)/m_; 
  alpha[4] = -(1.414213562373095*dApardt[3]*q_)/m_; 
  double alphaUp[8]; 
  alphaUp[0] = -(0.3535533905932737*(1.732050807568877*(2.0*dfac_v*m_*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*wv+Bmag[1]*dfac_x*((Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*dfac_v*q_-2.0*BdriftX[0]*m_)*wm-2.0*((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*q_+dfac_v*(((Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[3]+(Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*Phi[2])*dfac_y+((BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*Phi[3]+(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[1])*dfac_x)*q2)+dfac_v*(3.0*Bmag[1]*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_x*dfac_y*q_*wm+((BmagInv[1]*(3.0*Phi[1]*Apar[3]-3.0*Apar[1]*Phi[3])+BmagInv[0]*(3.0*Phi[1]*Apar[2]-3.0*Apar[1]*Phi[2]))*dfac_x*dfac_y+4.0*dApardtPrev[0])*q2)))/(dfac_v*m_*q_); 
  alphaUp[1] = -(0.07071067811865474*(1.732050807568877*(10.0*dfac_v*m_*(BdriftX[1]*Bmag[1]*dfac_x*wm+((BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*q_)*wv+Bmag[1]*dfac_x*(5.0*(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*dfac_v*q_-10.0*BdriftX[1]*m_)*wm-10.0*((BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*m_*q_+dfac_v*(((9.0*Apar[1]*BdriftY[1]+5.0*Apar[0]*BdriftY[0])*Phi[3]+5.0*(Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[2])*dfac_y+5.0*((BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*Phi[3]+(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*Phi[1])*dfac_x)*q2)+dfac_v*(15.0*Bmag[1]*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_x*dfac_y*q_*wm+((BmagInv[0]*(15.0*Phi[1]*Apar[3]-15.0*Apar[1]*Phi[3])+BmagInv[1]*(15.0*Phi[1]*Apar[2]-15.0*Apar[1]*Phi[2]))*dfac_x*dfac_y+20.0*dApardtPrev[1])*q2)))/(dfac_v*m_*q_); 
  alphaUp[2] = -(0.3535533905932737*(1.732050807568877*(dfac_v*(dfac_x*(2.0*BdriftX[0]*Phi[3]*m_*wv+Bmag[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*wm)+(((BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2])*Phi[3]+Phi[2]*(BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2]))*dfac_y+((Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[3]+Phi[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2]))*dfac_x)*q_)-2.0*BdriftX[0]*Phi[3]*dfac_x*m_)+dfac_v*(BmagInv[0]*(3.0*Apar[2]*Phi[3]-3.0*Phi[2]*Apar[3])*dfac_x*dfac_y+4.0*dApardtPrev[2])*q_))/(dfac_v*m_); 
  alphaUp[3] = -(0.3535533905932737*Bmag[1]*dfac_x*(dfac_v*(2.0*BdriftX[0]*m_*wv+(1.732050807568877*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_y+Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*q_)-2.0*BdriftX[0]*m_))/(dfac_m*dfac_v*m_*q_); 
  alphaUp[4] = -(0.07071067811865474*(1.732050807568877*(dfac_v*(dfac_x*(10.0*BdriftX[1]*Phi[3]*m_*wv+5.0*Bmag[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*wm)+(((9.0*BdriftY[1]*Apar[3]+5.0*BdriftY[0]*Apar[2])*Phi[3]+5.0*Phi[2]*(BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2]))*dfac_y+5.0*((Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*Phi[3]+Phi[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2]))*dfac_x)*q_)-10.0*BdriftX[1]*Phi[3]*dfac_x*m_)+dfac_v*(BmagInv[1]*(15.0*Apar[2]*Phi[3]-15.0*Phi[2]*Apar[3])*dfac_x*dfac_y+20.0*dApardtPrev[3])*q_))/(dfac_v*m_); 
  alphaUp[5] = -(0.3535533905932737*Bmag[1]*dfac_x*(dfac_v*(2.0*BdriftX[1]*m_*wv+(1.732050807568877*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_y+Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*q_)-2.0*BdriftX[1]*m_))/(dfac_m*dfac_v*m_*q_); 
  alphaUp[6] = -(0.3535533905932737*Bmag[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*dfac_x)/(dfac_m*m_); 
  alphaUp[7] = -(0.3535533905932737*Bmag[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*dfac_x)/(dfac_m*m_); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[11]+alpha[2]*fl[7]+alpha[1]*fl[6])+alpha[4]*fl[5]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[11]+alpha[4]*fl[7]+alpha[0]*fl[6])+alpha[2]*fl[5]+fl[2]*alpha[4]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[11]+alpha[0]*fl[7]+alpha[4]*fl[6])+alpha[1]*fl[5]+fl[1]*alpha[4]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.1767766952966368*(3.0*(alpha[4]*fl[11]+alpha[2]*fl[7]+alpha[1]*fl[6])+1.732050807568877*alpha[4]*fl[5]+3.0*alpha[0]*fl[3]+1.732050807568877*(alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+alpha[4]*fl[12]+1.732050807568877*alpha[0]*fl[10]+alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[11]+alpha[1]*fl[7]+alpha[2]*fl[6])+alpha[0]*fl[5]+(1.732050807568877*fl[3]+fl[0])*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*(alpha[2]*fl[11]+alpha[4]*fl[7]+alpha[0]*fl[6])+1.732050807568877*(alpha[2]*fl[5]+fl[2]*alpha[4])+3.0*alpha[1]*fl[3]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[1]*fl[11]+alpha[0]*fl[7]+alpha[4]*fl[6])+1.732050807568877*(alpha[1]*fl[5]+fl[1]*alpha[4])+3.0*alpha[2]*fl[3]+1.732050807568877*(alpha[0]*fl[2]+fl[0]*alpha[2]))*dfac_v; 
  incr[8] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+alpha[2]*fl[12]+1.732050807568877*alpha[1]*fl[10]+alpha[4]*fl[9]+alpha[0]*fl[8]+alpha[1]*fl[4])*dfac_v; 
  incr[9] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+alpha[1]*fl[12]+1.732050807568877*alpha[2]*fl[10]+alpha[0]*fl[9]+alpha[4]*fl[8]+alpha[2]*fl[4])*dfac_v; 
  incr[10] = -0.1767766952966368*(3.0*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+1.732050807568877*alpha[4]*fl[12]+3.0*alpha[0]*fl[10]+1.732050807568877*(alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4]))*dfac_v; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[0]*fl[11]+alpha[1]*fl[7]+alpha[2]*fl[6])+1.732050807568877*alpha[0]*fl[5]+3.0*fl[3]*alpha[4]+1.732050807568877*(fl[0]*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2]))*dfac_v; 
  incr[12] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+alpha[0]*fl[12]+1.732050807568877*alpha[4]*fl[10]+alpha[1]*fl[9]+alpha[2]*fl[8]+alpha[4]*fl[4])*dfac_v; 
  incr[13] = -0.1767766952966368*(3.0*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+1.732050807568877*alpha[2]*fl[12]+3.0*alpha[1]*fl[10]+1.732050807568877*(alpha[4]*fl[9]+alpha[0]*fl[8]+alpha[1]*fl[4]))*dfac_v; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+1.732050807568877*alpha[1]*fl[12]+3.0*alpha[2]*fl[10]+1.732050807568877*(alpha[0]*fl[9]+alpha[4]*fl[8]+alpha[2]*fl[4]))*dfac_v; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+1.732050807568877*alpha[0]*fl[12]+3.0*alpha[4]*fl[10]+1.732050807568877*(alpha[1]*fl[9]+alpha[2]*fl[8]+alpha[4]*fl[4]))*dfac_v; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[11]+alpha[2]*fr[7]+alpha[1]*fr[6])-1.0*alpha[4]*fr[5]+1.732050807568877*alpha[0]*fr[3]-1.0*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[11]+alpha[4]*fr[7]+alpha[0]*fr[6])-1.0*(alpha[2]*fr[5]+fr[2]*alpha[4])+1.732050807568877*alpha[1]*fr[3]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[11]+alpha[0]*fr[7]+alpha[4]*fr[6])-1.0*(alpha[1]*fr[5]+fr[1]*alpha[4])+1.732050807568877*alpha[2]*fr[3]-1.0*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[4]*fr[11]+alpha[2]*fr[7]+alpha[1]*fr[6])-1.732050807568877*alpha[4]*fr[5]+3.0*alpha[0]*fr[3]-1.732050807568877*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.0*alpha[4]*fr[12]+1.732050807568877*alpha[0]*fr[10]-1.0*(alpha[2]*fr[9]+alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[11]+alpha[1]*fr[7]+alpha[2]*fr[6])-1.0*alpha[0]*fr[5]+1.732050807568877*fr[3]*alpha[4]-1.0*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[2]*fr[11]+alpha[4]*fr[7]+alpha[0]*fr[6])-1.732050807568877*(alpha[2]*fr[5]+fr[2]*alpha[4])+3.0*alpha[1]*fr[3]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[1]*fr[11]+alpha[0]*fr[7]+alpha[4]*fr[6])-1.732050807568877*(alpha[1]*fr[5]+fr[1]*alpha[4])+3.0*alpha[2]*fr[3]-1.732050807568877*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.0*alpha[2]*fr[12]+1.732050807568877*alpha[1]*fr[10]-1.0*(alpha[4]*fr[9]+alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.0*alpha[1]*fr[12]+1.732050807568877*alpha[2]*fr[10]-1.0*(alpha[0]*fr[9]+alpha[4]*fr[8]+alpha[2]*fr[4]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.732050807568877*alpha[4]*fr[12]+3.0*alpha[0]*fr[10]-1.732050807568877*(alpha[2]*fr[9]+alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[0]*fr[11]+alpha[1]*fr[7]+alpha[2]*fr[6])-1.732050807568877*alpha[0]*fr[5]+3.0*fr[3]*alpha[4]-1.732050807568877*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.0*alpha[0]*fr[12]+1.732050807568877*alpha[4]*fr[10]-1.0*(alpha[1]*fr[9]+alpha[2]*fr[8]+alpha[4]*fr[4]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.732050807568877*alpha[2]*fr[12]+3.0*alpha[1]*fr[10]-1.732050807568877*(alpha[4]*fr[9]+alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.732050807568877*alpha[1]*fr[12]+3.0*alpha[2]*fr[10]-1.732050807568877*(alpha[0]*fr[9]+alpha[4]*fr[8]+alpha[2]*fr[4]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.732050807568877*alpha[0]*fr[12]+3.0*alpha[4]*fr[10]-1.732050807568877*(alpha[1]*fr[9]+alpha[2]*fr[8]+alpha[4]*fr[4]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[16];
double fupwindQuad[16];
double limQuad[16];
  if((-0.3535533905932737*alphaUp[7])+0.3535533905932737*(alphaUp[6]+alphaUp[5]+alphaUp[4])-0.3535533905932737*(alphaUp[3]+alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[0] = (-0.25*fl[15])+0.25*(fl[14]+fl[13])-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*fl[5]-0.25*fl[4]+0.25*fl[3]-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[0] = 0.25*fr[15]-0.25*(fr[14]+fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.25*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[7]+alphaUp[6])-0.3535533905932737*(alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2])+0.3535533905932737*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[1] = 0.25*(fl[15]+fl[14])-0.25*fl[13]+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*(fl[5]+fl[4])+0.25*fl[3]-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[1] = (-0.25*(fr[15]+fr[14]))+0.25*(fr[13]+fr[12]+fr[11]+fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5]+fr[4]+fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[7]-0.3535533905932737*alphaUp[6]+0.3535533905932737*alphaUp[5]-0.3535533905932737*(alphaUp[4]+alphaUp[3])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[2] = 0.25*fl[15]-0.25*fl[14]+0.25*(fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5]+fl[4])+0.25*(fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[2] = (-0.25*fr[15])+0.25*fr[14]-0.25*fr[13]+0.25*(fr[12]+fr[11]+fr[10])-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*(fr[5]+fr[4]+fr[3])+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if((-0.3535533905932737*(alphaUp[7]+alphaUp[6]+alphaUp[5]))+0.3535533905932737*alphaUp[4]-0.3535533905932737*alphaUp[3]+0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[3] = (-0.25*(fl[15]+fl[14]+fl[13]+fl[12]))+0.25*fl[11]-0.25*(fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5])-0.25*fl[4]+0.25*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[3] = 0.25*(fr[15]+fr[14]+fr[13])-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*fr[5]-0.25*(fr[4]+fr[3])+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[7])+0.3535533905932737*(alphaUp[6]+alphaUp[5]+alphaUp[4])-0.3535533905932737*(alphaUp[3]+alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[4] = 0.25*fl[15]-0.25*(fl[14]+fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])-0.25*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[4] = (-0.25*fr[15])+0.25*(fr[14]+fr[13])-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*fr[5]-0.25*fr[4]+0.25*fr[3]-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[7]+alphaUp[6])-0.3535533905932737*(alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2])+0.3535533905932737*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[5] = (-0.25*(fl[15]+fl[14]))+0.25*(fl[13]+fl[12]+fl[11]+fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5]+fl[4]+fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[5] = 0.25*(fr[15]+fr[14])-0.25*fr[13]+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*(fr[5]+fr[4])+0.25*fr[3]-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[7]-0.3535533905932737*alphaUp[6]+0.3535533905932737*alphaUp[5]-0.3535533905932737*(alphaUp[4]+alphaUp[3])+0.3535533905932737*alphaUp[2]-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[6] = (-0.25*fl[15])+0.25*fl[14]-0.25*fl[13]+0.25*(fl[12]+fl[11]+fl[10])-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*(fl[5]+fl[4]+fl[3])+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[6] = 0.25*fr[15]-0.25*fr[14]+0.25*(fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5]+fr[4])+0.25*(fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if((-0.3535533905932737*(alphaUp[7]+alphaUp[6]+alphaUp[5]))+0.3535533905932737*alphaUp[4]-0.3535533905932737*alphaUp[3]+0.3535533905932737*(alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[7] = 0.25*(fl[15]+fl[14]+fl[13])-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*fl[5]-0.25*(fl[4]+fl[3])+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[7] = (-0.25*(fr[15]+fr[14]+fr[13]+fr[12]))+0.25*fr[11]-0.25*(fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5])-0.25*fr[4]+0.25*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[7]-0.3535533905932737*(alphaUp[6]+alphaUp[5])+0.3535533905932737*(alphaUp[4]+alphaUp[3])-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[8] = 0.25*fl[15]-0.25*(fl[14]+fl[13])+0.25*(fl[12]+fl[11]+fl[10])-0.25*(fl[9]+fl[8]+fl[7]+fl[6])+0.25*(fl[5]+fl[4]+fl[3])-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[8] = (-0.25*fr[15])+0.25*(fr[14]+fr[13]+fr[12])-0.25*(fr[11]+fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5]+fr[4])-0.25*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  }
  if((-0.3535533905932737*(alphaUp[7]+alphaUp[6]))+0.3535533905932737*alphaUp[5]-0.3535533905932737*alphaUp[4]+0.3535533905932737*alphaUp[3]-0.3535533905932737*alphaUp[2]+0.3535533905932737*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[9] = (-0.25*(fl[15]+fl[14]))+0.25*fl[13]-0.25*(fl[12]+fl[11])+0.25*fl[10]-0.25*fl[9]+0.25*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*fl[5]+0.25*(fl[4]+fl[3])-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[9] = 0.25*(fr[15]+fr[14])-0.25*(fr[13]+fr[12])+0.25*fr[11]-0.25*(fr[10]+fr[9])+0.25*(fr[8]+fr[7])-0.25*(fr[6]+fr[5])+0.25*fr[4]-0.25*(fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[7])+0.3535533905932737*alphaUp[6]-0.3535533905932737*(alphaUp[5]+alphaUp[4])+0.3535533905932737*(alphaUp[3]+alphaUp[2])-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[10] = (-0.25*fl[15])+0.25*fl[14]-0.25*(fl[13]+fl[12]+fl[11])+0.25*(fl[10]+fl[9])-0.25*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5])+0.25*(fl[4]+fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[10] = 0.25*fr[15]-0.25*fr[14]+0.25*fr[13]-0.25*fr[12]+0.25*fr[11]-0.25*fr[10]+0.25*fr[9]-0.25*(fr[8]+fr[7])+0.25*fr[6]-0.25*fr[5]+0.25*fr[4]-0.25*fr[3]+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[7]+alphaUp[6]+alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[11] = 0.25*(fl[15]+fl[14]+fl[13]+fl[12]+fl[11]+fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5]+fl[4]+fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[11] = (-0.25*(fr[15]+fr[14]+fr[13]))+0.25*fr[12]-0.25*(fr[11]+fr[10])+0.25*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*(fr[5]+fr[4])-0.25*fr[3]+0.25*(fr[2]+fr[1]+fr[0]); 
  }
  if(0.3535533905932737*alphaUp[7]-0.3535533905932737*(alphaUp[6]+alphaUp[5])+0.3535533905932737*(alphaUp[4]+alphaUp[3])-0.3535533905932737*(alphaUp[2]+alphaUp[1])+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[12] = (-0.25*fl[15])+0.25*(fl[14]+fl[13]+fl[12])-0.25*(fl[11]+fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5]+fl[4])-0.25*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else {
  fupwindQuad[12] = 0.25*fr[15]-0.25*(fr[14]+fr[13])+0.25*(fr[12]+fr[11]+fr[10])-0.25*(fr[9]+fr[8]+fr[7]+fr[6])+0.25*(fr[5]+fr[4]+fr[3])-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  }
  if((-0.3535533905932737*(alphaUp[7]+alphaUp[6]))+0.3535533905932737*alphaUp[5]-0.3535533905932737*alphaUp[4]+0.3535533905932737*alphaUp[3]-0.3535533905932737*alphaUp[2]+0.3535533905932737*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[13] = 0.25*(fl[15]+fl[14])-0.25*(fl[13]+fl[12])+0.25*fl[11]-0.25*(fl[10]+fl[9])+0.25*(fl[8]+fl[7])-0.25*(fl[6]+fl[5])+0.25*fl[4]-0.25*(fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else {
  fupwindQuad[13] = (-0.25*(fr[15]+fr[14]))+0.25*fr[13]-0.25*(fr[12]+fr[11])+0.25*fr[10]-0.25*fr[9]+0.25*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*fr[5]+0.25*(fr[4]+fr[3])-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  }
  if((-0.3535533905932737*alphaUp[7])+0.3535533905932737*alphaUp[6]-0.3535533905932737*(alphaUp[5]+alphaUp[4])+0.3535533905932737*(alphaUp[3]+alphaUp[2])-0.3535533905932737*alphaUp[1]+0.3535533905932737*alphaUp[0] > 0) {
  fupwindQuad[14] = 0.25*fl[15]-0.25*fl[14]+0.25*fl[13]-0.25*fl[12]+0.25*fl[11]-0.25*fl[10]+0.25*fl[9]-0.25*(fl[8]+fl[7])+0.25*fl[6]-0.25*fl[5]+0.25*fl[4]-0.25*fl[3]+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else {
  fupwindQuad[14] = (-0.25*fr[15])+0.25*fr[14]-0.25*(fr[13]+fr[12]+fr[11])+0.25*(fr[10]+fr[9])-0.25*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5])+0.25*(fr[4]+fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  }
  if(0.3535533905932737*(alphaUp[7]+alphaUp[6]+alphaUp[5]+alphaUp[4]+alphaUp[3]+alphaUp[2]+alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[15] = (-0.25*(fl[15]+fl[14]+fl[13]))+0.25*fl[12]-0.25*(fl[11]+fl[10])+0.25*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*(fl[5]+fl[4])-0.25*fl[3]+0.25*(fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[15] = 0.25*(fr[15]+fr[14]+fr[13]+fr[12]+fr[11]+fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5]+fr[4]+fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8])+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[4] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[5] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  fupwind[6] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]-1.0*fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[7] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]-1.0*(fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[8] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*fupwindQuad[12]+fupwindQuad[11]-1.0*fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]-1.0*fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[9] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12])+fupwindQuad[11]+fupwindQuad[10]-1.0*(fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]-1.0*(fupwindQuad[3]+fupwindQuad[2])+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[10] = 0.25*(fupwindQuad[15]+fupwindQuad[14]+fupwindQuad[13]+fupwindQuad[12]-1.0*(fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]+fupwindQuad[8]+fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]+fupwindQuad[4])+fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[11] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*fupwindQuad[8]+fupwindQuad[7]-1.0*(fupwindQuad[6]+fupwindQuad[5])+fupwindQuad[4]-1.0*fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[12] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]+fupwindQuad[11]-1.0*(fupwindQuad[10]+fupwindQuad[9])+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]+fupwindQuad[5]-1.0*(fupwindQuad[4]+fupwindQuad[3])+fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[13] = 0.25*(fupwindQuad[15]-1.0*fupwindQuad[14]+fupwindQuad[13]-1.0*(fupwindQuad[12]+fupwindQuad[11])+fupwindQuad[10]-1.0*fupwindQuad[9]+fupwindQuad[8]-1.0*fupwindQuad[7]+fupwindQuad[6]-1.0*fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[14] = 0.25*(fupwindQuad[15]+fupwindQuad[14]-1.0*(fupwindQuad[13]+fupwindQuad[12]+fupwindQuad[11]+fupwindQuad[10])+fupwindQuad[9]+fupwindQuad[8]-1.0*(fupwindQuad[7]+fupwindQuad[6])+fupwindQuad[5]+fupwindQuad[4]+fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[15] = 0.25*(fupwindQuad[15]-1.0*(fupwindQuad[14]+fupwindQuad[13])+fupwindQuad[12]-1.0*fupwindQuad[11]+fupwindQuad[10]+fupwindQuad[9]-1.0*(fupwindQuad[8]+fupwindQuad[7])+fupwindQuad[6]+fupwindQuad[5]-1.0*fupwindQuad[4]+fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fupwind[11]+alpha[2]*fupwind[7]+alpha[1]*fupwind[6])-1.0*alpha[4]*fupwind[5]+1.732050807568877*alpha[0]*fupwind[3]-1.0*(alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[11]+alpha[4]*fupwind[7]+alpha[0]*fupwind[6])-1.0*(alpha[2]*fupwind[5]+fupwind[2]*alpha[4])+1.732050807568877*alpha[1]*fupwind[3]-1.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fupwind[11]+alpha[0]*fupwind[7]+alpha[4]*fupwind[6])-1.0*(alpha[1]*fupwind[5]+fupwind[1]*alpha[4])+1.732050807568877*alpha[2]*fupwind[3]-1.0*(alpha[0]*fupwind[2]+fupwind[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[4]*fupwind[11]+alpha[2]*fupwind[7]+alpha[1]*fupwind[6])-1.732050807568877*alpha[4]*fupwind[5]+3.0*alpha[0]*fupwind[3]-1.732050807568877*(alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fupwind[15]+alpha[2]*fupwind[14]+alpha[1]*fupwind[13])-1.0*alpha[4]*fupwind[12]+1.732050807568877*alpha[0]*fupwind[10]-1.0*(alpha[2]*fupwind[9]+alpha[1]*fupwind[8]+alpha[0]*fupwind[4]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fupwind[11]+alpha[1]*fupwind[7]+alpha[2]*fupwind[6])-1.0*alpha[0]*fupwind[5]+1.732050807568877*fupwind[3]*alpha[4]-1.0*(fupwind[0]*alpha[4]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[11]+alpha[4]*fupwind[7]+alpha[0]*fupwind[6])-1.732050807568877*(alpha[2]*fupwind[5]+fupwind[2]*alpha[4])+3.0*alpha[1]*fupwind[3]-1.732050807568877*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[1]*fupwind[11]+alpha[0]*fupwind[7]+alpha[4]*fupwind[6])-1.732050807568877*(alpha[1]*fupwind[5]+fupwind[1]*alpha[4])+3.0*alpha[2]*fupwind[3]-1.732050807568877*(alpha[0]*fupwind[2]+fupwind[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fupwind[15]+alpha[4]*fupwind[14]+alpha[0]*fupwind[13])-1.0*alpha[2]*fupwind[12]+1.732050807568877*alpha[1]*fupwind[10]-1.0*(alpha[4]*fupwind[9]+alpha[0]*fupwind[8]+alpha[1]*fupwind[4]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14]+alpha[4]*fupwind[13])-1.0*alpha[1]*fupwind[12]+1.732050807568877*alpha[2]*fupwind[10]-1.0*(alpha[0]*fupwind[9]+alpha[4]*fupwind[8]+alpha[2]*fupwind[4]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[4]*fupwind[15]+alpha[2]*fupwind[14]+alpha[1]*fupwind[13])-1.732050807568877*alpha[4]*fupwind[12]+3.0*alpha[0]*fupwind[10]-1.732050807568877*(alpha[2]*fupwind[9]+alpha[1]*fupwind[8]+alpha[0]*fupwind[4]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[0]*fupwind[11]+alpha[1]*fupwind[7]+alpha[2]*fupwind[6])-1.732050807568877*alpha[0]*fupwind[5]+3.0*fupwind[3]*alpha[4]-1.732050807568877*(fupwind[0]*alpha[4]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14]+alpha[2]*fupwind[13])-1.0*alpha[0]*fupwind[12]+1.732050807568877*alpha[4]*fupwind[10]-1.0*(alpha[1]*fupwind[9]+alpha[2]*fupwind[8]+alpha[4]*fupwind[4]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fupwind[15]+alpha[4]*fupwind[14]+alpha[0]*fupwind[13])-1.732050807568877*alpha[2]*fupwind[12]+3.0*alpha[1]*fupwind[10]-1.732050807568877*(alpha[4]*fupwind[9]+alpha[0]*fupwind[8]+alpha[1]*fupwind[4]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fupwind[15]+alpha[0]*fupwind[14]+alpha[4]*fupwind[13])-1.732050807568877*alpha[1]*fupwind[12]+3.0*alpha[2]*fupwind[10]-1.732050807568877*(alpha[0]*fupwind[9]+alpha[4]*fupwind[8]+alpha[2]*fupwind[4]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fupwind[15]+alpha[1]*fupwind[14]+alpha[2]*fupwind[13])-1.732050807568877*alpha[0]*fupwind[12]+3.0*alpha[4]*fupwind[10]-1.732050807568877*(alpha[1]*fupwind[9]+alpha[2]*fupwind[8]+alpha[4]*fupwind[4]))*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
