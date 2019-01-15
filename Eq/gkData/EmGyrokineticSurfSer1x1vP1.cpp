#include <GyrokineticModDecl.h> 
#define upwindType SURFAVG
double EmGyrokineticSurf1x1vSer_X_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*Gradpar[0]*wv; 

  double alpha[2]; 
  alpha[0] = Gradpar[0]*wv; 
  double alphaUp[2]; 
  alphaUp[0] = Gradpar[0]*wv; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.3535533905932737*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*(1.732050807568877*fl[3]+fl[2])*dfac_x; 
  incr[3] = -0.3535533905932737*alpha[0]*(3.0*fl[3]+1.732050807568877*fl[2])*dfac_x; 
  } else { 
  incr[0] = -0.3535533905932737*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.3535533905932737*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.3535533905932737*alpha[0]*(1.732050807568877*fr[3]-1.0*fr[2])*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*(3.0*fr[3]-1.732050807568877*fr[2])*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double limQuad[4];
  if(0.7071067811865475*alphaUp[0] > 0) {
  fupwindQuad[0] = 0.5*(fl[1]+fl[0])-0.5*(fl[3]+fl[2]); 
  } else {
  fupwindQuad[0] = 0.5*fr[3]-0.5*(fr[2]+fr[1])+0.5*fr[0]; 
  }
  if(0.7071067811865475*alphaUp[0] > 0) {
  fupwindQuad[1] = 0.5*fl[3]-0.5*(fl[2]+fl[1])+0.5*fl[0]; 
  } else {
  fupwindQuad[1] = 0.5*(fr[1]+fr[0])-0.5*(fr[3]+fr[2]); 
  }
  if(0.7071067811865475*alphaUp[0] > 0) {
  fupwindQuad[2] = 0.5*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[2] = (-0.5*fr[3])+0.5*fr[2]-0.5*fr[1]+0.5*fr[0]; 
  }
  if(0.7071067811865475*alphaUp[0] > 0) {
  fupwindQuad[3] = (-0.5*fl[3])+0.5*fl[2]-0.5*fl[1]+0.5*fl[0]; 
  } else {
  fupwindQuad[3] = 0.5*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*(fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.5*(fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  incr[0] = -0.3535533905932737*alpha[0]*(1.732050807568877*fupwind[1]-1.0*fupwind[0])*dfac_x; 
  incr[1] = 0.3535533905932737*alpha[0]*(3.0*fupwind[1]-1.732050807568877*fupwind[0])*dfac_x; 
  incr[2] = -0.3535533905932737*alpha[0]*(1.732050807568877*fupwind[3]-1.0*fupwind[2])*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*(3.0*fupwind[3]-1.732050807568877*fupwind[2])*dfac_x; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x1vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  double incrEmMod[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 

  double alpha[2]; 
  alpha[0] = -(1.224744871391589*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  double alphaUp[2]; 
  alphaUp[0] = -(0.7071067811865475*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 
  alphaUp[1] = -(1.0*dApardtPrev[1]*q_)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  emModR[0] += 0.5*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  emModR[1] += 0.5*(1.732050807568877*fl[3]+fl[1])*dfac_v; 
  emModR[2] += -0.5*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  emModR[3] += -0.5*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_v; 
  emModL[0] += -0.5*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  emModL[1] += -0.5*(1.732050807568877*fl[3]+fl[1])*dfac_v; 
  emModL[2] += -0.5*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  emModL[3] += -0.5*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_v; 
  incr[0] = 0.3535533905932737*alpha[0]*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  incr[1] = 0.3535533905932737*alpha[0]*(1.732050807568877*fl[3]+fl[1])*dfac_v; 
  incr[2] = -0.3535533905932737*alpha[0]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  incr[3] = -0.3535533905932737*alpha[0]*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_v; 
  } else { 
  emModR[0] += -0.5*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  emModR[1] += -0.5*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_v; 
  emModR[2] += 0.5*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  emModR[3] += 0.5*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_v; 
  emModL[0] += 0.5*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  emModL[1] += 0.5*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_v; 
  emModL[2] += 0.5*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  emModL[3] += 0.5*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_v; 
  incr[0] = -0.3535533905932737*alpha[0]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  incr[1] = -0.3535533905932737*alpha[0]*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_v; 
  incr[2] = 0.3535533905932737*alpha[0]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  incr[3] = 0.3535533905932737*alpha[0]*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double limQuad[4];
  if(0.7071067811865475*alphaUp[0]-0.7071067811865475*alphaUp[1] > 0) {
  fupwindQuad[0] = (-0.5*fl[3])+0.5*fl[2]-0.5*fl[1]+0.5*fl[0]; 
  } else {
  fupwindQuad[0] = 0.5*fr[3]-0.5*(fr[2]+fr[1])+0.5*fr[0]; 
  }
  if(0.7071067811865475*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[1] = 0.5*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[1] = 0.5*(fr[1]+fr[0])-0.5*(fr[3]+fr[2]); 
  }
  if(0.7071067811865475*alphaUp[0]-0.7071067811865475*alphaUp[1] > 0) {
  fupwindQuad[2] = 0.5*fl[3]-0.5*(fl[2]+fl[1])+0.5*fl[0]; 
  } else {
  fupwindQuad[2] = (-0.5*fr[3])+0.5*fr[2]-0.5*fr[1]+0.5*fr[0]; 
  }
  if(0.7071067811865475*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[3] = 0.5*(fl[1]+fl[0])-0.5*(fl[3]+fl[2]); 
  } else {
  fupwindQuad[3] = 0.5*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*(fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.5*(fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  emModR[0] += -0.5*(1.732050807568877*fupwind[2]-1.0*fupwind[0])*dfac_v; 
  emModR[1] += -0.5*(1.732050807568877*fupwind[3]-1.0*fupwind[1])*dfac_v; 
  emModR[2] += 0.5*(3.0*fupwind[2]-1.732050807568877*fupwind[0])*dfac_v; 
  emModR[3] += 0.5*(3.0*fupwind[3]-1.732050807568877*fupwind[1])*dfac_v; 
  emModL[0] += 0.5*(1.732050807568877*fupwind[2]-1.0*fupwind[0])*dfac_v; 
  emModL[1] += 0.5*(1.732050807568877*fupwind[3]-1.0*fupwind[1])*dfac_v; 
  emModL[2] += 0.5*(3.0*fupwind[2]-1.732050807568877*fupwind[0])*dfac_v; 
  emModL[3] += 0.5*(3.0*fupwind[3]-1.732050807568877*fupwind[1])*dfac_v; 
  incr[0] = -0.3535533905932737*alpha[0]*(1.732050807568877*fupwind[2]-1.0*fupwind[0])*dfac_v; 
  incr[1] = -0.3535533905932737*alpha[0]*(1.732050807568877*fupwind[3]-1.0*fupwind[1])*dfac_v; 
  incr[2] = 0.3535533905932737*alpha[0]*(3.0*fupwind[2]-1.732050807568877*fupwind[0])*dfac_v; 
  incr[3] = 0.3535533905932737*alpha[0]*(3.0*fupwind[3]-1.732050807568877*fupwind[1])*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x1vSerStep2_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 

  double alpha[2]; 
  alpha[0] = -(1.0*dApardt[0]*q_)/m_; 
  alpha[1] = -(1.0*dApardt[1]*q_)/m_; 
  double alphaUp[2]; 
  alphaUp[0] = -(0.7071067811865475*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 
  alphaUp[1] = -(1.0*dApardtPrev[1]*q_)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*(1.732050807568877*(alpha[1]*fl[3]+alpha[0]*fl[2])+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(1.732050807568877*(alpha[0]*fl[3]+alpha[1]*fl[2])+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.3535533905932737*(3.0*(alpha[1]*fl[3]+alpha[0]*fl[2])+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[3] = -0.3535533905932737*(3.0*(alpha[0]*fl[3]+alpha[1]*fl[2])+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  } else { 
  incr[0] = -0.3535533905932737*(1.732050807568877*(alpha[1]*fr[3]+alpha[0]*fr[2])-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.3535533905932737*(1.732050807568877*(alpha[0]*fr[3]+alpha[1]*fr[2])-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = 0.3535533905932737*(3.0*(alpha[1]*fr[3]+alpha[0]*fr[2])-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[3] = 0.3535533905932737*(3.0*(alpha[0]*fr[3]+alpha[1]*fr[2])-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double limQuad[4];
  if(0.7071067811865475*alphaUp[0]-0.7071067811865475*alphaUp[1] > 0) {
  fupwindQuad[0] = (-0.5*fl[3])+0.5*fl[2]-0.5*fl[1]+0.5*fl[0]; 
  } else {
  fupwindQuad[0] = 0.5*fr[3]-0.5*(fr[2]+fr[1])+0.5*fr[0]; 
  }
  if(0.7071067811865475*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[1] = 0.5*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[1] = 0.5*(fr[1]+fr[0])-0.5*(fr[3]+fr[2]); 
  }
  if(0.7071067811865475*alphaUp[0]-0.7071067811865475*alphaUp[1] > 0) {
  fupwindQuad[2] = 0.5*fl[3]-0.5*(fl[2]+fl[1])+0.5*fl[0]; 
  } else {
  fupwindQuad[2] = (-0.5*fr[3])+0.5*fr[2]-0.5*fr[1]+0.5*fr[0]; 
  }
  if(0.7071067811865475*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[3] = 0.5*(fl[1]+fl[0])-0.5*(fl[3]+fl[2]); 
  } else {
  fupwindQuad[3] = 0.5*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*(fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.5*(fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  incr[0] = -0.3535533905932737*(1.732050807568877*(alpha[1]*fupwind[3]+alpha[0]*fupwind[2])-1.0*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[1] = -0.3535533905932737*(1.732050807568877*(alpha[0]*fupwind[3]+alpha[1]*fupwind[2])-1.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 
  incr[2] = 0.3535533905932737*(3.0*(alpha[1]*fupwind[3]+alpha[0]*fupwind[2])-1.732050807568877*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[3] = 0.3535533905932737*(3.0*(alpha[0]*fupwind[3]+alpha[1]*fupwind[2])-1.732050807568877*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x1vSer_X_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.3535533905932737*(1.732050807568877*Gradpar[1]-1.0*Gradpar[0])*wv; 

  double alpha[2]; 
  alpha[0] = (Gradpar[0]-1.732050807568877*Gradpar[1])*wv; 
  double alphaUp[2]; 
  alphaUp[0] = (Gradpar[0]-1.732050807568877*Gradpar[1])*wv; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.3535533905932737*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*(1.732050807568877*fl[3]+fl[2])*dfac_x; 
  incr[3] = -0.3535533905932737*alpha[0]*(3.0*fl[3]+1.732050807568877*fl[2])*dfac_x; 
  } else { 
  incr[0] = -0.3535533905932737*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.3535533905932737*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.3535533905932737*alpha[0]*(1.732050807568877*fr[3]-1.0*fr[2])*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*(3.0*fr[3]-1.732050807568877*fr[2])*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double limQuad[4];
  if(0.7071067811865475*alphaUp[0] > 0) {
  fupwindQuad[0] = 0.5*(fl[1]+fl[0])-0.5*(fl[3]+fl[2]); 
  } else {
  fupwindQuad[0] = 0.5*fr[3]-0.5*(fr[2]+fr[1])+0.5*fr[0]; 
  }
  if(0.7071067811865475*alphaUp[0] > 0) {
  fupwindQuad[1] = 0.5*fl[3]-0.5*(fl[2]+fl[1])+0.5*fl[0]; 
  } else {
  fupwindQuad[1] = 0.5*(fr[1]+fr[0])-0.5*(fr[3]+fr[2]); 
  }
  if(0.7071067811865475*alphaUp[0] > 0) {
  fupwindQuad[2] = 0.5*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[2] = (-0.5*fr[3])+0.5*fr[2]-0.5*fr[1]+0.5*fr[0]; 
  }
  if(0.7071067811865475*alphaUp[0] > 0) {
  fupwindQuad[3] = (-0.5*fl[3])+0.5*fl[2]-0.5*fl[1]+0.5*fl[0]; 
  } else {
  fupwindQuad[3] = 0.5*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*(fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.5*(fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  incr[0] = -0.3535533905932737*alpha[0]*(1.732050807568877*fupwind[1]-1.0*fupwind[0])*dfac_x; 
  incr[1] = 0.3535533905932737*alpha[0]*(3.0*fupwind[1]-1.732050807568877*fupwind[0])*dfac_x; 
  incr[2] = -0.3535533905932737*alpha[0]*(1.732050807568877*fupwind[3]-1.0*fupwind[2])*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*(3.0*fupwind[3]-1.732050807568877*fupwind[2])*dfac_x; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x1vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  double incrEmMod[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 

  double alpha[2]; 
  alpha[0] = -(1.224744871391589*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(1.224744871391589*Gradpar[1]*Phi[1]*dfac_x*q_)/m_; 
  double alphaUp[2]; 
  alphaUp[0] = -(0.7071067811865475*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 
  alphaUp[1] = -(0.7071067811865475*(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[1])*q_)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  emModR[0] += 0.5*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  emModR[1] += 0.5*(1.732050807568877*fl[3]+fl[1])*dfac_v; 
  emModR[2] += -0.5*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  emModR[3] += -0.5*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_v; 
  emModL[0] += -0.5*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  emModL[1] += -0.5*(1.732050807568877*fl[3]+fl[1])*dfac_v; 
  emModL[2] += -0.5*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  emModL[3] += -0.5*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_v; 
  incr[0] = 0.3535533905932737*(1.732050807568877*(alpha[1]*fl[3]+alpha[0]*fl[2])+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(1.732050807568877*(alpha[0]*fl[3]+alpha[1]*fl[2])+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.3535533905932737*(3.0*(alpha[1]*fl[3]+alpha[0]*fl[2])+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[3] = -0.3535533905932737*(3.0*(alpha[0]*fl[3]+alpha[1]*fl[2])+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  } else { 
  emModR[0] += -0.5*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  emModR[1] += -0.5*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_v; 
  emModR[2] += 0.5*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  emModR[3] += 0.5*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_v; 
  emModL[0] += 0.5*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  emModL[1] += 0.5*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_v; 
  emModL[2] += 0.5*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  emModL[3] += 0.5*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_v; 
  incr[0] = -0.3535533905932737*(1.732050807568877*(alpha[1]*fr[3]+alpha[0]*fr[2])-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.3535533905932737*(1.732050807568877*(alpha[0]*fr[3]+alpha[1]*fr[2])-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = 0.3535533905932737*(3.0*(alpha[1]*fr[3]+alpha[0]*fr[2])-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[3] = 0.3535533905932737*(3.0*(alpha[0]*fr[3]+alpha[1]*fr[2])-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double limQuad[4];
  if(0.7071067811865475*alphaUp[0]-0.7071067811865475*alphaUp[1] > 0) {
  fupwindQuad[0] = (-0.5*fl[3])+0.5*fl[2]-0.5*fl[1]+0.5*fl[0]; 
  } else {
  fupwindQuad[0] = 0.5*fr[3]-0.5*(fr[2]+fr[1])+0.5*fr[0]; 
  }
  if(0.7071067811865475*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[1] = 0.5*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[1] = 0.5*(fr[1]+fr[0])-0.5*(fr[3]+fr[2]); 
  }
  if(0.7071067811865475*alphaUp[0]-0.7071067811865475*alphaUp[1] > 0) {
  fupwindQuad[2] = 0.5*fl[3]-0.5*(fl[2]+fl[1])+0.5*fl[0]; 
  } else {
  fupwindQuad[2] = (-0.5*fr[3])+0.5*fr[2]-0.5*fr[1]+0.5*fr[0]; 
  }
  if(0.7071067811865475*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[3] = 0.5*(fl[1]+fl[0])-0.5*(fl[3]+fl[2]); 
  } else {
  fupwindQuad[3] = 0.5*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*(fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.5*(fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  emModR[0] += -0.5*(1.732050807568877*fupwind[2]-1.0*fupwind[0])*dfac_v; 
  emModR[1] += -0.5*(1.732050807568877*fupwind[3]-1.0*fupwind[1])*dfac_v; 
  emModR[2] += 0.5*(3.0*fupwind[2]-1.732050807568877*fupwind[0])*dfac_v; 
  emModR[3] += 0.5*(3.0*fupwind[3]-1.732050807568877*fupwind[1])*dfac_v; 
  emModL[0] += 0.5*(1.732050807568877*fupwind[2]-1.0*fupwind[0])*dfac_v; 
  emModL[1] += 0.5*(1.732050807568877*fupwind[3]-1.0*fupwind[1])*dfac_v; 
  emModL[2] += 0.5*(3.0*fupwind[2]-1.732050807568877*fupwind[0])*dfac_v; 
  emModL[3] += 0.5*(3.0*fupwind[3]-1.732050807568877*fupwind[1])*dfac_v; 
  incr[0] = -0.3535533905932737*(1.732050807568877*(alpha[1]*fupwind[3]+alpha[0]*fupwind[2])-1.0*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[1] = -0.3535533905932737*(1.732050807568877*(alpha[0]*fupwind[3]+alpha[1]*fupwind[2])-1.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 
  incr[2] = 0.3535533905932737*(3.0*(alpha[1]*fupwind[3]+alpha[0]*fupwind[2])-1.732050807568877*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[3] = 0.3535533905932737*(3.0*(alpha[0]*fupwind[3]+alpha[1]*fupwind[2])-1.732050807568877*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x1vSerStep2_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 

  double alpha[2]; 
  alpha[0] = -(1.0*dApardt[0]*q_)/m_; 
  alpha[1] = -(1.0*dApardt[1]*q_)/m_; 
  double alphaUp[2]; 
  alphaUp[0] = -(0.7071067811865475*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 
  alphaUp[1] = -(0.7071067811865475*(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[1])*q_)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*(1.732050807568877*(alpha[1]*fl[3]+alpha[0]*fl[2])+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(1.732050807568877*(alpha[0]*fl[3]+alpha[1]*fl[2])+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.3535533905932737*(3.0*(alpha[1]*fl[3]+alpha[0]*fl[2])+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[3] = -0.3535533905932737*(3.0*(alpha[0]*fl[3]+alpha[1]*fl[2])+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  } else { 
  incr[0] = -0.3535533905932737*(1.732050807568877*(alpha[1]*fr[3]+alpha[0]*fr[2])-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.3535533905932737*(1.732050807568877*(alpha[0]*fr[3]+alpha[1]*fr[2])-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = 0.3535533905932737*(3.0*(alpha[1]*fr[3]+alpha[0]*fr[2])-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[3] = 0.3535533905932737*(3.0*(alpha[0]*fr[3]+alpha[1]*fr[2])-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double limQuad[4];
  if(0.7071067811865475*alphaUp[0]-0.7071067811865475*alphaUp[1] > 0) {
  fupwindQuad[0] = (-0.5*fl[3])+0.5*fl[2]-0.5*fl[1]+0.5*fl[0]; 
  } else {
  fupwindQuad[0] = 0.5*fr[3]-0.5*(fr[2]+fr[1])+0.5*fr[0]; 
  }
  if(0.7071067811865475*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[1] = 0.5*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[1] = 0.5*(fr[1]+fr[0])-0.5*(fr[3]+fr[2]); 
  }
  if(0.7071067811865475*alphaUp[0]-0.7071067811865475*alphaUp[1] > 0) {
  fupwindQuad[2] = 0.5*fl[3]-0.5*(fl[2]+fl[1])+0.5*fl[0]; 
  } else {
  fupwindQuad[2] = (-0.5*fr[3])+0.5*fr[2]-0.5*fr[1]+0.5*fr[0]; 
  }
  if(0.7071067811865475*(alphaUp[1]+alphaUp[0]) > 0) {
  fupwindQuad[3] = 0.5*(fl[1]+fl[0])-0.5*(fl[3]+fl[2]); 
  } else {
  fupwindQuad[3] = 0.5*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*(fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.5*(fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  incr[0] = -0.3535533905932737*(1.732050807568877*(alpha[1]*fupwind[3]+alpha[0]*fupwind[2])-1.0*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[1] = -0.3535533905932737*(1.732050807568877*(alpha[0]*fupwind[3]+alpha[1]*fupwind[2])-1.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 
  incr[2] = 0.3535533905932737*(3.0*(alpha[1]*fupwind[3]+alpha[0]*fupwind[2])-1.732050807568877*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0]))*dfac_v; 
  incr[3] = 0.3535533905932737*(3.0*(alpha[0]*fupwind[3]+alpha[1]*fupwind[2])-1.732050807568877*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  return std::abs(alpha0); 
} 
