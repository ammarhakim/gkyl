#include <GyrokineticModDecl.h> 
double EmGyrokineticSurf1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*Gradpar[0]*wv; 

  double alpha[4]; 
  alpha[0] = 1.414213562373095*Gradpar[0]*wv; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_x; 
  incr[4] = -0.25*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_x; 
  incr[5] = -0.25*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_x; 
  incr[6] = 0.25*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_x; 
  incr[7] = -0.25*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_x; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_x; 
  incr[4] = 0.25*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_x; 
  incr[5] = 0.25*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_x; 
  incr[6] = -0.25*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_x; 
  incr[7] = 0.25*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double alphaQuad;
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[0] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.3535533905932737*fr[6]+0.3535533905932737*fl[6]-0.6123724356957944*(fr[5]+fl[5]+fr[4]+fl[4])+0.3535533905932737*(fr[3]+fr[2])-0.3535533905932737*(fl[3]+fl[2])+0.6123724356957944*(fr[1]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)-0.6123724356957944*fr[7]+0.6123724356957944*fl[7]+0.3535533905932737*(fr[6]+fl[6])+0.6123724356957944*(fr[5]+fr[4])-0.6123724356957944*(fl[5]+fl[4])-0.3535533905932737*(fr[3]+fl[3]+fr[2]+fl[2])-0.6123724356957944*fr[1]+0.6123724356957944*fl[1]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.3535533905932737*fr[6]-0.3535533905932737*fl[6]-0.6123724356957944*(fr[5]+fl[5])+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*(fl[3]+fr[2])+0.3535533905932737*fl[2]+0.6123724356957944*(fr[1]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)+0.6123724356957944*fr[7]-0.6123724356957944*fl[7]-0.3535533905932737*(fr[6]+fl[6])+0.6123724356957944*fr[5]-0.6123724356957944*(fl[5]+fr[4])+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])+0.3535533905932737*(fr[2]+fl[2])-0.6123724356957944*fr[1]+0.6123724356957944*fl[1]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[2] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.3535533905932737*fr[6]-0.3535533905932737*fl[6]+0.6123724356957944*(fr[5]+fl[5])-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*(fl[3]+fr[2])-0.3535533905932737*fl[2]+0.6123724356957944*(fr[1]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)+0.6123724356957944*fr[7]-0.6123724356957944*fl[7]-0.3535533905932737*(fr[6]+fl[6])-0.6123724356957944*fr[5]+0.6123724356957944*(fl[5]+fr[4])-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.3535533905932737*(fr[2]+fl[2])-0.6123724356957944*fr[1]+0.6123724356957944*fl[1]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[3] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.3535533905932737*fr[6]+0.3535533905932737*fl[6]+0.6123724356957944*(fr[5]+fl[5]+fr[4]+fl[4])-0.3535533905932737*(fr[3]+fr[2])+0.3535533905932737*(fl[3]+fl[2])+0.6123724356957944*(fr[1]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)-0.6123724356957944*fr[7]+0.6123724356957944*fl[7]+0.3535533905932737*(fr[6]+fl[6])-0.6123724356957944*(fr[5]+fr[4])+0.6123724356957944*(fl[5]+fl[4])+0.3535533905932737*(fr[3]+fl[3]+fr[2]+fl[2])-0.6123724356957944*fr[1]+0.6123724356957944*fl[1]+0.3535533905932737*(fr[0]+fl[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*fupwindQuad[3]-0.5*fupwindQuad[2]+0.5*fupwindQuad[1]-0.5*fupwindQuad[0]; 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2])-0.5*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.5*fupwindQuad[3]-0.5*(fupwindQuad[2]+fupwindQuad[1])+0.5*fupwindQuad[0]; 
  incr[0] = 0.3535533905932737*alpha[0]*fupwind[0]*dfac_x; 
  incr[1] = -0.6123724356957944*alpha[0]*fupwind[0]*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*fupwind[1]*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*fupwind[2]*dfac_x; 
  incr[4] = -0.6123724356957944*alpha[0]*fupwind[1]*dfac_x; 
  incr[5] = -0.6123724356957944*alpha[0]*fupwind[2]*dfac_x; 
  incr[6] = 0.3535533905932737*alpha[0]*fupwind[3]*dfac_x; 
  incr[7] = -0.6123724356957944*alpha[0]*fupwind[3]*dfac_x; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  double incrEmMod[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 

  double alpha[4]; 
  alpha[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  double alphaUp[4]; 
  alphaUp[0] = -(1.0*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 
  alphaUp[1] = -(1.414213562373095*dApardtPrev[1]*q_)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  incr[1] = 0.25*alpha[0]*(1.732050807568877*fl[4]+fl[1])*dfac_v; 
  incr[2] = -0.25*alpha[0]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  incr[3] = 0.25*alpha[0]*(1.732050807568877*fl[6]+fl[3])*dfac_v; 
  incr[4] = -0.25*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[1])*dfac_v; 
  incr[5] = 0.25*alpha[0]*(1.732050807568877*fl[7]+fl[5])*dfac_v; 
  incr[6] = -0.25*alpha[0]*(3.0*fl[6]+1.732050807568877*fl[3])*dfac_v; 
  incr[7] = -0.25*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[5])*dfac_v; 
  incrEmMod[0] = 0.5*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  incrEmMod[1] = 0.5*(1.732050807568877*fl[4]+fl[1])*dfac_v; 
  incrEmMod[2] = -0.5*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  incrEmMod[3] = 0.5*(1.732050807568877*fl[6]+fl[3])*dfac_v; 
  incrEmMod[4] = -0.5*(3.0*fl[4]+1.732050807568877*fl[1])*dfac_v; 
  incrEmMod[5] = 0.5*(1.732050807568877*fl[7]+fl[5])*dfac_v; 
  incrEmMod[6] = -0.5*(3.0*fl[6]+1.732050807568877*fl[3])*dfac_v; 
  incrEmMod[7] = -0.5*(3.0*fl[7]+1.732050807568877*fl[5])*dfac_v; 
  } else { 
  incrEmMod[0] = -0.5*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  incrEmMod[1] = -0.5*(1.732050807568877*fr[4]-1.0*fr[1])*dfac_v; 
  incrEmMod[2] = 0.5*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  incrEmMod[3] = -0.5*(1.732050807568877*fr[6]-1.0*fr[3])*dfac_v; 
  incrEmMod[4] = 0.5*(3.0*fr[4]-1.732050807568877*fr[1])*dfac_v; 
  incrEmMod[5] = -0.5*(1.732050807568877*fr[7]-1.0*fr[5])*dfac_v; 
  incrEmMod[6] = 0.5*(3.0*fr[6]-1.732050807568877*fr[3])*dfac_v; 
  incrEmMod[7] = 0.5*(3.0*fr[7]-1.732050807568877*fr[5])*dfac_v; 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  incr[1] = -0.25*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[1])*dfac_v; 
  incr[2] = 0.25*alpha[0]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  incr[3] = -0.25*alpha[0]*(1.732050807568877*fr[6]-1.0*fr[3])*dfac_v; 
  incr[4] = 0.25*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[1])*dfac_v; 
  incr[5] = -0.25*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[5])*dfac_v; 
  incr[6] = 0.25*alpha[0]*(3.0*fr[6]-1.732050807568877*fr[3])*dfac_v; 
  incr[7] = 0.25*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[5])*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double alphaQuad;
  alphaQuad = 0.5*alphaUp[0]-0.5*alphaUp[1]; 
  fupwindQuad[0] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)-0.6123724356957944*fr[7]+0.6123724356957944*(fl[7]+fr[6])-0.6123724356957944*fl[6]+0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*(alphaUp[1]+alphaUp[0]); 
  fupwindQuad[1] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)+0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6])-0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaQuad = 0.5*alphaUp[0]-0.5*alphaUp[1]; 
  fupwindQuad[2] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)+0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*fl[6]-0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*(alphaUp[1]+alphaUp[0]); 
  fupwindQuad[3] = 0.5*((0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)-0.6123724356957944*(fr[7]+fr[6])+0.6123724356957944*(fl[7]+fl[6])+0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*fupwindQuad[3]-0.5*fupwindQuad[2]+0.5*fupwindQuad[1]-0.5*fupwindQuad[0]; 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2])-0.5*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.5*fupwindQuad[3]-0.5*(fupwindQuad[2]+fupwindQuad[1])+0.5*fupwindQuad[0]; 
  incrEmMod[0] = 0.7071067811865475*fupwind[0]*dfac_v; 
  incrEmMod[1] = 0.7071067811865475*fupwind[1]*dfac_v; 
  incrEmMod[2] = -1.224744871391589*fupwind[0]*dfac_v; 
  incrEmMod[3] = 0.7071067811865475*fupwind[2]*dfac_v; 
  incrEmMod[4] = -1.224744871391589*fupwind[1]*dfac_v; 
  incrEmMod[5] = 0.7071067811865475*fupwind[3]*dfac_v; 
  incrEmMod[6] = -1.224744871391589*fupwind[2]*dfac_v; 
  incrEmMod[7] = -1.224744871391589*fupwind[3]*dfac_v; 
  incr[0] = 0.3535533905932737*alpha[0]*fupwind[0]*dfac_v; 
  incr[1] = 0.3535533905932737*alpha[0]*fupwind[1]*dfac_v; 
  incr[2] = -0.6123724356957944*alpha[0]*fupwind[0]*dfac_v; 
  incr[3] = 0.3535533905932737*alpha[0]*fupwind[2]*dfac_v; 
  incr[4] = -0.6123724356957944*alpha[0]*fupwind[1]*dfac_v; 
  incr[5] = 0.3535533905932737*alpha[0]*fupwind[3]*dfac_v; 
  incr[6] = -0.6123724356957944*alpha[0]*fupwind[2]*dfac_v; 
  incr[7] = -0.6123724356957944*alpha[0]*fupwind[3]*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  emModR[0] += incrEmMod[0]; 
  emModR[1] += incrEmMod[1]; 
  emModR[2] += incrEmMod[2]; 
  emModR[3] += incrEmMod[3]; 
  emModR[4] += incrEmMod[4]; 
  emModR[5] += incrEmMod[5]; 
  emModR[6] += incrEmMod[6]; 
  emModR[7] += incrEmMod[7]; 
  emModL[0] += -1.0*incrEmMod[0]; 
  emModL[1] += -1.0*incrEmMod[1]; 
  emModL[2] += incrEmMod[2]; 
  emModL[3] += -1.0*incrEmMod[3]; 
  emModL[4] += incrEmMod[4]; 
  emModL[5] += -1.0*incrEmMod[5]; 
  emModL[6] += incrEmMod[6]; 
  emModL[7] += incrEmMod[7]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSerStep2_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 

  double alpha[4]; 
  alpha[0] = -(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 
  double alphaUp[4]; 
  alphaUp[0] = -(1.0*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_)/m_; 
  alphaUp[1] = -(1.414213562373095*dApardtPrev[1]*q_)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*(alpha[1]*fl[4]+alpha[0]*fl[2])+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*(alpha[0]*fl[4]+alpha[1]*fl[2])+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*(alpha[1]*fl[4]+alpha[0]*fl[2])+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[3] = 0.25*(1.732050807568877*(alpha[1]*fl[7]+alpha[0]*fl[6])+alpha[1]*fl[5]+alpha[0]*fl[3])*dfac_v; 
  incr[4] = -0.25*(3.0*(alpha[0]*fl[4]+alpha[1]*fl[2])+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[5] = 0.25*(1.732050807568877*(alpha[0]*fl[7]+alpha[1]*fl[6])+alpha[0]*fl[5]+alpha[1]*fl[3])*dfac_v; 
  incr[6] = -0.25*(3.0*(alpha[1]*fl[7]+alpha[0]*fl[6])+1.732050807568877*(alpha[1]*fl[5]+alpha[0]*fl[3]))*dfac_v; 
  incr[7] = -0.25*(3.0*(alpha[0]*fl[7]+alpha[1]*fl[6])+1.732050807568877*(alpha[0]*fl[5]+alpha[1]*fl[3]))*dfac_v; 
  } else { 
  incr[0] = -0.25*(1.732050807568877*(alpha[1]*fr[4]+alpha[0]*fr[2])-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*(alpha[0]*fr[4]+alpha[1]*fr[2])-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = 0.25*(3.0*(alpha[1]*fr[4]+alpha[0]*fr[2])-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[3] = -0.25*(1.732050807568877*(alpha[1]*fr[7]+alpha[0]*fr[6])-1.0*(alpha[1]*fr[5]+alpha[0]*fr[3]))*dfac_v; 
  incr[4] = 0.25*(3.0*(alpha[0]*fr[4]+alpha[1]*fr[2])-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[5] = -0.25*(1.732050807568877*(alpha[0]*fr[7]+alpha[1]*fr[6])-1.0*(alpha[0]*fr[5]+alpha[1]*fr[3]))*dfac_v; 
  incr[6] = 0.25*(3.0*(alpha[1]*fr[7]+alpha[0]*fr[6])-1.732050807568877*(alpha[1]*fr[5]+alpha[0]*fr[3]))*dfac_v; 
  incr[7] = 0.25*(3.0*(alpha[0]*fr[7]+alpha[1]*fr[6])-1.732050807568877*(alpha[0]*fr[5]+alpha[1]*fr[3]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double alphaQuad;
  alphaQuad = 0.5*alphaUp[0]-0.5*alphaUp[1]; 
  fupwindQuad[0] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)-0.6123724356957944*fr[7]+0.6123724356957944*(fl[7]+fr[6])-0.6123724356957944*fl[6]+0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*(alphaUp[1]+alphaUp[0]); 
  fupwindQuad[1] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)+0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6])-0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaQuad = 0.5*alphaUp[0]-0.5*alphaUp[1]; 
  fupwindQuad[2] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)+0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*fl[6]-0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*(alphaUp[1]+alphaUp[0]); 
  fupwindQuad[3] = 0.5*((0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)-0.6123724356957944*(fr[7]+fr[6])+0.6123724356957944*(fl[7]+fl[6])+0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*fupwindQuad[3]-0.5*fupwindQuad[2]+0.5*fupwindQuad[1]-0.5*fupwindQuad[0]; 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2])-0.5*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.5*fupwindQuad[3]-0.5*(fupwindQuad[2]+fupwindQuad[1])+0.5*fupwindQuad[0]; 
  incr[0] = 0.3535533905932737*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac_v; 
  incr[2] = -0.6123724356957944*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 
  incr[3] = 0.3535533905932737*(alpha[1]*fupwind[3]+alpha[0]*fupwind[2])*dfac_v; 
  incr[4] = -0.6123724356957944*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac_v; 
  incr[5] = 0.3535533905932737*(alpha[0]*fupwind[3]+alpha[1]*fupwind[2])*dfac_v; 
  incr[6] = -0.6123724356957944*(alpha[1]*fupwind[3]+alpha[0]*fupwind[2])*dfac_v; 
  incr[7] = -0.6123724356957944*(alpha[0]*fupwind[3]+alpha[1]*fupwind[2])*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.3535533905932737*(1.732050807568877*Gradpar[1]-1.0*Gradpar[0])*wv; 

  double alpha[4]; 
  alpha[0] = (1.414213562373095*Gradpar[0]-2.449489742783178*Gradpar[1])*wv; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_x; 
  incr[4] = -0.25*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_x; 
  incr[5] = -0.25*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_x; 
  incr[6] = 0.25*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_x; 
  incr[7] = -0.25*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_x; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_x; 
  incr[4] = 0.25*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_x; 
  incr[5] = 0.25*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_x; 
  incr[6] = -0.25*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_x; 
  incr[7] = 0.25*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double alphaQuad;
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[0] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.3535533905932737*fr[6]+0.3535533905932737*fl[6]-0.6123724356957944*(fr[5]+fl[5]+fr[4]+fl[4])+0.3535533905932737*(fr[3]+fr[2])-0.3535533905932737*(fl[3]+fl[2])+0.6123724356957944*(fr[1]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)-0.6123724356957944*fr[7]+0.6123724356957944*fl[7]+0.3535533905932737*(fr[6]+fl[6])+0.6123724356957944*(fr[5]+fr[4])-0.6123724356957944*(fl[5]+fl[4])-0.3535533905932737*(fr[3]+fl[3]+fr[2]+fl[2])-0.6123724356957944*fr[1]+0.6123724356957944*fl[1]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.3535533905932737*fr[6]-0.3535533905932737*fl[6]-0.6123724356957944*(fr[5]+fl[5])+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*(fl[3]+fr[2])+0.3535533905932737*fl[2]+0.6123724356957944*(fr[1]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)+0.6123724356957944*fr[7]-0.6123724356957944*fl[7]-0.3535533905932737*(fr[6]+fl[6])+0.6123724356957944*fr[5]-0.6123724356957944*(fl[5]+fr[4])+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])+0.3535533905932737*(fr[2]+fl[2])-0.6123724356957944*fr[1]+0.6123724356957944*fl[1]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[2] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.3535533905932737*fr[6]-0.3535533905932737*fl[6]+0.6123724356957944*(fr[5]+fl[5])-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*(fl[3]+fr[2])-0.3535533905932737*fl[2]+0.6123724356957944*(fr[1]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)+0.6123724356957944*fr[7]-0.6123724356957944*fl[7]-0.3535533905932737*(fr[6]+fl[6])-0.6123724356957944*fr[5]+0.6123724356957944*(fl[5]+fr[4])-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.3535533905932737*(fr[2]+fl[2])-0.6123724356957944*fr[1]+0.6123724356957944*fl[1]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[3] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.3535533905932737*fr[6]+0.3535533905932737*fl[6]+0.6123724356957944*(fr[5]+fl[5]+fr[4]+fl[4])-0.3535533905932737*(fr[3]+fr[2])+0.3535533905932737*(fl[3]+fl[2])+0.6123724356957944*(fr[1]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)-0.6123724356957944*fr[7]+0.6123724356957944*fl[7]+0.3535533905932737*(fr[6]+fl[6])-0.6123724356957944*(fr[5]+fr[4])+0.6123724356957944*(fl[5]+fl[4])+0.3535533905932737*(fr[3]+fl[3]+fr[2]+fl[2])-0.6123724356957944*fr[1]+0.6123724356957944*fl[1]+0.3535533905932737*(fr[0]+fl[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*fupwindQuad[3]-0.5*fupwindQuad[2]+0.5*fupwindQuad[1]-0.5*fupwindQuad[0]; 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2])-0.5*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.5*fupwindQuad[3]-0.5*(fupwindQuad[2]+fupwindQuad[1])+0.5*fupwindQuad[0]; 
  incr[0] = 0.3535533905932737*alpha[0]*fupwind[0]*dfac_x; 
  incr[1] = -0.6123724356957944*alpha[0]*fupwind[0]*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*fupwind[1]*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*fupwind[2]*dfac_x; 
  incr[4] = -0.6123724356957944*alpha[0]*fupwind[1]*dfac_x; 
  incr[5] = -0.6123724356957944*alpha[0]*fupwind[2]*dfac_x; 
  incr[6] = 0.3535533905932737*alpha[0]*fupwind[3]*dfac_x; 
  incr[7] = -0.6123724356957944*alpha[0]*fupwind[3]*dfac_x; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  double incrEmMod[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(1.732050807568877*Gradpar[0]*Bmag[1]*dfac_x*wm+(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_))/m_; 

  double alpha[4]; 
  alpha[0] = -(1.732050807568877*Gradpar[0]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))/m_; 
  alpha[1] = -(1.732050807568877*Gradpar[1]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))/m_; 
  alpha[2] = -(1.0*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[3] = -(1.0*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
  double alphaUp[4]; 
  alphaUp[0] = -(1.0*(1.732050807568877*Gradpar[0]*Bmag[1]*dfac_x*wm+(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_))/m_; 
  alphaUp[1] = -(1.0*(1.732050807568877*Bmag[1]*Gradpar[1]*dfac_x*wm+(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[1])*q_))/m_; 
  alphaUp[2] = -(1.0*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alphaUp[3] = -(1.0*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*(alpha[3]*fl[7]+alpha[2]*fl[6])+alpha[3]*fl[5]+1.732050807568877*alpha[1]*fl[4]+alpha[2]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*(alpha[2]*fl[7]+alpha[3]*fl[6])+alpha[2]*fl[5]+1.732050807568877*alpha[0]*fl[4]+alpha[3]*fl[3]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*(alpha[3]*fl[7]+alpha[2]*fl[6])+1.732050807568877*alpha[3]*fl[5]+3.0*alpha[1]*fl[4]+1.732050807568877*alpha[2]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[3] = 0.25*(1.732050807568877*(alpha[1]*fl[7]+alpha[0]*fl[6])+alpha[1]*fl[5]+1.732050807568877*alpha[3]*fl[4]+alpha[0]*fl[3]+fl[1]*alpha[3]+alpha[2]*(1.732050807568877*fl[2]+fl[0]))*dfac_v; 
  incr[4] = -0.25*(3.0*(alpha[2]*fl[7]+alpha[3]*fl[6])+1.732050807568877*alpha[2]*fl[5]+3.0*alpha[0]*fl[4]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[1]*fl[2]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[5] = 0.25*(1.732050807568877*(alpha[0]*fl[7]+alpha[1]*fl[6])+alpha[0]*fl[5]+1.732050807568877*alpha[2]*fl[4]+alpha[1]*fl[3]+(1.732050807568877*fl[2]+fl[0])*alpha[3]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.25*(3.0*(alpha[1]*fl[7]+alpha[0]*fl[6])+1.732050807568877*alpha[1]*fl[5]+3.0*alpha[3]*fl[4]+1.732050807568877*(alpha[0]*fl[3]+fl[1]*alpha[3])+alpha[2]*(3.0*fl[2]+1.732050807568877*fl[0]))*dfac_v; 
  incr[7] = -0.25*(3.0*(alpha[0]*fl[7]+alpha[1]*fl[6])+1.732050807568877*alpha[0]*fl[5]+3.0*alpha[2]*fl[4]+1.732050807568877*alpha[1]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*(fl[0]*alpha[3]+fl[1]*alpha[2]))*dfac_v; 
  incrEmMod[0] = 0.5*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  incrEmMod[1] = 0.5*(1.732050807568877*fl[4]+fl[1])*dfac_v; 
  incrEmMod[2] = -0.5*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  incrEmMod[3] = 0.5*(1.732050807568877*fl[6]+fl[3])*dfac_v; 
  incrEmMod[4] = -0.5*(3.0*fl[4]+1.732050807568877*fl[1])*dfac_v; 
  incrEmMod[5] = 0.5*(1.732050807568877*fl[7]+fl[5])*dfac_v; 
  incrEmMod[6] = -0.5*(3.0*fl[6]+1.732050807568877*fl[3])*dfac_v; 
  incrEmMod[7] = -0.5*(3.0*fl[7]+1.732050807568877*fl[5])*dfac_v; 
  } else { 
  incrEmMod[0] = -0.5*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  incrEmMod[1] = -0.5*(1.732050807568877*fr[4]-1.0*fr[1])*dfac_v; 
  incrEmMod[2] = 0.5*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  incrEmMod[3] = -0.5*(1.732050807568877*fr[6]-1.0*fr[3])*dfac_v; 
  incrEmMod[4] = 0.5*(3.0*fr[4]-1.732050807568877*fr[1])*dfac_v; 
  incrEmMod[5] = -0.5*(1.732050807568877*fr[7]-1.0*fr[5])*dfac_v; 
  incrEmMod[6] = 0.5*(3.0*fr[6]-1.732050807568877*fr[3])*dfac_v; 
  incrEmMod[7] = 0.5*(3.0*fr[7]-1.732050807568877*fr[5])*dfac_v; 
  incr[0] = -0.25*(1.732050807568877*(alpha[3]*fr[7]+alpha[2]*fr[6])-1.0*alpha[3]*fr[5]+1.732050807568877*alpha[1]*fr[4]-1.0*alpha[2]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*(alpha[2]*fr[7]+alpha[3]*fr[6])-1.0*alpha[2]*fr[5]+1.732050807568877*alpha[0]*fr[4]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[1]*fr[2]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = 0.25*(3.0*(alpha[3]*fr[7]+alpha[2]*fr[6])-1.732050807568877*alpha[3]*fr[5]+3.0*alpha[1]*fr[4]-1.732050807568877*alpha[2]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[3] = -0.25*(1.732050807568877*(alpha[1]*fr[7]+alpha[0]*fr[6])-1.0*alpha[1]*fr[5]+1.732050807568877*alpha[3]*fr[4]-1.0*(alpha[0]*fr[3]+fr[1]*alpha[3])+alpha[2]*(1.732050807568877*fr[2]-1.0*fr[0]))*dfac_v; 
  incr[4] = 0.25*(3.0*(alpha[2]*fr[7]+alpha[3]*fr[6])-1.732050807568877*alpha[2]*fr[5]+3.0*alpha[0]*fr[4]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[1]*fr[2]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[5] = -0.25*(1.732050807568877*(alpha[0]*fr[7]+alpha[1]*fr[6])-1.0*alpha[0]*fr[5]+1.732050807568877*alpha[2]*fr[4]-1.0*alpha[1]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*(fr[0]*alpha[3]+fr[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.25*(3.0*(alpha[1]*fr[7]+alpha[0]*fr[6])-1.732050807568877*alpha[1]*fr[5]+3.0*alpha[3]*fr[4]-1.732050807568877*(alpha[0]*fr[3]+fr[1]*alpha[3])+alpha[2]*(3.0*fr[2]-1.732050807568877*fr[0]))*dfac_v; 
  incr[7] = 0.25*(3.0*(alpha[0]*fr[7]+alpha[1]*fr[6])-1.732050807568877*alpha[0]*fr[5]+3.0*alpha[2]*fr[4]-1.732050807568877*alpha[1]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*(fr[0]*alpha[3]+fr[1]*alpha[2]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double alphaQuad;
  alphaQuad = 0.5*alphaUp[3]-0.5*(alphaUp[2]+alphaUp[1])+0.5*alphaUp[0]; 
  fupwindQuad[0] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)-0.6123724356957944*fr[7]+0.6123724356957944*(fl[7]+fr[6])-0.6123724356957944*fl[6]+0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*(alphaUp[1]+alphaUp[0])-0.5*(alphaUp[3]+alphaUp[2]); 
  fupwindQuad[1] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)+0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6])-0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaQuad = (-0.5*alphaUp[3])+0.5*alphaUp[2]-0.5*alphaUp[1]+0.5*alphaUp[0]; 
  fupwindQuad[2] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)+0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*fl[6]-0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*(alphaUp[3]+alphaUp[2]+alphaUp[1]+alphaUp[0]); 
  fupwindQuad[3] = 0.5*((0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)-0.6123724356957944*(fr[7]+fr[6])+0.6123724356957944*(fl[7]+fl[6])+0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*fupwindQuad[3]-0.5*fupwindQuad[2]+0.5*fupwindQuad[1]-0.5*fupwindQuad[0]; 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2])-0.5*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.5*fupwindQuad[3]-0.5*(fupwindQuad[2]+fupwindQuad[1])+0.5*fupwindQuad[0]; 
  incrEmMod[0] = 0.7071067811865475*fupwind[0]*dfac_v; 
  incrEmMod[1] = 0.7071067811865475*fupwind[1]*dfac_v; 
  incrEmMod[2] = -1.224744871391589*fupwind[0]*dfac_v; 
  incrEmMod[3] = 0.7071067811865475*fupwind[2]*dfac_v; 
  incrEmMod[4] = -1.224744871391589*fupwind[1]*dfac_v; 
  incrEmMod[5] = 0.7071067811865475*fupwind[3]*dfac_v; 
  incrEmMod[6] = -1.224744871391589*fupwind[2]*dfac_v; 
  incrEmMod[7] = -1.224744871391589*fupwind[3]*dfac_v; 
  incr[0] = 0.3535533905932737*(alpha[3]*fupwind[3]+alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(alpha[2]*fupwind[3]+fupwind[2]*alpha[3]+alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac_v; 
  incr[2] = -0.6123724356957944*(alpha[3]*fupwind[3]+alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 
  incr[3] = 0.3535533905932737*(alpha[1]*fupwind[3]+fupwind[1]*alpha[3]+alpha[0]*fupwind[2]+fupwind[0]*alpha[2])*dfac_v; 
  incr[4] = -0.6123724356957944*(alpha[2]*fupwind[3]+fupwind[2]*alpha[3]+alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac_v; 
  incr[5] = 0.3535533905932737*(alpha[0]*fupwind[3]+fupwind[0]*alpha[3]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2])*dfac_v; 
  incr[6] = -0.6123724356957944*(alpha[1]*fupwind[3]+fupwind[1]*alpha[3]+alpha[0]*fupwind[2]+fupwind[0]*alpha[2])*dfac_v; 
  incr[7] = -0.6123724356957944*(alpha[0]*fupwind[3]+fupwind[0]*alpha[3]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2])*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  emModR[0] += incrEmMod[0]; 
  emModR[1] += incrEmMod[1]; 
  emModR[2] += incrEmMod[2]; 
  emModR[3] += incrEmMod[3]; 
  emModR[4] += incrEmMod[4]; 
  emModR[5] += incrEmMod[5]; 
  emModR[6] += incrEmMod[6]; 
  emModR[7] += incrEmMod[7]; 
  emModL[0] += -1.0*incrEmMod[0]; 
  emModL[1] += -1.0*incrEmMod[1]; 
  emModL[2] += incrEmMod[2]; 
  emModL[3] += -1.0*incrEmMod[3]; 
  emModL[4] += incrEmMod[4]; 
  emModL[5] += -1.0*incrEmMod[5]; 
  emModL[6] += incrEmMod[6]; 
  emModL[7] += incrEmMod[7]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSerStep2_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fl, const double *fr, double *outl, double *outr, double *emModL, double *emModR) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(1.732050807568877*Gradpar[0]*Bmag[1]*dfac_x*wm+(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_))/m_; 

  double alpha[4]; 
  alpha[0] = -(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 
  double alphaUp[4]; 
  alphaUp[0] = -(1.0*(1.732050807568877*Gradpar[0]*Bmag[1]*dfac_x*wm+(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[0])*q_))/m_; 
  alphaUp[1] = -(1.0*(1.732050807568877*Bmag[1]*Gradpar[1]*dfac_x*wm+(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x+1.414213562373095*dApardtPrev[1])*q_))/m_; 
  alphaUp[2] = -(1.0*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alphaUp[3] = -(1.0*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*(alpha[1]*fl[4]+alpha[0]*fl[2])+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*(alpha[0]*fl[4]+alpha[1]*fl[2])+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*(alpha[1]*fl[4]+alpha[0]*fl[2])+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[3] = 0.25*(1.732050807568877*(alpha[1]*fl[7]+alpha[0]*fl[6])+alpha[1]*fl[5]+alpha[0]*fl[3])*dfac_v; 
  incr[4] = -0.25*(3.0*(alpha[0]*fl[4]+alpha[1]*fl[2])+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[5] = 0.25*(1.732050807568877*(alpha[0]*fl[7]+alpha[1]*fl[6])+alpha[0]*fl[5]+alpha[1]*fl[3])*dfac_v; 
  incr[6] = -0.25*(3.0*(alpha[1]*fl[7]+alpha[0]*fl[6])+1.732050807568877*(alpha[1]*fl[5]+alpha[0]*fl[3]))*dfac_v; 
  incr[7] = -0.25*(3.0*(alpha[0]*fl[7]+alpha[1]*fl[6])+1.732050807568877*(alpha[0]*fl[5]+alpha[1]*fl[3]))*dfac_v; 
  } else { 
  incr[0] = -0.25*(1.732050807568877*(alpha[1]*fr[4]+alpha[0]*fr[2])-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*(alpha[0]*fr[4]+alpha[1]*fr[2])-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = 0.25*(3.0*(alpha[1]*fr[4]+alpha[0]*fr[2])-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[3] = -0.25*(1.732050807568877*(alpha[1]*fr[7]+alpha[0]*fr[6])-1.0*(alpha[1]*fr[5]+alpha[0]*fr[3]))*dfac_v; 
  incr[4] = 0.25*(3.0*(alpha[0]*fr[4]+alpha[1]*fr[2])-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[5] = -0.25*(1.732050807568877*(alpha[0]*fr[7]+alpha[1]*fr[6])-1.0*(alpha[0]*fr[5]+alpha[1]*fr[3]))*dfac_v; 
  incr[6] = 0.25*(3.0*(alpha[1]*fr[7]+alpha[0]*fr[6])-1.732050807568877*(alpha[1]*fr[5]+alpha[0]*fr[3]))*dfac_v; 
  incr[7] = 0.25*(3.0*(alpha[0]*fr[7]+alpha[1]*fr[6])-1.732050807568877*(alpha[0]*fr[5]+alpha[1]*fr[3]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double alphaQuad;
  alphaQuad = 0.5*alphaUp[3]-0.5*(alphaUp[2]+alphaUp[1])+0.5*alphaUp[0]; 
  fupwindQuad[0] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)-0.6123724356957944*fr[7]+0.6123724356957944*(fl[7]+fr[6])-0.6123724356957944*fl[6]+0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*(alphaUp[1]+alphaUp[0])-0.5*(alphaUp[3]+alphaUp[2]); 
  fupwindQuad[1] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)+0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6])-0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaQuad = (-0.5*alphaUp[3])+0.5*alphaUp[2]-0.5*alphaUp[1]+0.5*alphaUp[0]; 
  fupwindQuad[2] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)+0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*fl[6]-0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*(alphaUp[3]+alphaUp[2]+alphaUp[1]+alphaUp[0]); 
  fupwindQuad[3] = 0.5*((0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)-0.6123724356957944*(fr[7]+fr[6])+0.6123724356957944*(fl[7]+fl[6])+0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*fupwindQuad[3]-0.5*fupwindQuad[2]+0.5*fupwindQuad[1]-0.5*fupwindQuad[0]; 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2])-0.5*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.5*fupwindQuad[3]-0.5*(fupwindQuad[2]+fupwindQuad[1])+0.5*fupwindQuad[0]; 
  incr[0] = 0.3535533905932737*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac_v; 
  incr[2] = -0.6123724356957944*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 
  incr[3] = 0.3535533905932737*(alpha[1]*fupwind[3]+alpha[0]*fupwind[2])*dfac_v; 
  incr[4] = -0.6123724356957944*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac_v; 
  incr[5] = 0.3535533905932737*(alpha[0]*fupwind[3]+alpha[1]*fupwind[2])*dfac_v; 
  incr[6] = -0.6123724356957944*(alpha[1]*fupwind[3]+alpha[0]*fupwind[2])*dfac_v; 
  incr[7] = -0.6123724356957944*(alpha[0]*fupwind[3]+alpha[1]*fupwind[2])*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
