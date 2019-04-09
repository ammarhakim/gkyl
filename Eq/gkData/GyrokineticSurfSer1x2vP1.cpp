#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *f0l, const double *f0r, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_z = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wz = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  double hamil[8]; 
  hamil[0] = (0.2357022603955158*(3.0*dfac_v2*(2.0*m_*wv2+2.828427124746191*(Bmag[0]*wm+Phi[0]*q_))+2.0*m_))/dfac_v2; 
  hamil[1] = 2.0*Phi[1]*q_; 
  hamil[2] = (1.632993161855453*m_*wv)/dfac_v; 
  hamil[3] = (1.154700538379252*Bmag[0])/dfac_m; 
  double BstarX_by_Bmag[8]; 
  double BstarY_by_Bmag[8]; 
  double BstarZ_by_Bmag[8]; 
  BstarZ_by_Bmag[0] = 2.0*Gradpar[0]; 

  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.1082531754730548*BstarZ_by_Bmag[0]*hamil[2]*dfac_v)/m_; 

  double alpha[4]; 
  alpha[0] = (0.4330127018922193*BstarZ_by_Bmag[0]*hamil[2]*dfac_v)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_z; 
  incr[1] = -0.25*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_z; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_z; 
  incr[3] = 0.25*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_z; 
  incr[4] = -0.25*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_z; 
  incr[5] = -0.25*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_z; 
  incr[6] = 0.25*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_z; 
  incr[7] = -0.25*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_z; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_z; 
  incr[1] = 0.25*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_z; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_z; 
  incr[3] = -0.25*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_z; 
  incr[4] = 0.25*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_z; 
  incr[5] = 0.25*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_z; 
  incr[6] = -0.25*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_z; 
  incr[7] = 0.25*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_z; 
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
  incr[0] = 0.3535533905932737*alpha[0]*fupwind[0]*dfac_z; 
  incr[1] = -0.6123724356957944*alpha[0]*fupwind[0]*dfac_z; 
  incr[2] = 0.3535533905932737*alpha[0]*fupwind[1]*dfac_z; 
  incr[3] = 0.3535533905932737*alpha[0]*fupwind[2]*dfac_z; 
  incr[4] = -0.6123724356957944*alpha[0]*fupwind[1]*dfac_z; 
  incr[5] = -0.6123724356957944*alpha[0]*fupwind[2]*dfac_z; 
  incr[6] = 0.3535533905932737*alpha[0]*fupwind[3]*dfac_z; 
  incr[7] = -0.6123724356957944*alpha[0]*fupwind[3]*dfac_z; 

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

  // ensure cancellation of zeroth order terms for f0=F_M 
#if upwindType == SURFAVG 
  alpha[0] = -(0.8660254037844386*Gradpar[0]*hamil[2]*dfac_v)/m_; 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*f0l[1]+f0l[0])*dfac_z; 
  incr[1] = -0.25*alpha[0]*(3.0*f0l[1]+1.732050807568877*f0l[0])*dfac_z; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*f0l[4]+f0l[2])*dfac_z; 
  incr[3] = 0.25*alpha[0]*(1.732050807568877*f0l[5]+f0l[3])*dfac_z; 
  incr[4] = -0.25*alpha[0]*(3.0*f0l[4]+1.732050807568877*f0l[2])*dfac_z; 
  incr[5] = -0.25*alpha[0]*(3.0*f0l[5]+1.732050807568877*f0l[3])*dfac_z; 
  incr[6] = 0.25*alpha[0]*(1.732050807568877*f0l[7]+f0l[6])*dfac_z; 
  incr[7] = -0.25*alpha[0]*(3.0*f0l[7]+1.732050807568877*f0l[6])*dfac_z; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*f0r[1]-1.0*f0r[0])*dfac_z; 
  incr[1] = 0.25*alpha[0]*(3.0*f0r[1]-1.732050807568877*f0r[0])*dfac_z; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*f0r[4]-1.0*f0r[2])*dfac_z; 
  incr[3] = -0.25*alpha[0]*(1.732050807568877*f0r[5]-1.0*f0r[3])*dfac_z; 
  incr[4] = 0.25*alpha[0]*(3.0*f0r[4]-1.732050807568877*f0r[2])*dfac_z; 
  incr[5] = 0.25*alpha[0]*(3.0*f0r[5]-1.732050807568877*f0r[3])*dfac_z; 
  incr[6] = -0.25*alpha[0]*(1.732050807568877*f0r[7]-1.0*f0r[6])*dfac_z; 
  incr[7] = 0.25*alpha[0]*(3.0*f0r[7]-1.732050807568877*f0r[6])*dfac_z; 
  }
#elif upwindType == QUAD 
//double fupwind[4];
//double fupwindQuad[4];
//double alphaQuad;
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[0] = 0.5*((0.6123724356957944*(f0r[7]+f0l[7])-0.3535533905932737*f0r[6]+0.3535533905932737*f0l[6]-0.6123724356957944*(f0r[5]+f0l[5]+f0r[4]+f0l[4])+0.3535533905932737*(f0r[3]+f0r[2])-0.3535533905932737*(f0l[3]+f0l[2])+0.6123724356957944*(f0r[1]+f0l[1])-0.3535533905932737*f0r[0]+0.3535533905932737*f0l[0])*sgn(alphaQuad)-0.6123724356957944*f0r[7]+0.6123724356957944*f0l[7]+0.3535533905932737*(f0r[6]+f0l[6])+0.6123724356957944*(f0r[5]+f0r[4])-0.6123724356957944*(f0l[5]+f0l[4])-0.3535533905932737*(f0r[3]+f0l[3]+f0r[2]+f0l[2])-0.6123724356957944*f0r[1]+0.6123724356957944*f0l[1]+0.3535533905932737*(f0r[0]+f0l[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-0.6123724356957944*(f0r[7]+f0l[7]))+0.3535533905932737*f0r[6]-0.3535533905932737*f0l[6]-0.6123724356957944*(f0r[5]+f0l[5])+0.6123724356957944*(f0r[4]+f0l[4])+0.3535533905932737*f0r[3]-0.3535533905932737*(f0l[3]+f0r[2])+0.3535533905932737*f0l[2]+0.6123724356957944*(f0r[1]+f0l[1])-0.3535533905932737*f0r[0]+0.3535533905932737*f0l[0])*sgn(alphaQuad)+0.6123724356957944*f0r[7]-0.6123724356957944*f0l[7]-0.3535533905932737*(f0r[6]+f0l[6])+0.6123724356957944*f0r[5]-0.6123724356957944*(f0l[5]+f0r[4])+0.6123724356957944*f0l[4]-0.3535533905932737*(f0r[3]+f0l[3])+0.3535533905932737*(f0r[2]+f0l[2])-0.6123724356957944*f0r[1]+0.6123724356957944*f0l[1]+0.3535533905932737*(f0r[0]+f0l[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[2] = 0.5*(((-0.6123724356957944*(f0r[7]+f0l[7]))+0.3535533905932737*f0r[6]-0.3535533905932737*f0l[6]+0.6123724356957944*(f0r[5]+f0l[5])-0.6123724356957944*(f0r[4]+f0l[4])-0.3535533905932737*f0r[3]+0.3535533905932737*(f0l[3]+f0r[2])-0.3535533905932737*f0l[2]+0.6123724356957944*(f0r[1]+f0l[1])-0.3535533905932737*f0r[0]+0.3535533905932737*f0l[0])*sgn(alphaQuad)+0.6123724356957944*f0r[7]-0.6123724356957944*f0l[7]-0.3535533905932737*(f0r[6]+f0l[6])-0.6123724356957944*f0r[5]+0.6123724356957944*(f0l[5]+f0r[4])-0.6123724356957944*f0l[4]+0.3535533905932737*(f0r[3]+f0l[3])-0.3535533905932737*(f0r[2]+f0l[2])-0.6123724356957944*f0r[1]+0.6123724356957944*f0l[1]+0.3535533905932737*(f0r[0]+f0l[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[3] = 0.5*((0.6123724356957944*(f0r[7]+f0l[7])-0.3535533905932737*f0r[6]+0.3535533905932737*f0l[6]+0.6123724356957944*(f0r[5]+f0l[5]+f0r[4]+f0l[4])-0.3535533905932737*(f0r[3]+f0r[2])+0.3535533905932737*(f0l[3]+f0l[2])+0.6123724356957944*(f0r[1]+f0l[1])-0.3535533905932737*f0r[0]+0.3535533905932737*f0l[0])*sgn(alphaQuad)-0.6123724356957944*f0r[7]+0.6123724356957944*f0l[7]+0.3535533905932737*(f0r[6]+f0l[6])-0.6123724356957944*(f0r[5]+f0r[4])+0.6123724356957944*(f0l[5]+f0l[4])+0.3535533905932737*(f0r[3]+f0l[3]+f0r[2]+f0l[2])-0.6123724356957944*f0r[1]+0.6123724356957944*f0l[1]+0.3535533905932737*(f0r[0]+f0l[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*fupwindQuad[3]-0.5*fupwindQuad[2]+0.5*fupwindQuad[1]-0.5*fupwindQuad[0]; 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2])-0.5*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.5*fupwindQuad[3]-0.5*(fupwindQuad[2]+fupwindQuad[1])+0.5*fupwindQuad[0]; 
  alpha[0] = -(0.8660254037844386*Gradpar[0]*hamil[2]*dfac_v)/m_; 
  incr[0] = 0.3535533905932737*alpha[0]*fupwind[0]*dfac_z; 
  incr[1] = -0.6123724356957944*alpha[0]*fupwind[0]*dfac_z; 
  incr[2] = 0.3535533905932737*alpha[0]*fupwind[1]*dfac_z; 
  incr[3] = 0.3535533905932737*alpha[0]*fupwind[2]*dfac_z; 
  incr[4] = -0.6123724356957944*alpha[0]*fupwind[1]*dfac_z; 
  incr[5] = -0.6123724356957944*alpha[0]*fupwind[2]*dfac_z; 
  incr[6] = 0.3535533905932737*alpha[0]*fupwind[3]*dfac_z; 
  incr[7] = -0.6123724356957944*alpha[0]*fupwind[3]*dfac_z; 

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
double GyrokineticSurf1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *f0l, const double *f0r, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_z = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wz = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  double hamil[8]; 
  hamil[0] = (0.2357022603955158*(3.0*dfac_v2*(2.0*m_*wv2+2.828427124746191*(Bmag[0]*wm+Phi[0]*q_))+2.0*m_))/dfac_v2; 
  hamil[1] = 2.0*Phi[1]*q_; 
  hamil[2] = (1.632993161855453*m_*wv)/dfac_v; 
  hamil[3] = (1.154700538379252*Bmag[0])/dfac_m; 
  double BstarX_by_Bmag[8]; 
  double BstarY_by_Bmag[8]; 
  double BstarZ_by_Bmag[8]; 
  BstarZ_by_Bmag[0] = 2.0*Gradpar[0]; 

  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.1082531754730548*BstarZ_by_Bmag[0]*hamil[1]*dfac_z)/m_; 

  double alpha[4]; 
  alpha[0] = -(0.4330127018922193*BstarZ_by_Bmag[0]*hamil[1]*dfac_z)/m_; 
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
  } else { 
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
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[0] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)-0.6123724356957944*fr[7]+0.6123724356957944*(fl[7]+fr[6])-0.6123724356957944*fl[6]+0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)+0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6])-0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[2] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)+0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*fl[6]-0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[3] = 0.5*((0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)-0.6123724356957944*(fr[7]+fr[6])+0.6123724356957944*(fl[7]+fl[6])+0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*fupwindQuad[3]-0.5*fupwindQuad[2]+0.5*fupwindQuad[1]-0.5*fupwindQuad[0]; 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2])-0.5*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.5*fupwindQuad[3]-0.5*(fupwindQuad[2]+fupwindQuad[1])+0.5*fupwindQuad[0]; 
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

  // ensure cancellation of zeroth order terms for f0=F_M 
  // alpha == 0, so nothing to do 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *f0l, const double *f0r, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_z = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wz = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  double hamil[8]; 
  hamil[0] = (0.2357022603955158*(3.0*dfac_v2*(2.0*m_*wv2+2.828427124746191*(Bmag[0]*wm+Phi[0]*q_))+2.0*m_))/dfac_v2; 
  hamil[1] = 2.0*(Bmag[1]*wm+Phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wv)/dfac_v; 
  hamil[3] = (1.154700538379252*Bmag[0])/dfac_m; 
  hamil[5] = (1.154700538379252*Bmag[1])/dfac_m; 
  double BstarX_by_Bmag[8]; 
  double BstarY_by_Bmag[8]; 
  double BstarZ_by_Bmag[8]; 
  BstarZ_by_Bmag[0] = 2.0*Gradpar[0]; 
  BstarZ_by_Bmag[1] = 2.0*Gradpar[1]; 

  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.0625*(3.0*BstarZ_by_Bmag[1]-1.732050807568877*BstarZ_by_Bmag[0])*hamil[2]*dfac_v)/m_; 

  double alpha[4]; 
  alpha[0] = -(0.25*(3.0*BstarZ_by_Bmag[1]-1.732050807568877*BstarZ_by_Bmag[0])*hamil[2]*dfac_v)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_z; 
  incr[1] = -0.25*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_z; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_z; 
  incr[3] = 0.25*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_z; 
  incr[4] = -0.25*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_z; 
  incr[5] = -0.25*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_z; 
  incr[6] = 0.25*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_z; 
  incr[7] = -0.25*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_z; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_z; 
  incr[1] = 0.25*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_z; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_z; 
  incr[3] = -0.25*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_z; 
  incr[4] = 0.25*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_z; 
  incr[5] = 0.25*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_z; 
  incr[6] = -0.25*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_z; 
  incr[7] = 0.25*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_z; 
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
  incr[0] = 0.3535533905932737*alpha[0]*fupwind[0]*dfac_z; 
  incr[1] = -0.6123724356957944*alpha[0]*fupwind[0]*dfac_z; 
  incr[2] = 0.3535533905932737*alpha[0]*fupwind[1]*dfac_z; 
  incr[3] = 0.3535533905932737*alpha[0]*fupwind[2]*dfac_z; 
  incr[4] = -0.6123724356957944*alpha[0]*fupwind[1]*dfac_z; 
  incr[5] = -0.6123724356957944*alpha[0]*fupwind[2]*dfac_z; 
  incr[6] = 0.3535533905932737*alpha[0]*fupwind[3]*dfac_z; 
  incr[7] = -0.6123724356957944*alpha[0]*fupwind[3]*dfac_z; 

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

  // ensure cancellation of zeroth order terms for f0=F_M 
#if upwindType == SURFAVG 
  alpha[0] = (0.5*(3.0*Gradpar[1]-1.732050807568877*Gradpar[0])*hamil[2]*dfac_v)/m_; 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*f0l[1]+f0l[0])*dfac_z; 
  incr[1] = -0.25*alpha[0]*(3.0*f0l[1]+1.732050807568877*f0l[0])*dfac_z; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*f0l[4]+f0l[2])*dfac_z; 
  incr[3] = 0.25*alpha[0]*(1.732050807568877*f0l[5]+f0l[3])*dfac_z; 
  incr[4] = -0.25*alpha[0]*(3.0*f0l[4]+1.732050807568877*f0l[2])*dfac_z; 
  incr[5] = -0.25*alpha[0]*(3.0*f0l[5]+1.732050807568877*f0l[3])*dfac_z; 
  incr[6] = 0.25*alpha[0]*(1.732050807568877*f0l[7]+f0l[6])*dfac_z; 
  incr[7] = -0.25*alpha[0]*(3.0*f0l[7]+1.732050807568877*f0l[6])*dfac_z; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*f0r[1]-1.0*f0r[0])*dfac_z; 
  incr[1] = 0.25*alpha[0]*(3.0*f0r[1]-1.732050807568877*f0r[0])*dfac_z; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*f0r[4]-1.0*f0r[2])*dfac_z; 
  incr[3] = -0.25*alpha[0]*(1.732050807568877*f0r[5]-1.0*f0r[3])*dfac_z; 
  incr[4] = 0.25*alpha[0]*(3.0*f0r[4]-1.732050807568877*f0r[2])*dfac_z; 
  incr[5] = 0.25*alpha[0]*(3.0*f0r[5]-1.732050807568877*f0r[3])*dfac_z; 
  incr[6] = -0.25*alpha[0]*(1.732050807568877*f0r[7]-1.0*f0r[6])*dfac_z; 
  incr[7] = 0.25*alpha[0]*(3.0*f0r[7]-1.732050807568877*f0r[6])*dfac_z; 
  }
#elif upwindType == QUAD 
//double fupwind[4];
//double fupwindQuad[4];
//double alphaQuad;
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[0] = 0.5*((0.6123724356957944*(f0r[7]+f0l[7])-0.3535533905932737*f0r[6]+0.3535533905932737*f0l[6]-0.6123724356957944*(f0r[5]+f0l[5]+f0r[4]+f0l[4])+0.3535533905932737*(f0r[3]+f0r[2])-0.3535533905932737*(f0l[3]+f0l[2])+0.6123724356957944*(f0r[1]+f0l[1])-0.3535533905932737*f0r[0]+0.3535533905932737*f0l[0])*sgn(alphaQuad)-0.6123724356957944*f0r[7]+0.6123724356957944*f0l[7]+0.3535533905932737*(f0r[6]+f0l[6])+0.6123724356957944*(f0r[5]+f0r[4])-0.6123724356957944*(f0l[5]+f0l[4])-0.3535533905932737*(f0r[3]+f0l[3]+f0r[2]+f0l[2])-0.6123724356957944*f0r[1]+0.6123724356957944*f0l[1]+0.3535533905932737*(f0r[0]+f0l[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-0.6123724356957944*(f0r[7]+f0l[7]))+0.3535533905932737*f0r[6]-0.3535533905932737*f0l[6]-0.6123724356957944*(f0r[5]+f0l[5])+0.6123724356957944*(f0r[4]+f0l[4])+0.3535533905932737*f0r[3]-0.3535533905932737*(f0l[3]+f0r[2])+0.3535533905932737*f0l[2]+0.6123724356957944*(f0r[1]+f0l[1])-0.3535533905932737*f0r[0]+0.3535533905932737*f0l[0])*sgn(alphaQuad)+0.6123724356957944*f0r[7]-0.6123724356957944*f0l[7]-0.3535533905932737*(f0r[6]+f0l[6])+0.6123724356957944*f0r[5]-0.6123724356957944*(f0l[5]+f0r[4])+0.6123724356957944*f0l[4]-0.3535533905932737*(f0r[3]+f0l[3])+0.3535533905932737*(f0r[2]+f0l[2])-0.6123724356957944*f0r[1]+0.6123724356957944*f0l[1]+0.3535533905932737*(f0r[0]+f0l[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[2] = 0.5*(((-0.6123724356957944*(f0r[7]+f0l[7]))+0.3535533905932737*f0r[6]-0.3535533905932737*f0l[6]+0.6123724356957944*(f0r[5]+f0l[5])-0.6123724356957944*(f0r[4]+f0l[4])-0.3535533905932737*f0r[3]+0.3535533905932737*(f0l[3]+f0r[2])-0.3535533905932737*f0l[2]+0.6123724356957944*(f0r[1]+f0l[1])-0.3535533905932737*f0r[0]+0.3535533905932737*f0l[0])*sgn(alphaQuad)+0.6123724356957944*f0r[7]-0.6123724356957944*f0l[7]-0.3535533905932737*(f0r[6]+f0l[6])-0.6123724356957944*f0r[5]+0.6123724356957944*(f0l[5]+f0r[4])-0.6123724356957944*f0l[4]+0.3535533905932737*(f0r[3]+f0l[3])-0.3535533905932737*(f0r[2]+f0l[2])-0.6123724356957944*f0r[1]+0.6123724356957944*f0l[1]+0.3535533905932737*(f0r[0]+f0l[0])); 
  alphaQuad = 0.5*alpha[0]; 
  fupwindQuad[3] = 0.5*((0.6123724356957944*(f0r[7]+f0l[7])-0.3535533905932737*f0r[6]+0.3535533905932737*f0l[6]+0.6123724356957944*(f0r[5]+f0l[5]+f0r[4]+f0l[4])-0.3535533905932737*(f0r[3]+f0r[2])+0.3535533905932737*(f0l[3]+f0l[2])+0.6123724356957944*(f0r[1]+f0l[1])-0.3535533905932737*f0r[0]+0.3535533905932737*f0l[0])*sgn(alphaQuad)-0.6123724356957944*f0r[7]+0.6123724356957944*f0l[7]+0.3535533905932737*(f0r[6]+f0l[6])-0.6123724356957944*(f0r[5]+f0r[4])+0.6123724356957944*(f0l[5]+f0l[4])+0.3535533905932737*(f0r[3]+f0l[3]+f0r[2]+f0l[2])-0.6123724356957944*f0r[1]+0.6123724356957944*f0l[1]+0.3535533905932737*(f0r[0]+f0l[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*fupwindQuad[3]-0.5*fupwindQuad[2]+0.5*fupwindQuad[1]-0.5*fupwindQuad[0]; 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2])-0.5*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.5*fupwindQuad[3]-0.5*(fupwindQuad[2]+fupwindQuad[1])+0.5*fupwindQuad[0]; 
  alpha[0] = (0.5*(3.0*Gradpar[1]-1.732050807568877*Gradpar[0])*hamil[2]*dfac_v)/m_; 
  incr[0] = 0.3535533905932737*alpha[0]*fupwind[0]*dfac_z; 
  incr[1] = -0.6123724356957944*alpha[0]*fupwind[0]*dfac_z; 
  incr[2] = 0.3535533905932737*alpha[0]*fupwind[1]*dfac_z; 
  incr[3] = 0.3535533905932737*alpha[0]*fupwind[2]*dfac_z; 
  incr[4] = -0.6123724356957944*alpha[0]*fupwind[1]*dfac_z; 
  incr[5] = -0.6123724356957944*alpha[0]*fupwind[2]*dfac_z; 
  incr[6] = 0.3535533905932737*alpha[0]*fupwind[3]*dfac_z; 
  incr[7] = -0.6123724356957944*alpha[0]*fupwind[3]*dfac_z; 

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
double GyrokineticSurf1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *f0l, const double *f0r, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_z = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wz = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  double hamil[8]; 
  hamil[0] = (0.2357022603955158*(3.0*dfac_v2*(2.0*m_*wv2+2.828427124746191*(Bmag[0]*wm+Phi[0]*q_))+2.0*m_))/dfac_v2; 
  hamil[1] = 2.0*(Bmag[1]*wm+Phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wv)/dfac_v; 
  hamil[3] = (1.154700538379252*Bmag[0])/dfac_m; 
  hamil[5] = (1.154700538379252*Bmag[1])/dfac_m; 
  double BstarX_by_Bmag[8]; 
  double BstarY_by_Bmag[8]; 
  double BstarZ_by_Bmag[8]; 
  BstarZ_by_Bmag[0] = 2.0*Gradpar[0]; 
  BstarZ_by_Bmag[1] = 2.0*Gradpar[1]; 

  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.1082531754730548*BstarZ_by_Bmag[0]*hamil[1]*dfac_z)/m_; 

  double alpha[4]; 
  alpha[0] = -(0.4330127018922193*BstarZ_by_Bmag[0]*hamil[1]*dfac_z)/m_; 
  alpha[1] = -(0.4330127018922193*BstarZ_by_Bmag[1]*hamil[1]*dfac_z)/m_; 
  alpha[2] = -(0.4330127018922193*BstarZ_by_Bmag[0]*hamil[5]*dfac_z)/m_; 
  alpha[3] = -(0.4330127018922193*BstarZ_by_Bmag[1]*hamil[5]*dfac_z)/m_; 
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
  } else { 
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
  alphaQuad = 0.5*alpha[3]-0.5*(alpha[2]+alpha[1])+0.5*alpha[0]; 
  fupwindQuad[0] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)-0.6123724356957944*fr[7]+0.6123724356957944*(fl[7]+fr[6])-0.6123724356957944*fl[6]+0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*(alpha[1]+alpha[0])-0.5*(alpha[3]+alpha[2]); 
  fupwindQuad[1] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)+0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6])-0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaQuad = (-0.5*alpha[3])+0.5*alpha[2]-0.5*alpha[1]+0.5*alpha[0]; 
  fupwindQuad[2] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaQuad)+0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*fl[6]-0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 0.5*(alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  fupwindQuad[3] = 0.5*((0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaQuad)-0.6123724356957944*(fr[7]+fr[6])+0.6123724356957944*(fl[7]+fl[6])+0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*fupwindQuad[3]-0.5*fupwindQuad[2]+0.5*fupwindQuad[1]-0.5*fupwindQuad[0]; 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2])-0.5*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.5*fupwindQuad[3]-0.5*(fupwindQuad[2]+fupwindQuad[1])+0.5*fupwindQuad[0]; 
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

  // ensure cancellation of zeroth order terms for f0=F_M 
#if upwindType == SURFAVG 
  alpha[0] = (1.732050807568877*Gradpar[0]*Bmag[1]*dfac_z*wm)/m_; 
  alpha[1] = (1.732050807568877*Bmag[1]*Gradpar[1]*dfac_z*wm)/m_; 
  alpha[2] = (Gradpar[0]*Bmag[1]*dfac_z)/(dfac_m*m_); 
  alpha[3] = (Bmag[1]*Gradpar[1]*dfac_z)/(dfac_m*m_); 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*(alpha[3]*f0l[7]+alpha[2]*f0l[6])+alpha[3]*f0l[5]+1.732050807568877*alpha[1]*f0l[4]+alpha[2]*f0l[3]+1.732050807568877*alpha[0]*f0l[2]+alpha[1]*f0l[1]+alpha[0]*f0l[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*(alpha[2]*f0l[7]+alpha[3]*f0l[6])+alpha[2]*f0l[5]+1.732050807568877*alpha[0]*f0l[4]+alpha[3]*f0l[3]+1.732050807568877*alpha[1]*f0l[2]+alpha[0]*f0l[1]+f0l[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*(alpha[3]*f0l[7]+alpha[2]*f0l[6])+1.732050807568877*alpha[3]*f0l[5]+3.0*alpha[1]*f0l[4]+1.732050807568877*alpha[2]*f0l[3]+3.0*alpha[0]*f0l[2]+1.732050807568877*(alpha[1]*f0l[1]+alpha[0]*f0l[0]))*dfac_v; 
  incr[3] = 0.25*(1.732050807568877*(alpha[1]*f0l[7]+alpha[0]*f0l[6])+alpha[1]*f0l[5]+1.732050807568877*alpha[3]*f0l[4]+alpha[0]*f0l[3]+f0l[1]*alpha[3]+alpha[2]*(1.732050807568877*f0l[2]+f0l[0]))*dfac_v; 
  incr[4] = -0.25*(3.0*(alpha[2]*f0l[7]+alpha[3]*f0l[6])+1.732050807568877*alpha[2]*f0l[5]+3.0*alpha[0]*f0l[4]+1.732050807568877*alpha[3]*f0l[3]+3.0*alpha[1]*f0l[2]+1.732050807568877*(alpha[0]*f0l[1]+f0l[0]*alpha[1]))*dfac_v; 
  incr[5] = 0.25*(1.732050807568877*(alpha[0]*f0l[7]+alpha[1]*f0l[6])+alpha[0]*f0l[5]+1.732050807568877*alpha[2]*f0l[4]+alpha[1]*f0l[3]+(1.732050807568877*f0l[2]+f0l[0])*alpha[3]+f0l[1]*alpha[2])*dfac_v; 
  incr[6] = -0.25*(3.0*(alpha[1]*f0l[7]+alpha[0]*f0l[6])+1.732050807568877*alpha[1]*f0l[5]+3.0*alpha[3]*f0l[4]+1.732050807568877*(alpha[0]*f0l[3]+f0l[1]*alpha[3])+alpha[2]*(3.0*f0l[2]+1.732050807568877*f0l[0]))*dfac_v; 
  incr[7] = -0.25*(3.0*(alpha[0]*f0l[7]+alpha[1]*f0l[6])+1.732050807568877*alpha[0]*f0l[5]+3.0*alpha[2]*f0l[4]+1.732050807568877*alpha[1]*f0l[3]+3.0*f0l[2]*alpha[3]+1.732050807568877*(f0l[0]*alpha[3]+f0l[1]*alpha[2]))*dfac_v; 
  } else { 
  incr[0] = -0.25*(1.732050807568877*(alpha[3]*f0r[7]+alpha[2]*f0r[6])-1.0*alpha[3]*f0r[5]+1.732050807568877*alpha[1]*f0r[4]-1.0*alpha[2]*f0r[3]+1.732050807568877*alpha[0]*f0r[2]-1.0*(alpha[1]*f0r[1]+alpha[0]*f0r[0]))*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*(alpha[2]*f0r[7]+alpha[3]*f0r[6])-1.0*alpha[2]*f0r[5]+1.732050807568877*alpha[0]*f0r[4]-1.0*alpha[3]*f0r[3]+1.732050807568877*alpha[1]*f0r[2]-1.0*(alpha[0]*f0r[1]+f0r[0]*alpha[1]))*dfac_v; 
  incr[2] = 0.25*(3.0*(alpha[3]*f0r[7]+alpha[2]*f0r[6])-1.732050807568877*alpha[3]*f0r[5]+3.0*alpha[1]*f0r[4]-1.732050807568877*alpha[2]*f0r[3]+3.0*alpha[0]*f0r[2]-1.732050807568877*(alpha[1]*f0r[1]+alpha[0]*f0r[0]))*dfac_v; 
  incr[3] = -0.25*(1.732050807568877*(alpha[1]*f0r[7]+alpha[0]*f0r[6])-1.0*alpha[1]*f0r[5]+1.732050807568877*alpha[3]*f0r[4]-1.0*(alpha[0]*f0r[3]+f0r[1]*alpha[3])+alpha[2]*(1.732050807568877*f0r[2]-1.0*f0r[0]))*dfac_v; 
  incr[4] = 0.25*(3.0*(alpha[2]*f0r[7]+alpha[3]*f0r[6])-1.732050807568877*alpha[2]*f0r[5]+3.0*alpha[0]*f0r[4]-1.732050807568877*alpha[3]*f0r[3]+3.0*alpha[1]*f0r[2]-1.732050807568877*(alpha[0]*f0r[1]+f0r[0]*alpha[1]))*dfac_v; 
  incr[5] = -0.25*(1.732050807568877*(alpha[0]*f0r[7]+alpha[1]*f0r[6])-1.0*alpha[0]*f0r[5]+1.732050807568877*alpha[2]*f0r[4]-1.0*alpha[1]*f0r[3]+1.732050807568877*f0r[2]*alpha[3]-1.0*(f0r[0]*alpha[3]+f0r[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.25*(3.0*(alpha[1]*f0r[7]+alpha[0]*f0r[6])-1.732050807568877*alpha[1]*f0r[5]+3.0*alpha[3]*f0r[4]-1.732050807568877*(alpha[0]*f0r[3]+f0r[1]*alpha[3])+alpha[2]*(3.0*f0r[2]-1.732050807568877*f0r[0]))*dfac_v; 
  incr[7] = 0.25*(3.0*(alpha[0]*f0r[7]+alpha[1]*f0r[6])-1.732050807568877*alpha[0]*f0r[5]+3.0*alpha[2]*f0r[4]-1.732050807568877*alpha[1]*f0r[3]+3.0*f0r[2]*alpha[3]-1.732050807568877*(f0r[0]*alpha[3]+f0r[1]*alpha[2]))*dfac_v; 
  }
#elif upwindType == QUAD 
//double fupwind[4];
//double fupwindQuad[4];
//double alphaQuad;
  alphaQuad = 0.5*alpha[3]-0.5*(alpha[2]+alpha[1])+0.5*alpha[0]; 
  fupwindQuad[0] = 0.5*((0.6123724356957944*(f0r[7]+f0l[7])-0.6123724356957944*(f0r[6]+f0l[6])-0.3535533905932737*f0r[5]+0.3535533905932737*f0l[5]-0.6123724356957944*(f0r[4]+f0l[4])+0.3535533905932737*f0r[3]-0.3535533905932737*f0l[3]+0.6123724356957944*(f0r[2]+f0l[2])+0.3535533905932737*f0r[1]-0.3535533905932737*(f0l[1]+f0r[0])+0.3535533905932737*f0l[0])*sgn(alphaQuad)-0.6123724356957944*f0r[7]+0.6123724356957944*(f0l[7]+f0r[6])-0.6123724356957944*f0l[6]+0.3535533905932737*(f0r[5]+f0l[5])+0.6123724356957944*f0r[4]-0.6123724356957944*f0l[4]-0.3535533905932737*(f0r[3]+f0l[3])-0.6123724356957944*f0r[2]+0.6123724356957944*f0l[2]-0.3535533905932737*(f0r[1]+f0l[1])+0.3535533905932737*(f0r[0]+f0l[0])); 
  alphaQuad = 0.5*(alpha[1]+alpha[0])-0.5*(alpha[3]+alpha[2]); 
  fupwindQuad[1] = 0.5*(((-0.6123724356957944*(f0r[7]+f0l[7]+f0r[6]+f0l[6]))+0.3535533905932737*f0r[5]-0.3535533905932737*f0l[5]+0.6123724356957944*(f0r[4]+f0l[4])+0.3535533905932737*f0r[3]-0.3535533905932737*f0l[3]+0.6123724356957944*(f0r[2]+f0l[2])-0.3535533905932737*(f0r[1]+f0r[0])+0.3535533905932737*(f0l[1]+f0l[0]))*sgn(alphaQuad)+0.6123724356957944*(f0r[7]+f0r[6])-0.6123724356957944*(f0l[7]+f0l[6])-0.3535533905932737*(f0r[5]+f0l[5])-0.6123724356957944*f0r[4]+0.6123724356957944*f0l[4]-0.3535533905932737*(f0r[3]+f0l[3])-0.6123724356957944*f0r[2]+0.6123724356957944*f0l[2]+0.3535533905932737*(f0r[1]+f0l[1]+f0r[0]+f0l[0])); 
  alphaQuad = (-0.5*alpha[3])+0.5*alpha[2]-0.5*alpha[1]+0.5*alpha[0]; 
  fupwindQuad[2] = 0.5*(((-0.6123724356957944*(f0r[7]+f0l[7]))+0.6123724356957944*(f0r[6]+f0l[6])+0.3535533905932737*f0r[5]-0.3535533905932737*f0l[5]-0.6123724356957944*(f0r[4]+f0l[4])-0.3535533905932737*f0r[3]+0.3535533905932737*f0l[3]+0.6123724356957944*(f0r[2]+f0l[2])+0.3535533905932737*f0r[1]-0.3535533905932737*(f0l[1]+f0r[0])+0.3535533905932737*f0l[0])*sgn(alphaQuad)+0.6123724356957944*f0r[7]-0.6123724356957944*(f0l[7]+f0r[6])+0.6123724356957944*f0l[6]-0.3535533905932737*(f0r[5]+f0l[5])+0.6123724356957944*f0r[4]-0.6123724356957944*f0l[4]+0.3535533905932737*(f0r[3]+f0l[3])-0.6123724356957944*f0r[2]+0.6123724356957944*f0l[2]-0.3535533905932737*(f0r[1]+f0l[1])+0.3535533905932737*(f0r[0]+f0l[0])); 
  alphaQuad = 0.5*(alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  fupwindQuad[3] = 0.5*((0.6123724356957944*(f0r[7]+f0l[7]+f0r[6]+f0l[6])-0.3535533905932737*f0r[5]+0.3535533905932737*f0l[5]+0.6123724356957944*(f0r[4]+f0l[4])-0.3535533905932737*f0r[3]+0.3535533905932737*f0l[3]+0.6123724356957944*(f0r[2]+f0l[2])-0.3535533905932737*(f0r[1]+f0r[0])+0.3535533905932737*(f0l[1]+f0l[0]))*sgn(alphaQuad)-0.6123724356957944*(f0r[7]+f0r[6])+0.6123724356957944*(f0l[7]+f0l[6])+0.3535533905932737*(f0r[5]+f0l[5])-0.6123724356957944*f0r[4]+0.6123724356957944*f0l[4]+0.3535533905932737*(f0r[3]+f0l[3])-0.6123724356957944*f0r[2]+0.6123724356957944*f0l[2]+0.3535533905932737*(f0r[1]+f0l[1]+f0r[0]+f0l[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*fupwindQuad[3]-0.5*fupwindQuad[2]+0.5*fupwindQuad[1]-0.5*fupwindQuad[0]; 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2])-0.5*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.5*fupwindQuad[3]-0.5*(fupwindQuad[2]+fupwindQuad[1])+0.5*fupwindQuad[0]; 
  alpha[0] = (1.732050807568877*Gradpar[0]*Bmag[1]*dfac_z*wm)/m_; 
  alpha[1] = (1.732050807568877*Bmag[1]*Gradpar[1]*dfac_z*wm)/m_; 
  alpha[2] = (Gradpar[0]*Bmag[1]*dfac_z)/(dfac_m*m_); 
  alpha[3] = (Bmag[1]*Gradpar[1]*dfac_z)/(dfac_m*m_); 
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
return std::abs(alpha0); 
} 
