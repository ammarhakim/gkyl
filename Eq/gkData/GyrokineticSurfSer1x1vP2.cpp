#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x1vSer_X_P2_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*Gradpar[0]*wv; 

  double alpha[3]; 
  alpha[0] = Gradpar[0]*wv; 
  alpha[1] = (0.5773502691896258*Gradpar[0])/dfac_v; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*alpha[0]*(2.23606797749979*fl[4]+1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.3535533905932737*alpha[0]*(3.872983346207417*fl[4]+3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.1178511301977579*alpha[0]*(6.708203932499369*fl[6]+5.196152422706631*fl[3]+3.0*fl[2])*dfac_x; 
  incr[3] = -0.3535533905932737*alpha[0]*(3.872983346207417*fl[6]+3.0*fl[3]+1.732050807568877*fl[2])*dfac_x; 
  incr[4] = 0.3535533905932737*alpha[0]*(5.0*fl[4]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*dfac_x; 
  incr[5] = 0.07071067811865474*alpha[0]*(8.660254037844387*fl[7]+5.0*fl[5])*dfac_x; 
  incr[6] = 0.1178511301977579*alpha[0]*(15.0*fl[6]+11.61895003862225*fl[3]+6.708203932499369*fl[2])*dfac_x; 
  incr[7] = -0.07071067811865474*alpha[0]*(15.0*fl[7]+8.660254037844387*fl[5])*dfac_x; 
  } else { 
  incr[0] = 0.3535533905932737*alpha[0]*(2.23606797749979*fr[4]-1.732050807568877*fr[1]+fr[0])*dfac_x; 
  incr[1] = -0.3535533905932737*alpha[0]*(3.872983346207417*fr[4]-3.0*fr[1]+1.732050807568877*fr[0])*dfac_x; 
  incr[2] = 0.1178511301977579*alpha[0]*(6.708203932499369*fr[6]-5.196152422706631*fr[3]+3.0*fr[2])*dfac_x; 
  incr[3] = -0.3535533905932737*alpha[0]*(3.872983346207417*fr[6]-3.0*fr[3]+1.732050807568877*fr[2])*dfac_x; 
  incr[4] = 0.3535533905932737*alpha[0]*(5.0*fr[4]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*dfac_x; 
  incr[5] = -0.07071067811865474*alpha[0]*(8.660254037844387*fr[7]-5.0*fr[5])*dfac_x; 
  incr[6] = 0.1178511301977579*alpha[0]*(15.0*fr[6]-11.61895003862225*fr[3]+6.708203932499369*fr[2])*dfac_x; 
  incr[7] = 0.07071067811865474*alpha[0]*(15.0*fr[7]-8.660254037844387*fr[5])*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[3];
double fupwindQuad[3];
double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[0] = 0.5*((0.7745966692414833*(fr[7]+fl[7])+1.5*fr[6]-1.5*fl[6]-0.4472135954999579*fr[5]+0.4472135954999579*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]-1.161895003862225*(fr[3]+fl[3])+0.6708203932499369*fr[2]-0.6708203932499369*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.7745966692414833*fr[7]+0.7745966692414833*fl[7]-1.5*(fr[6]+fl[6])+0.4472135954999579*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])+1.161895003862225*fr[3]-1.161895003862225*fl[3]-0.6708203932499369*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-0.9682458365518543*(fr[7]+fl[7]))+0.5590169943749475*fr[5]-0.5590169943749475*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)+0.9682458365518543*fr[7]-0.9682458365518543*fl[7]-0.5590169943749475*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[2] = 0.5*((0.7745966692414833*(fr[7]+fl[7])-1.5*fr[6]+1.5*fl[6]-0.4472135954999579*fr[5]+0.4472135954999579*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]+1.161895003862225*(fr[3]+fl[3])-0.6708203932499369*fr[2]+0.6708203932499369*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.7745966692414833*fr[7]+0.7745966692414833*fl[7]+1.5*(fr[6]+fl[6])+0.4472135954999579*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])-1.161895003862225*fr[3]+1.161895003862225*fl[3]+0.6708203932499369*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  fupwind[0] = 0.3928371006591931*fupwindQuad[2]+0.6285393610547092*fupwindQuad[1]+0.3928371006591931*fupwindQuad[0]; 
  fupwind[1] = 0.52704627669473*fupwindQuad[2]-0.52704627669473*fupwindQuad[0]; 
  fupwind[2] = 0.3513641844631533*fupwindQuad[2]-0.7027283689263066*fupwindQuad[1]+0.3513641844631533*fupwindQuad[0]; 
  incr[0] = 0.5*alpha[0]*fupwind[0]*dfac_x; 
  incr[1] = -0.8660254037844386*alpha[0]*fupwind[0]*dfac_x; 
  incr[2] = 0.5*alpha[0]*fupwind[1]*dfac_x; 
  incr[3] = -0.8660254037844386*alpha[0]*fupwind[1]*dfac_x; 
  incr[4] = 1.118033988749895*alpha[0]*fupwind[0]*dfac_x; 
  incr[5] = 0.5*alpha[0]*fupwind[2]*dfac_x; 
  incr[6] = 1.118033988749895*alpha[0]*fupwind[1]*dfac_x; 
  incr[7] = -0.8660254037844387*alpha[0]*fupwind[2]*dfac_x; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.4330127018922193*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[3]; 
  alpha[0] = -(1.224744871391589*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(2.738612787525831*Gradpar[0]*Phi[2]*dfac_x*q_)/m_; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[3]; 
  double favg[3]; 
  favg[0] = 0.7071067811865475*(2.23606797749979*(fr[5]+fl[5])+1.732050807568877*(fl[2]-1.0*fr[2])+fr[0]+fl[0]); 
  favg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[7]+fl[7])+3.0*(fl[3]-1.0*fr[3]))+3.0*(fr[1]+fl[1])); 
  favg[2] = -0.1414213562373095*(8.660254037844387*fr[6]-1.0*(8.660254037844387*fl[6]+5.0*(fr[4]+fl[4]))); 

  Ghat[0] = -0.25*((3.16227766016838*fr[5]-3.16227766016838*fl[5]-2.449489742783178*(fr[2]+fl[2])+1.414213562373095*fr[0]-1.414213562373095*fl[0])*amax-1.414213562373095*alpha[0]*favg[0]); 
  Ghat[1] = -0.08333333333333333*((9.48683298050514*fr[7]-9.48683298050514*fl[7]-7.348469228349534*(fr[3]+fl[3])+4.242640687119286*fr[1]-4.242640687119286*fl[1])*amax-4.242640687119286*alpha[0]*favg[1]); 
  Ghat[2] = 0.05*((12.24744871391589*(fr[6]+fl[6])-7.071067811865476*fr[4]+7.071067811865476*fl[4])*amax+7.071067811865476*alpha[0]*favg[2]); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_v; 
  incr[1] = 0.7071067811865475*Ghat[1]*dfac_v; 
  incr[2] = -1.224744871391589*Ghat[0]*dfac_v; 
  incr[3] = -1.224744871391589*Ghat[1]*dfac_v; 
  incr[4] = 0.7071067811865475*Ghat[2]*dfac_v; 
  incr[5] = 1.58113883008419*Ghat[0]*dfac_v; 
  incr[6] = -1.224744871391589*Ghat[2]*dfac_v; 
  incr[7] = 1.58113883008419*Ghat[1]*dfac_v; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_X_P2_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*(2.23606797749979*Gradpar[2]-1.732050807568877*Gradpar[1]+Gradpar[0])*wv; 

  double alpha[3]; 
  alpha[0] = (2.23606797749979*Gradpar[2]-1.732050807568877*Gradpar[1]+Gradpar[0])*wv; 
  alpha[1] = (0.3333333333333333*(3.872983346207417*Gradpar[2]-3.0*Gradpar[1]+1.732050807568877*Gradpar[0]))/dfac_v; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*alpha[0]*(2.23606797749979*fl[4]+1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.3535533905932737*alpha[0]*(3.872983346207417*fl[4]+3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.1178511301977579*alpha[0]*(6.708203932499369*fl[6]+5.196152422706631*fl[3]+3.0*fl[2])*dfac_x; 
  incr[3] = -0.3535533905932737*alpha[0]*(3.872983346207417*fl[6]+3.0*fl[3]+1.732050807568877*fl[2])*dfac_x; 
  incr[4] = 0.3535533905932737*alpha[0]*(5.0*fl[4]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*dfac_x; 
  incr[5] = 0.07071067811865474*alpha[0]*(8.660254037844387*fl[7]+5.0*fl[5])*dfac_x; 
  incr[6] = 0.1178511301977579*alpha[0]*(15.0*fl[6]+11.61895003862225*fl[3]+6.708203932499369*fl[2])*dfac_x; 
  incr[7] = -0.07071067811865474*alpha[0]*(15.0*fl[7]+8.660254037844387*fl[5])*dfac_x; 
  } else { 
  incr[0] = 0.3535533905932737*alpha[0]*(2.23606797749979*fr[4]-1.732050807568877*fr[1]+fr[0])*dfac_x; 
  incr[1] = -0.3535533905932737*alpha[0]*(3.872983346207417*fr[4]-3.0*fr[1]+1.732050807568877*fr[0])*dfac_x; 
  incr[2] = 0.1178511301977579*alpha[0]*(6.708203932499369*fr[6]-5.196152422706631*fr[3]+3.0*fr[2])*dfac_x; 
  incr[3] = -0.3535533905932737*alpha[0]*(3.872983346207417*fr[6]-3.0*fr[3]+1.732050807568877*fr[2])*dfac_x; 
  incr[4] = 0.3535533905932737*alpha[0]*(5.0*fr[4]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*dfac_x; 
  incr[5] = -0.07071067811865474*alpha[0]*(8.660254037844387*fr[7]-5.0*fr[5])*dfac_x; 
  incr[6] = 0.1178511301977579*alpha[0]*(15.0*fr[6]-11.61895003862225*fr[3]+6.708203932499369*fr[2])*dfac_x; 
  incr[7] = 0.07071067811865474*alpha[0]*(15.0*fr[7]-8.660254037844387*fr[5])*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[3];
double fupwindQuad[3];
double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[0] = 0.5*((0.7745966692414833*(fr[7]+fl[7])+1.5*fr[6]-1.5*fl[6]-0.4472135954999579*fr[5]+0.4472135954999579*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]-1.161895003862225*(fr[3]+fl[3])+0.6708203932499369*fr[2]-0.6708203932499369*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.7745966692414833*fr[7]+0.7745966692414833*fl[7]-1.5*(fr[6]+fl[6])+0.4472135954999579*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])+1.161895003862225*fr[3]-1.161895003862225*fl[3]-0.6708203932499369*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-0.9682458365518543*(fr[7]+fl[7]))+0.5590169943749475*fr[5]-0.5590169943749475*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)+0.9682458365518543*fr[7]-0.9682458365518543*fl[7]-0.5590169943749475*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[2] = 0.5*((0.7745966692414833*(fr[7]+fl[7])-1.5*fr[6]+1.5*fl[6]-0.4472135954999579*fr[5]+0.4472135954999579*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]+1.161895003862225*(fr[3]+fl[3])-0.6708203932499369*fr[2]+0.6708203932499369*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.7745966692414833*fr[7]+0.7745966692414833*fl[7]+1.5*(fr[6]+fl[6])+0.4472135954999579*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])-1.161895003862225*fr[3]+1.161895003862225*fl[3]+0.6708203932499369*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  fupwind[0] = 0.3928371006591931*fupwindQuad[2]+0.6285393610547092*fupwindQuad[1]+0.3928371006591931*fupwindQuad[0]; 
  fupwind[1] = 0.52704627669473*fupwindQuad[2]-0.52704627669473*fupwindQuad[0]; 
  fupwind[2] = 0.3513641844631533*fupwindQuad[2]-0.7027283689263066*fupwindQuad[1]+0.3513641844631533*fupwindQuad[0]; 
  incr[0] = 0.5*alpha[0]*fupwind[0]*dfac_x; 
  incr[1] = -0.8660254037844386*alpha[0]*fupwind[0]*dfac_x; 
  incr[2] = 0.5*alpha[0]*fupwind[1]*dfac_x; 
  incr[3] = -0.8660254037844386*alpha[0]*fupwind[1]*dfac_x; 
  incr[4] = 1.118033988749895*alpha[0]*fupwind[0]*dfac_x; 
  incr[5] = 0.5*alpha[0]*fupwind[2]*dfac_x; 
  incr[6] = 1.118033988749895*alpha[0]*fupwind[1]*dfac_x; 
  incr[7] = -0.8660254037844387*alpha[0]*fupwind[2]*dfac_x; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.4330127018922193*(2.23606797749979*Gradpar[1]*Phi[2]+Gradpar[0]*Phi[1])*dfac_x*q_)/m_; 

  double alpha[3]; 
  alpha[0] = -(0.7071067811865475*(3.872983346207417*Gradpar[1]*Phi[2]+1.732050807568877*Gradpar[0]*Phi[1])*dfac_x*q_)/m_; 
  alpha[1] = -(0.7071067811865475*((3.464101615137754*Gradpar[2]+3.872983346207417*Gradpar[0])*Phi[2]+1.732050807568877*Gradpar[1]*Phi[1])*dfac_x*q_)/m_; 
  alpha[2] = -(0.7071067811865475*(3.464101615137754*Gradpar[1]*Phi[2]+1.732050807568877*Phi[1]*Gradpar[2])*dfac_x*q_)/m_; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[3]; 
  double favg[3]; 
  favg[0] = 0.7071067811865475*(2.23606797749979*(fr[5]+fl[5])+1.732050807568877*(fl[2]-1.0*fr[2])+fr[0]+fl[0]); 
  favg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[7]+fl[7])+3.0*(fl[3]-1.0*fr[3]))+3.0*(fr[1]+fl[1])); 
  favg[2] = -0.1414213562373095*(8.660254037844387*fr[6]-1.0*(8.660254037844387*fl[6]+5.0*(fr[4]+fl[4]))); 

  Ghat[0] = -0.125*((6.324555320336761*fr[5]-6.324555320336761*fl[5]-4.898979485566357*(fr[2]+fl[2])+2.828427124746191*fr[0]-2.828427124746191*fl[0])*amax+favg[0]*(3.16227766016838*alpha[2]-2.828427124746191*alpha[0])); 
  Ghat[1] = -0.04166666666666666*((18.97366596101028*fr[7]-18.97366596101028*fl[7]-14.69693845669907*(fr[3]+fl[3])+8.485281374238571*fr[1]-8.485281374238571*fl[1])*amax+favg[1]*(9.48683298050514*alpha[2]-8.485281374238571*alpha[0])); 
  Ghat[2] = 0.025*((24.49489742783179*(fr[6]+fl[6])-14.14213562373095*fr[4]+14.14213562373095*fl[4])*amax+(14.14213562373095*alpha[0]-15.8113883008419*alpha[2])*favg[2]); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_v; 
  incr[1] = 0.7071067811865475*Ghat[1]*dfac_v; 
  incr[2] = -1.224744871391589*Ghat[0]*dfac_v; 
  incr[3] = -1.224744871391589*Ghat[1]*dfac_v; 
  incr[4] = 0.7071067811865475*Ghat[2]*dfac_v; 
  incr[5] = 1.58113883008419*Ghat[0]*dfac_v; 
  incr[6] = -1.224744871391589*Ghat[2]*dfac_v; 
  incr[7] = 1.58113883008419*Ghat[1]*dfac_v; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
return std::abs(alpha0); 
} 
