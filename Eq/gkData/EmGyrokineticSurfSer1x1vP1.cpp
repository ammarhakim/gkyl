#include <GyrokineticModDecl.h> 
double EmGyrokineticSurf1x1vSer_X_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[5]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*Gradpar[0]*wv; 

  double alpha[3]; 
  alpha[0] = Gradpar[0]*wv; 
  alpha[1] = (0.5773502691896258*Gradpar[0])/dfac_v; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*(alpha[1]*(1.732050807568877*fl[3]+fl[2])+alpha[0]*(1.732050807568877*fl[1]+fl[0]))*dfac_x; 
  incr[1] = -0.3535533905932737*(alpha[1]*(3.0*fl[3]+1.732050807568877*fl[2])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[2] = 0.07071067811865474*(4.47213595499958*alpha[1]*fl[4]+alpha[0]*(8.660254037844386*fl[3]+5.0*fl[2])+alpha[1]*(8.660254037844386*fl[1]+5.0*fl[0]))*dfac_x; 
  incr[3] = -0.07071067811865474*(7.745966692414834*alpha[1]*fl[4]+alpha[0]*(15.0*fl[3]+8.660254037844386*fl[2])+alpha[1]*(15.0*fl[1]+8.660254037844386*fl[0]))*dfac_x; 
  incr[4] = 0.07071067811865474*(5.0*alpha[0]*fl[4]+alpha[1]*(7.745966692414834*fl[3]+4.47213595499958*fl[2]))*dfac_x; 
  } else { 
  incr[0] = -0.3535533905932737*(alpha[1]*(1.732050807568877*fr[3]-1.0*fr[2])+alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac_x; 
  incr[1] = 0.3535533905932737*(alpha[1]*(3.0*fr[3]-1.732050807568877*fr[2])+alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac_x; 
  incr[2] = 0.07071067811865474*(4.47213595499958*alpha[1]*fr[4]+alpha[0]*(5.0*fr[2]-8.660254037844386*fr[3])+alpha[1]*(5.0*fr[0]-8.660254037844386*fr[1]))*dfac_x; 
  incr[3] = -0.07071067811865474*(7.745966692414834*alpha[1]*fr[4]+alpha[0]*(8.660254037844386*fr[2]-15.0*fr[3])+alpha[1]*(8.660254037844386*fr[0]-15.0*fr[1]))*dfac_x; 
  incr[4] = 0.07071067811865474*(5.0*alpha[0]*fr[4]+alpha[1]*(4.47213595499958*fr[2]-7.745966692414834*fr[3]))*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[3];
double fupwindQuad[3];
double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]-0.7071067811865475*alpha[1]; 
  fupwindQuad[0] = 0.5*(((-0.8660254037844386*(fr[3]+fl[3]))+0.5*fr[2]-0.5*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)+0.8660254037844386*fr[3]-0.8660254037844386*fl[3]-0.5*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*(alpha[1]+alpha[0]); 
  fupwindQuad[1] = 0.5*((0.8660254037844386*(fr[3]+fl[3])-0.5*fr[2]+0.5*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.8660254037844386*fr[3]+0.8660254037844386*fl[3]+0.5*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[2] = 0.5*((0.5590169943749475*fr[4]-0.5590169943749475*fl[4]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.5590169943749475*(fr[4]+fl[4])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  fupwind[0] = 0.7071067811865477*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.7071067811865477*fupwindQuad[1]-0.7071067811865477*fupwindQuad[0]; 
  fupwind[2] = 0.6324555320336761*(fupwindQuad[1]+fupwindQuad[0])-1.264911064067352*fupwindQuad[2]; 
  incr[0] = 0.5*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_x; 
  incr[1] = -0.8660254037844386*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_x; 
  incr[2] = 0.1*(4.47213595499958*alpha[1]*fupwind[2]+5.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_x; 
  incr[3] = -0.1732050807568877*(4.47213595499958*alpha[1]*fupwind[2]+5.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_x; 
  incr[4] = 0.1*(5.0*alpha[0]*fupwind[2]+4.47213595499958*alpha[1]*fupwind[1])*dfac_x; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x1vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[5]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardt[0])*q_)/m_; 

  double alpha[2]; 
  alpha[0] = -(0.7071067811865475*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardt[0])*q_)/m_; 
  alpha[1] = -(1.0*dApardt[1]*q_)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*(2.23606797749979*alpha[0]*fl[4]+1.732050807568877*(alpha[1]*fl[3]+alpha[0]*fl[2])+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(2.23606797749979*alpha[1]*fl[4]+1.732050807568877*(alpha[0]*fl[3]+alpha[1]*fl[2])+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.3535533905932737*(3.872983346207417*alpha[0]*fl[4]+3.0*(alpha[1]*fl[3]+alpha[0]*fl[2])+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[3] = -0.3535533905932737*(3.872983346207417*alpha[1]*fl[4]+3.0*(alpha[0]*fl[3]+alpha[1]*fl[2])+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[4] = 0.3535533905932737*(5.0*alpha[0]*fl[4]+3.872983346207417*(alpha[1]*fl[3]+alpha[0]*fl[2])+2.23606797749979*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  } else { 
  incr[0] = 0.3535533905932737*(2.23606797749979*alpha[0]*fr[4]-1.732050807568877*(alpha[1]*fr[3]+alpha[0]*fr[2])+alpha[1]*fr[1]+alpha[0]*fr[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(2.23606797749979*alpha[1]*fr[4]-1.732050807568877*(alpha[0]*fr[3]+alpha[1]*fr[2])+alpha[0]*fr[1]+fr[0]*alpha[1])*dfac_v; 
  incr[2] = -0.3535533905932737*(3.872983346207417*alpha[0]*fr[4]-3.0*(alpha[1]*fr[3]+alpha[0]*fr[2])+1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[3] = -0.3535533905932737*(3.872983346207417*alpha[1]*fr[4]-3.0*(alpha[0]*fr[3]+alpha[1]*fr[2])+1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[4] = 0.3535533905932737*(5.0*alpha[0]*fr[4]-3.872983346207417*(alpha[1]*fr[3]+alpha[0]*fr[2])+2.23606797749979*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[2];
double fupwindQuad[2];
double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]-0.7071067811865475*alpha[1]; 
  fupwindQuad[0] = 0.5*(((-1.118033988749895*fr[4])+1.118033988749895*fl[4]-0.8660254037844386*(fr[3]+fl[3])+0.8660254037844386*(fr[2]+fl[2])+0.5*fr[1]-0.5*(fl[1]+fr[0])+0.5*fl[0])*sgn(alphaQuad)+1.118033988749895*(fr[4]+fl[4])+0.8660254037844386*fr[3]-0.8660254037844386*(fl[3]+fr[2])+0.8660254037844386*fl[2]-0.5*(fr[1]+fl[1])+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*(alpha[1]+alpha[0]); 
  fupwindQuad[1] = 0.5*(((-1.118033988749895*fr[4])+1.118033988749895*fl[4]+0.8660254037844386*(fr[3]+fl[3]+fr[2]+fl[2])-0.5*(fr[1]+fr[0])+0.5*(fl[1]+fl[0]))*sgn(alphaQuad)+1.118033988749895*(fr[4]+fl[4])-0.8660254037844386*(fr[3]+fr[2])+0.8660254037844386*(fl[3]+fl[2])+0.5*(fr[1]+fl[1]+fr[0]+fl[0])); 
  fupwind[0] = 0.7071067811865476*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.7071067811865476*fupwindQuad[1]-0.7071067811865476*fupwindQuad[0]; 
  incr[0] = 0.5*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 
  incr[1] = 0.5*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac_v; 
  incr[2] = -0.8660254037844386*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 
  incr[3] = -0.8660254037844386*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac_v; 
  incr[4] = 1.118033988749895*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x1vSer_X_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[5]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.3535533905932737*(1.732050807568877*Gradpar[1]-1.0*Gradpar[0])*wv; 

  double alpha[3]; 
  alpha[0] = (Gradpar[0]-1.732050807568877*Gradpar[1])*wv; 
  alpha[1] = -(0.3333333333333333*(3.0*Gradpar[1]-1.732050807568877*Gradpar[0]))/dfac_v; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*(alpha[1]*(1.732050807568877*fl[3]+fl[2])+alpha[0]*(1.732050807568877*fl[1]+fl[0]))*dfac_x; 
  incr[1] = -0.3535533905932737*(alpha[1]*(3.0*fl[3]+1.732050807568877*fl[2])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[2] = 0.07071067811865474*(4.47213595499958*alpha[1]*fl[4]+alpha[0]*(8.660254037844386*fl[3]+5.0*fl[2])+alpha[1]*(8.660254037844386*fl[1]+5.0*fl[0]))*dfac_x; 
  incr[3] = -0.07071067811865474*(7.745966692414834*alpha[1]*fl[4]+alpha[0]*(15.0*fl[3]+8.660254037844386*fl[2])+alpha[1]*(15.0*fl[1]+8.660254037844386*fl[0]))*dfac_x; 
  incr[4] = 0.07071067811865474*(5.0*alpha[0]*fl[4]+alpha[1]*(7.745966692414834*fl[3]+4.47213595499958*fl[2]))*dfac_x; 
  } else { 
  incr[0] = -0.3535533905932737*(alpha[1]*(1.732050807568877*fr[3]-1.0*fr[2])+alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac_x; 
  incr[1] = 0.3535533905932737*(alpha[1]*(3.0*fr[3]-1.732050807568877*fr[2])+alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac_x; 
  incr[2] = 0.07071067811865474*(4.47213595499958*alpha[1]*fr[4]+alpha[0]*(5.0*fr[2]-8.660254037844386*fr[3])+alpha[1]*(5.0*fr[0]-8.660254037844386*fr[1]))*dfac_x; 
  incr[3] = -0.07071067811865474*(7.745966692414834*alpha[1]*fr[4]+alpha[0]*(8.660254037844386*fr[2]-15.0*fr[3])+alpha[1]*(8.660254037844386*fr[0]-15.0*fr[1]))*dfac_x; 
  incr[4] = 0.07071067811865474*(5.0*alpha[0]*fr[4]+alpha[1]*(4.47213595499958*fr[2]-7.745966692414834*fr[3]))*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[3];
double fupwindQuad[3];
double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]-0.7071067811865475*alpha[1]; 
  fupwindQuad[0] = 0.5*(((-0.8660254037844386*(fr[3]+fl[3]))+0.5*fr[2]-0.5*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)+0.8660254037844386*fr[3]-0.8660254037844386*fl[3]-0.5*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*(alpha[1]+alpha[0]); 
  fupwindQuad[1] = 0.5*((0.8660254037844386*(fr[3]+fl[3])-0.5*fr[2]+0.5*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.8660254037844386*fr[3]+0.8660254037844386*fl[3]+0.5*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[2] = 0.5*((0.5590169943749475*fr[4]-0.5590169943749475*fl[4]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.5590169943749475*(fr[4]+fl[4])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  fupwind[0] = 0.7071067811865477*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.7071067811865477*fupwindQuad[1]-0.7071067811865477*fupwindQuad[0]; 
  fupwind[2] = 0.6324555320336761*(fupwindQuad[1]+fupwindQuad[0])-1.264911064067352*fupwindQuad[2]; 
  incr[0] = 0.5*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_x; 
  incr[1] = -0.8660254037844386*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_x; 
  incr[2] = 0.1*(4.47213595499958*alpha[1]*fupwind[2]+5.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_x; 
  incr[3] = -0.1732050807568877*(4.47213595499958*alpha[1]*fupwind[2]+5.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_x; 
  incr[4] = 0.1*(5.0*alpha[0]*fupwind[2]+4.47213595499958*alpha[1]*fupwind[1])*dfac_x; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x1vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[5]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardt[0])*q_)/m_; 

  double alpha[2]; 
  alpha[0] = -(0.7071067811865475*(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x+1.414213562373095*dApardt[0])*q_)/m_; 
  alpha[1] = -(0.7071067811865475*(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x+1.414213562373095*dApardt[1])*q_)/m_; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*(2.23606797749979*alpha[0]*fl[4]+1.732050807568877*(alpha[1]*fl[3]+alpha[0]*fl[2])+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(2.23606797749979*alpha[1]*fl[4]+1.732050807568877*(alpha[0]*fl[3]+alpha[1]*fl[2])+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.3535533905932737*(3.872983346207417*alpha[0]*fl[4]+3.0*(alpha[1]*fl[3]+alpha[0]*fl[2])+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[3] = -0.3535533905932737*(3.872983346207417*alpha[1]*fl[4]+3.0*(alpha[0]*fl[3]+alpha[1]*fl[2])+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[4] = 0.3535533905932737*(5.0*alpha[0]*fl[4]+3.872983346207417*(alpha[1]*fl[3]+alpha[0]*fl[2])+2.23606797749979*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  } else { 
  incr[0] = 0.3535533905932737*(2.23606797749979*alpha[0]*fr[4]-1.732050807568877*(alpha[1]*fr[3]+alpha[0]*fr[2])+alpha[1]*fr[1]+alpha[0]*fr[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(2.23606797749979*alpha[1]*fr[4]-1.732050807568877*(alpha[0]*fr[3]+alpha[1]*fr[2])+alpha[0]*fr[1]+fr[0]*alpha[1])*dfac_v; 
  incr[2] = -0.3535533905932737*(3.872983346207417*alpha[0]*fr[4]-3.0*(alpha[1]*fr[3]+alpha[0]*fr[2])+1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[3] = -0.3535533905932737*(3.872983346207417*alpha[1]*fr[4]-3.0*(alpha[0]*fr[3]+alpha[1]*fr[2])+1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[4] = 0.3535533905932737*(5.0*alpha[0]*fr[4]-3.872983346207417*(alpha[1]*fr[3]+alpha[0]*fr[2])+2.23606797749979*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  }
#elif upwindType == QUAD 
double fupwind[2];
double fupwindQuad[2];
double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]-0.7071067811865475*alpha[1]; 
  fupwindQuad[0] = 0.5*(((-1.118033988749895*fr[4])+1.118033988749895*fl[4]-0.8660254037844386*(fr[3]+fl[3])+0.8660254037844386*(fr[2]+fl[2])+0.5*fr[1]-0.5*(fl[1]+fr[0])+0.5*fl[0])*sgn(alphaQuad)+1.118033988749895*(fr[4]+fl[4])+0.8660254037844386*fr[3]-0.8660254037844386*(fl[3]+fr[2])+0.8660254037844386*fl[2]-0.5*(fr[1]+fl[1])+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*(alpha[1]+alpha[0]); 
  fupwindQuad[1] = 0.5*(((-1.118033988749895*fr[4])+1.118033988749895*fl[4]+0.8660254037844386*(fr[3]+fl[3]+fr[2]+fl[2])-0.5*(fr[1]+fr[0])+0.5*(fl[1]+fl[0]))*sgn(alphaQuad)+1.118033988749895*(fr[4]+fl[4])-0.8660254037844386*(fr[3]+fr[2])+0.8660254037844386*(fl[3]+fl[2])+0.5*(fr[1]+fl[1]+fr[0]+fl[0])); 
  fupwind[0] = 0.7071067811865476*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.7071067811865476*fupwindQuad[1]-0.7071067811865476*fupwindQuad[0]; 
  incr[0] = 0.5*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 
  incr[1] = 0.5*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac_v; 
  incr[2] = -0.8660254037844386*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 
  incr[3] = -0.8660254037844386*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac_v; 
  incr[4] = 1.118033988749895*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_v; 

#endif 
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  return std::abs(alpha0); 
} 
