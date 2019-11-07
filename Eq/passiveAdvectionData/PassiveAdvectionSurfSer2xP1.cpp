#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionSurf2xSer_X1_P1(const double *cflFreqCtrlL, const double *cflFreqCtrlR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  const double *v1 = &fr[4]; 
  const double *v2 = &fr[8]; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.25*(1.732050807568877*v1[1]-1.0*v1[0]); 

  double alpha[2]; 
  alpha[0] = -0.7071067811865475*(1.732050807568877*v1[1]-1.0*v1[0]); 
  alpha[1] = -0.7071067811865475*(1.732050807568877*v1[3]-1.0*v1[2]); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*(alpha[1]*(1.732050807568877*fl[3]+fl[2])+alpha[0]*(1.732050807568877*fl[1]+fl[0]))*dfac1; 
  incr[1] = -0.3535533905932737*(alpha[1]*(3.0*fl[3]+1.732050807568877*fl[2])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac1; 
  incr[2] = 0.3535533905932737*(alpha[0]*(1.732050807568877*fl[3]+fl[2])+alpha[1]*(1.732050807568877*fl[1]+fl[0]))*dfac1; 
  incr[3] = -0.3535533905932737*(alpha[0]*(3.0*fl[3]+1.732050807568877*fl[2])+alpha[1]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac1; 
  } else { 
  incr[0] = -0.3535533905932737*(alpha[1]*(1.732050807568877*fr[3]-1.0*fr[2])+alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac1; 
  incr[1] = 0.3535533905932737*(alpha[1]*(3.0*fr[3]-1.732050807568877*fr[2])+alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac1; 
  incr[2] = -0.3535533905932737*(alpha[0]*(1.732050807568877*fr[3]-1.0*fr[2])+alpha[1]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac1; 
  incr[3] = 0.3535533905932737*(alpha[0]*(3.0*fr[3]-1.732050807568877*fr[2])+alpha[1]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac1; 
  }
#elif upwindType == QUAD 
double fupwind[2];
double fupwindQuad[2];
double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]-1.224744871391589*alpha[1]; 
  fupwindQuad[0] = 0.5*(((-1.5*(fr[3]+fl[3]))+0.8660254037844386*(fr[2]+fr[1])-0.8660254037844386*fl[2]+0.8660254037844386*fl[1]-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)+1.5*fr[3]-1.5*fl[3]-0.8660254037844386*(fr[2]+fl[2]+fr[1])+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 1.224744871391589*alpha[1]+0.7071067811865475*alpha[0]; 
  fupwindQuad[1] = 0.5*((1.5*(fr[3]+fl[3])-0.8660254037844386*fr[2]+0.8660254037844386*(fl[2]+fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-1.5*fr[3]+1.5*fl[3]+0.8660254037844386*(fr[2]+fl[2]+fl[1])-0.8660254037844386*fr[1]+0.5*(fr[0]+fl[0])); 
  fupwind[0] = 0.7071067811865476*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.408248290463863*fupwindQuad[1]-0.408248290463863*fupwindQuad[0]; 
  incr[0] = 0.5*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac1; 
  incr[1] = -0.8660254037844386*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac1; 
  incr[2] = 0.5*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac1; 
  incr[3] = -0.8660254037844386*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac1; 

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
double PassiveAdvectionSurf2xSer_X2_P1(const double *cflFreqCtrlL, const double *cflFreqCtrlR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  const double *v1 = &fr[4]; 
  const double *v2 = &fr[8]; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.25*(1.732050807568877*v2[2]-1.0*v2[0]); 

  double alpha[2]; 
  alpha[0] = -0.7071067811865475*(1.732050807568877*v2[2]-1.0*v2[0]); 
  alpha[1] = -0.7071067811865475*(1.732050807568877*v2[3]-1.0*v2[1]); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*(1.732050807568877*(alpha[1]*fl[3]+alpha[0]*fl[2])+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac2; 
  incr[1] = 0.3535533905932737*(1.732050807568877*(alpha[0]*fl[3]+alpha[1]*fl[2])+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac2; 
  incr[2] = -0.3535533905932737*(3.0*(alpha[1]*fl[3]+alpha[0]*fl[2])+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac2; 
  incr[3] = -0.3535533905932737*(3.0*(alpha[0]*fl[3]+alpha[1]*fl[2])+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac2; 
  } else { 
  incr[0] = -0.3535533905932737*(1.732050807568877*(alpha[1]*fr[3]+alpha[0]*fr[2])-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac2; 
  incr[1] = -0.3535533905932737*(1.732050807568877*(alpha[0]*fr[3]+alpha[1]*fr[2])-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac2; 
  incr[2] = 0.3535533905932737*(3.0*(alpha[1]*fr[3]+alpha[0]*fr[2])-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac2; 
  incr[3] = 0.3535533905932737*(3.0*(alpha[0]*fr[3]+alpha[1]*fr[2])-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac2; 
  }
#elif upwindType == QUAD 
double fupwind[2];
double fupwindQuad[2];
double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]-1.224744871391589*alpha[1]; 
  fupwindQuad[0] = 0.5*(((-1.5*(fr[3]+fl[3]))+0.8660254037844386*(fr[2]+fl[2]+fr[1])-0.8660254037844386*fl[1]-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)+1.5*fr[3]-1.5*fl[3]-0.8660254037844386*(fr[2]+fr[1])+0.8660254037844386*fl[2]-0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 1.224744871391589*alpha[1]+0.7071067811865475*alpha[0]; 
  fupwindQuad[1] = 0.5*((1.5*(fr[3]+fl[3])+0.8660254037844386*(fr[2]+fl[2]+fl[1])-0.8660254037844386*fr[1]-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-1.5*fr[3]+1.5*fl[3]-0.8660254037844386*fr[2]+0.8660254037844386*(fl[2]+fr[1]+fl[1])+0.5*(fr[0]+fl[0])); 
  fupwind[0] = 0.7071067811865476*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.408248290463863*fupwindQuad[1]-0.408248290463863*fupwindQuad[0]; 
  incr[0] = 0.5*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac2; 
  incr[1] = 0.5*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac2; 
  incr[2] = -0.8660254037844386*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac2; 
  incr[3] = -0.8660254037844386*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac2; 

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
