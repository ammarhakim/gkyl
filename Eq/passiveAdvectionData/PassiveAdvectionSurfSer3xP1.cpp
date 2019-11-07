#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionSurf3xSer_X1_P1(const double *cflFreqCtrlL, const double *cflFreqCtrlR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  double dfac3 = 2.0/dxv[2]; 
  double w3 = w[2]; 
  const double *v1 = &fr[8]; 
  const double *v2 = &fr[16]; 
  const double *v3 = &fr[24]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.1767766952966368*(1.732050807568877*v1[1]-1.0*v1[0]); 

  double alpha[4]; 
  alpha[0] = -0.5*(2.449489742783178*v1[1]-1.414213562373095*v1[0]); 
  alpha[1] = -0.5*(2.449489742783178*v1[4]-1.414213562373095*v1[2]); 
  alpha[2] = -0.5*(2.449489742783178*v1[5]-1.414213562373095*v1[3]); 
  alpha[3] = -0.5*(2.449489742783178*v1[7]-1.414213562373095*v1[6]); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*(alpha[3]*(1.732050807568877*fl[7]+fl[6])+1.732050807568877*(alpha[2]*fl[5]+alpha[1]*fl[4])+alpha[2]*fl[3]+alpha[1]*fl[2]+alpha[0]*(1.732050807568877*fl[1]+fl[0]))*dfac1; 
  incr[1] = -0.25*(alpha[3]*(3.0*fl[7]+1.732050807568877*fl[6])+3.0*(alpha[2]*fl[5]+alpha[1]*fl[4])+1.732050807568877*(alpha[2]*fl[3]+alpha[1]*fl[2])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac1; 
  incr[2] = 0.25*(alpha[2]*(1.732050807568877*fl[7]+fl[6])+1.732050807568877*(alpha[3]*fl[5]+alpha[0]*fl[4])+alpha[3]*fl[3]+alpha[0]*fl[2]+alpha[1]*(1.732050807568877*fl[1]+fl[0]))*dfac1; 
  incr[3] = 0.25*(alpha[1]*(1.732050807568877*fl[7]+fl[6])+1.732050807568877*(alpha[0]*fl[5]+alpha[3]*fl[4])+alpha[0]*fl[3]+fl[2]*alpha[3]+(1.732050807568877*fl[1]+fl[0])*alpha[2])*dfac1; 
  incr[4] = -0.25*(alpha[2]*(3.0*fl[7]+1.732050807568877*fl[6])+3.0*(alpha[3]*fl[5]+alpha[0]*fl[4])+1.732050807568877*(alpha[3]*fl[3]+alpha[0]*fl[2])+alpha[1]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac1; 
  incr[5] = -0.25*(alpha[1]*(3.0*fl[7]+1.732050807568877*fl[6])+3.0*(alpha[0]*fl[5]+alpha[3]*fl[4])+1.732050807568877*(alpha[0]*fl[3]+fl[2]*alpha[3])+(3.0*fl[1]+1.732050807568877*fl[0])*alpha[2])*dfac1; 
  incr[6] = 0.25*(alpha[0]*(1.732050807568877*fl[7]+fl[6])+1.732050807568877*(alpha[1]*fl[5]+alpha[2]*fl[4])+alpha[1]*fl[3]+(1.732050807568877*fl[1]+fl[0])*alpha[3]+alpha[2]*fl[2])*dfac1; 
  incr[7] = -0.25*(alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])+3.0*(alpha[1]*fl[5]+alpha[2]*fl[4])+1.732050807568877*alpha[1]*fl[3]+3.0*fl[1]*alpha[3]+1.732050807568877*(fl[0]*alpha[3]+alpha[2]*fl[2]))*dfac1; 
  } else { 
  incr[0] = -0.25*(alpha[3]*(1.732050807568877*fr[7]-1.0*fr[6])+1.732050807568877*(alpha[2]*fr[5]+alpha[1]*fr[4])-1.0*(alpha[2]*fr[3]+alpha[1]*fr[2])+alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac1; 
  incr[1] = 0.25*(alpha[3]*(3.0*fr[7]-1.732050807568877*fr[6])+3.0*(alpha[2]*fr[5]+alpha[1]*fr[4])-1.732050807568877*(alpha[2]*fr[3]+alpha[1]*fr[2])+alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac1; 
  incr[2] = -0.25*(alpha[2]*(1.732050807568877*fr[7]-1.0*fr[6])+1.732050807568877*(alpha[3]*fr[5]+alpha[0]*fr[4])-1.0*(alpha[3]*fr[3]+alpha[0]*fr[2])+alpha[1]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac1; 
  incr[3] = -0.25*(alpha[1]*(1.732050807568877*fr[7]-1.0*fr[6])+1.732050807568877*(alpha[0]*fr[5]+alpha[3]*fr[4])-1.0*(alpha[0]*fr[3]+fr[2]*alpha[3])+(1.732050807568877*fr[1]-1.0*fr[0])*alpha[2])*dfac1; 
  incr[4] = 0.25*(alpha[2]*(3.0*fr[7]-1.732050807568877*fr[6])+3.0*(alpha[3]*fr[5]+alpha[0]*fr[4])-1.732050807568877*(alpha[3]*fr[3]+alpha[0]*fr[2])+alpha[1]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac1; 
  incr[5] = 0.25*(alpha[1]*(3.0*fr[7]-1.732050807568877*fr[6])+3.0*(alpha[0]*fr[5]+alpha[3]*fr[4])-1.732050807568877*(alpha[0]*fr[3]+fr[2]*alpha[3])+(3.0*fr[1]-1.732050807568877*fr[0])*alpha[2])*dfac1; 
  incr[6] = -0.25*(alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])+1.732050807568877*(alpha[1]*fr[5]+alpha[2]*fr[4])-1.0*alpha[1]*fr[3]+1.732050807568877*fr[1]*alpha[3]-1.0*(fr[0]*alpha[3]+alpha[2]*fr[2]))*dfac1; 
  incr[7] = 0.25*(alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])+3.0*(alpha[1]*fr[5]+alpha[2]*fr[4])-1.732050807568877*alpha[1]*fr[3]+3.0*fr[1]*alpha[3]-1.732050807568877*(fr[0]*alpha[3]+alpha[2]*fr[2]))*dfac1; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double alphaQuad;
  alphaQuad = 1.5*alpha[3]-0.8660254037844386*(alpha[2]+alpha[1])+0.5*alpha[0]; 
  fupwindQuad[0] = 0.5*((1.837117307087383*(fr[7]+fl[7])-1.060660171779821*(fr[6]+fr[5]+fr[4])+1.060660171779821*fl[6]-1.060660171779821*(fl[5]+fl[4])+0.6123724356957944*(fr[3]+fr[2]+fr[1])-0.6123724356957944*(fl[3]+fl[2])+0.6123724356957944*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)-1.837117307087383*fr[7]+1.837117307087383*fl[7]+1.060660171779821*(fr[6]+fl[6]+fr[5]+fr[4])-1.060660171779821*(fl[5]+fl[4])-0.6123724356957944*(fr[3]+fl[3]+fr[2]+fl[2]+fr[1])+0.6123724356957944*fl[1]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = (-1.5*alpha[3])-0.8660254037844386*alpha[2]+0.8660254037844386*alpha[1]+0.5*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-1.837117307087383*(fr[7]+fl[7]))+1.060660171779821*fr[6]-1.060660171779821*(fl[6]+fr[5]+fl[5])+1.060660171779821*(fr[4]+fl[4])+0.6123724356957944*fr[3]-0.6123724356957944*(fl[3]+fr[2])+0.6123724356957944*(fl[2]+fr[1]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)+1.837117307087383*fr[7]-1.837117307087383*fl[7]-1.060660171779821*(fr[6]+fl[6]+fl[5])+1.060660171779821*fr[5]-1.060660171779821*fr[4]+1.060660171779821*fl[4]-0.6123724356957944*(fr[3]+fl[3])+0.6123724356957944*(fr[2]+fl[2]+fl[1])-0.6123724356957944*fr[1]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = (-1.5*alpha[3])+0.8660254037844386*alpha[2]-0.8660254037844386*alpha[1]+0.5*alpha[0]; 
  fupwindQuad[2] = 0.5*(((-1.837117307087383*(fr[7]+fl[7]))+1.060660171779821*(fr[6]+fr[5])-1.060660171779821*fl[6]+1.060660171779821*fl[5]-1.060660171779821*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*(fl[3]+fr[2]+fr[1])-0.6123724356957944*fl[2]+0.6123724356957944*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)+1.837117307087383*fr[7]-1.837117307087383*fl[7]-1.060660171779821*(fr[6]+fl[6]+fr[5])+1.060660171779821*(fl[5]+fr[4])-1.060660171779821*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.6123724356957944*(fr[2]+fl[2]+fr[1])+0.6123724356957944*fl[1]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 1.5*alpha[3]+0.8660254037844386*(alpha[2]+alpha[1])+0.5*alpha[0]; 
  fupwindQuad[3] = 0.5*((1.837117307087383*(fr[7]+fl[7])-1.060660171779821*fr[6]+1.060660171779821*(fl[6]+fr[5]+fl[5]+fr[4]+fl[4])-0.6123724356957944*(fr[3]+fr[2])+0.6123724356957944*(fl[3]+fl[2]+fr[1]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)-1.837117307087383*fr[7]+1.837117307087383*fl[7]+1.060660171779821*(fr[6]+fl[6]+fl[5]+fl[4])-1.060660171779821*(fr[5]+fr[4])+0.6123724356957944*(fr[3]+fl[3]+fr[2]+fl[2]+fl[1])-0.6123724356957944*fr[1]+0.3535533905932737*(fr[0]+fl[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.2886751345948129*fupwindQuad[3]-0.2886751345948129*fupwindQuad[2]+0.2886751345948129*fupwindQuad[1]-0.2886751345948129*fupwindQuad[0]; 
  fupwind[2] = 0.2886751345948129*(fupwindQuad[3]+fupwindQuad[2])-0.2886751345948129*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.1666666666666666*fupwindQuad[3]-0.1666666666666666*(fupwindQuad[2]+fupwindQuad[1])+0.1666666666666666*fupwindQuad[0]; 
  incr[0] = 0.3535533905932737*(alpha[3]*fupwind[3]+alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac1; 
  incr[1] = -0.6123724356957944*(alpha[3]*fupwind[3]+alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac1; 
  incr[2] = 0.3535533905932737*(alpha[2]*fupwind[3]+fupwind[2]*alpha[3]+alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac1; 
  incr[3] = 0.3535533905932737*(alpha[1]*fupwind[3]+fupwind[1]*alpha[3]+alpha[0]*fupwind[2]+fupwind[0]*alpha[2])*dfac1; 
  incr[4] = -0.6123724356957944*(alpha[2]*fupwind[3]+fupwind[2]*alpha[3]+alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac1; 
  incr[5] = -0.6123724356957944*(alpha[1]*fupwind[3]+fupwind[1]*alpha[3]+alpha[0]*fupwind[2]+fupwind[0]*alpha[2])*dfac1; 
  incr[6] = 0.3535533905932737*(alpha[0]*fupwind[3]+fupwind[0]*alpha[3]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2])*dfac1; 
  incr[7] = -0.6123724356957944*(alpha[0]*fupwind[3]+fupwind[0]*alpha[3]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2])*dfac1; 

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
double PassiveAdvectionSurf3xSer_X2_P1(const double *cflFreqCtrlL, const double *cflFreqCtrlR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  double dfac3 = 2.0/dxv[2]; 
  double w3 = w[2]; 
  const double *v1 = &fr[8]; 
  const double *v2 = &fr[16]; 
  const double *v3 = &fr[24]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.1767766952966368*(1.732050807568877*v2[2]-1.0*v2[0]); 

  double alpha[4]; 
  alpha[0] = -0.5*(2.449489742783178*v2[2]-1.414213562373095*v2[0]); 
  alpha[1] = -0.5*(2.449489742783178*v2[4]-1.414213562373095*v2[1]); 
  alpha[2] = -0.5*(2.449489742783178*v2[6]-1.414213562373095*v2[3]); 
  alpha[3] = -0.5*(2.449489742783178*v2[7]-1.414213562373095*v2[5]); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*(alpha[3]*fl[7]+alpha[2]*fl[6])+alpha[3]*fl[5]+1.732050807568877*alpha[1]*fl[4]+alpha[2]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac2; 
  incr[1] = 0.25*(1.732050807568877*(alpha[2]*fl[7]+alpha[3]*fl[6])+alpha[2]*fl[5]+1.732050807568877*alpha[0]*fl[4]+alpha[3]*fl[3]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac2; 
  incr[2] = -0.25*(3.0*(alpha[3]*fl[7]+alpha[2]*fl[6])+1.732050807568877*alpha[3]*fl[5]+3.0*alpha[1]*fl[4]+1.732050807568877*alpha[2]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac2; 
  incr[3] = 0.25*(1.732050807568877*(alpha[1]*fl[7]+alpha[0]*fl[6])+alpha[1]*fl[5]+1.732050807568877*alpha[3]*fl[4]+alpha[0]*fl[3]+fl[1]*alpha[3]+alpha[2]*(1.732050807568877*fl[2]+fl[0]))*dfac2; 
  incr[4] = -0.25*(3.0*(alpha[2]*fl[7]+alpha[3]*fl[6])+1.732050807568877*alpha[2]*fl[5]+3.0*alpha[0]*fl[4]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[1]*fl[2]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac2; 
  incr[5] = 0.25*(1.732050807568877*(alpha[0]*fl[7]+alpha[1]*fl[6])+alpha[0]*fl[5]+1.732050807568877*alpha[2]*fl[4]+alpha[1]*fl[3]+(1.732050807568877*fl[2]+fl[0])*alpha[3]+fl[1]*alpha[2])*dfac2; 
  incr[6] = -0.25*(3.0*(alpha[1]*fl[7]+alpha[0]*fl[6])+1.732050807568877*alpha[1]*fl[5]+3.0*alpha[3]*fl[4]+1.732050807568877*(alpha[0]*fl[3]+fl[1]*alpha[3])+alpha[2]*(3.0*fl[2]+1.732050807568877*fl[0]))*dfac2; 
  incr[7] = -0.25*(3.0*(alpha[0]*fl[7]+alpha[1]*fl[6])+1.732050807568877*alpha[0]*fl[5]+3.0*alpha[2]*fl[4]+1.732050807568877*alpha[1]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*(fl[0]*alpha[3]+fl[1]*alpha[2]))*dfac2; 
  } else { 
  incr[0] = -0.25*(1.732050807568877*(alpha[3]*fr[7]+alpha[2]*fr[6])-1.0*alpha[3]*fr[5]+1.732050807568877*alpha[1]*fr[4]-1.0*alpha[2]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac2; 
  incr[1] = -0.25*(1.732050807568877*(alpha[2]*fr[7]+alpha[3]*fr[6])-1.0*alpha[2]*fr[5]+1.732050807568877*alpha[0]*fr[4]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[1]*fr[2]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac2; 
  incr[2] = 0.25*(3.0*(alpha[3]*fr[7]+alpha[2]*fr[6])-1.732050807568877*alpha[3]*fr[5]+3.0*alpha[1]*fr[4]-1.732050807568877*alpha[2]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac2; 
  incr[3] = -0.25*(1.732050807568877*(alpha[1]*fr[7]+alpha[0]*fr[6])-1.0*alpha[1]*fr[5]+1.732050807568877*alpha[3]*fr[4]-1.0*(alpha[0]*fr[3]+fr[1]*alpha[3])+alpha[2]*(1.732050807568877*fr[2]-1.0*fr[0]))*dfac2; 
  incr[4] = 0.25*(3.0*(alpha[2]*fr[7]+alpha[3]*fr[6])-1.732050807568877*alpha[2]*fr[5]+3.0*alpha[0]*fr[4]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[1]*fr[2]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac2; 
  incr[5] = -0.25*(1.732050807568877*(alpha[0]*fr[7]+alpha[1]*fr[6])-1.0*alpha[0]*fr[5]+1.732050807568877*alpha[2]*fr[4]-1.0*alpha[1]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*(fr[0]*alpha[3]+fr[1]*alpha[2]))*dfac2; 
  incr[6] = 0.25*(3.0*(alpha[1]*fr[7]+alpha[0]*fr[6])-1.732050807568877*alpha[1]*fr[5]+3.0*alpha[3]*fr[4]-1.732050807568877*(alpha[0]*fr[3]+fr[1]*alpha[3])+alpha[2]*(3.0*fr[2]-1.732050807568877*fr[0]))*dfac2; 
  incr[7] = 0.25*(3.0*(alpha[0]*fr[7]+alpha[1]*fr[6])-1.732050807568877*alpha[0]*fr[5]+3.0*alpha[2]*fr[4]-1.732050807568877*alpha[1]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*(fr[0]*alpha[3]+fr[1]*alpha[2]))*dfac2; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double alphaQuad;
  alphaQuad = 1.5*alpha[3]-0.8660254037844386*(alpha[2]+alpha[1])+0.5*alpha[0]; 
  fupwindQuad[0] = 0.5*((1.837117307087383*(fr[7]+fl[7])-1.060660171779821*(fr[6]+fl[6]+fr[5]+fr[4])+1.060660171779821*fl[5]-1.060660171779821*fl[4]+0.6123724356957944*(fr[3]+fr[2]+fr[1])-0.6123724356957944*fl[3]+0.6123724356957944*fl[2]-0.6123724356957944*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)-1.837117307087383*fr[7]+1.837117307087383*fl[7]+1.060660171779821*(fr[6]+fr[5]+fr[4])-1.060660171779821*fl[6]+1.060660171779821*fl[5]-1.060660171779821*fl[4]-0.6123724356957944*(fr[3]+fl[3]+fr[2]+fr[1])+0.6123724356957944*fl[2]-0.6123724356957944*fl[1]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = (-1.5*alpha[3])-0.8660254037844386*alpha[2]+0.8660254037844386*alpha[1]+0.5*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-1.837117307087383*(fr[7]+fl[7]))-1.060660171779821*(fr[6]+fl[6]+fl[5])+1.060660171779821*(fr[5]+fr[4]+fl[4])+0.6123724356957944*(fr[3]+fr[2])-0.6123724356957944*fl[3]+0.6123724356957944*(fl[2]+fl[1])-0.6123724356957944*fr[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)+1.837117307087383*fr[7]-1.837117307087383*fl[7]+1.060660171779821*fr[6]-1.060660171779821*(fl[6]+fr[5]+fl[5]+fr[4])+1.060660171779821*fl[4]-0.6123724356957944*(fr[3]+fl[3]+fr[2])+0.6123724356957944*(fl[2]+fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = (-1.5*alpha[3])+0.8660254037844386*alpha[2]-0.8660254037844386*alpha[1]+0.5*alpha[0]; 
  fupwindQuad[2] = 0.5*(((-1.837117307087383*(fr[7]+fl[7]))+1.060660171779821*(fr[6]+fl[6]+fr[5])-1.060660171779821*(fl[5]+fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*(fl[3]+fr[2]+fl[2]+fr[1])-0.6123724356957944*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)+1.837117307087383*fr[7]-1.837117307087383*fl[7]-1.060660171779821*(fr[6]+fr[5])+1.060660171779821*fl[6]-1.060660171779821*(fl[5]+fl[4])+1.060660171779821*fr[4]+0.6123724356957944*(fr[3]+fl[3]+fl[2])-0.6123724356957944*(fr[2]+fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 1.5*alpha[3]+0.8660254037844386*(alpha[2]+alpha[1])+0.5*alpha[0]; 
  fupwindQuad[3] = 0.5*((1.837117307087383*(fr[7]+fl[7])+1.060660171779821*(fr[6]+fl[6]+fl[5]+fl[4])-1.060660171779821*fr[5]+1.060660171779821*fr[4]-0.6123724356957944*fr[3]+0.6123724356957944*(fl[3]+fr[2]+fl[2]+fl[1])-0.6123724356957944*fr[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)-1.837117307087383*fr[7]+1.837117307087383*fl[7]-1.060660171779821*fr[6]+1.060660171779821*(fl[6]+fr[5]+fl[5]+fl[4])-1.060660171779821*fr[4]+0.6123724356957944*(fr[3]+fl[3]+fl[2]+fl[1])-0.6123724356957944*fr[2]+0.6123724356957944*fr[1]+0.3535533905932737*(fr[0]+fl[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.2886751345948129*fupwindQuad[3]-0.2886751345948129*fupwindQuad[2]+0.2886751345948129*fupwindQuad[1]-0.2886751345948129*fupwindQuad[0]; 
  fupwind[2] = 0.2886751345948129*(fupwindQuad[3]+fupwindQuad[2])-0.2886751345948129*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.1666666666666666*fupwindQuad[3]-0.1666666666666666*(fupwindQuad[2]+fupwindQuad[1])+0.1666666666666666*fupwindQuad[0]; 
  incr[0] = 0.3535533905932737*(alpha[3]*fupwind[3]+alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac2; 
  incr[1] = 0.3535533905932737*(alpha[2]*fupwind[3]+fupwind[2]*alpha[3]+alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac2; 
  incr[2] = -0.6123724356957944*(alpha[3]*fupwind[3]+alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac2; 
  incr[3] = 0.3535533905932737*(alpha[1]*fupwind[3]+fupwind[1]*alpha[3]+alpha[0]*fupwind[2]+fupwind[0]*alpha[2])*dfac2; 
  incr[4] = -0.6123724356957944*(alpha[2]*fupwind[3]+fupwind[2]*alpha[3]+alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac2; 
  incr[5] = 0.3535533905932737*(alpha[0]*fupwind[3]+fupwind[0]*alpha[3]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2])*dfac2; 
  incr[6] = -0.6123724356957944*(alpha[1]*fupwind[3]+fupwind[1]*alpha[3]+alpha[0]*fupwind[2]+fupwind[0]*alpha[2])*dfac2; 
  incr[7] = -0.6123724356957944*(alpha[0]*fupwind[3]+fupwind[0]*alpha[3]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2])*dfac2; 

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
double PassiveAdvectionSurf3xSer_X3_P1(const double *cflFreqCtrlL, const double *cflFreqCtrlR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  double dfac3 = 2.0/dxv[2]; 
  double w3 = w[2]; 
  const double *v1 = &fr[8]; 
  const double *v2 = &fr[16]; 
  const double *v3 = &fr[24]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.1767766952966368*(1.732050807568877*v3[3]-1.0*v3[0]); 

  double alpha[4]; 
  alpha[0] = -0.5*(2.449489742783178*v3[3]-1.414213562373095*v3[0]); 
  alpha[1] = -0.5*(2.449489742783178*v3[5]-1.414213562373095*v3[1]); 
  alpha[2] = -0.5*(2.449489742783178*v3[6]-1.414213562373095*v3[2]); 
  alpha[3] = -0.5*(2.449489742783178*v3[7]-1.414213562373095*v3[4]); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*(alpha[3]*fl[7]+alpha[2]*fl[6]+alpha[1]*fl[5])+alpha[3]*fl[4]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac3; 
  incr[1] = 0.25*(1.732050807568877*(alpha[2]*fl[7]+alpha[3]*fl[6]+alpha[0]*fl[5])+alpha[2]*fl[4]+1.732050807568877*alpha[1]*fl[3]+fl[2]*alpha[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac3; 
  incr[2] = 0.25*(1.732050807568877*(alpha[1]*fl[7]+alpha[0]*fl[6]+alpha[3]*fl[5])+alpha[1]*fl[4]+1.732050807568877*alpha[2]*fl[3]+fl[1]*alpha[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac3; 
  incr[3] = -0.25*(3.0*(alpha[3]*fl[7]+alpha[2]*fl[6]+alpha[1]*fl[5])+1.732050807568877*alpha[3]*fl[4]+3.0*alpha[0]*fl[3]+1.732050807568877*(alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac3; 
  incr[4] = 0.25*(1.732050807568877*(alpha[0]*fl[7]+alpha[1]*fl[6]+alpha[2]*fl[5])+alpha[0]*fl[4]+alpha[3]*(1.732050807568877*fl[3]+fl[0])+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac3; 
  incr[5] = -0.25*(3.0*(alpha[2]*fl[7]+alpha[3]*fl[6]+alpha[0]*fl[5])+1.732050807568877*alpha[2]*fl[4]+3.0*alpha[1]*fl[3]+1.732050807568877*(fl[2]*alpha[3]+alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac3; 
  incr[6] = -0.25*(3.0*(alpha[1]*fl[7]+alpha[0]*fl[6]+alpha[3]*fl[5])+1.732050807568877*alpha[1]*fl[4]+3.0*alpha[2]*fl[3]+1.732050807568877*(fl[1]*alpha[3]+alpha[0]*fl[2]+fl[0]*alpha[2]))*dfac3; 
  incr[7] = -0.25*(3.0*(alpha[0]*fl[7]+alpha[1]*fl[6]+alpha[2]*fl[5])+1.732050807568877*alpha[0]*fl[4]+3.0*alpha[3]*fl[3]+1.732050807568877*(fl[0]*alpha[3]+alpha[1]*fl[2]+fl[1]*alpha[2]))*dfac3; 
  } else { 
  incr[0] = -0.25*(1.732050807568877*(alpha[3]*fr[7]+alpha[2]*fr[6]+alpha[1]*fr[5])-1.0*alpha[3]*fr[4]+1.732050807568877*alpha[0]*fr[3]-1.0*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac3; 
  incr[1] = -0.25*(1.732050807568877*(alpha[2]*fr[7]+alpha[3]*fr[6]+alpha[0]*fr[5])-1.0*alpha[2]*fr[4]+1.732050807568877*alpha[1]*fr[3]-1.0*(fr[2]*alpha[3]+alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac3; 
  incr[2] = -0.25*(1.732050807568877*(alpha[1]*fr[7]+alpha[0]*fr[6]+alpha[3]*fr[5])-1.0*alpha[1]*fr[4]+1.732050807568877*alpha[2]*fr[3]-1.0*(fr[1]*alpha[3]+alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac3; 
  incr[3] = 0.25*(3.0*(alpha[3]*fr[7]+alpha[2]*fr[6]+alpha[1]*fr[5])-1.732050807568877*alpha[3]*fr[4]+3.0*alpha[0]*fr[3]-1.732050807568877*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac3; 
  incr[4] = -0.25*(1.732050807568877*(alpha[0]*fr[7]+alpha[1]*fr[6]+alpha[2]*fr[5])-1.0*alpha[0]*fr[4]+1.732050807568877*alpha[3]*fr[3]-1.0*(fr[0]*alpha[3]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac3; 
  incr[5] = 0.25*(3.0*(alpha[2]*fr[7]+alpha[3]*fr[6]+alpha[0]*fr[5])-1.732050807568877*alpha[2]*fr[4]+3.0*alpha[1]*fr[3]-1.732050807568877*(fr[2]*alpha[3]+alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac3; 
  incr[6] = 0.25*(3.0*(alpha[1]*fr[7]+alpha[0]*fr[6]+alpha[3]*fr[5])-1.732050807568877*alpha[1]*fr[4]+3.0*alpha[2]*fr[3]-1.732050807568877*(fr[1]*alpha[3]+alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac3; 
  incr[7] = 0.25*(3.0*(alpha[0]*fr[7]+alpha[1]*fr[6]+alpha[2]*fr[5])-1.732050807568877*alpha[0]*fr[4]+3.0*alpha[3]*fr[3]-1.732050807568877*(fr[0]*alpha[3]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac3; 
  }
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
double alphaQuad;
  alphaQuad = 1.5*alpha[3]-0.8660254037844386*(alpha[2]+alpha[1])+0.5*alpha[0]; 
  fupwindQuad[0] = 0.5*((1.837117307087383*(fr[7]+fl[7])-1.060660171779821*(fr[6]+fl[6]+fr[5]+fl[5]+fr[4])+1.060660171779821*fl[4]+0.6123724356957944*(fr[3]+fl[3]+fr[2]+fr[1])-0.6123724356957944*(fl[2]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)-1.837117307087383*fr[7]+1.837117307087383*fl[7]+1.060660171779821*(fr[6]+fr[5]+fr[4])-1.060660171779821*(fl[6]+fl[5])+1.060660171779821*fl[4]-0.6123724356957944*(fr[3]+fr[2]+fr[1])+0.6123724356957944*fl[3]-0.6123724356957944*(fl[2]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = (-1.5*alpha[3])-0.8660254037844386*alpha[2]+0.8660254037844386*alpha[1]+0.5*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-1.837117307087383*(fr[7]+fl[7]))-1.060660171779821*(fr[6]+fl[6])+1.060660171779821*(fr[5]+fl[5]+fr[4])-1.060660171779821*fl[4]+0.6123724356957944*(fr[3]+fl[3]+fr[2])-0.6123724356957944*(fl[2]+fr[1])+0.6123724356957944*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)+1.837117307087383*fr[7]-1.837117307087383*fl[7]+1.060660171779821*fr[6]-1.060660171779821*(fl[6]+fr[5]+fr[4])+1.060660171779821*fl[5]-1.060660171779821*fl[4]-0.6123724356957944*(fr[3]+fr[2])+0.6123724356957944*fl[3]-0.6123724356957944*fl[2]+0.6123724356957944*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = (-1.5*alpha[3])+0.8660254037844386*alpha[2]-0.8660254037844386*alpha[1]+0.5*alpha[0]; 
  fupwindQuad[2] = 0.5*(((-1.837117307087383*(fr[7]+fl[7]))+1.060660171779821*(fr[6]+fl[6])-1.060660171779821*(fr[5]+fl[5]+fl[4])+1.060660171779821*fr[4]+0.6123724356957944*(fr[3]+fl[3]+fl[2])-0.6123724356957944*fr[2]+0.6123724356957944*fr[1]-0.6123724356957944*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)+1.837117307087383*fr[7]-1.837117307087383*fl[7]-1.060660171779821*fr[6]+1.060660171779821*(fl[6]+fr[5])-1.060660171779821*(fl[5]+fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*(fl[3]+fr[2]+fl[2])-0.6123724356957944*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaQuad = 1.5*alpha[3]+0.8660254037844386*(alpha[2]+alpha[1])+0.5*alpha[0]; 
  fupwindQuad[3] = 0.5*((1.837117307087383*(fr[7]+fl[7])+1.060660171779821*(fr[6]+fl[6]+fr[5]+fl[5]+fl[4])-1.060660171779821*fr[4]+0.6123724356957944*(fr[3]+fl[3]+fl[2]+fl[1])-0.6123724356957944*(fr[2]+fr[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaQuad)-1.837117307087383*fr[7]+1.837117307087383*fl[7]-1.060660171779821*(fr[6]+fr[5])+1.060660171779821*(fl[6]+fl[5]+fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*(fl[3]+fr[2]+fl[2]+fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.2886751345948129*fupwindQuad[3]-0.2886751345948129*fupwindQuad[2]+0.2886751345948129*fupwindQuad[1]-0.2886751345948129*fupwindQuad[0]; 
  fupwind[2] = 0.2886751345948129*(fupwindQuad[3]+fupwindQuad[2])-0.2886751345948129*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[3] = 0.1666666666666666*fupwindQuad[3]-0.1666666666666666*(fupwindQuad[2]+fupwindQuad[1])+0.1666666666666666*fupwindQuad[0]; 
  incr[0] = 0.3535533905932737*(alpha[3]*fupwind[3]+alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac3; 
  incr[1] = 0.3535533905932737*(alpha[2]*fupwind[3]+fupwind[2]*alpha[3]+alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac3; 
  incr[2] = 0.3535533905932737*(alpha[1]*fupwind[3]+fupwind[1]*alpha[3]+alpha[0]*fupwind[2]+fupwind[0]*alpha[2])*dfac3; 
  incr[3] = -0.6123724356957944*(alpha[3]*fupwind[3]+alpha[2]*fupwind[2]+alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac3; 
  incr[4] = 0.3535533905932737*(alpha[0]*fupwind[3]+fupwind[0]*alpha[3]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2])*dfac3; 
  incr[5] = -0.6123724356957944*(alpha[2]*fupwind[3]+fupwind[2]*alpha[3]+alpha[0]*fupwind[1]+fupwind[0]*alpha[1])*dfac3; 
  incr[6] = -0.6123724356957944*(alpha[1]*fupwind[3]+fupwind[1]*alpha[3]+alpha[0]*fupwind[2]+fupwind[0]*alpha[2])*dfac3; 
  incr[7] = -0.6123724356957944*(alpha[0]*fupwind[3]+fupwind[0]*alpha[3]+alpha[1]*fupwind[2]+fupwind[1]*alpha[2])*dfac3; 

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
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
