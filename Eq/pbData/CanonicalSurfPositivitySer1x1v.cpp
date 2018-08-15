#include <CanonicalModDecl.h> 
#include <cmath> 
double CanonicalSurfPositivity1x1vSer_X_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[0]; 
  double incr[4]; 

  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.0*(3.0*H[3]-1.732050807568877*H[2]))/dxv1; 

  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[2] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv0); 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -(1.0*(5.196152422706631*fhatALVal[1]*H[3]+3.0*fhatALVal[0]*H[3]-3.0*fhatALVal[1]*H[2]-1.732050807568877*fhatALVal[0]*H[2]))/(dxv0*dxv1); 
  incr[1] = (3.0*(3.0*fhatALVal[1]*H[3]+1.732050807568877*fhatALVal[0]*H[3]-1.732050807568877*fhatALVal[1]*H[2]-1.0*fhatALVal[0]*H[2]))/(dxv0*dxv1); 
  incr[2] = -(1.0*(5.196152422706631*H[3]*fhatALVal[3]-3.0*H[2]*fhatALVal[3]+3.0*fhatALVal[2]*H[3]-1.732050807568877*H[2]*fhatALVal[2]))/(dxv0*dxv1); 
  incr[3] = (3.0*(3.0*H[3]*fhatALVal[3]-1.732050807568877*H[2]*fhatALVal[3]+1.732050807568877*fhatALVal[2]*H[3]-1.0*H[2]*fhatALVal[2]))/(dxv0*dxv1); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[2] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv0); 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = (5.196152422706631*fhatALVal[1]*H[3]-3.0*fhatALVal[0]*H[3]-3.0*fhatALVal[1]*H[2]+1.732050807568877*fhatALVal[0]*H[2])/(dxv0*dxv1); 
  incr[1] = -(3.0*(3.0*fhatALVal[1]*H[3]-1.732050807568877*fhatALVal[0]*H[3]-1.732050807568877*fhatALVal[1]*H[2]+fhatALVal[0]*H[2]))/(dxv0*dxv1); 
  incr[2] = (5.196152422706631*H[3]*fhatALVal[3]-3.0*H[2]*fhatALVal[3]-3.0*fhatALVal[2]*H[3]+1.732050807568877*H[2]*fhatALVal[2])/(dxv0*dxv1); 
  incr[3] = -(3.0*(3.0*H[3]*fhatALVal[3]-1.732050807568877*H[2]*fhatALVal[3]-1.732050807568877*fhatALVal[2]*H[3]+H[2]*fhatALVal[2]))/(dxv0*dxv1); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double CanonicalSurfPositivity1x1vSer_VX_P1(const double dt, const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[1]; 
  double incr[4]; 

  // surface-averaged phase velocity in this direction 
  double alpha0 = (3.0*H[3]-1.732050807568877*H[1])/dxv0; 

  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv1); 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = (5.196152422706631*fhatALVal[2]*H[3]+3.0*fhatALVal[0]*H[3]-3.0*H[1]*fhatALVal[2]-1.732050807568877*fhatALVal[0]*H[1])/(dxv0*dxv1); 
  incr[1] = (5.196152422706631*H[3]*fhatALVal[3]-3.0*H[1]*fhatALVal[3]+3.0*fhatALVal[1]*H[3]-1.732050807568877*H[1]*fhatALVal[1])/(dxv0*dxv1); 
  incr[2] = -(3.0*(3.0*fhatALVal[2]*H[3]+1.732050807568877*fhatALVal[0]*H[3]-1.732050807568877*H[1]*fhatALVal[2]-1.0*fhatALVal[0]*H[1]))/(dxv0*dxv1); 
  incr[3] = -(3.0*(3.0*H[3]*fhatALVal[3]-1.732050807568877*H[1]*fhatALVal[3]+1.732050807568877*fhatALVal[1]*H[3]-1.0*H[1]*fhatALVal[1]))/(dxv0*dxv1); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv1); 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -(1.0*(5.196152422706631*fhatALVal[2]*H[3]-3.0*fhatALVal[0]*H[3]-3.0*H[1]*fhatALVal[2]+1.732050807568877*fhatALVal[0]*H[1]))/(dxv0*dxv1); 
  incr[1] = -(1.0*(5.196152422706631*H[3]*fhatALVal[3]-3.0*H[1]*fhatALVal[3]-3.0*fhatALVal[1]*H[3]+1.732050807568877*H[1]*fhatALVal[1]))/(dxv0*dxv1); 
  incr[2] = (3.0*(3.0*fhatALVal[2]*H[3]-1.732050807568877*fhatALVal[0]*H[3]-1.732050807568877*H[1]*fhatALVal[2]+fhatALVal[0]*H[1]))/(dxv0*dxv1); 
  incr[3] = (3.0*(3.0*H[3]*fhatALVal[3]-1.732050807568877*H[1]*fhatALVal[3]-1.732050807568877*fhatALVal[1]*H[3]+H[1]*fhatALVal[1]))/(dxv0*dxv1); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
