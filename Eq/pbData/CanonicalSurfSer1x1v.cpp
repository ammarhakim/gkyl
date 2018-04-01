#include <CanonicalModDecl.h> 
#include <cmath> 
double CanonicalSurf1x1vSer_X_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
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
  incr[0] = -(1.0*(5.196152422706631*fl[1]*H[3]+3.0*fl[0]*H[3]-3.0*fl[1]*H[2]-1.732050807568877*fl[0]*H[2]))/(dxv0*dxv1); 
  incr[1] = (3.0*(3.0*fl[1]*H[3]+1.732050807568877*fl[0]*H[3]-1.732050807568877*fl[1]*H[2]-1.0*fl[0]*H[2]))/(dxv0*dxv1); 
  incr[2] = -(1.0*(5.196152422706631*H[3]*fl[3]-3.0*H[2]*fl[3]+3.0*fl[2]*H[3]-1.732050807568877*H[2]*fl[2]))/(dxv0*dxv1); 
  incr[3] = (3.0*(3.0*H[3]*fl[3]-1.732050807568877*H[2]*fl[3]+1.732050807568877*fl[2]*H[3]-1.0*H[2]*fl[2]))/(dxv0*dxv1); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = (5.196152422706631*fr[1]*H[3]-3.0*fr[0]*H[3]-3.0*fr[1]*H[2]+1.732050807568877*fr[0]*H[2])/(dxv0*dxv1); 
  incr[1] = -(3.0*(3.0*fr[1]*H[3]-1.732050807568877*fr[0]*H[3]-1.732050807568877*fr[1]*H[2]+fr[0]*H[2]))/(dxv0*dxv1); 
  incr[2] = (5.196152422706631*H[3]*fr[3]-3.0*H[2]*fr[3]-3.0*fr[2]*H[3]+1.732050807568877*H[2]*fr[2])/(dxv0*dxv1); 
  incr[3] = -(3.0*(3.0*H[3]*fr[3]-1.732050807568877*H[2]*fr[3]-1.732050807568877*fr[2]*H[3]+H[2]*fr[2]))/(dxv0*dxv1); 

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
double CanonicalSurf1x1vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
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
  incr[0] = (5.196152422706631*fl[2]*H[3]+3.0*fl[0]*H[3]-3.0*H[1]*fl[2]-1.732050807568877*fl[0]*H[1])/(dxv0*dxv1); 
  incr[1] = (5.196152422706631*H[3]*fl[3]-3.0*H[1]*fl[3]+3.0*fl[1]*H[3]-1.732050807568877*H[1]*fl[1])/(dxv0*dxv1); 
  incr[2] = -(3.0*(3.0*fl[2]*H[3]+1.732050807568877*fl[0]*H[3]-1.732050807568877*H[1]*fl[2]-1.0*fl[0]*H[1]))/(dxv0*dxv1); 
  incr[3] = -(3.0*(3.0*H[3]*fl[3]-1.732050807568877*H[1]*fl[3]+1.732050807568877*fl[1]*H[3]-1.0*H[1]*fl[1]))/(dxv0*dxv1); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -(1.0*(5.196152422706631*fr[2]*H[3]-3.0*fr[0]*H[3]-3.0*H[1]*fr[2]+1.732050807568877*fr[0]*H[1]))/(dxv0*dxv1); 
  incr[1] = -(1.0*(5.196152422706631*H[3]*fr[3]-3.0*H[1]*fr[3]-3.0*fr[1]*H[3]+1.732050807568877*H[1]*fr[1]))/(dxv0*dxv1); 
  incr[2] = (3.0*(3.0*fr[2]*H[3]-1.732050807568877*fr[0]*H[3]-1.732050807568877*H[1]*fr[2]+fr[0]*H[1]))/(dxv0*dxv1); 
  incr[3] = (3.0*(3.0*H[3]*fr[3]-1.732050807568877*H[1]*fr[3]-1.732050807568877*fr[1]*H[3]+H[1]*fr[1]))/(dxv0*dxv1); 

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
