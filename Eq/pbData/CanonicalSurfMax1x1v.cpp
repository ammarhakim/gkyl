#include <CanonicalModDecl.h> 
void CanonicalSurf1x1vMax_X_P1(const double *w, const double *dxv, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[0]; 
  double incr[3]; 

  // surface-averaged phase velocity in this direction 
  alpha0 = (1.732050807568877*H[2])/dxv1; 

  if (alpha0>0) { 
  incr[0] = 0.25*(3.0*fl[1]+1.732050807568877*fl[0])*H[2]; 
  incr[1] = -0.75*(1.732050807568877*fl[1]+fl[0])*H[2]; 
  incr[2] = 0.4330127018922193*H[2]*fl[2]; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  } else { 
  incr[0] = -0.25*(3.0*fr[1]-1.732050807568877*fr[0])*H[2]; 
  incr[1] = 0.75*(1.732050807568877*fr[1]-1.0*fr[0])*H[2]; 
  incr[2] = 0.4330127018922193*H[2]*fr[2]; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  } 
} 
void CanonicalSurf1x1vMax_X_P2(const double *w, const double *dxv, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[0]; 
  double incr[6]; 

  // surface-averaged phase velocity in this direction 
  alpha0 = -(1.0*(3.0*H[3]-1.732050807568877*H[2]))/dxv1; 

  if (alpha0>0) { 
  incr[0] = 0.25*(6.708203932499369*fl[3]*H[5]+3.872983346207417*fl[2]*H[5]-6.708203932499369*H[3]*fl[4]+3.872983346207417*H[2]*fl[4]-5.196152422706631*fl[1]*H[3]-3.0*fl[0]*H[3]+3.0*fl[1]*H[2]+1.732050807568877*fl[0]*H[2]); 
  incr[1] = -0.75*(3.872983346207417*fl[3]*H[5]+2.23606797749979*fl[2]*H[5]-3.872983346207417*H[3]*fl[4]+2.23606797749979*H[2]*fl[4]-3.0*fl[1]*H[3]-1.732050807568877*fl[0]*H[3]+1.732050807568877*fl[1]*H[2]+fl[0]*H[2]); 
  incr[2] = 0.25*(3.464101615137754*H[5]*fl[5]+8.660254037844386*fl[4]*H[5]+6.708203932499369*fl[1]*H[5]+3.872983346207417*fl[0]*H[5]-5.196152422706631*H[3]*fl[3]+3.0*H[2]*fl[3]-3.0*fl[2]*H[3]+1.732050807568877*H[2]*fl[2]); 
  incr[3] = -0.75*(2.0*H[5]*fl[5]+5.0*fl[4]*H[5]+3.872983346207417*fl[1]*H[5]+2.23606797749979*fl[0]*H[5]-3.0*H[3]*fl[3]+1.732050807568877*H[2]*fl[3]-1.732050807568877*fl[2]*H[3]+H[2]*fl[2]); 
  incr[4] = 0.25*(15.0*fl[3]*H[5]+8.660254037844386*fl[2]*H[5]-15.0*H[3]*fl[4]+8.660254037844386*H[2]*fl[4]-11.61895003862225*fl[1]*H[3]-6.708203932499369*fl[0]*H[3]+6.708203932499369*fl[1]*H[2]+3.872983346207417*fl[0]*H[2]); 
  incr[5] = -0.25*(3.0*H[3]*fl[5]-1.732050807568877*H[2]*fl[5]-6.0*fl[3]*H[5]-3.464101615137754*fl[2]*H[5]); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  } else { 
  incr[0] = -0.25*(6.708203932499369*fr[3]*H[5]-3.872983346207417*fr[2]*H[5]+6.708203932499369*H[3]*fr[4]-3.872983346207417*H[2]*fr[4]-5.196152422706631*fr[1]*H[3]+3.0*fr[0]*H[3]+3.0*fr[1]*H[2]-1.732050807568877*fr[0]*H[2]); 
  incr[1] = 0.75*(3.872983346207417*fr[3]*H[5]-2.23606797749979*fr[2]*H[5]+3.872983346207417*H[3]*fr[4]-2.23606797749979*H[2]*fr[4]-3.0*fr[1]*H[3]+1.732050807568877*fr[0]*H[3]+1.732050807568877*fr[1]*H[2]-1.0*fr[0]*H[2]); 
  incr[2] = 0.25*(3.464101615137754*H[5]*fr[5]+8.660254037844386*fr[4]*H[5]-6.708203932499369*fr[1]*H[5]+3.872983346207417*fr[0]*H[5]+5.196152422706631*H[3]*fr[3]-3.0*H[2]*fr[3]-3.0*fr[2]*H[3]+1.732050807568877*H[2]*fr[2]); 
  incr[3] = -0.75*(2.0*H[5]*fr[5]+5.0*fr[4]*H[5]-3.872983346207417*fr[1]*H[5]+2.23606797749979*fr[0]*H[5]+3.0*H[3]*fr[3]-1.732050807568877*H[2]*fr[3]-1.732050807568877*fr[2]*H[3]+H[2]*fr[2]); 
  incr[4] = -0.25*(15.0*fr[3]*H[5]-8.660254037844386*fr[2]*H[5]+15.0*H[3]*fr[4]-8.660254037844386*H[2]*fr[4]-11.61895003862225*fr[1]*H[3]+6.708203932499369*fr[0]*H[3]+6.708203932499369*fr[1]*H[2]-3.872983346207417*fr[0]*H[2]); 
  incr[5] = -0.25*(3.0*H[3]*fr[5]-1.732050807568877*H[2]*fr[5]+6.0*fr[3]*H[5]-3.464101615137754*fr[2]*H[5]); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  } 
} 
void CanonicalSurf1x1vMax_VX_P1(const double *w, const double *dxv, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[1]; 
  double incr[3]; 

  // surface-averaged phase velocity in this direction 
  alpha0 = -(1.732050807568877*H[1])/dxv0; 

  if (alpha0>0) { 
  incr[0] = -0.25*H[1]*(3.0*fl[2]+1.732050807568877*fl[0]); 
  incr[1] = -0.4330127018922193*H[1]*fl[1]; 
  incr[2] = 0.75*H[1]*(1.732050807568877*fl[2]+fl[0]); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  } else { 
  incr[0] = 0.25*H[1]*(3.0*fr[2]-1.732050807568877*fr[0]); 
  incr[1] = -0.4330127018922193*H[1]*fr[1]; 
  incr[2] = -0.75*H[1]*(1.732050807568877*fr[2]-1.0*fr[0]); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  } 
} 
void CanonicalSurf1x1vMax_VX_P2(const double *w, const double *dxv, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[1]; 
  double incr[6]; 

  // surface-averaged phase velocity in this direction 
  alpha0 = (3.0*H[3]-1.732050807568877*H[1])/dxv0; 

  if (alpha0>0) { 
  incr[0] = 0.25*(6.708203932499369*H[3]*fl[5]-3.872983346207417*H[1]*fl[5]-6.708203932499369*fl[3]*H[4]-3.872983346207417*fl[1]*H[4]+5.196152422706631*fl[2]*H[3]+3.0*fl[0]*H[3]-3.0*H[1]*fl[2]-1.732050807568877*fl[0]*H[1]); 
  incr[1] = -0.25*(8.660254037844386*H[4]*fl[5]+3.464101615137754*H[4]*fl[4]+6.708203932499369*fl[2]*H[4]+3.872983346207417*fl[0]*H[4]-5.196152422706631*H[3]*fl[3]+3.0*H[1]*fl[3]-3.0*fl[1]*H[3]+1.732050807568877*H[1]*fl[1]); 
  incr[2] = -0.75*(3.872983346207417*H[3]*fl[5]-2.23606797749979*H[1]*fl[5]-3.872983346207417*fl[3]*H[4]-2.23606797749979*fl[1]*H[4]+3.0*fl[2]*H[3]+1.732050807568877*fl[0]*H[3]-1.732050807568877*H[1]*fl[2]-1.0*fl[0]*H[1]); 
  incr[3] = 0.75*(5.0*H[4]*fl[5]+2.0*H[4]*fl[4]+3.872983346207417*fl[2]*H[4]+2.23606797749979*fl[0]*H[4]-3.0*H[3]*fl[3]+1.732050807568877*H[1]*fl[3]-1.732050807568877*fl[1]*H[3]+H[1]*fl[1]); 
  incr[4] = 0.25*(3.0*H[3]*fl[4]-1.732050807568877*H[1]*fl[4]-6.0*fl[3]*H[4]-3.464101615137754*fl[1]*H[4]); 
  incr[5] = 0.25*(15.0*H[3]*fl[5]-8.660254037844386*H[1]*fl[5]-15.0*fl[3]*H[4]-8.660254037844386*fl[1]*H[4]+11.61895003862225*fl[2]*H[3]+6.708203932499369*fl[0]*H[3]-6.708203932499369*H[1]*fl[2]-3.872983346207417*fl[0]*H[1]); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  } else { 
  incr[0] = 0.25*(6.708203932499369*H[3]*fr[5]-3.872983346207417*H[1]*fr[5]+6.708203932499369*fr[3]*H[4]-3.872983346207417*fr[1]*H[4]-5.196152422706631*fr[2]*H[3]+3.0*fr[0]*H[3]+3.0*H[1]*fr[2]-1.732050807568877*fr[0]*H[1]); 
  incr[1] = -0.25*(8.660254037844386*H[4]*fr[5]+3.464101615137754*H[4]*fr[4]-6.708203932499369*fr[2]*H[4]+3.872983346207417*fr[0]*H[4]+5.196152422706631*H[3]*fr[3]-3.0*H[1]*fr[3]-3.0*fr[1]*H[3]+1.732050807568877*H[1]*fr[1]); 
  incr[2] = -0.75*(3.872983346207417*H[3]*fr[5]-2.23606797749979*H[1]*fr[5]+3.872983346207417*fr[3]*H[4]-2.23606797749979*fr[1]*H[4]-3.0*fr[2]*H[3]+1.732050807568877*fr[0]*H[3]+1.732050807568877*H[1]*fr[2]-1.0*fr[0]*H[1]); 
  incr[3] = 0.75*(5.0*H[4]*fr[5]+2.0*H[4]*fr[4]-3.872983346207417*fr[2]*H[4]+2.23606797749979*fr[0]*H[4]+3.0*H[3]*fr[3]-1.732050807568877*H[1]*fr[3]-1.732050807568877*fr[1]*H[3]+H[1]*fr[1]); 
  incr[4] = 0.25*(3.0*H[3]*fr[4]-1.732050807568877*H[1]*fr[4]+6.0*fr[3]*H[4]-3.464101615137754*fr[1]*H[4]); 
  incr[5] = 0.25*(15.0*H[3]*fr[5]-8.660254037844386*H[1]*fr[5]+15.0*fr[3]*H[4]-8.660254037844386*fr[1]*H[4]-11.61895003862225*fr[2]*H[3]+6.708203932499369*fr[0]*H[3]+6.708203932499369*H[1]*fr[2]-3.872983346207417*fr[0]*H[1]); 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  } 
} 
