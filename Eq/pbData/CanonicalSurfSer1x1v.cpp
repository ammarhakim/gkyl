#include <CanonicalModDecl.h> 
double CanonicalSurf1x1vSer_X_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. H: Hamiltonian. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(alpha0) for use in determining amax in cfl and global lax flux 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[0]; 
  double incr[4]; 

  double Ghat[4]; 

  double favg[4]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = -1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  double fjump[4]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(-1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  // cell-average phase velocity in this direction 
  double alpha0 = -(1.0*(3.0*H[3]-1.732050807568877*H[2]))/dxv1; 

  Ghat[0] = (1.5*favg[1]*H[2]+0.8660254037844386*favg[0]*H[2])/dxv1-0.8660254037844386*fjump[1]-0.5*fjump[0]; 
  Ghat[1] = (1.5*favg[1]*H[3]+0.8660254037844386*favg[0]*H[3])/dxv1; 
  Ghat[2] = (1.5*H[2]*favg[3]+0.8660254037844386*H[2]*favg[2])/dxv1-0.8660254037844386*fjump[3]-0.5*fjump[2]; 
  Ghat[3] = (1.5*H[3]*favg[3]+0.8660254037844386*favg[2]*H[3])/dxv1; 

  outr[0] += (1.732050807568877*Ghat[1])/dxv0+Ghat[0]/dxv0; 
  outr[1] += (-(3.0*Ghat[1])/dxv0)-(1.732050807568877*Ghat[0])/dxv0; 
  outr[2] += (1.732050807568877*Ghat[3])/dxv0+Ghat[2]/dxv0; 
  outr[3] += (-(3.0*Ghat[3])/dxv0)-(1.732050807568877*Ghat[2])/dxv0; 

  outl[0] += (-(1.732050807568877*Ghat[1])/dxv0)-(1.0*Ghat[0])/dxv0; 
  outl[1] += (-(3.0*Ghat[1])/dxv0)-(1.732050807568877*Ghat[0])/dxv0; 
  outl[2] += (-(1.732050807568877*Ghat[3])/dxv0)-(1.0*Ghat[2])/dxv0; 
  outl[3] += (-(3.0*Ghat[3])/dxv0)-(1.732050807568877*Ghat[2])/dxv0; 
return std::abs(alpha0); 
} 
double CanonicalSurf1x1vSer_X_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. H: Hamiltonian. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(alpha0) for use in determining amax in cfl and global lax flux 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[0]; 
  double incr[8]; 

  double Ghat[8]; 

  double favg[8]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = -1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  double fjump[8]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(-1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(1*fr[6]-fl[6]); 
  fjump[7] = amax*(-1*fr[7]-fl[7]); 
  // cell-average phase velocity in this direction 
  double alpha0 = (3.872983346207417*H[6]-3.0*H[3]+1.732050807568877*H[2])/dxv1; 

  Ghat[0] = (4.330127018922194*H[5]*favg[6]+3.354101966249685*favg[3]*H[5]+1.936491673103709*favg[2]*H[5]+1.936491673103709*H[2]*favg[4]+1.5*favg[1]*H[2]+0.8660254037844386*favg[0]*H[2])/dxv1-1.118033988749895*fjump[4]-0.8660254037844386*fjump[1]-0.5*fjump[0]; 
  Ghat[1] = (4.330127018922193*favg[6]*H[7]+3.354101966249684*favg[3]*H[7]+1.936491673103709*favg[2]*H[7]+1.936491673103709*H[3]*favg[4]+1.5*favg[1]*H[3]+0.8660254037844386*favg[0]*H[3])/dxv1; 
  Ghat[2] = (3.0*H[5]*favg[7]+1.936491673103709*H[2]*favg[6]+1.732050807568877*H[5]*favg[5]+4.330127018922193*favg[4]*H[5]+3.354101966249685*favg[1]*H[5]+1.936491673103709*favg[0]*H[5]+1.5*H[2]*favg[3]+0.8660254037844386*H[2]*favg[2])/dxv1-1.118033988749895*fjump[6]-0.8660254037844386*fjump[3]-0.5*fjump[2]; 
  Ghat[3] = (3.0*H[7]*favg[7]+1.732050807568877*favg[5]*H[7]+4.330127018922194*favg[4]*H[7]+3.354101966249684*favg[1]*H[7]+1.936491673103709*favg[0]*H[7]+1.936491673103709*H[3]*favg[6]+1.5*H[3]*favg[3]+0.8660254037844386*favg[2]*H[3])/dxv1; 
  Ghat[4] = (1.936491673103709*favg[4]*H[6]+1.5*favg[1]*H[6]+0.8660254037844387*favg[0]*H[6])/dxv1; 
  Ghat[5] = (1.5*H[2]*favg[7]+3.872983346207417*H[5]*favg[6]+0.8660254037844386*H[2]*favg[5]+3.0*favg[3]*H[5]+1.732050807568877*favg[2]*H[5])/dxv1-0.8660254037844387*fjump[7]-0.5*fjump[5]; 
  Ghat[6] = (1.936491673103709*H[6]*favg[6]+1.5*favg[3]*H[6]+0.8660254037844386*favg[2]*H[6])/dxv1; 
  Ghat[7] = (1.5*H[3]*favg[7]+3.872983346207417*favg[6]*H[7]+3.0*favg[3]*H[7]+1.732050807568877*favg[2]*H[7]+0.8660254037844387*H[3]*favg[5])/dxv1; 

  outr[0] += (2.23606797749979*Ghat[4])/dxv0+(1.732050807568877*Ghat[1])/dxv0+Ghat[0]/dxv0; 
  outr[1] += (-(3.872983346207417*Ghat[4])/dxv0)-(3.0*Ghat[1])/dxv0-(1.732050807568877*Ghat[0])/dxv0; 
  outr[2] += (2.23606797749979*Ghat[6])/dxv0+(1.732050807568877*Ghat[3])/dxv0+Ghat[2]/dxv0; 
  outr[3] += (-(3.872983346207417*Ghat[6])/dxv0)-(3.0*Ghat[3])/dxv0-(1.732050807568877*Ghat[2])/dxv0; 
  outr[4] += (5.0*Ghat[4])/dxv0+(3.872983346207417*Ghat[1])/dxv0+(2.23606797749979*Ghat[0])/dxv0; 
  outr[5] += (1.732050807568877*Ghat[7])/dxv0+Ghat[5]/dxv0; 
  outr[6] += (5.0*Ghat[6])/dxv0+(3.872983346207417*Ghat[3])/dxv0+(2.23606797749979*Ghat[2])/dxv0; 
  outr[7] += (-(3.0*Ghat[7])/dxv0)-(1.732050807568877*Ghat[5])/dxv0; 

  outl[0] += (-(2.23606797749979*Ghat[4])/dxv0)-(1.732050807568877*Ghat[1])/dxv0-(1.0*Ghat[0])/dxv0; 
  outl[1] += (-(3.872983346207417*Ghat[4])/dxv0)-(3.0*Ghat[1])/dxv0-(1.732050807568877*Ghat[0])/dxv0; 
  outl[2] += (-(2.23606797749979*Ghat[6])/dxv0)-(1.732050807568877*Ghat[3])/dxv0-(1.0*Ghat[2])/dxv0; 
  outl[3] += (-(3.872983346207417*Ghat[6])/dxv0)-(3.0*Ghat[3])/dxv0-(1.732050807568877*Ghat[2])/dxv0; 
  outl[4] += (-(5.0*Ghat[4])/dxv0)-(3.872983346207417*Ghat[1])/dxv0-(2.23606797749979*Ghat[0])/dxv0; 
  outl[5] += (-(1.732050807568877*Ghat[7])/dxv0)-(1.0*Ghat[5])/dxv0; 
  outl[6] += (-(5.0*Ghat[6])/dxv0)-(3.872983346207417*Ghat[3])/dxv0-(2.23606797749979*Ghat[2])/dxv0; 
  outl[7] += (-(3.0*Ghat[7])/dxv0)-(1.732050807568877*Ghat[5])/dxv0; 
return std::abs(alpha0); 
} 
double CanonicalSurf1x1vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. H: Hamiltonian. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(alpha0) for use in determining amax in cfl and global lax flux 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[1]; 
  double incr[4]; 

  double Ghat[4]; 

  double favg[4]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  double fjump[4]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  // cell-average phase velocity in this direction 
  double alpha0 = (3.0*H[3]-1.732050807568877*H[1])/dxv0; 

  Ghat[0] = ((-1.5*H[1]*favg[2])-0.8660254037844386*favg[0]*H[1])/dxv0-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = ((-1.5*H[1]*favg[3])-0.8660254037844386*H[1]*favg[1])/dxv0-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[2] = ((-1.5*favg[2]*H[3])-0.8660254037844386*favg[0]*H[3])/dxv0; 
  Ghat[3] = ((-1.5*H[3]*favg[3])-0.8660254037844386*favg[1]*H[3])/dxv0; 

  outr[0] += (1.732050807568877*Ghat[2])/dxv1+Ghat[0]/dxv1; 
  outr[1] += (1.732050807568877*Ghat[3])/dxv1+Ghat[1]/dxv1; 
  outr[2] += (-(3.0*Ghat[2])/dxv1)-(1.732050807568877*Ghat[0])/dxv1; 
  outr[3] += (-(3.0*Ghat[3])/dxv1)-(1.732050807568877*Ghat[1])/dxv1; 

  outl[0] += (-(1.732050807568877*Ghat[2])/dxv1)-(1.0*Ghat[0])/dxv1; 
  outl[1] += (-(1.732050807568877*Ghat[3])/dxv1)-(1.0*Ghat[1])/dxv1; 
  outl[2] += (-(3.0*Ghat[2])/dxv1)-(1.732050807568877*Ghat[0])/dxv1; 
  outl[3] += (-(3.0*Ghat[3])/dxv1)-(1.732050807568877*Ghat[1])/dxv1; 
return std::abs(alpha0); 
} 
double CanonicalSurf1x1vSer_VX_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. H: Hamiltonian. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(alpha0) for use in determining amax in cfl and global lax flux 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[1]; 
  double incr[8]; 

  double Ghat[8]; 

  double favg[8]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  double fjump[8]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  // cell-average phase velocity in this direction 
  double alpha0 = -(1.0*(3.872983346207417*H[7]-3.0*H[3]+1.732050807568877*H[1]))/dxv0; 

  Ghat[0] = ((-4.330127018922194*H[4]*favg[7])-1.936491673103709*H[1]*favg[5]-3.354101966249685*favg[3]*H[4]-1.936491673103709*favg[1]*H[4]-1.5*H[1]*favg[2]-0.8660254037844386*favg[0]*H[1])/dxv0-1.118033988749895*fjump[5]-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = ((-1.936491673103709*H[1]*favg[7])-3.0*H[4]*favg[6]-4.330127018922193*H[4]*favg[5]-1.732050807568877*H[4]*favg[4]-3.354101966249685*favg[2]*H[4]-1.936491673103709*favg[0]*H[4]-1.5*H[1]*favg[3]-0.8660254037844386*H[1]*favg[1])/dxv0-1.118033988749895*fjump[7]-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[2] = ((-4.330127018922193*H[6]*favg[7])-3.354101966249684*favg[3]*H[6]-1.936491673103709*favg[1]*H[6]-1.936491673103709*H[3]*favg[5]-1.5*favg[2]*H[3]-0.8660254037844386*favg[0]*H[3])/dxv0; 
  Ghat[3] = ((-1.936491673103709*H[3]*favg[7])-3.0*H[6]*favg[6]-4.330127018922194*favg[5]*H[6]-1.732050807568877*favg[4]*H[6]-3.354101966249684*favg[2]*H[6]-1.936491673103709*favg[0]*H[6]-1.5*H[3]*favg[3]-0.8660254037844386*favg[1]*H[3])/dxv0; 
  Ghat[4] = ((-3.872983346207417*H[4]*favg[7])-1.5*H[1]*favg[6]-0.8660254037844386*H[1]*favg[4]-3.0*favg[3]*H[4]-1.732050807568877*favg[1]*H[4])/dxv0-0.8660254037844387*fjump[6]-0.5*fjump[4]; 
  Ghat[5] = ((-1.936491673103709*favg[5]*H[7])-1.5*favg[2]*H[7]-0.8660254037844387*favg[0]*H[7])/dxv0; 
  Ghat[6] = ((-3.872983346207417*H[6]*favg[7])-1.5*H[3]*favg[6]-3.0*favg[3]*H[6]-1.732050807568877*favg[1]*H[6]-0.8660254037844387*H[3]*favg[4])/dxv0; 
  Ghat[7] = ((-1.936491673103709*H[7]*favg[7])-1.5*favg[3]*H[7]-0.8660254037844386*favg[1]*H[7])/dxv0; 

  outr[0] += (2.23606797749979*Ghat[5])/dxv1+(1.732050807568877*Ghat[2])/dxv1+Ghat[0]/dxv1; 
  outr[1] += (2.23606797749979*Ghat[7])/dxv1+(1.732050807568877*Ghat[3])/dxv1+Ghat[1]/dxv1; 
  outr[2] += (-(3.872983346207417*Ghat[5])/dxv1)-(3.0*Ghat[2])/dxv1-(1.732050807568877*Ghat[0])/dxv1; 
  outr[3] += (-(3.872983346207417*Ghat[7])/dxv1)-(3.0*Ghat[3])/dxv1-(1.732050807568877*Ghat[1])/dxv1; 
  outr[4] += (1.732050807568877*Ghat[6])/dxv1+Ghat[4]/dxv1; 
  outr[5] += (5.0*Ghat[5])/dxv1+(3.872983346207417*Ghat[2])/dxv1+(2.23606797749979*Ghat[0])/dxv1; 
  outr[6] += (-(3.0*Ghat[6])/dxv1)-(1.732050807568877*Ghat[4])/dxv1; 
  outr[7] += (5.0*Ghat[7])/dxv1+(3.872983346207417*Ghat[3])/dxv1+(2.23606797749979*Ghat[1])/dxv1; 

  outl[0] += (-(2.23606797749979*Ghat[5])/dxv1)-(1.732050807568877*Ghat[2])/dxv1-(1.0*Ghat[0])/dxv1; 
  outl[1] += (-(2.23606797749979*Ghat[7])/dxv1)-(1.732050807568877*Ghat[3])/dxv1-(1.0*Ghat[1])/dxv1; 
  outl[2] += (-(3.872983346207417*Ghat[5])/dxv1)-(3.0*Ghat[2])/dxv1-(1.732050807568877*Ghat[0])/dxv1; 
  outl[3] += (-(3.872983346207417*Ghat[7])/dxv1)-(3.0*Ghat[3])/dxv1-(1.732050807568877*Ghat[1])/dxv1; 
  outl[4] += (-(1.732050807568877*Ghat[6])/dxv1)-(1.0*Ghat[4])/dxv1; 
  outl[5] += (-(5.0*Ghat[5])/dxv1)-(3.872983346207417*Ghat[2])/dxv1-(2.23606797749979*Ghat[0])/dxv1; 
  outl[6] += (-(3.0*Ghat[6])/dxv1)-(1.732050807568877*Ghat[4])/dxv1; 
  outl[7] += (-(5.0*Ghat[7])/dxv1)-(3.872983346207417*Ghat[3])/dxv1-(2.23606797749979*Ghat[1])/dxv1; 
return std::abs(alpha0); 
} 
