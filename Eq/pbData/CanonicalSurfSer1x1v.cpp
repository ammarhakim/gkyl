#include <CanonicalModDecl.h> 
double CanonicalSurf1x1vSer_X_P1(const double *w, const double *dxv, const double amax_in, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[0]; 
  double incr[4]; 

  double Ghat[4]; 

  // surface-averaged phase velocity in this direction 
  double alpha0 = std::abs(-(1.0*(3.0*H[3]-1.732050807568877*H[2]))/dxv1); 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = alpha0; 
  else 
    amax = amax_in; 

  Ghat[0] = (2.598076211353316*fr[1]*H[3])/dxv1-(2.598076211353316*fl[1]*H[3])/dxv1-(1.5*fr[0]*H[3])/dxv1-(1.5*fl[0]*H[3])/dxv1-(1.5*fr[1]*H[2])/dxv1+(1.5*fl[1]*H[2])/dxv1+(0.8660254037844386*fr[0]*H[2])/dxv1+(0.8660254037844386*fl[0]*H[2])/dxv1+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = (2.598076211353316*H[3]*fr[3])/dxv1-(1.5*H[2]*fr[3])/dxv1-(2.598076211353316*H[3]*fl[3])/dxv1+(1.5*H[2]*fl[3])/dxv1-(1.5*fr[2]*H[3])/dxv1-(1.5*fl[2]*H[3])/dxv1+(0.8660254037844386*H[2]*fr[2])/dxv1+(0.8660254037844386*H[2]*fl[2])/dxv1+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 

  incr[0] = Ghat[0]/dxv0; 
  incr[1] = -(1.732050807568877*Ghat[0])/dxv0; 
  incr[2] = Ghat[2]/dxv0; 
  incr[3] = -(1.732050807568877*Ghat[2])/dxv0; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
return alpha0; 
} 
double CanonicalSurf1x1vSer_X_P2(const double *w, const double *dxv, const double amax_in, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[0]; 
  double incr[8]; 

  double Ghat[8]; 

  // surface-averaged phase velocity in this direction 
  double alpha0 = std::abs((3.872983346207417*H[6]-3.0*H[3]+1.732050807568877*H[2])/dxv1); 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = alpha0; 
  else 
    amax = amax_in; 

  Ghat[0] = (-(7.5*fr[6]*H[7])/dxv1)-(7.5*fl[6]*H[7])/dxv1+(5.809475019311126*fr[3]*H[7])/dxv1-(5.809475019311126*fl[3]*H[7])/dxv1-(3.354101966249684*fr[2]*H[7])/dxv1-(3.354101966249684*fl[2]*H[7])/dxv1+(4.330127018922194*H[5]*fr[6])/dxv1+(4.330127018922194*H[5]*fl[6])/dxv1+(4.330127018922194*fr[4]*H[6])/dxv1+(4.330127018922194*fl[4]*H[6])/dxv1-(3.354101966249684*fr[1]*H[6])/dxv1+(3.354101966249684*fl[1]*H[6])/dxv1+(1.936491673103709*fr[0]*H[6])/dxv1+(1.936491673103709*fl[0]*H[6])/dxv1-(3.354101966249685*fr[3]*H[5])/dxv1+(3.354101966249685*fl[3]*H[5])/dxv1+(1.936491673103709*fr[2]*H[5])/dxv1+(1.936491673103709*fl[2]*H[5])/dxv1-(3.354101966249685*H[3]*fr[4])/dxv1+(1.936491673103709*H[2]*fr[4])/dxv1-(3.354101966249685*H[3]*fl[4])/dxv1+(1.936491673103709*H[2]*fl[4])/dxv1+(2.598076211353316*fr[1]*H[3])/dxv1-(2.598076211353316*fl[1]*H[3])/dxv1-(1.5*fr[0]*H[3])/dxv1-(1.5*fl[0]*H[3])/dxv1-(1.5*fr[1]*H[2])/dxv1+(1.5*fl[1]*H[2])/dxv1+(0.8660254037844386*fr[0]*H[2])/dxv1+(0.8660254037844386*fl[0]*H[2])/dxv1-1.118033988749895*fr[4]*amax+1.118033988749895*fl[4]*amax+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = (5.196152422706631*H[7]*fr[7])/dxv1-(3.0*H[5]*fr[7])/dxv1-(5.196152422706631*H[7]*fl[7])/dxv1+(3.0*H[5]*fl[7])/dxv1-(3.0*fr[5]*H[7])/dxv1-(3.0*fl[5]*H[7])/dxv1-(7.500000000000001*fr[4]*H[7])/dxv1-(7.500000000000001*fl[4]*H[7])/dxv1+(5.809475019311126*fr[1]*H[7])/dxv1-(5.809475019311126*fl[1]*H[7])/dxv1-(3.354101966249684*fr[0]*H[7])/dxv1-(3.354101966249684*fl[0]*H[7])/dxv1+(4.330127018922193*H[6]*fr[6])/dxv1-(3.354101966249684*H[3]*fr[6])/dxv1+(1.936491673103709*H[2]*fr[6])/dxv1+(4.330127018922193*H[6]*fl[6])/dxv1-(3.354101966249684*H[3]*fl[6])/dxv1+(1.936491673103709*H[2]*fl[6])/dxv1-(3.354101966249684*fr[3]*H[6])/dxv1+(3.354101966249684*fl[3]*H[6])/dxv1+(1.936491673103709*fr[2]*H[6])/dxv1+(1.936491673103709*fl[2]*H[6])/dxv1+(1.732050807568877*H[5]*fr[5])/dxv1+(1.732050807568877*H[5]*fl[5])/dxv1+(4.330127018922193*fr[4]*H[5])/dxv1+(4.330127018922193*fl[4]*H[5])/dxv1-(3.354101966249685*fr[1]*H[5])/dxv1+(3.354101966249685*fl[1]*H[5])/dxv1+(1.936491673103709*fr[0]*H[5])/dxv1+(1.936491673103709*fl[0]*H[5])/dxv1+(2.598076211353316*H[3]*fr[3])/dxv1-(1.5*H[2]*fr[3])/dxv1-(2.598076211353316*H[3]*fl[3])/dxv1+(1.5*H[2]*fl[3])/dxv1-(1.5*fr[2]*H[3])/dxv1-(1.5*fl[2]*H[3])/dxv1+(0.8660254037844386*H[2]*fr[2])/dxv1+(0.8660254037844386*H[2]*fl[2])/dxv1-1.118033988749895*fr[6]*amax+1.118033988749895*fl[6]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[5] = (-(3.354101966249685*H[6]*fr[7])/dxv1)+(2.598076211353316*H[3]*fr[7])/dxv1-(1.5*H[2]*fr[7])/dxv1+(3.354101966249685*H[6]*fl[7])/dxv1-(2.598076211353316*H[3]*fl[7])/dxv1+(1.5*H[2]*fl[7])/dxv1-(6.708203932499369*fr[6]*H[7])/dxv1-(6.708203932499369*fl[6]*H[7])/dxv1+(5.196152422706632*fr[3]*H[7])/dxv1-(5.196152422706632*fl[3]*H[7])/dxv1-(3.0*fr[2]*H[7])/dxv1-(3.0*fl[2]*H[7])/dxv1+(3.872983346207417*H[5]*fr[6])/dxv1+(3.872983346207417*H[5]*fl[6])/dxv1+(1.936491673103709*fr[5]*H[6])/dxv1+(1.936491673103709*fl[5]*H[6])/dxv1-(1.5*H[3]*fr[5])/dxv1+(0.8660254037844386*H[2]*fr[5])/dxv1-(1.5*H[3]*fl[5])/dxv1+(0.8660254037844386*H[2]*fl[5])/dxv1-(3.0*fr[3]*H[5])/dxv1+(3.0*fl[3]*H[5])/dxv1+(1.732050807568877*fr[2]*H[5])/dxv1+(1.732050807568877*fl[2]*H[5])/dxv1+0.8660254037844387*fr[7]*amax+0.8660254037844387*fl[7]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 

  incr[0] = Ghat[0]/dxv0; 
  incr[1] = -(1.732050807568877*Ghat[0])/dxv0; 
  incr[2] = Ghat[2]/dxv0; 
  incr[3] = -(1.732050807568877*Ghat[2])/dxv0; 
  incr[4] = (2.23606797749979*Ghat[0])/dxv0; 
  incr[5] = Ghat[5]/dxv0; 
  incr[6] = (2.23606797749979*Ghat[2])/dxv0; 
  incr[7] = -(1.732050807568877*Ghat[5])/dxv0; 

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
return alpha0; 
} 
double CanonicalSurf1x1vSer_VX_P1(const double *w, const double *dxv, const double amax_in, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[1]; 
  double incr[4]; 

  double Ghat[4]; 

  // surface-averaged phase velocity in this direction 
  double alpha0 = std::abs((3.0*H[3]-1.732050807568877*H[1])/dxv0); 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = alpha0; 
  else 
    amax = amax_in; 

  Ghat[0] = (-(2.598076211353316*fr[2]*H[3])/dxv0)+(2.598076211353316*fl[2]*H[3])/dxv0+(1.5*fr[0]*H[3])/dxv0+(1.5*fl[0]*H[3])/dxv0+(1.5*H[1]*fr[2])/dxv0-(1.5*H[1]*fl[2])/dxv0-(0.8660254037844386*fr[0]*H[1])/dxv0-(0.8660254037844386*fl[0]*H[1])/dxv0+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(2.598076211353316*H[3]*fr[3])/dxv0)+(1.5*H[1]*fr[3])/dxv0+(2.598076211353316*H[3]*fl[3])/dxv0-(1.5*H[1]*fl[3])/dxv0+(1.5*fr[1]*H[3])/dxv0+(1.5*fl[1]*H[3])/dxv0-(0.8660254037844386*H[1]*fr[1])/dxv0-(0.8660254037844386*H[1]*fl[1])/dxv0+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 

  incr[0] = Ghat[0]/dxv1; 
  incr[1] = Ghat[1]/dxv1; 
  incr[2] = -(1.732050807568877*Ghat[0])/dxv1; 
  incr[3] = -(1.732050807568877*Ghat[1])/dxv1; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
return alpha0; 
} 
double CanonicalSurf1x1vSer_VX_P2(const double *w, const double *dxv, const double amax_in, const double *H, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[1]; 
  double incr[8]; 

  double Ghat[8]; 

  // surface-averaged phase velocity in this direction 
  double alpha0 = std::abs(-(1.0*(3.872983346207417*H[7]-3.0*H[3]+1.732050807568877*H[1]))/dxv0); 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = alpha0; 
  else 
    amax = amax_in; 

  Ghat[0] = (7.5*H[6]*fr[7])/dxv0-(4.330127018922194*H[4]*fr[7])/dxv0+(7.5*H[6]*fl[7])/dxv0-(4.330127018922194*H[4]*fl[7])/dxv0-(4.330127018922194*fr[5]*H[7])/dxv0-(4.330127018922194*fl[5]*H[7])/dxv0+(3.354101966249684*fr[2]*H[7])/dxv0-(3.354101966249684*fl[2]*H[7])/dxv0-(1.936491673103709*fr[0]*H[7])/dxv0-(1.936491673103709*fl[0]*H[7])/dxv0-(5.809475019311126*fr[3]*H[6])/dxv0+(5.809475019311126*fl[3]*H[6])/dxv0+(3.354101966249684*fr[1]*H[6])/dxv0+(3.354101966249684*fl[1]*H[6])/dxv0+(3.354101966249685*H[3]*fr[5])/dxv0-(1.936491673103709*H[1]*fr[5])/dxv0+(3.354101966249685*H[3]*fl[5])/dxv0-(1.936491673103709*H[1]*fl[5])/dxv0+(3.354101966249685*fr[3]*H[4])/dxv0-(3.354101966249685*fl[3]*H[4])/dxv0-(1.936491673103709*fr[1]*H[4])/dxv0-(1.936491673103709*fl[1]*H[4])/dxv0-(2.598076211353316*fr[2]*H[3])/dxv0+(2.598076211353316*fl[2]*H[3])/dxv0+(1.5*fr[0]*H[3])/dxv0+(1.5*fl[0]*H[3])/dxv0+(1.5*H[1]*fr[2])/dxv0-(1.5*H[1]*fl[2])/dxv0-(0.8660254037844386*fr[0]*H[1])/dxv0-(0.8660254037844386*fl[0]*H[1])/dxv0-1.118033988749895*fr[5]*amax+1.118033988749895*fl[5]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(4.330127018922193*H[7]*fr[7])/dxv0)+(3.354101966249684*H[3]*fr[7])/dxv0-(1.936491673103709*H[1]*fr[7])/dxv0-(4.330127018922193*H[7]*fl[7])/dxv0+(3.354101966249684*H[3]*fl[7])/dxv0-(1.936491673103709*H[1]*fl[7])/dxv0+(3.354101966249684*fr[3]*H[7])/dxv0-(3.354101966249684*fl[3]*H[7])/dxv0-(1.936491673103709*fr[1]*H[7])/dxv0-(1.936491673103709*fl[1]*H[7])/dxv0-(5.196152422706631*H[6]*fr[6])/dxv0+(3.0*H[4]*fr[6])/dxv0+(5.196152422706631*H[6]*fl[6])/dxv0-(3.0*H[4]*fl[6])/dxv0+(7.500000000000001*fr[5]*H[6])/dxv0+(7.500000000000001*fl[5]*H[6])/dxv0+(3.0*fr[4]*H[6])/dxv0+(3.0*fl[4]*H[6])/dxv0-(5.809475019311126*fr[2]*H[6])/dxv0+(5.809475019311126*fl[2]*H[6])/dxv0+(3.354101966249684*fr[0]*H[6])/dxv0+(3.354101966249684*fl[0]*H[6])/dxv0-(4.330127018922193*H[4]*fr[5])/dxv0-(4.330127018922193*H[4]*fl[5])/dxv0-(1.732050807568877*H[4]*fr[4])/dxv0-(1.732050807568877*H[4]*fl[4])/dxv0+(3.354101966249685*fr[2]*H[4])/dxv0-(3.354101966249685*fl[2]*H[4])/dxv0-(1.936491673103709*fr[0]*H[4])/dxv0-(1.936491673103709*fl[0]*H[4])/dxv0-(2.598076211353316*H[3]*fr[3])/dxv0+(1.5*H[1]*fr[3])/dxv0+(2.598076211353316*H[3]*fl[3])/dxv0-(1.5*H[1]*fl[3])/dxv0+(1.5*fr[1]*H[3])/dxv0+(1.5*fl[1]*H[3])/dxv0-(0.8660254037844386*H[1]*fr[1])/dxv0-(0.8660254037844386*H[1]*fl[1])/dxv0-1.118033988749895*fr[7]*amax+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[4] = (6.708203932499369*H[6]*fr[7])/dxv0-(3.872983346207417*H[4]*fr[7])/dxv0+(6.708203932499369*H[6]*fl[7])/dxv0-(3.872983346207417*H[4]*fl[7])/dxv0+(3.354101966249685*fr[6]*H[7])/dxv0-(3.354101966249685*fl[6]*H[7])/dxv0-(1.936491673103709*fr[4]*H[7])/dxv0-(1.936491673103709*fl[4]*H[7])/dxv0-(2.598076211353316*H[3]*fr[6])/dxv0+(1.5*H[1]*fr[6])/dxv0+(2.598076211353316*H[3]*fl[6])/dxv0-(1.5*H[1]*fl[6])/dxv0-(5.196152422706632*fr[3]*H[6])/dxv0+(5.196152422706632*fl[3]*H[6])/dxv0+(3.0*fr[1]*H[6])/dxv0+(3.0*fl[1]*H[6])/dxv0+(1.5*H[3]*fr[4])/dxv0-(0.8660254037844386*H[1]*fr[4])/dxv0+(1.5*H[3]*fl[4])/dxv0-(0.8660254037844386*H[1]*fl[4])/dxv0+(3.0*fr[3]*H[4])/dxv0-(3.0*fl[3]*H[4])/dxv0-(1.732050807568877*fr[1]*H[4])/dxv0-(1.732050807568877*fl[1]*H[4])/dxv0+0.8660254037844387*fr[6]*amax+0.8660254037844387*fl[6]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax; 

  incr[0] = Ghat[0]/dxv1; 
  incr[1] = Ghat[1]/dxv1; 
  incr[2] = -(1.732050807568877*Ghat[0])/dxv1; 
  incr[3] = -(1.732050807568877*Ghat[1])/dxv1; 
  incr[4] = Ghat[4]/dxv1; 
  incr[5] = (2.23606797749979*Ghat[0])/dxv1; 
  incr[6] = -(1.732050807568877*Ghat[4])/dxv1; 
  incr[7] = (2.23606797749979*Ghat[1])/dxv1; 

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
return alpha0; 
} 
