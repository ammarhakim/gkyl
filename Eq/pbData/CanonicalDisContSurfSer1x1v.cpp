#include <CanonicalModDecl.h> 
void CanonicalDisContCorrectionSurf1x1vSer_X_P1(const double *w, const double *dxv, const double amax, const double *Hl, const double* Hr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[0]; 
  double incr[4]; 

  // check if hamiltonian is continuous in this direction, and if so exit 
  double tol = 1e-10; 
  if (
       std::abs(1.732050807568877*Hr[1]+1.732050807568877*Hl[1]-1.0*Hr[0]+Hl[0]) < tol && 
       std::abs(1.732050807568877*Hr[3]+1.732050807568877*Hl[3]-1.0*Hr[2]+Hl[2]) < tol && 
       true 
     ) 
  { 
    return; 
  } 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.0*(3.0*Hr[3]-1.732050807568877*Hr[2]))/dxv1; 

  if (alpha0>0) { 
  outr[2] += (5.196152422706631*Hl[3]*fl[3]+3.0*Hl[2]*fl[3]+3.0*fl[2]*Hl[3]+1.732050807568877*Hl[2]*fl[2]+5.196152422706631*Hl[1]*fl[1]+3.0*Hl[0]*fl[1]+3.0*fl[0]*Hl[1]+1.732050807568877*Hl[0]*fl[0])/(dxv0*dxv1); 
  outr[3] += -(3.0*(3.0*Hl[3]*fl[3]+1.732050807568877*Hl[2]*fl[3]+1.732050807568877*fl[2]*Hl[3]+Hl[2]*fl[2]+3.0*Hl[1]*fl[1]+1.732050807568877*Hl[0]*fl[1]+1.732050807568877*fl[0]*Hl[1]+Hl[0]*fl[0]))/(dxv0*dxv1); 

  outl[2] += ((5.196152422706631*Hr[3]-3.0*Hr[2])*fl[3]+3.0*fl[2]*Hr[3]-1.732050807568877*Hr[2]*fl[2]+(5.196152422706631*Hr[1]-3.0*Hr[0])*fl[1]+3.0*fl[0]*Hr[1]-1.732050807568877*Hr[0]*fl[0])/(dxv0*dxv1); 
  outl[3] += ((9.0*Hr[3]-5.196152422706631*Hr[2])*fl[3]+5.196152422706631*fl[2]*Hr[3]-3.0*Hr[2]*fl[2]+(9.0*Hr[1]-5.196152422706631*Hr[0])*fl[1]+5.196152422706631*fl[0]*Hr[1]-3.0*Hr[0]*fl[0])/(dxv0*dxv1); 
  } else { 
  outr[2] += -(1.0*(5.196152422706631*Hl[3]*fr[3]+3.0*Hl[2]*fr[3]-3.0*fr[2]*Hl[3]-1.732050807568877*Hl[2]*fr[2]+5.196152422706631*Hl[1]*fr[1]+3.0*Hl[0]*fr[1]-3.0*fr[0]*Hl[1]-1.732050807568877*Hl[0]*fr[0]))/(dxv0*dxv1); 
  outr[3] += (3.0*(3.0*Hl[3]*fr[3]+1.732050807568877*Hl[2]*fr[3]-1.732050807568877*fr[2]*Hl[3]-1.0*Hl[2]*fr[2]+3.0*Hl[1]*fr[1]+1.732050807568877*Hl[0]*fr[1]-1.732050807568877*fr[0]*Hl[1]-1.0*Hl[0]*fr[0]))/(dxv0*dxv1); 

  outl[2] += -(1.0*((5.196152422706631*Hr[3]-3.0*Hr[2])*fr[3]-3.0*fr[2]*Hr[3]+1.732050807568877*Hr[2]*fr[2]+(5.196152422706631*Hr[1]-3.0*Hr[0])*fr[1]-3.0*fr[0]*Hr[1]+1.732050807568877*Hr[0]*fr[0]))/(dxv0*dxv1); 
  outl[3] += -(1.0*((9.0*Hr[3]-5.196152422706631*Hr[2])*fr[3]-5.196152422706631*fr[2]*Hr[3]+3.0*Hr[2]*fr[2]+(9.0*Hr[1]-5.196152422706631*Hr[0])*fr[1]-5.196152422706631*fr[0]*Hr[1]+3.0*Hr[0]*fr[0]))/(dxv0*dxv1); 
  } 
} 
void CanonicalDisContCorrectionSurf1x1vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *Hl, const double* Hr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. H: Hamiltonian, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double wd = w[1]; 
  double incr[4]; 

  // check if hamiltonian is continuous in this direction, and if so exit 
  double tol = 1e-10; 
  if (
       std::abs(1.732050807568877*Hr[2]+1.732050807568877*Hl[2]-1.0*Hr[0]+Hl[0]) < tol && 
       std::abs(1.732050807568877*Hr[3]+1.732050807568877*Hl[3]-1.0*Hr[1]+Hl[1]) < tol && 
       true 
     ) 
  { 
    return; 
  } 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (3.0*Hr[3]-1.732050807568877*Hr[1])/dxv0; 

  if (alpha0>0) { 
  outr[1] += -(1.0*(5.196152422706631*Hl[3]*fl[3]+3.0*Hl[1]*fl[3]+3.0*fl[1]*Hl[3]+5.196152422706631*Hl[2]*fl[2]+3.0*Hl[0]*fl[2]+3.0*fl[0]*Hl[2]+1.732050807568877*Hl[1]*fl[1]+1.732050807568877*Hl[0]*fl[0]))/(dxv0*dxv1); 
  outr[3] += (3.0*(3.0*Hl[3]*fl[3]+1.732050807568877*Hl[1]*fl[3]+1.732050807568877*fl[1]*Hl[3]+3.0*Hl[2]*fl[2]+1.732050807568877*Hl[0]*fl[2]+1.732050807568877*fl[0]*Hl[2]+Hl[1]*fl[1]+Hl[0]*fl[0]))/(dxv0*dxv1); 

  outl[1] += -(1.0*((5.196152422706631*Hr[3]-3.0*Hr[1])*fl[3]+3.0*fl[1]*Hr[3]+(5.196152422706631*Hr[2]-3.0*Hr[0])*fl[2]+3.0*fl[0]*Hr[2]-1.732050807568877*Hr[1]*fl[1]-1.732050807568877*Hr[0]*fl[0]))/(dxv0*dxv1); 
  outl[3] += -(1.0*((9.0*Hr[3]-5.196152422706631*Hr[1])*fl[3]+5.196152422706631*fl[1]*Hr[3]+(9.0*Hr[2]-5.196152422706631*Hr[0])*fl[2]+5.196152422706631*fl[0]*Hr[2]-3.0*Hr[1]*fl[1]-3.0*Hr[0]*fl[0]))/(dxv0*dxv1); 
  } else { 
  outr[1] += (5.196152422706631*Hl[3]*fr[3]+3.0*Hl[1]*fr[3]-3.0*fr[1]*Hl[3]+5.196152422706631*Hl[2]*fr[2]+3.0*Hl[0]*fr[2]-3.0*fr[0]*Hl[2]-1.732050807568877*Hl[1]*fr[1]-1.732050807568877*Hl[0]*fr[0])/(dxv0*dxv1); 
  outr[3] += -(3.0*(3.0*Hl[3]*fr[3]+1.732050807568877*Hl[1]*fr[3]-1.732050807568877*fr[1]*Hl[3]+3.0*Hl[2]*fr[2]+1.732050807568877*Hl[0]*fr[2]-1.732050807568877*fr[0]*Hl[2]-1.0*Hl[1]*fr[1]-1.0*Hl[0]*fr[0]))/(dxv0*dxv1); 

  outl[1] += ((5.196152422706631*Hr[3]-3.0*Hr[1])*fr[3]-3.0*fr[1]*Hr[3]+(5.196152422706631*Hr[2]-3.0*Hr[0])*fr[2]-3.0*fr[0]*Hr[2]+1.732050807568877*Hr[1]*fr[1]+1.732050807568877*Hr[0]*fr[0])/(dxv0*dxv1); 
  outl[3] += ((9.0*Hr[3]-5.196152422706631*Hr[1])*fr[3]-5.196152422706631*fr[1]*Hr[3]+(9.0*Hr[2]-5.196152422706631*Hr[0])*fr[2]-5.196152422706631*fr[0]*Hr[2]+3.0*Hr[1]*fr[1]+3.0*Hr[0]*fr[0])/(dxv0*dxv1); 
  } 
} 
