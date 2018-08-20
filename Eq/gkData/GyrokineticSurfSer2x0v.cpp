#include <GyrokineticModDecl.h> 
double GyrokineticSurf2x0vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.25*(3.0*BmagInv[0]*Phi[3]-1.732050807568877*BmagInv[0]*Phi[2])*dfac_y; 

  double alpha[4]; 
  alpha[0] = 1.5*BmagInv[0]*Phi[3]*dfac_y-0.8660254037844386*BmagInv[0]*Phi[2]*dfac_y; 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fl[3]+fl[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(3.0*fl[3]+1.732050807568877*fl[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fr[3]-1.0*fr[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(3.0*fr[3]-1.732050807568877*fr[2])*dfac_x; 

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
double GyrokineticSurf2x0vSer_X_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.25*(3.872983346207417*BmagInv[0]*Phi[6]-3.0*BmagInv[0]*Phi[3]+1.732050807568877*BmagInv[0]*Phi[2])*dfac_y; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  Ghat[0] = 1.875*BmagInv[0]*fr[6]*Phi[7]*dfac_y+1.875*BmagInv[0]*fl[6]*Phi[7]*dfac_y-1.452368754827781*BmagInv[0]*fr[3]*Phi[7]*dfac_y+1.452368754827781*BmagInv[0]*fl[3]*Phi[7]*dfac_y+0.8385254915624211*BmagInv[0]*fr[2]*Phi[7]*dfac_y+0.8385254915624211*BmagInv[0]*fl[2]*Phi[7]*dfac_y-1.082531754730548*BmagInv[0]*Phi[5]*fr[6]*dfac_y-1.082531754730548*BmagInv[0]*Phi[5]*fl[6]*dfac_y-1.082531754730548*BmagInv[0]*fr[4]*Phi[6]*dfac_y-1.082531754730548*BmagInv[0]*fl[4]*Phi[6]*dfac_y+0.8385254915624211*BmagInv[0]*fr[1]*Phi[6]*dfac_y-0.8385254915624211*BmagInv[0]*fl[1]*Phi[6]*dfac_y-0.4841229182759271*BmagInv[0]*fr[0]*Phi[6]*dfac_y-0.4841229182759271*BmagInv[0]*fl[0]*Phi[6]*dfac_y+0.8385254915624212*BmagInv[0]*fr[3]*Phi[5]*dfac_y-0.8385254915624212*BmagInv[0]*fl[3]*Phi[5]*dfac_y-0.4841229182759271*BmagInv[0]*fr[2]*Phi[5]*dfac_y-0.4841229182759271*BmagInv[0]*fl[2]*Phi[5]*dfac_y+0.8385254915624212*BmagInv[0]*Phi[3]*fr[4]*dfac_y-0.4841229182759271*BmagInv[0]*Phi[2]*fr[4]*dfac_y+0.8385254915624212*BmagInv[0]*Phi[3]*fl[4]*dfac_y-0.4841229182759271*BmagInv[0]*Phi[2]*fl[4]*dfac_y-0.6495190528383289*BmagInv[0]*fr[1]*Phi[3]*dfac_y+0.6495190528383289*BmagInv[0]*fl[1]*Phi[3]*dfac_y+0.375*BmagInv[0]*fr[0]*Phi[3]*dfac_y+0.375*BmagInv[0]*fl[0]*Phi[3]*dfac_y+0.375*BmagInv[0]*fr[1]*Phi[2]*dfac_y-0.375*BmagInv[0]*fl[1]*Phi[2]*dfac_y-0.2165063509461096*BmagInv[0]*fr[0]*Phi[2]*dfac_y-0.2165063509461096*BmagInv[0]*fl[0]*Phi[2]*dfac_y-1.118033988749895*fr[4]*amax+1.118033988749895*fl[4]*amax+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = (-1.299038105676658*BmagInv[0]*Phi[7]*fr[7]*dfac_y)+0.75*BmagInv[0]*Phi[5]*fr[7]*dfac_y+1.299038105676658*BmagInv[0]*Phi[7]*fl[7]*dfac_y-0.75*BmagInv[0]*Phi[5]*fl[7]*dfac_y+0.75*BmagInv[0]*fr[5]*Phi[7]*dfac_y+0.75*BmagInv[0]*fl[5]*Phi[7]*dfac_y+1.875*BmagInv[0]*fr[4]*Phi[7]*dfac_y+1.875*BmagInv[0]*fl[4]*Phi[7]*dfac_y-1.452368754827781*BmagInv[0]*fr[1]*Phi[7]*dfac_y+1.452368754827781*BmagInv[0]*fl[1]*Phi[7]*dfac_y+0.8385254915624211*BmagInv[0]*fr[0]*Phi[7]*dfac_y+0.8385254915624211*BmagInv[0]*fl[0]*Phi[7]*dfac_y-1.082531754730548*BmagInv[0]*Phi[6]*fr[6]*dfac_y+0.8385254915624211*BmagInv[0]*Phi[3]*fr[6]*dfac_y-0.4841229182759271*BmagInv[0]*Phi[2]*fr[6]*dfac_y-1.082531754730548*BmagInv[0]*Phi[6]*fl[6]*dfac_y+0.8385254915624211*BmagInv[0]*Phi[3]*fl[6]*dfac_y-0.4841229182759271*BmagInv[0]*Phi[2]*fl[6]*dfac_y+0.8385254915624211*BmagInv[0]*fr[3]*Phi[6]*dfac_y-0.8385254915624211*BmagInv[0]*fl[3]*Phi[6]*dfac_y-0.4841229182759271*BmagInv[0]*fr[2]*Phi[6]*dfac_y-0.4841229182759271*BmagInv[0]*fl[2]*Phi[6]*dfac_y-0.4330127018922193*BmagInv[0]*Phi[5]*fr[5]*dfac_y-0.4330127018922193*BmagInv[0]*Phi[5]*fl[5]*dfac_y-1.082531754730548*BmagInv[0]*fr[4]*Phi[5]*dfac_y-1.082531754730548*BmagInv[0]*fl[4]*Phi[5]*dfac_y+0.8385254915624212*BmagInv[0]*fr[1]*Phi[5]*dfac_y-0.8385254915624212*BmagInv[0]*fl[1]*Phi[5]*dfac_y-0.4841229182759271*BmagInv[0]*fr[0]*Phi[5]*dfac_y-0.4841229182759271*BmagInv[0]*fl[0]*Phi[5]*dfac_y-0.6495190528383289*BmagInv[0]*Phi[3]*fr[3]*dfac_y+0.375*BmagInv[0]*Phi[2]*fr[3]*dfac_y+0.6495190528383289*BmagInv[0]*Phi[3]*fl[3]*dfac_y-0.375*BmagInv[0]*Phi[2]*fl[3]*dfac_y+0.375*BmagInv[0]*fr[2]*Phi[3]*dfac_y+0.375*BmagInv[0]*fl[2]*Phi[3]*dfac_y-0.2165063509461096*BmagInv[0]*Phi[2]*fr[2]*dfac_y-0.2165063509461096*BmagInv[0]*Phi[2]*fl[2]*dfac_y-1.118033988749895*fr[6]*amax+1.118033988749895*fl[6]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[5] = 0.8385254915624212*BmagInv[0]*Phi[6]*fr[7]*dfac_y-0.649519052838329*BmagInv[0]*Phi[3]*fr[7]*dfac_y+0.375*BmagInv[0]*Phi[2]*fr[7]*dfac_y-0.8385254915624212*BmagInv[0]*Phi[6]*fl[7]*dfac_y+0.649519052838329*BmagInv[0]*Phi[3]*fl[7]*dfac_y-0.375*BmagInv[0]*Phi[2]*fl[7]*dfac_y+1.677050983124842*BmagInv[0]*fr[6]*Phi[7]*dfac_y+1.677050983124842*BmagInv[0]*fl[6]*Phi[7]*dfac_y-1.299038105676658*BmagInv[0]*fr[3]*Phi[7]*dfac_y+1.299038105676658*BmagInv[0]*fl[3]*Phi[7]*dfac_y+0.75*BmagInv[0]*fr[2]*Phi[7]*dfac_y+0.75*BmagInv[0]*fl[2]*Phi[7]*dfac_y-0.9682458365518543*BmagInv[0]*Phi[5]*fr[6]*dfac_y-0.9682458365518543*BmagInv[0]*Phi[5]*fl[6]*dfac_y-0.4841229182759271*BmagInv[0]*fr[5]*Phi[6]*dfac_y-0.4841229182759271*BmagInv[0]*fl[5]*Phi[6]*dfac_y+0.375*BmagInv[0]*Phi[3]*fr[5]*dfac_y-0.2165063509461096*BmagInv[0]*Phi[2]*fr[5]*dfac_y+0.375*BmagInv[0]*Phi[3]*fl[5]*dfac_y-0.2165063509461096*BmagInv[0]*Phi[2]*fl[5]*dfac_y+0.75*BmagInv[0]*fr[3]*Phi[5]*dfac_y-0.75*BmagInv[0]*fl[3]*Phi[5]*dfac_y-0.4330127018922193*BmagInv[0]*fr[2]*Phi[5]*dfac_y-0.4330127018922193*BmagInv[0]*fl[2]*Phi[5]*dfac_y+0.8660254037844387*fr[7]*amax+0.8660254037844387*fl[7]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[4] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[5] = 0.5*Ghat[5]*dfac_x; 
  incr[6] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[7] = -0.8660254037844387*Ghat[5]*dfac_x; 

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
double GyrokineticSurf2x0vSer_Y_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.25*(3.0*BmagInv[0]*Phi[3]-1.732050807568877*BmagInv[0]*Phi[1])*dfac_x; 

  double alpha[4]; 
  alpha[0] = 0.8660254037844386*BmagInv[0]*Phi[1]*dfac_x-1.5*BmagInv[0]*Phi[3]*dfac_x; 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[2]+fl[0])*dfac_y; 
  incr[1] = 0.25*alpha[0]*(1.732050807568877*fl[3]+fl[1])*dfac_y; 
  incr[2] = -0.25*alpha[0]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_y; 
  incr[3] = -0.25*alpha[0]*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_y; 
  incr[1] = -0.25*alpha[0]*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_y; 
  incr[2] = 0.25*alpha[0]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_y; 
  incr[3] = 0.25*alpha[0]*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_y; 

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
double GyrokineticSurf2x0vSer_Y_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.25*(3.872983346207417*BmagInv[0]*Phi[7]-3.0*BmagInv[0]*Phi[3]+1.732050807568877*BmagInv[0]*Phi[1])*dfac_x; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  Ghat[0] = (-1.875*BmagInv[0]*Phi[6]*fr[7]*dfac_x)+1.082531754730548*BmagInv[0]*Phi[4]*fr[7]*dfac_x-1.875*BmagInv[0]*Phi[6]*fl[7]*dfac_x+1.082531754730548*BmagInv[0]*Phi[4]*fl[7]*dfac_x+1.082531754730548*BmagInv[0]*fr[5]*Phi[7]*dfac_x+1.082531754730548*BmagInv[0]*fl[5]*Phi[7]*dfac_x-0.8385254915624211*BmagInv[0]*fr[2]*Phi[7]*dfac_x+0.8385254915624211*BmagInv[0]*fl[2]*Phi[7]*dfac_x+0.4841229182759271*BmagInv[0]*fr[0]*Phi[7]*dfac_x+0.4841229182759271*BmagInv[0]*fl[0]*Phi[7]*dfac_x+1.452368754827781*BmagInv[0]*fr[3]*Phi[6]*dfac_x-1.452368754827781*BmagInv[0]*fl[3]*Phi[6]*dfac_x-0.8385254915624211*BmagInv[0]*fr[1]*Phi[6]*dfac_x-0.8385254915624211*BmagInv[0]*fl[1]*Phi[6]*dfac_x-0.8385254915624212*BmagInv[0]*Phi[3]*fr[5]*dfac_x+0.4841229182759271*BmagInv[0]*Phi[1]*fr[5]*dfac_x-0.8385254915624212*BmagInv[0]*Phi[3]*fl[5]*dfac_x+0.4841229182759271*BmagInv[0]*Phi[1]*fl[5]*dfac_x-0.8385254915624212*BmagInv[0]*fr[3]*Phi[4]*dfac_x+0.8385254915624212*BmagInv[0]*fl[3]*Phi[4]*dfac_x+0.4841229182759271*BmagInv[0]*fr[1]*Phi[4]*dfac_x+0.4841229182759271*BmagInv[0]*fl[1]*Phi[4]*dfac_x+0.6495190528383289*BmagInv[0]*fr[2]*Phi[3]*dfac_x-0.6495190528383289*BmagInv[0]*fl[2]*Phi[3]*dfac_x-0.375*BmagInv[0]*fr[0]*Phi[3]*dfac_x-0.375*BmagInv[0]*fl[0]*Phi[3]*dfac_x-0.375*BmagInv[0]*Phi[1]*fr[2]*dfac_x+0.375*BmagInv[0]*Phi[1]*fl[2]*dfac_x+0.2165063509461096*BmagInv[0]*fr[0]*Phi[1]*dfac_x+0.2165063509461096*BmagInv[0]*fl[0]*Phi[1]*dfac_x-1.118033988749895*fr[5]*amax+1.118033988749895*fl[5]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = 1.082531754730548*BmagInv[0]*Phi[7]*fr[7]*dfac_x-0.8385254915624211*BmagInv[0]*Phi[3]*fr[7]*dfac_x+0.4841229182759271*BmagInv[0]*Phi[1]*fr[7]*dfac_x+1.082531754730548*BmagInv[0]*Phi[7]*fl[7]*dfac_x-0.8385254915624211*BmagInv[0]*Phi[3]*fl[7]*dfac_x+0.4841229182759271*BmagInv[0]*Phi[1]*fl[7]*dfac_x-0.8385254915624211*BmagInv[0]*fr[3]*Phi[7]*dfac_x+0.8385254915624211*BmagInv[0]*fl[3]*Phi[7]*dfac_x+0.4841229182759271*BmagInv[0]*fr[1]*Phi[7]*dfac_x+0.4841229182759271*BmagInv[0]*fl[1]*Phi[7]*dfac_x+1.299038105676658*BmagInv[0]*Phi[6]*fr[6]*dfac_x-0.75*BmagInv[0]*Phi[4]*fr[6]*dfac_x-1.299038105676658*BmagInv[0]*Phi[6]*fl[6]*dfac_x+0.75*BmagInv[0]*Phi[4]*fl[6]*dfac_x-1.875*BmagInv[0]*fr[5]*Phi[6]*dfac_x-1.875*BmagInv[0]*fl[5]*Phi[6]*dfac_x-0.75*BmagInv[0]*fr[4]*Phi[6]*dfac_x-0.75*BmagInv[0]*fl[4]*Phi[6]*dfac_x+1.452368754827781*BmagInv[0]*fr[2]*Phi[6]*dfac_x-1.452368754827781*BmagInv[0]*fl[2]*Phi[6]*dfac_x-0.8385254915624211*BmagInv[0]*fr[0]*Phi[6]*dfac_x-0.8385254915624211*BmagInv[0]*fl[0]*Phi[6]*dfac_x+1.082531754730548*BmagInv[0]*Phi[4]*fr[5]*dfac_x+1.082531754730548*BmagInv[0]*Phi[4]*fl[5]*dfac_x+0.4330127018922193*BmagInv[0]*Phi[4]*fr[4]*dfac_x+0.4330127018922193*BmagInv[0]*Phi[4]*fl[4]*dfac_x-0.8385254915624212*BmagInv[0]*fr[2]*Phi[4]*dfac_x+0.8385254915624212*BmagInv[0]*fl[2]*Phi[4]*dfac_x+0.4841229182759271*BmagInv[0]*fr[0]*Phi[4]*dfac_x+0.4841229182759271*BmagInv[0]*fl[0]*Phi[4]*dfac_x+0.6495190528383289*BmagInv[0]*Phi[3]*fr[3]*dfac_x-0.375*BmagInv[0]*Phi[1]*fr[3]*dfac_x-0.6495190528383289*BmagInv[0]*Phi[3]*fl[3]*dfac_x+0.375*BmagInv[0]*Phi[1]*fl[3]*dfac_x-0.375*BmagInv[0]*fr[1]*Phi[3]*dfac_x-0.375*BmagInv[0]*fl[1]*Phi[3]*dfac_x+0.2165063509461096*BmagInv[0]*Phi[1]*fr[1]*dfac_x+0.2165063509461096*BmagInv[0]*Phi[1]*fl[1]*dfac_x-1.118033988749895*fr[7]*amax+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[4] = (-1.677050983124842*BmagInv[0]*Phi[6]*fr[7]*dfac_x)+0.9682458365518543*BmagInv[0]*Phi[4]*fr[7]*dfac_x-1.677050983124842*BmagInv[0]*Phi[6]*fl[7]*dfac_x+0.9682458365518543*BmagInv[0]*Phi[4]*fl[7]*dfac_x-0.8385254915624212*BmagInv[0]*fr[6]*Phi[7]*dfac_x+0.8385254915624212*BmagInv[0]*fl[6]*Phi[7]*dfac_x+0.4841229182759271*BmagInv[0]*fr[4]*Phi[7]*dfac_x+0.4841229182759271*BmagInv[0]*fl[4]*Phi[7]*dfac_x+0.649519052838329*BmagInv[0]*Phi[3]*fr[6]*dfac_x-0.375*BmagInv[0]*Phi[1]*fr[6]*dfac_x-0.649519052838329*BmagInv[0]*Phi[3]*fl[6]*dfac_x+0.375*BmagInv[0]*Phi[1]*fl[6]*dfac_x+1.299038105676658*BmagInv[0]*fr[3]*Phi[6]*dfac_x-1.299038105676658*BmagInv[0]*fl[3]*Phi[6]*dfac_x-0.75*BmagInv[0]*fr[1]*Phi[6]*dfac_x-0.75*BmagInv[0]*fl[1]*Phi[6]*dfac_x-0.375*BmagInv[0]*Phi[3]*fr[4]*dfac_x+0.2165063509461096*BmagInv[0]*Phi[1]*fr[4]*dfac_x-0.375*BmagInv[0]*Phi[3]*fl[4]*dfac_x+0.2165063509461096*BmagInv[0]*Phi[1]*fl[4]*dfac_x-0.75*BmagInv[0]*fr[3]*Phi[4]*dfac_x+0.75*BmagInv[0]*fl[3]*Phi[4]*dfac_x+0.4330127018922193*BmagInv[0]*fr[1]*Phi[4]*dfac_x+0.4330127018922193*BmagInv[0]*fl[1]*Phi[4]*dfac_x+0.8660254037844387*fr[6]*amax+0.8660254037844387*fl[6]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_y; 
  incr[1] = 0.5*Ghat[1]*dfac_y; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_y; 
  incr[3] = -0.8660254037844386*Ghat[1]*dfac_y; 
  incr[4] = 0.5*Ghat[4]*dfac_y; 
  incr[5] = 1.118033988749895*Ghat[0]*dfac_y; 
  incr[6] = -0.8660254037844387*Ghat[4]*dfac_y; 
  incr[7] = 1.118033988749895*Ghat[1]*dfac_y; 

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
