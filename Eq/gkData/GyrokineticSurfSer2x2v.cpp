#include <GyrokineticModDecl.h> 
double GyrokineticSurf2x2vSer_X_P1(const double mcByq, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double Ghat[16]; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = std::abs((3.0*Phi[3]-1.732050807568877*Phi[2])*dfac_y*mcByq); 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = alpha0; 
  else 
    amax = amax_in; 

  Ghat[0] = (-1.299038105676658*fr[1]*Phi[3]*dfac_y*mcByq)+1.299038105676658*fl[1]*Phi[3]*dfac_y*mcByq+0.75*fr[0]*Phi[3]*dfac_y*mcByq+0.75*fl[0]*Phi[3]*dfac_y*mcByq+0.75*fr[1]*Phi[2]*dfac_y*mcByq-0.75*fl[1]*Phi[2]*dfac_y*mcByq-0.4330127018922193*fr[0]*Phi[2]*dfac_y*mcByq-0.4330127018922193*fl[0]*Phi[2]*dfac_y*mcByq+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = (-1.299038105676658*Phi[3]*fr[5]*dfac_y*mcByq)+0.75*Phi[2]*fr[5]*dfac_y*mcByq+1.299038105676658*Phi[3]*fl[5]*dfac_y*mcByq-0.75*Phi[2]*fl[5]*dfac_y*mcByq+0.75*fr[2]*Phi[3]*dfac_y*mcByq+0.75*fl[2]*Phi[3]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fr[2]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fl[2]*dfac_y*mcByq+0.8660254037844386*fr[5]*amax+0.8660254037844386*fl[5]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[3] = (-1.299038105676658*Phi[3]*fr[6]*dfac_y*mcByq)+0.75*Phi[2]*fr[6]*dfac_y*mcByq+1.299038105676658*Phi[3]*fl[6]*dfac_y*mcByq-0.75*Phi[2]*fl[6]*dfac_y*mcByq+0.75*Phi[3]*fr[3]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fr[3]*dfac_y*mcByq+0.75*Phi[3]*fl[3]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fl[3]*dfac_y*mcByq+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[4] = (-1.299038105676658*Phi[3]*fr[8]*dfac_y*mcByq)+0.75*Phi[2]*fr[8]*dfac_y*mcByq+1.299038105676658*Phi[3]*fl[8]*dfac_y*mcByq-0.75*Phi[2]*fl[8]*dfac_y*mcByq+0.75*Phi[3]*fr[4]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fr[4]*dfac_y*mcByq+0.75*Phi[3]*fl[4]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fl[4]*dfac_y*mcByq+0.8660254037844386*fr[8]*amax+0.8660254037844386*fl[8]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax; 
  Ghat[7] = (-1.299038105676658*Phi[3]*fr[11]*dfac_y*mcByq)+0.75*Phi[2]*fr[11]*dfac_y*mcByq+1.299038105676658*Phi[3]*fl[11]*dfac_y*mcByq-0.75*Phi[2]*fl[11]*dfac_y*mcByq+0.75*Phi[3]*fr[7]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fr[7]*dfac_y*mcByq+0.75*Phi[3]*fl[7]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fl[7]*dfac_y*mcByq+0.8660254037844386*fr[11]*amax+0.8660254037844386*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax; 
  Ghat[9] = (-1.299038105676658*Phi[3]*fr[12]*dfac_y*mcByq)+0.75*Phi[2]*fr[12]*dfac_y*mcByq+1.299038105676658*Phi[3]*fl[12]*dfac_y*mcByq-0.75*Phi[2]*fl[12]*dfac_y*mcByq+0.75*Phi[3]*fr[9]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fr[9]*dfac_y*mcByq+0.75*Phi[3]*fl[9]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fl[9]*dfac_y*mcByq+0.8660254037844386*fr[12]*amax+0.8660254037844386*fl[12]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[10] = (-1.299038105676658*Phi[3]*fr[13]*dfac_y*mcByq)+0.75*Phi[2]*fr[13]*dfac_y*mcByq+1.299038105676658*Phi[3]*fl[13]*dfac_y*mcByq-0.75*Phi[2]*fl[13]*dfac_y*mcByq+0.75*Phi[3]*fr[10]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fr[10]*dfac_y*mcByq+0.75*Phi[3]*fl[10]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fl[10]*dfac_y*mcByq+0.8660254037844386*fr[13]*amax+0.8660254037844386*fl[13]*amax-0.5*fr[10]*amax+0.5*fl[10]*amax; 
  Ghat[14] = (-1.299038105676658*Phi[3]*fr[15]*dfac_y*mcByq)+0.75*Phi[2]*fr[15]*dfac_y*mcByq+1.299038105676658*Phi[3]*fl[15]*dfac_y*mcByq-0.75*Phi[2]*fl[15]*dfac_y*mcByq+0.75*Phi[3]*fr[14]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fr[14]*dfac_y*mcByq+0.75*Phi[3]*fl[14]*dfac_y*mcByq-0.4330127018922193*Phi[2]*fl[14]*dfac_y*mcByq+0.8660254037844386*fr[15]*amax+0.8660254037844386*fl[15]*amax-0.5*fr[14]*amax+0.5*fl[14]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = 0.5*Ghat[3]*dfac_x; 
  incr[4] = 0.5*Ghat[4]*dfac_x; 
  incr[5] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[6] = -0.8660254037844386*Ghat[3]*dfac_x; 
  incr[7] = 0.5*Ghat[7]*dfac_x; 
  incr[8] = -0.8660254037844386*Ghat[4]*dfac_x; 
  incr[9] = 0.5*Ghat[9]*dfac_x; 
  incr[10] = 0.5*Ghat[10]*dfac_x; 
  incr[11] = -0.8660254037844386*Ghat[7]*dfac_x; 
  incr[12] = -0.8660254037844386*Ghat[9]*dfac_x; 
  incr[13] = -0.8660254037844386*Ghat[10]*dfac_x; 
  incr[14] = 0.5*Ghat[14]*dfac_x; 
  incr[15] = -0.8660254037844386*Ghat[14]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
return alpha0; 
} 
double GyrokineticSurf2x2vSer_Y_P1(const double mcByq, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double Ghat[16]; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = std::abs(BcurvY[0]*mcByq*wv2+1.732050807568877*Bmag[1]*dfac_x*mcByq*wm+(1.732050807568877*Phi[1]-3.0*Phi[3])*dfac_x*mcByq); 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = alpha0; 
  else 
    amax = amax_in; 

  Ghat[0] = (-0.4330127018922193*BcurvY[1]*fr[5]*mcByq*wv2)+0.4330127018922193*BcurvY[1]*fl[5]*mcByq*wv2-0.4330127018922193*BcurvY[0]*fr[2]*mcByq*wv2+0.4330127018922193*BcurvY[0]*fl[2]*mcByq*wv2+0.25*BcurvY[1]*fr[1]*mcByq*wv2+0.25*BcurvY[1]*fl[1]*mcByq*wv2+0.25*BcurvY[0]*fr[0]*mcByq*wv2+0.25*BcurvY[0]*fl[0]*mcByq*wv2-(0.25*BcurvY[1]*fr[11]*mcByq*wv)/dfac_v+(0.25*BcurvY[1]*fl[11]*mcByq*wv)/dfac_v-(0.25*BcurvY[0]*fr[7]*mcByq*wv)/dfac_v+(0.25*BcurvY[0]*fl[7]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fr[6]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fl[6]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fr[3]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fl[3]*mcByq*wv)/dfac_v-0.75*Bmag[1]*fr[2]*dfac_x*mcByq*wm+0.75*Bmag[1]*fl[2]*dfac_x*mcByq*wm+0.4330127018922193*fr[0]*Bmag[1]*dfac_x*mcByq*wm+0.4330127018922193*fl[0]*Bmag[1]*dfac_x*mcByq*wm-(0.4330127018922193*Bmag[1]*fr[9]*dfac_x*mcByq)/dfac_m+(0.4330127018922193*Bmag[1]*fl[9]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fr[4]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fl[4]*dfac_x*mcByq)/dfac_m+1.299038105676658*fr[2]*Phi[3]*dfac_x*mcByq-1.299038105676658*fl[2]*Phi[3]*dfac_x*mcByq-0.75*fr[0]*Phi[3]*dfac_x*mcByq-0.75*fl[0]*Phi[3]*dfac_x*mcByq-0.75*Phi[1]*fr[2]*dfac_x*mcByq+0.75*Phi[1]*fl[2]*dfac_x*mcByq+0.4330127018922193*fr[0]*Phi[1]*dfac_x*mcByq+0.4330127018922193*fl[0]*Phi[1]*dfac_x*mcByq+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-0.4330127018922193*BcurvY[0]*fr[5]*mcByq*wv2)+0.4330127018922193*BcurvY[0]*fl[5]*mcByq*wv2-0.4330127018922193*BcurvY[1]*fr[2]*mcByq*wv2+0.4330127018922193*BcurvY[1]*fl[2]*mcByq*wv2+0.25*BcurvY[0]*fr[1]*mcByq*wv2+0.25*BcurvY[0]*fl[1]*mcByq*wv2+0.25*fr[0]*BcurvY[1]*mcByq*wv2+0.25*fl[0]*BcurvY[1]*mcByq*wv2-(0.25*BcurvY[0]*fr[11]*mcByq*wv)/dfac_v+(0.25*BcurvY[0]*fl[11]*mcByq*wv)/dfac_v-(0.25*BcurvY[1]*fr[7]*mcByq*wv)/dfac_v+(0.25*BcurvY[1]*fl[7]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fr[6]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fl[6]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fr[3]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fl[3]*mcByq*wv)/dfac_v-0.75*Bmag[1]*fr[5]*dfac_x*mcByq*wm+0.75*Bmag[1]*fl[5]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fr[1]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fl[1]*dfac_x*mcByq*wm-(0.4330127018922193*Bmag[1]*fr[12]*dfac_x*mcByq)/dfac_m+(0.4330127018922193*Bmag[1]*fl[12]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fr[8]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fl[8]*dfac_x*mcByq)/dfac_m+1.299038105676658*Phi[3]*fr[5]*dfac_x*mcByq-0.75*Phi[1]*fr[5]*dfac_x*mcByq-1.299038105676658*Phi[3]*fl[5]*dfac_x*mcByq+0.75*Phi[1]*fl[5]*dfac_x*mcByq-0.75*fr[1]*Phi[3]*dfac_x*mcByq-0.75*fl[1]*Phi[3]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fr[1]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fl[1]*dfac_x*mcByq+0.8660254037844386*fr[5]*amax+0.8660254037844386*fl[5]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[3] = (-0.4330127018922193*BcurvY[1]*fr[11]*mcByq*wv2)+0.4330127018922193*BcurvY[1]*fl[11]*mcByq*wv2-0.4330127018922193*BcurvY[0]*fr[7]*mcByq*wv2+0.4330127018922193*BcurvY[0]*fl[7]*mcByq*wv2+0.25*BcurvY[1]*fr[6]*mcByq*wv2+0.25*BcurvY[1]*fl[6]*mcByq*wv2+0.25*BcurvY[0]*fr[3]*mcByq*wv2+0.25*BcurvY[0]*fl[3]*mcByq*wv2-(0.25*BcurvY[1]*fr[5]*mcByq*wv)/dfac_v+(0.25*BcurvY[1]*fl[5]*mcByq*wv)/dfac_v-(0.25*BcurvY[0]*fr[2]*mcByq*wv)/dfac_v+(0.25*BcurvY[0]*fl[2]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fr[1]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fl[1]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fr[0]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fl[0]*mcByq*wv)/dfac_v-0.75*Bmag[1]*fr[7]*dfac_x*mcByq*wm+0.75*Bmag[1]*fl[7]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fr[3]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fl[3]*dfac_x*mcByq*wm-(0.4330127018922193*Bmag[1]*fr[14]*dfac_x*mcByq)/dfac_m+(0.4330127018922193*Bmag[1]*fl[14]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fr[10]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fl[10]*dfac_x*mcByq)/dfac_m+1.299038105676658*Phi[3]*fr[7]*dfac_x*mcByq-0.75*Phi[1]*fr[7]*dfac_x*mcByq-1.299038105676658*Phi[3]*fl[7]*dfac_x*mcByq+0.75*Phi[1]*fl[7]*dfac_x*mcByq-0.75*Phi[3]*fr[3]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fr[3]*dfac_x*mcByq-0.75*Phi[3]*fl[3]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fl[3]*dfac_x*mcByq+0.8660254037844386*fr[7]*amax+0.8660254037844386*fl[7]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[4] = (-0.4330127018922193*BcurvY[1]*fr[12]*mcByq*wv2)+0.4330127018922193*BcurvY[1]*fl[12]*mcByq*wv2-0.4330127018922193*BcurvY[0]*fr[9]*mcByq*wv2+0.4330127018922193*BcurvY[0]*fl[9]*mcByq*wv2+0.25*BcurvY[1]*fr[8]*mcByq*wv2+0.25*BcurvY[1]*fl[8]*mcByq*wv2+0.25*BcurvY[0]*fr[4]*mcByq*wv2+0.25*BcurvY[0]*fl[4]*mcByq*wv2-(0.25*BcurvY[1]*fr[15]*mcByq*wv)/dfac_v+(0.25*BcurvY[1]*fl[15]*mcByq*wv)/dfac_v-(0.25*BcurvY[0]*fr[14]*mcByq*wv)/dfac_v+(0.25*BcurvY[0]*fl[14]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fr[13]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fl[13]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fr[10]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fl[10]*mcByq*wv)/dfac_v-0.75*Bmag[1]*fr[9]*dfac_x*mcByq*wm+0.75*Bmag[1]*fl[9]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fr[4]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fl[4]*dfac_x*mcByq*wm-(0.4330127018922193*Bmag[1]*fr[2]*dfac_x*mcByq)/dfac_m+(0.4330127018922193*Bmag[1]*fl[2]*dfac_x*mcByq)/dfac_m+(0.25*fr[0]*Bmag[1]*dfac_x*mcByq)/dfac_m+(0.25*fl[0]*Bmag[1]*dfac_x*mcByq)/dfac_m+1.299038105676658*Phi[3]*fr[9]*dfac_x*mcByq-0.75*Phi[1]*fr[9]*dfac_x*mcByq-1.299038105676658*Phi[3]*fl[9]*dfac_x*mcByq+0.75*Phi[1]*fl[9]*dfac_x*mcByq-0.75*Phi[3]*fr[4]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fr[4]*dfac_x*mcByq-0.75*Phi[3]*fl[4]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fl[4]*dfac_x*mcByq+0.8660254037844386*fr[9]*amax+0.8660254037844386*fl[9]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax; 
  Ghat[6] = (-0.4330127018922193*BcurvY[0]*fr[11]*mcByq*wv2)+0.4330127018922193*BcurvY[0]*fl[11]*mcByq*wv2-0.4330127018922193*BcurvY[1]*fr[7]*mcByq*wv2+0.4330127018922193*BcurvY[1]*fl[7]*mcByq*wv2+0.25*BcurvY[0]*fr[6]*mcByq*wv2+0.25*BcurvY[0]*fl[6]*mcByq*wv2+0.25*BcurvY[1]*fr[3]*mcByq*wv2+0.25*BcurvY[1]*fl[3]*mcByq*wv2-(0.25*BcurvY[0]*fr[5]*mcByq*wv)/dfac_v+(0.25*BcurvY[0]*fl[5]*mcByq*wv)/dfac_v-(0.25*BcurvY[1]*fr[2]*mcByq*wv)/dfac_v+(0.25*BcurvY[1]*fl[2]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fr[1]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fl[1]*mcByq*wv)/dfac_v+(0.1443375672974065*fr[0]*BcurvY[1]*mcByq*wv)/dfac_v+(0.1443375672974065*fl[0]*BcurvY[1]*mcByq*wv)/dfac_v-0.75*Bmag[1]*fr[11]*dfac_x*mcByq*wm+0.75*Bmag[1]*fl[11]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fr[6]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fl[6]*dfac_x*mcByq*wm-(0.4330127018922193*Bmag[1]*fr[15]*dfac_x*mcByq)/dfac_m+(0.4330127018922193*Bmag[1]*fl[15]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fr[13]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fl[13]*dfac_x*mcByq)/dfac_m+1.299038105676658*Phi[3]*fr[11]*dfac_x*mcByq-0.75*Phi[1]*fr[11]*dfac_x*mcByq-1.299038105676658*Phi[3]*fl[11]*dfac_x*mcByq+0.75*Phi[1]*fl[11]*dfac_x*mcByq-0.75*Phi[3]*fr[6]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fr[6]*dfac_x*mcByq-0.75*Phi[3]*fl[6]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fl[6]*dfac_x*mcByq+0.8660254037844386*fr[11]*amax+0.8660254037844386*fl[11]*amax-0.5*fr[6]*amax+0.5*fl[6]*amax; 
  Ghat[8] = (-0.4330127018922193*BcurvY[0]*fr[12]*mcByq*wv2)+0.4330127018922193*BcurvY[0]*fl[12]*mcByq*wv2-0.4330127018922193*BcurvY[1]*fr[9]*mcByq*wv2+0.4330127018922193*BcurvY[1]*fl[9]*mcByq*wv2+0.25*BcurvY[0]*fr[8]*mcByq*wv2+0.25*BcurvY[0]*fl[8]*mcByq*wv2+0.25*BcurvY[1]*fr[4]*mcByq*wv2+0.25*BcurvY[1]*fl[4]*mcByq*wv2-(0.25*BcurvY[0]*fr[15]*mcByq*wv)/dfac_v+(0.25*BcurvY[0]*fl[15]*mcByq*wv)/dfac_v-(0.25*BcurvY[1]*fr[14]*mcByq*wv)/dfac_v+(0.25*BcurvY[1]*fl[14]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fr[13]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fl[13]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fr[10]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fl[10]*mcByq*wv)/dfac_v-0.75*Bmag[1]*fr[12]*dfac_x*mcByq*wm+0.75*Bmag[1]*fl[12]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fr[8]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fl[8]*dfac_x*mcByq*wm-(0.4330127018922193*Bmag[1]*fr[5]*dfac_x*mcByq)/dfac_m+(0.4330127018922193*Bmag[1]*fl[5]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fr[1]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fl[1]*dfac_x*mcByq)/dfac_m+1.299038105676658*Phi[3]*fr[12]*dfac_x*mcByq-0.75*Phi[1]*fr[12]*dfac_x*mcByq-1.299038105676658*Phi[3]*fl[12]*dfac_x*mcByq+0.75*Phi[1]*fl[12]*dfac_x*mcByq-0.75*Phi[3]*fr[8]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fr[8]*dfac_x*mcByq-0.75*Phi[3]*fl[8]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fl[8]*dfac_x*mcByq+0.8660254037844386*fr[12]*amax+0.8660254037844386*fl[12]*amax-0.5*fr[8]*amax+0.5*fl[8]*amax; 
  Ghat[10] = (-0.4330127018922193*BcurvY[1]*fr[15]*mcByq*wv2)+0.4330127018922193*BcurvY[1]*fl[15]*mcByq*wv2-0.4330127018922193*BcurvY[0]*fr[14]*mcByq*wv2+0.4330127018922193*BcurvY[0]*fl[14]*mcByq*wv2+0.25*BcurvY[1]*fr[13]*mcByq*wv2+0.25*BcurvY[1]*fl[13]*mcByq*wv2+0.25*BcurvY[0]*fr[10]*mcByq*wv2+0.25*BcurvY[0]*fl[10]*mcByq*wv2-(0.25*BcurvY[1]*fr[12]*mcByq*wv)/dfac_v+(0.25*BcurvY[1]*fl[12]*mcByq*wv)/dfac_v-(0.25*BcurvY[0]*fr[9]*mcByq*wv)/dfac_v+(0.25*BcurvY[0]*fl[9]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fr[8]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fl[8]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fr[4]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fl[4]*mcByq*wv)/dfac_v-0.75*Bmag[1]*fr[14]*dfac_x*mcByq*wm+0.75*Bmag[1]*fl[14]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fr[10]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fl[10]*dfac_x*mcByq*wm-(0.4330127018922193*Bmag[1]*fr[7]*dfac_x*mcByq)/dfac_m+(0.4330127018922193*Bmag[1]*fl[7]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fr[3]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fl[3]*dfac_x*mcByq)/dfac_m+1.299038105676658*Phi[3]*fr[14]*dfac_x*mcByq-0.75*Phi[1]*fr[14]*dfac_x*mcByq-1.299038105676658*Phi[3]*fl[14]*dfac_x*mcByq+0.75*Phi[1]*fl[14]*dfac_x*mcByq-0.75*Phi[3]*fr[10]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fr[10]*dfac_x*mcByq-0.75*Phi[3]*fl[10]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fl[10]*dfac_x*mcByq+0.8660254037844386*fr[14]*amax+0.8660254037844386*fl[14]*amax-0.5*fr[10]*amax+0.5*fl[10]*amax; 
  Ghat[13] = (-0.4330127018922193*BcurvY[0]*fr[15]*mcByq*wv2)+0.4330127018922193*BcurvY[0]*fl[15]*mcByq*wv2-0.4330127018922193*BcurvY[1]*fr[14]*mcByq*wv2+0.4330127018922193*BcurvY[1]*fl[14]*mcByq*wv2+0.25*BcurvY[0]*fr[13]*mcByq*wv2+0.25*BcurvY[0]*fl[13]*mcByq*wv2+0.25*BcurvY[1]*fr[10]*mcByq*wv2+0.25*BcurvY[1]*fl[10]*mcByq*wv2-(0.25*BcurvY[0]*fr[12]*mcByq*wv)/dfac_v+(0.25*BcurvY[0]*fl[12]*mcByq*wv)/dfac_v-(0.25*BcurvY[1]*fr[9]*mcByq*wv)/dfac_v+(0.25*BcurvY[1]*fl[9]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fr[8]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[0]*fl[8]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fr[4]*mcByq*wv)/dfac_v+(0.1443375672974065*BcurvY[1]*fl[4]*mcByq*wv)/dfac_v-0.75*Bmag[1]*fr[15]*dfac_x*mcByq*wm+0.75*Bmag[1]*fl[15]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fr[13]*dfac_x*mcByq*wm+0.4330127018922193*Bmag[1]*fl[13]*dfac_x*mcByq*wm-(0.4330127018922193*Bmag[1]*fr[11]*dfac_x*mcByq)/dfac_m+(0.4330127018922193*Bmag[1]*fl[11]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fr[6]*dfac_x*mcByq)/dfac_m+(0.25*Bmag[1]*fl[6]*dfac_x*mcByq)/dfac_m+1.299038105676658*Phi[3]*fr[15]*dfac_x*mcByq-0.75*Phi[1]*fr[15]*dfac_x*mcByq-1.299038105676658*Phi[3]*fl[15]*dfac_x*mcByq+0.75*Phi[1]*fl[15]*dfac_x*mcByq-0.75*Phi[3]*fr[13]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fr[13]*dfac_x*mcByq-0.75*Phi[3]*fl[13]*dfac_x*mcByq+0.4330127018922193*Phi[1]*fl[13]*dfac_x*mcByq+0.8660254037844386*fr[15]*amax+0.8660254037844386*fl[15]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_y; 
  incr[1] = 0.5*Ghat[1]*dfac_y; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_y; 
  incr[3] = 0.5*Ghat[3]*dfac_y; 
  incr[4] = 0.5*Ghat[4]*dfac_y; 
  incr[5] = -0.8660254037844386*Ghat[1]*dfac_y; 
  incr[6] = 0.5*Ghat[6]*dfac_y; 
  incr[7] = -0.8660254037844386*Ghat[3]*dfac_y; 
  incr[8] = 0.5*Ghat[8]*dfac_y; 
  incr[9] = -0.8660254037844386*Ghat[4]*dfac_y; 
  incr[10] = 0.5*Ghat[10]*dfac_y; 
  incr[11] = -0.8660254037844386*Ghat[6]*dfac_y; 
  incr[12] = -0.8660254037844386*Ghat[8]*dfac_y; 
  incr[13] = 0.5*Ghat[13]*dfac_y; 
  incr[14] = -0.8660254037844386*Ghat[10]*dfac_y; 
  incr[15] = -0.8660254037844386*Ghat[13]*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
return alpha0; 
} 
double GyrokineticSurf2x2vSer_Vpar_P1(const double mcByq, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double Ghat[16]; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = std::abs(-(0.5*((1.732050807568877*BcurvY[1]*Phi[3]+1.732050807568877*BcurvY[0]*Phi[2])*dfac_v*dfac_y*mcByq*wv+((-1.732050807568877*BcurvY[1]*Phi[3])-1.732050807568877*BcurvY[0]*Phi[2])*dfac_y*mcByq))/dfac_v); 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = alpha0; 
  else 
    amax = amax_in; 

  Ghat[0] = 0.375*BcurvY[0]*Phi[3]*fr[6]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[2]*fr[6]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[3]*fl[6]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[2]*fl[6]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[3]*fr[3]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[2]*fr[3]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[3]*fl[3]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[2]*fl[3]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*fr[1]*Phi[3]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*fl[1]*Phi[3]*dfac_y*mcByq*wv-0.2165063509461096*fr[0]*BcurvY[1]*Phi[3]*dfac_y*mcByq*wv-0.2165063509461096*fl[0]*BcurvY[1]*Phi[3]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*fr[1]*Phi[2]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*fl[1]*Phi[2]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*fr[0]*Phi[2]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*fl[0]*Phi[2]*dfac_y*mcByq*wv-(0.375*BcurvY[0]*Phi[3]*fr[6]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[2]*fr[6]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[3]*fl[6]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[2]*fl[6]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[3]*fr[3]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[2]*fr[3]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[3]*fl[3]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[2]*fl[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*fr[1]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*fl[1]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*fr[0]*BcurvY[1]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*fl[0]*BcurvY[1]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*fr[1]*Phi[2]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*fl[1]*Phi[2]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*fr[0]*Phi[2]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*fl[0]*Phi[2]*dfac_y*mcByq)/dfac_v+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = 0.675*BcurvY[1]*Phi[3]*fr[6]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[2]*fr[6]*dfac_y*mcByq*wv-0.675*BcurvY[1]*Phi[3]*fl[6]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[2]*fl[6]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[3]*fr[3]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[2]*fr[3]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[3]*fl[3]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[2]*fl[3]*dfac_y*mcByq*wv-0.3897114317029973*BcurvY[1]*fr[1]*Phi[3]*dfac_y*mcByq*wv-0.3897114317029973*BcurvY[1]*fl[1]*Phi[3]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*fr[0]*Phi[3]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*fl[0]*Phi[3]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*fr[1]*Phi[2]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*fl[1]*Phi[2]*dfac_y*mcByq*wv-0.2165063509461096*fr[0]*BcurvY[1]*Phi[2]*dfac_y*mcByq*wv-0.2165063509461096*fl[0]*BcurvY[1]*Phi[2]*dfac_y*mcByq*wv-(0.675*BcurvY[1]*Phi[3]*fr[6]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[2]*fr[6]*dfac_y*mcByq)/dfac_v+(0.675*BcurvY[1]*Phi[3]*fl[6]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[2]*fl[6]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[3]*fr[3]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[2]*fr[3]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[3]*fl[3]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[2]*fl[3]*dfac_y*mcByq)/dfac_v+(0.3897114317029973*BcurvY[1]*fr[1]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.3897114317029973*BcurvY[1]*fl[1]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*fr[0]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*fl[0]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*fr[1]*Phi[2]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*fl[1]*Phi[2]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*fr[0]*BcurvY[1]*Phi[2]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*fl[0]*BcurvY[1]*Phi[2]*dfac_y*mcByq)/dfac_v+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[2] = 0.375*BcurvY[0]*Phi[3]*fr[11]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[2]*fr[11]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[3]*fl[11]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[2]*fl[11]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[3]*fr[7]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[2]*fr[7]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[3]*fl[7]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[2]*fl[7]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[3]*fr[5]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fr[5]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[3]*fl[5]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fl[5]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*fr[2]*Phi[3]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*fl[2]*Phi[3]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fr[2]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fl[2]*dfac_y*mcByq*wv-(0.375*BcurvY[0]*Phi[3]*fr[11]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[2]*fr[11]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[3]*fl[11]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[2]*fl[11]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[3]*fr[7]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[2]*fr[7]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[3]*fl[7]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[2]*fl[7]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[3]*fr[5]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fr[5]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[3]*fl[5]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fl[5]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*fr[2]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*fl[2]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fr[2]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fl[2]*dfac_y*mcByq)/dfac_v+0.8660254037844386*fr[7]*amax+0.8660254037844386*fl[7]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[4] = 0.375*BcurvY[0]*Phi[3]*fr[13]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[2]*fr[13]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[3]*fl[13]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[2]*fl[13]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[3]*fr[10]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[2]*fr[10]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[3]*fl[10]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[2]*fl[10]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[3]*fr[8]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fr[8]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[3]*fl[8]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fl[8]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[3]*fr[4]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fr[4]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[3]*fl[4]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fl[4]*dfac_y*mcByq*wv-(0.375*BcurvY[0]*Phi[3]*fr[13]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[2]*fr[13]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[3]*fl[13]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[2]*fl[13]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[3]*fr[10]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[2]*fr[10]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[3]*fl[10]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[2]*fl[10]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[3]*fr[8]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fr[8]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[3]*fl[8]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fl[8]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[3]*fr[4]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fr[4]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[3]*fl[4]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fl[4]*dfac_y*mcByq)/dfac_v+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax; 
  Ghat[5] = 0.675*BcurvY[1]*Phi[3]*fr[11]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[2]*fr[11]*dfac_y*mcByq*wv-0.675*BcurvY[1]*Phi[3]*fl[11]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[2]*fl[11]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[3]*fr[7]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[2]*fr[7]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[3]*fl[7]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[2]*fl[7]*dfac_y*mcByq*wv-0.3897114317029973*BcurvY[1]*Phi[3]*fr[5]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fr[5]*dfac_y*mcByq*wv-0.3897114317029973*BcurvY[1]*Phi[3]*fl[5]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fl[5]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*fr[2]*Phi[3]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*fl[2]*Phi[3]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fr[2]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fl[2]*dfac_y*mcByq*wv-(0.675*BcurvY[1]*Phi[3]*fr[11]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[2]*fr[11]*dfac_y*mcByq)/dfac_v+(0.675*BcurvY[1]*Phi[3]*fl[11]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[2]*fl[11]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[3]*fr[7]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[2]*fr[7]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[3]*fl[7]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[2]*fl[7]*dfac_y*mcByq)/dfac_v+(0.3897114317029973*BcurvY[1]*Phi[3]*fr[5]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fr[5]*dfac_y*mcByq)/dfac_v+(0.3897114317029973*BcurvY[1]*Phi[3]*fl[5]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fl[5]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*fr[2]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*fl[2]*Phi[3]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fr[2]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fl[2]*dfac_y*mcByq)/dfac_v+0.8660254037844386*fr[11]*amax+0.8660254037844386*fl[11]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 
  Ghat[8] = 0.675*BcurvY[1]*Phi[3]*fr[13]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[2]*fr[13]*dfac_y*mcByq*wv-0.675*BcurvY[1]*Phi[3]*fl[13]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[2]*fl[13]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[3]*fr[10]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[2]*fr[10]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[3]*fl[10]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[2]*fl[10]*dfac_y*mcByq*wv-0.3897114317029973*BcurvY[1]*Phi[3]*fr[8]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fr[8]*dfac_y*mcByq*wv-0.3897114317029973*BcurvY[1]*Phi[3]*fl[8]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fl[8]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[3]*fr[4]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fr[4]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[3]*fl[4]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fl[4]*dfac_y*mcByq*wv-(0.675*BcurvY[1]*Phi[3]*fr[13]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[2]*fr[13]*dfac_y*mcByq)/dfac_v+(0.675*BcurvY[1]*Phi[3]*fl[13]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[2]*fl[13]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[3]*fr[10]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[2]*fr[10]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[3]*fl[10]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[2]*fl[10]*dfac_y*mcByq)/dfac_v+(0.3897114317029973*BcurvY[1]*Phi[3]*fr[8]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fr[8]*dfac_y*mcByq)/dfac_v+(0.3897114317029973*BcurvY[1]*Phi[3]*fl[8]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fl[8]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[3]*fr[4]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fr[4]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[3]*fl[4]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fl[4]*dfac_y*mcByq)/dfac_v+0.8660254037844386*fr[13]*amax+0.8660254037844386*fl[13]*amax-0.5*fr[8]*amax+0.5*fl[8]*amax; 
  Ghat[9] = 0.375*BcurvY[0]*Phi[3]*fr[15]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[2]*fr[15]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[3]*fl[15]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[2]*fl[15]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[3]*fr[14]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[2]*fr[14]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[3]*fl[14]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[2]*fl[14]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[3]*fr[12]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fr[12]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[3]*fl[12]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fl[12]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[3]*fr[9]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fr[9]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[3]*fl[9]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fl[9]*dfac_y*mcByq*wv-(0.375*BcurvY[0]*Phi[3]*fr[15]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[2]*fr[15]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[3]*fl[15]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[2]*fl[15]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[3]*fr[14]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[2]*fr[14]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[3]*fl[14]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[2]*fl[14]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[3]*fr[12]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fr[12]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[3]*fl[12]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fl[12]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[3]*fr[9]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fr[9]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[3]*fl[9]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fl[9]*dfac_y*mcByq)/dfac_v+0.8660254037844386*fr[14]*amax+0.8660254037844386*fl[14]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[12] = 0.675*BcurvY[1]*Phi[3]*fr[15]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[2]*fr[15]*dfac_y*mcByq*wv-0.675*BcurvY[1]*Phi[3]*fl[15]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[2]*fl[15]*dfac_y*mcByq*wv+0.375*BcurvY[0]*Phi[3]*fr[14]*dfac_y*mcByq*wv+0.375*BcurvY[1]*Phi[2]*fr[14]*dfac_y*mcByq*wv-0.375*BcurvY[0]*Phi[3]*fl[14]*dfac_y*mcByq*wv-0.375*BcurvY[1]*Phi[2]*fl[14]*dfac_y*mcByq*wv-0.3897114317029973*BcurvY[1]*Phi[3]*fr[12]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fr[12]*dfac_y*mcByq*wv-0.3897114317029973*BcurvY[1]*Phi[3]*fl[12]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[2]*fl[12]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[3]*fr[9]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fr[9]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[0]*Phi[3]*fl[9]*dfac_y*mcByq*wv-0.2165063509461096*BcurvY[1]*Phi[2]*fl[9]*dfac_y*mcByq*wv-(0.675*BcurvY[1]*Phi[3]*fr[15]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[2]*fr[15]*dfac_y*mcByq)/dfac_v+(0.675*BcurvY[1]*Phi[3]*fl[15]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[2]*fl[15]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[0]*Phi[3]*fr[14]*dfac_y*mcByq)/dfac_v-(0.375*BcurvY[1]*Phi[2]*fr[14]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[0]*Phi[3]*fl[14]*dfac_y*mcByq)/dfac_v+(0.375*BcurvY[1]*Phi[2]*fl[14]*dfac_y*mcByq)/dfac_v+(0.3897114317029973*BcurvY[1]*Phi[3]*fr[12]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fr[12]*dfac_y*mcByq)/dfac_v+(0.3897114317029973*BcurvY[1]*Phi[3]*fl[12]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[2]*fl[12]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[3]*fr[9]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fr[9]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[0]*Phi[3]*fl[9]*dfac_y*mcByq)/dfac_v+(0.2165063509461096*BcurvY[1]*Phi[2]*fl[9]*dfac_y*mcByq)/dfac_v+0.8660254037844386*fr[15]*amax+0.8660254037844386*fl[15]*amax-0.5*fr[12]*amax+0.5*fl[12]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = 0.5*Ghat[2]*dfac_v; 
  incr[3] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[4] = 0.5*Ghat[4]*dfac_v; 
  incr[5] = 0.5*Ghat[5]*dfac_v; 
  incr[6] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[7] = -0.8660254037844386*Ghat[2]*dfac_v; 
  incr[8] = 0.5*Ghat[8]*dfac_v; 
  incr[9] = 0.5*Ghat[9]*dfac_v; 
  incr[10] = -0.8660254037844386*Ghat[4]*dfac_v; 
  incr[11] = -0.8660254037844386*Ghat[5]*dfac_v; 
  incr[12] = 0.5*Ghat[12]*dfac_v; 
  incr[13] = -0.8660254037844386*Ghat[8]*dfac_v; 
  incr[14] = -0.8660254037844386*Ghat[9]*dfac_v; 
  incr[15] = -0.8660254037844386*Ghat[12]*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
return alpha0; 
} 
