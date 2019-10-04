#include <IncompEulerModDecl.h> 
double IncompEulerSurf2xSer_X_P2(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // H/f: Input Hamiltonian/distribution function.
  // out: Incremented output.

  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[8]; 
  // Surface-averaged phase velocity in this direction.
  double alpha0 = 0.25*(3.872983346207417*Phi[6]-3.0*Phi[3]+1.732050807568877*Phi[2])*dfac_y; 

  double alpha[3]; 
  alpha[0] = 0.7071067811865475*(3.872983346207417*Phi[6]-3.0*Phi[3]+1.732050807568877*Phi[2])*dfac_y; 
  alpha[1] = -0.7071067811865475*(6.708203932499369*Phi[7]-3.872983346207417*Phi[5])*dfac_y; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[3]; 
  double favg[3]; 
  favg[0] = 0.7071067811865475*(2.23606797749979*(fr[4]+fl[4])+1.732050807568877*(fl[1]-1.0*fr[1])+fr[0]+fl[0]); 
  favg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[6]+fl[6])+3.0*(fl[3]-1.0*fr[3]))+3.0*(fr[2]+fl[2])); 
  favg[2] = -0.1414213562373095*(8.660254037844387*fr[7]-1.0*(8.660254037844387*fl[7]+5.0*(fr[5]+fl[5]))); 

  Ghat[0] = -0.25*((3.16227766016838*fr[4]-3.16227766016838*fl[4]-2.449489742783178*(fr[1]+fl[1])+1.414213562373095*fr[0]-1.414213562373095*fl[0])*amax-1.414213562373095*(alpha[1]*favg[1]+alpha[0]*favg[0])); 
  Ghat[1] = -0.01666666666666667*((47.43416490252569*fr[6]-47.43416490252569*fl[6]-36.74234614174767*(fr[3]+fl[3])+21.21320343559643*fr[2]-21.21320343559643*fl[2])*amax-18.97366596101028*alpha[1]*favg[2]-21.21320343559643*(alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.05*((12.24744871391589*(fr[7]+fl[7])-7.071067811865476*fr[5]+7.071067811865476*fl[5])*amax+7.071067811865476*alpha[0]*favg[2]+6.324555320336761*alpha[1]*favg[1]); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_x; 
  incr[1] = -1.224744871391589*Ghat[0]*dfac_x; 
  incr[2] = 0.7071067811865475*Ghat[1]*dfac_x; 
  incr[3] = -1.224744871391589*Ghat[1]*dfac_x; 
  incr[4] = 1.58113883008419*Ghat[0]*dfac_x; 
  incr[5] = 0.7071067811865475*Ghat[2]*dfac_x; 
  incr[6] = 1.58113883008419*Ghat[1]*dfac_x; 
  incr[7] = -1.224744871391589*Ghat[2]*dfac_x; 

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
double IncompEulerSurf2xSer_Y_P2(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // H/f: Input Hamiltonian/distribution function.
  // out: Incremented output.

  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[8]; 
  // Surface-averaged phase velocity in this direction.
  double alpha0 = -0.25*(3.872983346207417*Phi[7]-3.0*Phi[3]+1.732050807568877*Phi[1])*dfac_x; 

  double alpha[3]; 
  alpha[0] = -0.7071067811865475*(3.872983346207417*Phi[7]-3.0*Phi[3]+1.732050807568877*Phi[1])*dfac_x; 
  alpha[1] = 0.7071067811865475*(6.708203932499369*Phi[6]-3.872983346207417*Phi[4])*dfac_x; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[3]; 
  double favg[3]; 
  favg[0] = 0.7071067811865475*(2.23606797749979*(fr[5]+fl[5])+1.732050807568877*(fl[2]-1.0*fr[2])+fr[0]+fl[0]); 
  favg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[7]+fl[7])+3.0*(fl[3]-1.0*fr[3]))+3.0*(fr[1]+fl[1])); 
  favg[2] = -0.1414213562373095*(8.660254037844387*fr[6]-1.0*(8.660254037844387*fl[6]+5.0*(fr[4]+fl[4]))); 

  Ghat[0] = -0.25*((3.16227766016838*fr[5]-3.16227766016838*fl[5]-2.449489742783178*(fr[2]+fl[2])+1.414213562373095*fr[0]-1.414213562373095*fl[0])*amax-1.414213562373095*(alpha[1]*favg[1]+alpha[0]*favg[0])); 
  Ghat[1] = -0.01666666666666667*((47.43416490252569*fr[7]-47.43416490252569*fl[7]-36.74234614174767*(fr[3]+fl[3])+21.21320343559643*fr[1]-21.21320343559643*fl[1])*amax-18.97366596101028*alpha[1]*favg[2]-21.21320343559643*(alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.05*((12.24744871391589*(fr[6]+fl[6])-7.071067811865476*fr[4]+7.071067811865476*fl[4])*amax+7.071067811865476*alpha[0]*favg[2]+6.324555320336761*alpha[1]*favg[1]); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_y; 
  incr[1] = 0.7071067811865475*Ghat[1]*dfac_y; 
  incr[2] = -1.224744871391589*Ghat[0]*dfac_y; 
  incr[3] = -1.224744871391589*Ghat[1]*dfac_y; 
  incr[4] = 0.7071067811865475*Ghat[2]*dfac_y; 
  incr[5] = 1.58113883008419*Ghat[0]*dfac_y; 
  incr[6] = -1.224744871391589*Ghat[2]*dfac_y; 
  incr[7] = 1.58113883008419*Ghat[1]*dfac_y; 

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
