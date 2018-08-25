#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x1vSer_X_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.5*wv; 

  double alpha[3]; 
  alpha[0] = 1.414213562373095*wv; 
  alpha[1] = 0.8164965809277261/dfac_v; 
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
  Ghat[0] = 0.3535533905932737*(alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[4]-1.0*(2.23606797749979*fl[4]+1.732050807568877*(fr[1]+fl[1]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.07071067811865474*(4.47213595499958*alpha[1]*favg[2]+5.0*(alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[6]-1.0*(3.872983346207417*fl[6]+3.0*(fr[3]+fl[3])))+3.0*(fr[2]-1.0*fl[2]))*amax; 
  Ghat[2] = 0.07071067811865474*(8.660254037844387*(fr[7]+fl[7])+5.0*(fl[5]-1.0*fr[5]))*amax+0.07071067811865474*(5.0*alpha[0]*favg[2]+4.47213595499958*alpha[1]*favg[1]); 
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
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.6123724356957944*Phi[1]*dfac_x*q_)/m_; 

  double alpha[3]; 
  alpha[0] = -(1.732050807568877*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(3.872983346207417*Phi[2]*dfac_x*q_)/m_; 
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
  Ghat[0] = 0.3535533905932737*(alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[5]-1.0*(2.23606797749979*fl[5]+1.732050807568877*(fr[2]+fl[2]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.07071067811865474*(4.47213595499958*alpha[1]*favg[2]+5.0*(alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[7]-1.0*(3.872983346207417*fl[7]+3.0*(fr[3]+fl[3])))+3.0*(fr[1]-1.0*fl[1]))*amax; 
  Ghat[2] = 0.07071067811865474*(8.660254037844387*(fr[6]+fl[6])+5.0*(fl[4]-1.0*fr[4]))*amax+0.07071067811865474*(5.0*alpha[0]*favg[2]+4.47213595499958*alpha[1]*favg[1]); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_v; 
  incr[1] = 0.7071067811865475*Ghat[1]*dfac_v; 
  incr[2] = -1.224744871391589*Ghat[0]*dfac_v; 
  incr[3] = -1.224744871391589*Ghat[1]*dfac_v; 
  incr[4] = 0.7071067811865475*Ghat[2]*dfac_v; 
  incr[5] = 1.58113883008419*Ghat[0]*dfac_v; 
  incr[6] = -1.224744871391589*Ghat[2]*dfac_v; 
  incr[7] = 1.58113883008419*Ghat[1]*dfac_v; 

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
double GyrokineticSurf1x1vSer_X_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.5*wv; 

  double alpha[3]; 
  alpha[0] = 1.414213562373095*wv; 
  alpha[1] = 0.8164965809277261/dfac_v; 
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
  Ghat[0] = 0.3535533905932737*(alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[4]-1.0*(2.23606797749979*fl[4]+1.732050807568877*(fr[1]+fl[1]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.07071067811865474*(4.47213595499958*alpha[1]*favg[2]+5.0*(alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[6]-1.0*(3.872983346207417*fl[6]+3.0*(fr[3]+fl[3])))+3.0*(fr[2]-1.0*fl[2]))*amax; 
  Ghat[2] = 0.07071067811865474*(8.660254037844387*(fr[7]+fl[7])+5.0*(fl[5]-1.0*fr[5]))*amax+0.07071067811865474*(5.0*alpha[0]*favg[2]+4.47213595499958*alpha[1]*favg[1]); 
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
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.6123724356957944*Phi[1]*dfac_x*q_)/m_; 

  double alpha[3]; 
  alpha[0] = -(1.732050807568877*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(3.872983346207417*Phi[2]*dfac_x*q_)/m_; 
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
  Ghat[0] = 0.3535533905932737*(alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[5]-1.0*(2.23606797749979*fl[5]+1.732050807568877*(fr[2]+fl[2]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.07071067811865474*(4.47213595499958*alpha[1]*favg[2]+5.0*(alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[7]-1.0*(3.872983346207417*fl[7]+3.0*(fr[3]+fl[3])))+3.0*(fr[1]-1.0*fl[1]))*amax; 
  Ghat[2] = 0.07071067811865474*(8.660254037844387*(fr[6]+fl[6])+5.0*(fl[4]-1.0*fr[4]))*amax+0.07071067811865474*(5.0*alpha[0]*favg[2]+4.47213595499958*alpha[1]*favg[1]); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_v; 
  incr[1] = 0.7071067811865475*Ghat[1]*dfac_v; 
  incr[2] = -1.224744871391589*Ghat[0]*dfac_v; 
  incr[3] = -1.224744871391589*Ghat[1]*dfac_v; 
  incr[4] = 0.7071067811865475*Ghat[2]*dfac_v; 
  incr[5] = 1.58113883008419*Ghat[0]*dfac_v; 
  incr[6] = -1.224744871391589*Ghat[2]*dfac_v; 
  incr[7] = 1.58113883008419*Ghat[1]*dfac_v; 

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
