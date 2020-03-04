#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity2x0vSer_X_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *cflRateByDirL, const double *cflRateByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.125*BmagInv[0]*(1.732050807568877*Phi[2]-3.0*Phi[3])*dfac_y; 

  double alpha[2]; 
  alpha[0] = -0.3535533905932737*BmagInv[0]*(1.732050807568877*Phi[2]-3.0*Phi[3])*dfac_y; 
  double rCtrlL[2], rCtrlR[2];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrlL[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(1.732050807568877*(2.0*EPSILON+fl[0])-1.0*fl[2]); 
  rCtrlL[1] = (1.732050807568877*fl[3]+3.0*fl[1])/(1.732050807568877*(2.0*EPSILON+fl[0])+fl[2]); 
  rCtrlR[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(1.732050807568877*(2.0*EPSILON+fr[0])-1.0*fr[2]); 
  rCtrlR[1] = (1.732050807568877*fr[3]+3.0*fr[1])/(1.732050807568877*(2.0*EPSILON+fr[0])+fr[2]); 
  double fCtrlL[2], fCtrlR[2];  // fCtrl = anti-limited f evaluated at each control node on x surface 
  // control node [y] = [-1/3] 
  fCtrlL[0] = 0.2886751345948129*(1.732050807568877*fl[0]-1.0*fl[2])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = 0.2886751345948129*(1.732050807568877*fr[0]-1.0*fr[2])*limTheta(rCtrlR[0],-1.0); 
  // control node [y] = [1/3] 
  fCtrlL[1] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rCtrlR[1],-1.0); 
  double fL_AL[2], fR_AL[2];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.7071067811865475*(fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 1.224744871391589*(fCtrlL[1]-1.0*fCtrlL[0]); 
  fR_AL[0] = 0.7071067811865475*(fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 1.224744871391589*(fCtrlR[1]-1.0*fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[2], fhatAL[2]; 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = -0.5*(1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[0] = -0.5*(1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = 0.7071067811865476*(fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[1] = 0.7071067811865476*(fR_AL[1]+fR_AL[0]); 
  } 
  fhatAL[0] = 0.7071067811865475*(fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.7071067811865475*(fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double fluxFracL, fluxFracR, limfac, flim = 0.;
  double GhatCtrl[2];
  fluxFracL = cflRateByDirL[0] != -10. ? 0.5 : cflRateByDirL[1]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] != -10. ? 0.5 : cflRateByDirR[1]/cflRateByDirR[0]; 
  // control node [y] = [-1/3] 
  GhatCtrl[0] = alpha[0]*(0.5*fhatAL[0]-0.2886751345948129*fhatAL[1]); 
  if(GhatCtrl[0] > EPSILON) {
    flim = std::max(0., -0.1666666666666667*(fl[3]+1.732050807568877*(fl[2]-1.0*fl[1])-3.0*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[0]/dtApprox/dfac_x)); 
  } else if(GhatCtrl[0] < -EPSILON) {
    flim = std::max(0., -0.1666666666666667*((-1.0*fr[3])+1.732050807568877*(fr[2]+fr[1])-3.0*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[0]/dtApprox/dfac_x)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[0] *= limfac; 
  }
  // control node [y] = [1/3] 
  GhatCtrl[1] = alpha[0]*(0.2886751345948129*fhatAL[1]+0.5*fhatAL[0]); 
  if(GhatCtrl[1] > EPSILON) {
    flim = std::max(0., 0.1666666666666667*(fl[3]+1.732050807568877*(fl[2]+fl[1])+3.0*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[1]/dtApprox/dfac_x)); 
  } else if(GhatCtrl[1] < -EPSILON) {
    flim = std::max(0., 0.1666666666666667*((-1.0*fr[3])+1.732050807568877*(fr[2]-1.0*fr[1])+3.0*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[1]/dtApprox/dfac_x)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[1] *= limfac; 
  }

  incr[0] = 0.5*(GhatCtrl[1]+GhatCtrl[0])*dfac_x; 
  incr[1] = -0.8660254037844386*(GhatCtrl[1]+GhatCtrl[0])*dfac_x; 
  incr[2] = 0.8660254037844386*(GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_x; 
  incr[3] = -1.5*(GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x0vSer_Y_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *cflRateByDirL, const double *cflRateByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.125*BmagInv[0]*(1.732050807568877*Phi[1]-3.0*Phi[3])*dfac_x; 

  double alpha[2]; 
  alpha[0] = 0.3535533905932737*BmagInv[0]*(1.732050807568877*Phi[1]-3.0*Phi[3])*dfac_x; 
  double rCtrlL[2], rCtrlR[2];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrlL[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(1.732050807568877*(2.0*EPSILON+fl[0])-1.0*fl[1]); 
  rCtrlL[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(1.732050807568877*(2.0*EPSILON+fl[0])+fl[1]); 
  rCtrlR[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(1.732050807568877*(2.0*EPSILON+fr[0])-1.0*fr[1]); 
  rCtrlR[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(1.732050807568877*(2.0*EPSILON+fr[0])+fr[1]); 
  double fCtrlL[2], fCtrlR[2];  // fCtrl = anti-limited f evaluated at each control node on y surface 
  // control node [x] = [-1/3] 
  fCtrlL[0] = 0.2886751345948129*(1.732050807568877*fl[0]-1.0*fl[1])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = 0.2886751345948129*(1.732050807568877*fr[0]-1.0*fr[1])*limTheta(rCtrlR[0],-1.0); 
  // control node [x] = [1/3] 
  fCtrlL[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rCtrlR[1],-1.0); 
  double fL_AL[2], fR_AL[2];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.7071067811865475*(fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 1.224744871391589*(fCtrlL[1]-1.0*fCtrlL[0]); 
  fR_AL[0] = 0.7071067811865475*(fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 1.224744871391589*(fCtrlR[1]-1.0*fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[2], fhatAL[2]; 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = -0.5*(1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[0] = -0.5*(1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = 0.7071067811865476*(fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[1] = 0.7071067811865476*(fR_AL[1]+fR_AL[0]); 
  } 
  fhatAL[0] = 0.7071067811865475*(fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.7071067811865475*(fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double fluxFracL, fluxFracR, limfac, flim = 0.;
  double GhatCtrl[2];
  fluxFracL = cflRateByDirL[0] != -10. ? 0.5 : cflRateByDirL[2]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] != -10. ? 0.5 : cflRateByDirR[2]/cflRateByDirR[0]; 
  // control node [x] = [-1/3] 
  GhatCtrl[0] = alpha[0]*(0.5*fhatAL[0]-0.2886751345948129*fhatAL[1]); 
  if(GhatCtrl[0] > EPSILON) {
    flim = std::max(0., 0.1666666666666667*((-1.0*fl[3])+1.732050807568877*(fl[2]-1.0*fl[1])+3.0*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[0]/dtApprox/dfac_y)); 
  } else if(GhatCtrl[0] < -EPSILON) {
    flim = std::max(0., -0.1666666666666667*((-1.0*fr[3])+1.732050807568877*(fr[2]+fr[1])-3.0*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[0]/dtApprox/dfac_y)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[0] *= limfac; 
  }
  // control node [x] = [1/3] 
  GhatCtrl[1] = alpha[0]*(0.2886751345948129*fhatAL[1]+0.5*fhatAL[0]); 
  if(GhatCtrl[1] > EPSILON) {
    flim = std::max(0., 0.1666666666666667*(fl[3]+1.732050807568877*(fl[2]+fl[1])+3.0*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[1]/dtApprox/dfac_y)); 
  } else if(GhatCtrl[1] < -EPSILON) {
    flim = std::max(0., -0.1666666666666667*(fr[3]+1.732050807568877*(fr[2]-1.0*fr[1])-3.0*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[1]/dtApprox/dfac_y)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[1] *= limfac; 
  }

  incr[0] = 0.5*(GhatCtrl[1]+GhatCtrl[0])*dfac_y; 
  incr[1] = 0.8660254037844386*(GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_y; 
  incr[2] = -0.8660254037844386*(GhatCtrl[1]+GhatCtrl[0])*dfac_y; 
  incr[3] = -1.5*(GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  return std::abs(alpha0); 
} 
