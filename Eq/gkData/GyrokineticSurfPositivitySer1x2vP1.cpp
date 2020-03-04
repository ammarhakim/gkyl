#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *cflRateByDirL, const double *cflRateByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*Gradpar[0]*wv; 

  double alpha[4]; 
  alpha[0] = 1.414213562373095*Gradpar[0]*wv; 
  double rCtrlL[4], rCtrlR[4];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrlL[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rCtrlL[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rCtrlL[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rCtrlL[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rCtrlR[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rCtrlR[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rCtrlR[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rCtrlR[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  double fCtrlL[4], fCtrlR[4];  // fCtrl = anti-limited f evaluated at each control node on x surface 
  // control node [vx,vy] = [-1/3,-1/3] 
  fCtrlL[0] = 0.04811252243246882*(2.449489742783178*fl[6]-4.242640687119286*(fl[3]+fl[2])+7.348469228349534*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = 0.04811252243246882*(2.449489742783178*fr[6]-4.242640687119286*(fr[3]+fr[2])+7.348469228349534*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [vx,vy] = [1/3,-1/3] 
  fCtrlL[1] = -0.04811252243246882*(2.449489742783178*fl[6]+4.242640687119286*fl[3]-4.242640687119286*fl[2]-7.348469228349534*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = -0.04811252243246882*(2.449489742783178*fr[6]+4.242640687119286*fr[3]-4.242640687119286*fr[2]-7.348469228349534*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [vx,vy] = [-1/3,1/3] 
  fCtrlL[2] = -0.04811252243246882*(2.449489742783178*fl[6]-4.242640687119286*fl[3]+4.242640687119286*fl[2]-7.348469228349534*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = -0.04811252243246882*(2.449489742783178*fr[6]-4.242640687119286*fr[3]+4.242640687119286*fr[2]-7.348469228349534*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [vx,vy] = [1/3,1/3] 
  fCtrlL[3] = 0.04811252243246882*(2.449489742783178*fl[6]+4.242640687119286*(fl[3]+fl[2])+7.348469228349534*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = 0.04811252243246882*(2.449489742783178*fr[6]+4.242640687119286*(fr[3]+fr[2])+7.348469228349534*fr[0])*limTheta(rCtrlR[3],-1.0); 
  double fL_AL[4], fR_AL[4];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.5*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.8660254037844386*(fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.8660254037844386*(fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 1.5*(fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fR_AL[0] = 0.5*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.8660254037844386*(fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.8660254037844386*(fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 1.5*(fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[4], fhatAL[4]; 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = 0.5*(fL_AL[3]-1.0*(fL_AL[2]+fL_AL[1])+fL_AL[0]); 
  } else {
  fhatALQuad[0] = 0.5*(fR_AL[3]-1.0*(fR_AL[2]+fR_AL[1])+fR_AL[0]); 
  } 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = -0.5*(fL_AL[3]+fL_AL[2]-1.0*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = -0.5*(fR_AL[3]+fR_AL[2]-1.0*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = -0.5*(fL_AL[3]-1.0*fL_AL[2]+fL_AL[1]-1.0*fL_AL[0]); 
  } else {
  fhatALQuad[2] = -0.5*(fR_AL[3]-1.0*fR_AL[2]+fR_AL[1]-1.0*fR_AL[0]); 
  } 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[3] = 0.5*(fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[3] = 0.5*(fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
  } 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  double fluxFracL, fluxFracR, limfac, flim = 0.;
  double GhatCtrl[4];
  fluxFracL = cflRateByDirL[0] != -10. ? 0.3333333333333333 : cflRateByDirL[1]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] != -10. ? 0.3333333333333333 : cflRateByDirR[1]/cflRateByDirR[0]; 
  // control node [vpar,mu] = [-1/3,-1/3] 
  GhatCtrl[0] = alpha[0]*(0.08333333333333333*fhatAL[3]-0.1443375672974065*(fhatAL[2]+fhatAL[1])+0.25*fhatAL[0]); 
  if(GhatCtrl[0] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]+12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[0]/dtApprox/dfac_x)); 
  } else if(GhatCtrl[0] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[0]/dtApprox/dfac_x)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[0] *= limfac; 
  }
  // control node [vpar,mu] = [1/3,-1/3] 
  GhatCtrl[1] = alpha[0]*((-0.08333333333333333*fhatAL[3])-0.1443375672974065*fhatAL[2]+0.1443375672974065*fhatAL[1]+0.25*fhatAL[0]); 
  if(GhatCtrl[1] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[1]/dtApprox/dfac_x)); 
  } else if(GhatCtrl[1] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*fr[5]-4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]+12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[1]/dtApprox/dfac_x)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[1] *= limfac; 
  }
  // control node [vpar,mu] = [-1/3,1/3] 
  GhatCtrl[2] = alpha[0]*((-0.08333333333333333*fhatAL[3])+0.1443375672974065*fhatAL[2]-0.1443375672974065*fhatAL[1]+0.25*fhatAL[0]); 
  if(GhatCtrl[2] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*fl[5]+4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]-12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[2]/dtApprox/dfac_x)); 
  } else if(GhatCtrl[2] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[2]/dtApprox/dfac_x)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[2] *= limfac; 
  }
  // control node [vpar,mu] = [1/3,1/3] 
  GhatCtrl[3] = alpha[0]*(0.08333333333333333*fhatAL[3]+0.1443375672974065*(fhatAL[2]+fhatAL[1])+0.25*fhatAL[0]); 
  if(GhatCtrl[3] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[3]/dtApprox/dfac_x)); 
  } else if(GhatCtrl[3] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]-12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[3]/dtApprox/dfac_x)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[3] *= limfac; 
  }

  incr[0] = 0.3535533905932737*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])*dfac_x; 
  incr[1] = -0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])*dfac_x; 
  incr[2] = 0.6123724356957944*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_x; 
  incr[3] = 0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0]))*dfac_x; 
  incr[4] = -1.060660171779821*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_x; 
  incr[5] = -1.060660171779821*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0]))*dfac_x; 
  incr[6] = 1.060660171779821*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0])*dfac_x; 
  incr[7] = -1.837117307087383*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0])*dfac_x; 

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
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *cflRateByDirL, const double *cflRateByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.4330127018922193*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[4]; 
  alpha[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  double rCtrlL[4], rCtrlR[4];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrlL[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rCtrlL[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlR[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rCtrlR[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  double fCtrlL[4], fCtrlR[4];  // fCtrl = anti-limited f evaluated at each control node on vx surface 
  // control node [x,vy] = [-1/3,-1/3] 
  fCtrlL[0] = 0.04811252243246882*(2.449489742783178*fl[5]-4.242640687119286*(fl[3]+fl[1])+7.348469228349534*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = 0.04811252243246882*(2.449489742783178*fr[5]-4.242640687119286*(fr[3]+fr[1])+7.348469228349534*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,vy] = [1/3,-1/3] 
  fCtrlL[1] = -0.04811252243246882*(2.449489742783178*fl[5]+4.242640687119286*fl[3]-4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = -0.04811252243246882*(2.449489742783178*fr[5]+4.242640687119286*fr[3]-4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,vy] = [-1/3,1/3] 
  fCtrlL[2] = -0.04811252243246882*(2.449489742783178*fl[5]-4.242640687119286*fl[3]+4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = -0.04811252243246882*(2.449489742783178*fr[5]-4.242640687119286*fr[3]+4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,vy] = [1/3,1/3] 
  fCtrlL[3] = 0.04811252243246882*(2.449489742783178*fl[5]+4.242640687119286*(fl[3]+fl[1])+7.348469228349534*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = 0.04811252243246882*(2.449489742783178*fr[5]+4.242640687119286*(fr[3]+fr[1])+7.348469228349534*fr[0])*limTheta(rCtrlR[3],-1.0); 
  double fL_AL[4], fR_AL[4];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.5*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.8660254037844386*(fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.8660254037844386*(fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 1.5*(fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fR_AL[0] = 0.5*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.8660254037844386*(fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.8660254037844386*(fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 1.5*(fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[4], fhatAL[4]; 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = 0.5*(fL_AL[3]-1.0*(fL_AL[2]+fL_AL[1])+fL_AL[0]); 
  } else {
  fhatALQuad[0] = 0.5*(fR_AL[3]-1.0*(fR_AL[2]+fR_AL[1])+fR_AL[0]); 
  } 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = -0.5*(fL_AL[3]+fL_AL[2]-1.0*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = -0.5*(fR_AL[3]+fR_AL[2]-1.0*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = -0.5*(fL_AL[3]-1.0*fL_AL[2]+fL_AL[1]-1.0*fL_AL[0]); 
  } else {
  fhatALQuad[2] = -0.5*(fR_AL[3]-1.0*fR_AL[2]+fR_AL[1]-1.0*fR_AL[0]); 
  } 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[3] = 0.5*(fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[3] = 0.5*(fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
  } 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  double fluxFracL, fluxFracR, limfac, flim = 0.;
  double GhatCtrl[4];
  fluxFracL = cflRateByDirL[0] != -10. ? 0.3333333333333333 : cflRateByDirL[2]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] != -10. ? 0.3333333333333333 : cflRateByDirR[2]/cflRateByDirR[0]; 
  // control node [x,mu] = [-1/3,-1/3] 
  GhatCtrl[0] = alpha[0]*(0.08333333333333333*fhatAL[3]-0.1443375672974065*(fhatAL[2]+fhatAL[1])+0.25*fhatAL[0]); 
  if(GhatCtrl[0] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*fl[5]-4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]+12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[0]/dtApprox/dfac_v)); 
  } else if(GhatCtrl[0] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[0]/dtApprox/dfac_v)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[0] *= limfac; 
  }
  // control node [x,mu] = [1/3,-1/3] 
  GhatCtrl[1] = alpha[0]*((-0.08333333333333333*fhatAL[3])-0.1443375672974065*fhatAL[2]+0.1443375672974065*fhatAL[1]+0.25*fhatAL[0]); 
  if(GhatCtrl[1] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[1]/dtApprox/dfac_v)); 
  } else if(GhatCtrl[1] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[1]/dtApprox/dfac_v)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[1] *= limfac; 
  }
  // control node [x,mu] = [-1/3,1/3] 
  GhatCtrl[2] = alpha[0]*((-0.08333333333333333*fhatAL[3])+0.1443375672974065*fhatAL[2]-0.1443375672974065*fhatAL[1]+0.25*fhatAL[0]); 
  if(GhatCtrl[2] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[2]/dtApprox/dfac_v)); 
  } else if(GhatCtrl[2] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[2]/dtApprox/dfac_v)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[2] *= limfac; 
  }
  // control node [x,mu] = [1/3,1/3] 
  GhatCtrl[3] = alpha[0]*(0.08333333333333333*fhatAL[3]+0.1443375672974065*(fhatAL[2]+fhatAL[1])+0.25*fhatAL[0]); 
  if(GhatCtrl[3] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[3]/dtApprox/dfac_v)); 
  } else if(GhatCtrl[3] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*fr[5]+4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]-12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[3]/dtApprox/dfac_v)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[3] *= limfac; 
  }

  incr[0] = 0.3535533905932737*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])*dfac_v; 
  incr[1] = 0.6123724356957944*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_v; 
  incr[2] = -0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])*dfac_v; 
  incr[3] = 0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0]))*dfac_v; 
  incr[4] = -1.060660171779821*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_v; 
  incr[5] = 1.060660171779821*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0])*dfac_v; 
  incr[6] = -1.060660171779821*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0]))*dfac_v; 
  incr[7] = -1.837117307087383*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0])*dfac_v; 

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
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *cflRateByDirL, const double *cflRateByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.3535533905932737*(1.732050807568877*Gradpar[1]-1.0*Gradpar[0])*wv; 

  double alpha[4]; 
  alpha[0] = (1.414213562373095*Gradpar[0]-2.449489742783178*Gradpar[1])*wv; 
  double rCtrlL[4], rCtrlR[4];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrlL[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rCtrlL[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rCtrlL[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rCtrlL[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rCtrlR[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rCtrlR[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rCtrlR[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rCtrlR[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  double fCtrlL[4], fCtrlR[4];  // fCtrl = anti-limited f evaluated at each control node on x surface 
  // control node [vx,vy] = [-1/3,-1/3] 
  fCtrlL[0] = 0.04811252243246882*(2.449489742783178*fl[6]-4.242640687119286*(fl[3]+fl[2])+7.348469228349534*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = 0.04811252243246882*(2.449489742783178*fr[6]-4.242640687119286*(fr[3]+fr[2])+7.348469228349534*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [vx,vy] = [1/3,-1/3] 
  fCtrlL[1] = -0.04811252243246882*(2.449489742783178*fl[6]+4.242640687119286*fl[3]-4.242640687119286*fl[2]-7.348469228349534*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = -0.04811252243246882*(2.449489742783178*fr[6]+4.242640687119286*fr[3]-4.242640687119286*fr[2]-7.348469228349534*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [vx,vy] = [-1/3,1/3] 
  fCtrlL[2] = -0.04811252243246882*(2.449489742783178*fl[6]-4.242640687119286*fl[3]+4.242640687119286*fl[2]-7.348469228349534*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = -0.04811252243246882*(2.449489742783178*fr[6]-4.242640687119286*fr[3]+4.242640687119286*fr[2]-7.348469228349534*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [vx,vy] = [1/3,1/3] 
  fCtrlL[3] = 0.04811252243246882*(2.449489742783178*fl[6]+4.242640687119286*(fl[3]+fl[2])+7.348469228349534*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = 0.04811252243246882*(2.449489742783178*fr[6]+4.242640687119286*(fr[3]+fr[2])+7.348469228349534*fr[0])*limTheta(rCtrlR[3],-1.0); 
  double fL_AL[4], fR_AL[4];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.5*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.8660254037844386*(fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.8660254037844386*(fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 1.5*(fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fR_AL[0] = 0.5*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.8660254037844386*(fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.8660254037844386*(fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 1.5*(fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[4], fhatAL[4]; 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = 0.5*(fL_AL[3]-1.0*(fL_AL[2]+fL_AL[1])+fL_AL[0]); 
  } else {
  fhatALQuad[0] = 0.5*(fR_AL[3]-1.0*(fR_AL[2]+fR_AL[1])+fR_AL[0]); 
  } 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = -0.5*(fL_AL[3]+fL_AL[2]-1.0*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = -0.5*(fR_AL[3]+fR_AL[2]-1.0*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = -0.5*(fL_AL[3]-1.0*fL_AL[2]+fL_AL[1]-1.0*fL_AL[0]); 
  } else {
  fhatALQuad[2] = -0.5*(fR_AL[3]-1.0*fR_AL[2]+fR_AL[1]-1.0*fR_AL[0]); 
  } 
  alphaQuad = 0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[3] = 0.5*(fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[3] = 0.5*(fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
  } 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  double fluxFracL, fluxFracR, limfac, flim = 0.;
  double GhatCtrl[4];
  fluxFracL = cflRateByDirL[0] != -10. ? 0.3333333333333333 : cflRateByDirL[1]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] != -10. ? 0.3333333333333333 : cflRateByDirR[1]/cflRateByDirR[0]; 
  // control node [vpar,mu] = [-1/3,-1/3] 
  GhatCtrl[0] = alpha[0]*(0.08333333333333333*fhatAL[3]-0.1443375672974065*(fhatAL[2]+fhatAL[1])+0.25*fhatAL[0]); 
  if(GhatCtrl[0] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]+12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[0]/dtApprox/dfac_x)); 
  } else if(GhatCtrl[0] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[0]/dtApprox/dfac_x)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[0] *= limfac; 
  }
  // control node [vpar,mu] = [1/3,-1/3] 
  GhatCtrl[1] = alpha[0]*((-0.08333333333333333*fhatAL[3])-0.1443375672974065*fhatAL[2]+0.1443375672974065*fhatAL[1]+0.25*fhatAL[0]); 
  if(GhatCtrl[1] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[1]/dtApprox/dfac_x)); 
  } else if(GhatCtrl[1] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*fr[5]-4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]+12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[1]/dtApprox/dfac_x)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[1] *= limfac; 
  }
  // control node [vpar,mu] = [-1/3,1/3] 
  GhatCtrl[2] = alpha[0]*((-0.08333333333333333*fhatAL[3])+0.1443375672974065*fhatAL[2]-0.1443375672974065*fhatAL[1]+0.25*fhatAL[0]); 
  if(GhatCtrl[2] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*fl[5]+4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]-12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[2]/dtApprox/dfac_x)); 
  } else if(GhatCtrl[2] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[2]/dtApprox/dfac_x)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[2] *= limfac; 
  }
  // control node [vpar,mu] = [1/3,1/3] 
  GhatCtrl[3] = alpha[0]*(0.08333333333333333*fhatAL[3]+0.1443375672974065*(fhatAL[2]+fhatAL[1])+0.25*fhatAL[0]); 
  if(GhatCtrl[3] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[3]/dtApprox/dfac_x)); 
  } else if(GhatCtrl[3] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]-12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[3]/dtApprox/dfac_x)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[3] *= limfac; 
  }

  incr[0] = 0.3535533905932737*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])*dfac_x; 
  incr[1] = -0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])*dfac_x; 
  incr[2] = 0.6123724356957944*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_x; 
  incr[3] = 0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0]))*dfac_x; 
  incr[4] = -1.060660171779821*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_x; 
  incr[5] = -1.060660171779821*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0]))*dfac_x; 
  incr[6] = 1.060660171779821*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0])*dfac_x; 
  incr[7] = -1.837117307087383*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0])*dfac_x; 

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
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *cflRateByDirL, const double *cflRateByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.4330127018922193*Gradpar[0]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))/m_; 

  double alpha[4]; 
  alpha[0] = -(1.732050807568877*Gradpar[0]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))/m_; 
  alpha[1] = -(1.732050807568877*Gradpar[1]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))/m_; 
  alpha[2] = -(1.0*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[3] = -(1.0*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
  double rCtrlL[4], rCtrlR[4];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrlL[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rCtrlL[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlR[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rCtrlR[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  double fCtrlL[4], fCtrlR[4];  // fCtrl = anti-limited f evaluated at each control node on vx surface 
  // control node [x,vy] = [-1/3,-1/3] 
  fCtrlL[0] = 0.04811252243246882*(2.449489742783178*fl[5]-4.242640687119286*(fl[3]+fl[1])+7.348469228349534*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = 0.04811252243246882*(2.449489742783178*fr[5]-4.242640687119286*(fr[3]+fr[1])+7.348469228349534*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,vy] = [1/3,-1/3] 
  fCtrlL[1] = -0.04811252243246882*(2.449489742783178*fl[5]+4.242640687119286*fl[3]-4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = -0.04811252243246882*(2.449489742783178*fr[5]+4.242640687119286*fr[3]-4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,vy] = [-1/3,1/3] 
  fCtrlL[2] = -0.04811252243246882*(2.449489742783178*fl[5]-4.242640687119286*fl[3]+4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = -0.04811252243246882*(2.449489742783178*fr[5]-4.242640687119286*fr[3]+4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,vy] = [1/3,1/3] 
  fCtrlL[3] = 0.04811252243246882*(2.449489742783178*fl[5]+4.242640687119286*(fl[3]+fl[1])+7.348469228349534*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = 0.04811252243246882*(2.449489742783178*fr[5]+4.242640687119286*(fr[3]+fr[1])+7.348469228349534*fr[0])*limTheta(rCtrlR[3],-1.0); 
  double fL_AL[4], fR_AL[4];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.5*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.8660254037844386*(fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.8660254037844386*(fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 1.5*(fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fR_AL[0] = 0.5*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.8660254037844386*(fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.8660254037844386*(fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 1.5*(fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[4], fhatAL[4]; 
  alphaQuad = 0.5*alpha[3]-0.5*(alpha[2]+alpha[1])+0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = 0.5*(fL_AL[3]-1.0*(fL_AL[2]+fL_AL[1])+fL_AL[0]); 
  } else {
  fhatALQuad[0] = 0.5*(fR_AL[3]-1.0*(fR_AL[2]+fR_AL[1])+fR_AL[0]); 
  } 
  alphaQuad = 0.5*(alpha[1]+alpha[0])-0.5*(alpha[3]+alpha[2]); 
  if(alphaQuad > 0) {
  fhatALQuad[1] = -0.5*(fL_AL[3]+fL_AL[2]-1.0*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = -0.5*(fR_AL[3]+fR_AL[2]-1.0*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = (-0.5*alpha[3])+0.5*alpha[2]-0.5*alpha[1]+0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = -0.5*(fL_AL[3]-1.0*fL_AL[2]+fL_AL[1]-1.0*fL_AL[0]); 
  } else {
  fhatALQuad[2] = -0.5*(fR_AL[3]-1.0*fR_AL[2]+fR_AL[1]-1.0*fR_AL[0]); 
  } 
  alphaQuad = 0.5*(alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = 0.5*(fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[3] = 0.5*(fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
  } 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  double fluxFracL, fluxFracR, limfac, flim = 0.;
  double GhatCtrl[4];
  fluxFracL = cflRateByDirL[0] != -10. ? 0.3333333333333333 : cflRateByDirL[2]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] != -10. ? 0.3333333333333333 : cflRateByDirR[2]/cflRateByDirR[0]; 
  // control node [x,mu] = [-1/3,-1/3] 
  GhatCtrl[0] = 0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.1443375672974065*((alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]); 
  if(GhatCtrl[0] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*fl[5]-4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]+12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[0]/dtApprox/dfac_v)); 
  } else if(GhatCtrl[0] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[0]/dtApprox/dfac_v)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[0] *= limfac; 
  }
  // control node [x,mu] = [1/3,-1/3] 
  GhatCtrl[1] = 0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.1443375672974065*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])-0.1443375672974065*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]); 
  if(GhatCtrl[1] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[1]/dtApprox/dfac_v)); 
  } else if(GhatCtrl[1] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[1]/dtApprox/dfac_v)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[1] *= limfac; 
  }
  // control node [x,mu] = [-1/3,1/3] 
  GhatCtrl[2] = 0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.1443375672974065*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.1443375672974065*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]); 
  if(GhatCtrl[2] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[2]/dtApprox/dfac_v)); 
  } else if(GhatCtrl[2] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[2]/dtApprox/dfac_v)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[2] *= limfac; 
  }
  // control node [x,mu] = [1/3,1/3] 
  GhatCtrl[3] = 0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.1443375672974065*((alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]); 
  if(GhatCtrl[3] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0])); 
    limfac = std::min(1., std::abs(fluxFracL*flim/GhatCtrl[3]/dtApprox/dfac_v)); 
  } else if(GhatCtrl[3] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*fr[5]+4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]-12.72792206135786*fr[0])); 
    limfac = std::min(1., std::abs(fluxFracR*flim/GhatCtrl[3]/dtApprox/dfac_v)); 
  } else limfac = 0.; 
  if(limfac != 1.) { 
    GhatCtrl[3] *= limfac; 
  }

  incr[0] = 0.3535533905932737*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])*dfac_v; 
  incr[1] = 0.6123724356957944*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_v; 
  incr[2] = -0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])*dfac_v; 
  incr[3] = 0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0]))*dfac_v; 
  incr[4] = -1.060660171779821*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0])*dfac_v; 
  incr[5] = 1.060660171779821*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0])*dfac_v; 
  incr[6] = -1.060660171779821*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0]))*dfac_v; 
  incr[7] = -1.837117307087383*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0])*dfac_v; 

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
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
