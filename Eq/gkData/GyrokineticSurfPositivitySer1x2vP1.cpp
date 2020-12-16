#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, 
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
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_x; 
  incr[4] = -0.25*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_x; 
  incr[5] = -0.25*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_x; 
  incr[6] = 0.25*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_x; 
  incr[7] = -0.25*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_x; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_x; 
  incr[4] = 0.25*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_x; 
  incr[5] = 0.25*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_x; 
  incr[6] = -0.25*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_x; 
  incr[7] = 0.25*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_x; 
  }
#elif upwindType == QUAD 
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
 
  incr[0] = 0.3535533905932737*alpha[0]*fhatAL[0]*dfac_x; 
  incr[1] = -0.6123724356957944*alpha[0]*fhatAL[0]*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*fhatAL[1]*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*fhatAL[2]*dfac_x; 
  incr[4] = -0.6123724356957944*alpha[0]*fhatAL[1]*dfac_x; 
  incr[5] = -0.6123724356957944*alpha[0]*fhatAL[2]*dfac_x; 
  incr[6] = 0.3535533905932737*alpha[0]*fhatAL[3]*dfac_x; 
  incr[7] = -0.6123724356957944*alpha[0]*fhatAL[3]*dfac_x; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[8], outrPos[8]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.3333333333333333 : positivityWeightByDirL[1]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.3333333333333333 : positivityWeightByDirR[1]/positivityWeightByDirR[0]; 
  outlPos[0] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outlPos[1] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outlPos[2] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5])-4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outlPos[3] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*incr[5]-4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outlPos[4] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*incr[5]+4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outlPos[5] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5])+4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outlPos[6] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outlPos[7] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[0] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outrPos[1] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outrPos[2] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*incr[5]-4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outrPos[3] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5])-4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outrPos[4] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5])+4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outrPos[5] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*incr[5]+4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[6] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[7] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])+12.72792206135786*incr[0]); 
  if(outlPos[1] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]+12.72792206135786*fl[0]))/dtApprox/outlPos[1]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[1] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[1] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[3] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[2] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*fr[5]-4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPos[2]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[2] *= limFac; 
  outrPos[3] *= limFac; 
  outrPos[2] *= limFac; 
  if(outlPos[5] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*fl[5]+4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPos[5]); 
  else limFac = 1.0; 
  if(outrPos[4] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0]))/dtApprox/outrPos[4]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[5] *= limFac; 
  outlPos[4] *= limFac; 
  outrPos[5] *= limFac; 
  outrPos[4] *= limFac; 
  if(outlPos[7] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[6] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]-12.72792206135786*fr[0]))/dtApprox/outrPos[6]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[6] *= limFac; 
  outrPos[7] *= limFac; 
  outrPos[6] *= limFac; 
  outr[0] += 0.3535533905932737*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.3535533905932737*(1.732050807568877*outrPos[7]-1.732050807568877*outrPos[6]+1.732050807568877*outrPos[5]-1.732050807568877*outrPos[4]+1.732050807568877*outrPos[3]-1.732050807568877*outrPos[2]+1.732050807568877*outrPos[1]-1.732050807568877*outrPos[0]); 
  outr[2] += 0.3535533905932737*(1.732050807568877*(outrPos[7]+outrPos[6])-1.732050807568877*(outrPos[5]+outrPos[4])+1.732050807568877*(outrPos[3]+outrPos[2])-1.732050807568877*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.3535533905932737*(1.732050807568877*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])-1.732050807568877*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[4] += 0.3535533905932737*(3.0*outrPos[7]-3.0*(outrPos[6]+outrPos[5])+3.0*(outrPos[4]+outrPos[3])-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 
  outr[5] += 0.3535533905932737*(3.0*outrPos[7]-3.0*outrPos[6]+3.0*outrPos[5]-3.0*(outrPos[4]+outrPos[3])+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[6] += 0.3535533905932737*(3.0*(outrPos[7]+outrPos[6])-3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[7] += 0.3535533905932737*(5.196152422706631*outrPos[7]-5.196152422706631*(outrPos[6]+outrPos[5])+5.196152422706631*outrPos[4]-5.196152422706631*outrPos[3]+5.196152422706631*(outrPos[2]+outrPos[1])-5.196152422706631*outrPos[0]); 

  outl[0] += 0.3535533905932737*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.3535533905932737*(1.732050807568877*outlPos[7]-1.732050807568877*outlPos[6]+1.732050807568877*outlPos[5]-1.732050807568877*outlPos[4]+1.732050807568877*outlPos[3]-1.732050807568877*outlPos[2]+1.732050807568877*outlPos[1]-1.732050807568877*outlPos[0]); 
  outl[2] += 0.3535533905932737*(1.732050807568877*(outlPos[7]+outlPos[6])-1.732050807568877*(outlPos[5]+outlPos[4])+1.732050807568877*(outlPos[3]+outlPos[2])-1.732050807568877*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.3535533905932737*(1.732050807568877*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])-1.732050807568877*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[4] += 0.3535533905932737*(3.0*outlPos[7]-3.0*(outlPos[6]+outlPos[5])+3.0*(outlPos[4]+outlPos[3])-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  outl[5] += 0.3535533905932737*(3.0*outlPos[7]-3.0*outlPos[6]+3.0*outlPos[5]-3.0*(outlPos[4]+outlPos[3])+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[6] += 0.3535533905932737*(3.0*(outlPos[7]+outlPos[6])-3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[7] += 0.3535533905932737*(5.196152422706631*outlPos[7]-5.196152422706631*(outlPos[6]+outlPos[5])+5.196152422706631*outlPos[4]-5.196152422706631*outlPos[3]+5.196152422706631*(outlPos[2]+outlPos[1])-5.196152422706631*outlPos[0]); 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, 
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
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  incr[1] = 0.25*alpha[0]*(1.732050807568877*fl[4]+fl[1])*dfac_v; 
  incr[2] = -0.25*alpha[0]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  incr[3] = 0.25*alpha[0]*(1.732050807568877*fl[6]+fl[3])*dfac_v; 
  incr[4] = -0.25*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[1])*dfac_v; 
  incr[5] = 0.25*alpha[0]*(1.732050807568877*fl[7]+fl[5])*dfac_v; 
  incr[6] = -0.25*alpha[0]*(3.0*fl[6]+1.732050807568877*fl[3])*dfac_v; 
  incr[7] = -0.25*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[5])*dfac_v; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  incr[1] = -0.25*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[1])*dfac_v; 
  incr[2] = 0.25*alpha[0]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  incr[3] = -0.25*alpha[0]*(1.732050807568877*fr[6]-1.0*fr[3])*dfac_v; 
  incr[4] = 0.25*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[1])*dfac_v; 
  incr[5] = -0.25*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[5])*dfac_v; 
  incr[6] = 0.25*alpha[0]*(3.0*fr[6]-1.732050807568877*fr[3])*dfac_v; 
  incr[7] = 0.25*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[5])*dfac_v; 
  }
#elif upwindType == QUAD 
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
 
  incr[0] = 0.3535533905932737*alpha[0]*fhatAL[0]*dfac_v; 
  incr[1] = 0.3535533905932737*alpha[0]*fhatAL[1]*dfac_v; 
  incr[2] = -0.6123724356957944*alpha[0]*fhatAL[0]*dfac_v; 
  incr[3] = 0.3535533905932737*alpha[0]*fhatAL[2]*dfac_v; 
  incr[4] = -0.6123724356957944*alpha[0]*fhatAL[1]*dfac_v; 
  incr[5] = 0.3535533905932737*alpha[0]*fhatAL[3]*dfac_v; 
  incr[6] = -0.6123724356957944*alpha[0]*fhatAL[2]*dfac_v; 
  incr[7] = -0.6123724356957944*alpha[0]*fhatAL[3]*dfac_v; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[8], outrPos[8]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.3333333333333333 : positivityWeightByDirL[2]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.3333333333333333 : positivityWeightByDirR[2]/positivityWeightByDirR[0]; 
  outlPos[0] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*incr[5]-4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outlPos[1] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5])-4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outlPos[2] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outlPos[3] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outlPos[4] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outlPos[5] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outlPos[6] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5])+4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outlPos[7] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*incr[5]+4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[0] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outrPos[1] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outrPos[2] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*incr[5]-4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outrPos[3] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5])-4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outrPos[4] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5])+4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outrPos[5] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*incr[5]+4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[6] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[7] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])+12.72792206135786*incr[0]); 
  if(outlPos[2] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*fl[5]-4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]+12.72792206135786*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[2] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[3] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[1] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[1] *= limFac; 
  outrPos[3] *= limFac; 
  outrPos[1] *= limFac; 
  if(outlPos[6] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPos[6]); 
  else limFac = 1.0; 
  if(outrPos[4] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0]))/dtApprox/outrPos[4]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[6] *= limFac; 
  outlPos[4] *= limFac; 
  outrPos[6] *= limFac; 
  outrPos[4] *= limFac; 
  if(outlPos[7] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[5] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*fr[5]+4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]-12.72792206135786*fr[0]))/dtApprox/outrPos[5]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[5] *= limFac; 
  outrPos[7] *= limFac; 
  outrPos[5] *= limFac; 
  outr[0] += 0.3535533905932737*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.3535533905932737*(1.732050807568877*outrPos[7]-1.732050807568877*outrPos[6]+1.732050807568877*outrPos[5]-1.732050807568877*outrPos[4]+1.732050807568877*outrPos[3]-1.732050807568877*outrPos[2]+1.732050807568877*outrPos[1]-1.732050807568877*outrPos[0]); 
  outr[2] += 0.3535533905932737*(1.732050807568877*(outrPos[7]+outrPos[6])-1.732050807568877*(outrPos[5]+outrPos[4])+1.732050807568877*(outrPos[3]+outrPos[2])-1.732050807568877*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.3535533905932737*(1.732050807568877*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])-1.732050807568877*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[4] += 0.3535533905932737*(3.0*outrPos[7]-3.0*(outrPos[6]+outrPos[5])+3.0*(outrPos[4]+outrPos[3])-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 
  outr[5] += 0.3535533905932737*(3.0*outrPos[7]-3.0*outrPos[6]+3.0*outrPos[5]-3.0*(outrPos[4]+outrPos[3])+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[6] += 0.3535533905932737*(3.0*(outrPos[7]+outrPos[6])-3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[7] += 0.3535533905932737*(5.196152422706631*outrPos[7]-5.196152422706631*(outrPos[6]+outrPos[5])+5.196152422706631*outrPos[4]-5.196152422706631*outrPos[3]+5.196152422706631*(outrPos[2]+outrPos[1])-5.196152422706631*outrPos[0]); 

  outl[0] += 0.3535533905932737*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.3535533905932737*(1.732050807568877*outlPos[7]-1.732050807568877*outlPos[6]+1.732050807568877*outlPos[5]-1.732050807568877*outlPos[4]+1.732050807568877*outlPos[3]-1.732050807568877*outlPos[2]+1.732050807568877*outlPos[1]-1.732050807568877*outlPos[0]); 
  outl[2] += 0.3535533905932737*(1.732050807568877*(outlPos[7]+outlPos[6])-1.732050807568877*(outlPos[5]+outlPos[4])+1.732050807568877*(outlPos[3]+outlPos[2])-1.732050807568877*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.3535533905932737*(1.732050807568877*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])-1.732050807568877*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[4] += 0.3535533905932737*(3.0*outlPos[7]-3.0*(outlPos[6]+outlPos[5])+3.0*(outlPos[4]+outlPos[3])-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  outl[5] += 0.3535533905932737*(3.0*outlPos[7]-3.0*outlPos[6]+3.0*outlPos[5]-3.0*(outlPos[4]+outlPos[3])+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[6] += 0.3535533905932737*(3.0*(outlPos[7]+outlPos[6])-3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[7] += 0.3535533905932737*(5.196152422706631*outlPos[7]-5.196152422706631*(outlPos[6]+outlPos[5])+5.196152422706631*outlPos[4]-5.196152422706631*outlPos[3]+5.196152422706631*(outlPos[2]+outlPos[1])-5.196152422706631*outlPos[0]); 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, 
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
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_x; 
  incr[4] = -0.25*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_x; 
  incr[5] = -0.25*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_x; 
  incr[6] = 0.25*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_x; 
  incr[7] = -0.25*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_x; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_x; 
  incr[4] = 0.25*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_x; 
  incr[5] = 0.25*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_x; 
  incr[6] = -0.25*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_x; 
  incr[7] = 0.25*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_x; 
  }
#elif upwindType == QUAD 
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
 
  incr[0] = 0.3535533905932737*alpha[0]*fhatAL[0]*dfac_x; 
  incr[1] = -0.6123724356957944*alpha[0]*fhatAL[0]*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*fhatAL[1]*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*fhatAL[2]*dfac_x; 
  incr[4] = -0.6123724356957944*alpha[0]*fhatAL[1]*dfac_x; 
  incr[5] = -0.6123724356957944*alpha[0]*fhatAL[2]*dfac_x; 
  incr[6] = 0.3535533905932737*alpha[0]*fhatAL[3]*dfac_x; 
  incr[7] = -0.6123724356957944*alpha[0]*fhatAL[3]*dfac_x; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[8], outrPos[8]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.3333333333333333 : positivityWeightByDirL[1]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.3333333333333333 : positivityWeightByDirR[1]/positivityWeightByDirR[0]; 
  outlPos[0] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outlPos[1] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outlPos[2] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5])-4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outlPos[3] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*incr[5]-4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outlPos[4] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*incr[5]+4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outlPos[5] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5])+4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outlPos[6] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outlPos[7] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[0] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outrPos[1] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outrPos[2] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*incr[5]-4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outrPos[3] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5])-4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outrPos[4] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5])+4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outrPos[5] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*incr[5]+4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[6] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[7] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])+12.72792206135786*incr[0]); 
  if(outlPos[1] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]+12.72792206135786*fl[0]))/dtApprox/outlPos[1]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[1] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[1] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[3] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[2] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*fr[5]-4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPos[2]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[2] *= limFac; 
  outrPos[3] *= limFac; 
  outrPos[2] *= limFac; 
  if(outlPos[5] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*fl[5]+4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPos[5]); 
  else limFac = 1.0; 
  if(outrPos[4] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0]))/dtApprox/outrPos[4]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[5] *= limFac; 
  outlPos[4] *= limFac; 
  outrPos[5] *= limFac; 
  outrPos[4] *= limFac; 
  if(outlPos[7] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[6] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]-12.72792206135786*fr[0]))/dtApprox/outrPos[6]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[6] *= limFac; 
  outrPos[7] *= limFac; 
  outrPos[6] *= limFac; 
  outr[0] += 0.3535533905932737*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.3535533905932737*(1.732050807568877*outrPos[7]-1.732050807568877*outrPos[6]+1.732050807568877*outrPos[5]-1.732050807568877*outrPos[4]+1.732050807568877*outrPos[3]-1.732050807568877*outrPos[2]+1.732050807568877*outrPos[1]-1.732050807568877*outrPos[0]); 
  outr[2] += 0.3535533905932737*(1.732050807568877*(outrPos[7]+outrPos[6])-1.732050807568877*(outrPos[5]+outrPos[4])+1.732050807568877*(outrPos[3]+outrPos[2])-1.732050807568877*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.3535533905932737*(1.732050807568877*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])-1.732050807568877*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[4] += 0.3535533905932737*(3.0*outrPos[7]-3.0*(outrPos[6]+outrPos[5])+3.0*(outrPos[4]+outrPos[3])-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 
  outr[5] += 0.3535533905932737*(3.0*outrPos[7]-3.0*outrPos[6]+3.0*outrPos[5]-3.0*(outrPos[4]+outrPos[3])+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[6] += 0.3535533905932737*(3.0*(outrPos[7]+outrPos[6])-3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[7] += 0.3535533905932737*(5.196152422706631*outrPos[7]-5.196152422706631*(outrPos[6]+outrPos[5])+5.196152422706631*outrPos[4]-5.196152422706631*outrPos[3]+5.196152422706631*(outrPos[2]+outrPos[1])-5.196152422706631*outrPos[0]); 

  outl[0] += 0.3535533905932737*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.3535533905932737*(1.732050807568877*outlPos[7]-1.732050807568877*outlPos[6]+1.732050807568877*outlPos[5]-1.732050807568877*outlPos[4]+1.732050807568877*outlPos[3]-1.732050807568877*outlPos[2]+1.732050807568877*outlPos[1]-1.732050807568877*outlPos[0]); 
  outl[2] += 0.3535533905932737*(1.732050807568877*(outlPos[7]+outlPos[6])-1.732050807568877*(outlPos[5]+outlPos[4])+1.732050807568877*(outlPos[3]+outlPos[2])-1.732050807568877*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.3535533905932737*(1.732050807568877*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])-1.732050807568877*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[4] += 0.3535533905932737*(3.0*outlPos[7]-3.0*(outlPos[6]+outlPos[5])+3.0*(outlPos[4]+outlPos[3])-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  outl[5] += 0.3535533905932737*(3.0*outlPos[7]-3.0*outlPos[6]+3.0*outlPos[5]-3.0*(outlPos[4]+outlPos[3])+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[6] += 0.3535533905932737*(3.0*(outlPos[7]+outlPos[6])-3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[7] += 0.3535533905932737*(5.196152422706631*outlPos[7]-5.196152422706631*(outlPos[6]+outlPos[5])+5.196152422706631*outlPos[4]-5.196152422706631*outlPos[3]+5.196152422706631*(outlPos[2]+outlPos[1])-5.196152422706631*outlPos[0]); 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, 
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
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*(alpha[3]*fl[7]+alpha[2]*fl[6])+alpha[3]*fl[5]+1.732050807568877*alpha[1]*fl[4]+alpha[2]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*(alpha[2]*fl[7]+alpha[3]*fl[6])+alpha[2]*fl[5]+1.732050807568877*alpha[0]*fl[4]+alpha[3]*fl[3]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*(alpha[3]*fl[7]+alpha[2]*fl[6])+1.732050807568877*alpha[3]*fl[5]+3.0*alpha[1]*fl[4]+1.732050807568877*alpha[2]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[3] = 0.25*(1.732050807568877*(alpha[1]*fl[7]+alpha[0]*fl[6])+alpha[1]*fl[5]+1.732050807568877*alpha[3]*fl[4]+alpha[0]*fl[3]+fl[1]*alpha[3]+alpha[2]*(1.732050807568877*fl[2]+fl[0]))*dfac_v; 
  incr[4] = -0.25*(3.0*(alpha[2]*fl[7]+alpha[3]*fl[6])+1.732050807568877*alpha[2]*fl[5]+3.0*alpha[0]*fl[4]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[1]*fl[2]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[5] = 0.25*(1.732050807568877*(alpha[0]*fl[7]+alpha[1]*fl[6])+alpha[0]*fl[5]+1.732050807568877*alpha[2]*fl[4]+alpha[1]*fl[3]+(1.732050807568877*fl[2]+fl[0])*alpha[3]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.25*(3.0*(alpha[1]*fl[7]+alpha[0]*fl[6])+1.732050807568877*alpha[1]*fl[5]+3.0*alpha[3]*fl[4]+1.732050807568877*(alpha[0]*fl[3]+fl[1]*alpha[3])+alpha[2]*(3.0*fl[2]+1.732050807568877*fl[0]))*dfac_v; 
  incr[7] = -0.25*(3.0*(alpha[0]*fl[7]+alpha[1]*fl[6])+1.732050807568877*alpha[0]*fl[5]+3.0*alpha[2]*fl[4]+1.732050807568877*alpha[1]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*(fl[0]*alpha[3]+fl[1]*alpha[2]))*dfac_v; 
  } else { 
  incr[0] = -0.25*(1.732050807568877*(alpha[3]*fr[7]+alpha[2]*fr[6])-1.0*alpha[3]*fr[5]+1.732050807568877*alpha[1]*fr[4]-1.0*alpha[2]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*(alpha[2]*fr[7]+alpha[3]*fr[6])-1.0*alpha[2]*fr[5]+1.732050807568877*alpha[0]*fr[4]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[1]*fr[2]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = 0.25*(3.0*(alpha[3]*fr[7]+alpha[2]*fr[6])-1.732050807568877*alpha[3]*fr[5]+3.0*alpha[1]*fr[4]-1.732050807568877*alpha[2]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[3] = -0.25*(1.732050807568877*(alpha[1]*fr[7]+alpha[0]*fr[6])-1.0*alpha[1]*fr[5]+1.732050807568877*alpha[3]*fr[4]-1.0*(alpha[0]*fr[3]+fr[1]*alpha[3])+alpha[2]*(1.732050807568877*fr[2]-1.0*fr[0]))*dfac_v; 
  incr[4] = 0.25*(3.0*(alpha[2]*fr[7]+alpha[3]*fr[6])-1.732050807568877*alpha[2]*fr[5]+3.0*alpha[0]*fr[4]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[1]*fr[2]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[5] = -0.25*(1.732050807568877*(alpha[0]*fr[7]+alpha[1]*fr[6])-1.0*alpha[0]*fr[5]+1.732050807568877*alpha[2]*fr[4]-1.0*alpha[1]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*(fr[0]*alpha[3]+fr[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.25*(3.0*(alpha[1]*fr[7]+alpha[0]*fr[6])-1.732050807568877*alpha[1]*fr[5]+3.0*alpha[3]*fr[4]-1.732050807568877*(alpha[0]*fr[3]+fr[1]*alpha[3])+alpha[2]*(3.0*fr[2]-1.732050807568877*fr[0]))*dfac_v; 
  incr[7] = 0.25*(3.0*(alpha[0]*fr[7]+alpha[1]*fr[6])-1.732050807568877*alpha[0]*fr[5]+3.0*alpha[2]*fr[4]-1.732050807568877*alpha[1]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*(fr[0]*alpha[3]+fr[1]*alpha[2]))*dfac_v; 
  }
#elif upwindType == QUAD 
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
 
  incr[0] = 0.3535533905932737*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[2] = -0.6123724356957944*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[3] = 0.3535533905932737*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[4] = -0.6123724356957944*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[5] = 0.3535533905932737*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 
  incr[6] = -0.6123724356957944*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[7] = -0.6123724356957944*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[8], outrPos[8]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.3333333333333333 : positivityWeightByDirL[2]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.3333333333333333 : positivityWeightByDirR[2]/positivityWeightByDirR[0]; 
  outlPos[0] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*incr[5]-4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outlPos[1] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5])-4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outlPos[2] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outlPos[3] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outlPos[4] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outlPos[5] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outlPos[6] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5])+4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outlPos[7] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*incr[5]+4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[0] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outrPos[1] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outrPos[2] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*incr[5]-4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]+12.72792206135786*incr[0]); 
  outrPos[3] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5])-4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])-12.72792206135786*incr[0]); 
  outrPos[4] = 0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*(incr[6]+incr[5])+4.242640687119286*incr[4]+7.348469228349534*incr[3]-7.348469228349534*(incr[2]+incr[1])+12.72792206135786*incr[0]); 
  outrPos[5] = -0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*incr[6]-4.242640687119286*incr[5]+4.242640687119286*incr[4]-7.348469228349534*incr[3]+7.348469228349534*incr[2]-7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[6] = -0.02777777777777778*(2.449489742783178*incr[7]-4.242640687119286*incr[6]+4.242640687119286*(incr[5]+incr[4])-7.348469228349534*(incr[3]+incr[2])+7.348469228349534*incr[1]-12.72792206135786*incr[0]); 
  outrPos[7] = 0.02777777777777778*(2.449489742783178*incr[7]+4.242640687119286*(incr[6]+incr[5]+incr[4])+7.348469228349534*(incr[3]+incr[2]+incr[1])+12.72792206135786*incr[0]); 
  if(outlPos[2] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*fl[5]-4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]+12.72792206135786*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[2] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[3] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[1] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[1] *= limFac; 
  outrPos[3] *= limFac; 
  outrPos[1] *= limFac; 
  if(outlPos[6] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPos[6]); 
  else limFac = 1.0; 
  if(outrPos[4] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0]))/dtApprox/outrPos[4]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[6] *= limFac; 
  outlPos[4] *= limFac; 
  outrPos[6] *= limFac; 
  outrPos[4] *= limFac; 
  if(outlPos[7] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[5] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*fr[5]+4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]-12.72792206135786*fr[0]))/dtApprox/outrPos[5]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[5] *= limFac; 
  outrPos[7] *= limFac; 
  outrPos[5] *= limFac; 
  outr[0] += 0.3535533905932737*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.3535533905932737*(1.732050807568877*outrPos[7]-1.732050807568877*outrPos[6]+1.732050807568877*outrPos[5]-1.732050807568877*outrPos[4]+1.732050807568877*outrPos[3]-1.732050807568877*outrPos[2]+1.732050807568877*outrPos[1]-1.732050807568877*outrPos[0]); 
  outr[2] += 0.3535533905932737*(1.732050807568877*(outrPos[7]+outrPos[6])-1.732050807568877*(outrPos[5]+outrPos[4])+1.732050807568877*(outrPos[3]+outrPos[2])-1.732050807568877*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.3535533905932737*(1.732050807568877*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])-1.732050807568877*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[4] += 0.3535533905932737*(3.0*outrPos[7]-3.0*(outrPos[6]+outrPos[5])+3.0*(outrPos[4]+outrPos[3])-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 
  outr[5] += 0.3535533905932737*(3.0*outrPos[7]-3.0*outrPos[6]+3.0*outrPos[5]-3.0*(outrPos[4]+outrPos[3])+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[6] += 0.3535533905932737*(3.0*(outrPos[7]+outrPos[6])-3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[7] += 0.3535533905932737*(5.196152422706631*outrPos[7]-5.196152422706631*(outrPos[6]+outrPos[5])+5.196152422706631*outrPos[4]-5.196152422706631*outrPos[3]+5.196152422706631*(outrPos[2]+outrPos[1])-5.196152422706631*outrPos[0]); 

  outl[0] += 0.3535533905932737*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.3535533905932737*(1.732050807568877*outlPos[7]-1.732050807568877*outlPos[6]+1.732050807568877*outlPos[5]-1.732050807568877*outlPos[4]+1.732050807568877*outlPos[3]-1.732050807568877*outlPos[2]+1.732050807568877*outlPos[1]-1.732050807568877*outlPos[0]); 
  outl[2] += 0.3535533905932737*(1.732050807568877*(outlPos[7]+outlPos[6])-1.732050807568877*(outlPos[5]+outlPos[4])+1.732050807568877*(outlPos[3]+outlPos[2])-1.732050807568877*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.3535533905932737*(1.732050807568877*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])-1.732050807568877*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[4] += 0.3535533905932737*(3.0*outlPos[7]-3.0*(outlPos[6]+outlPos[5])+3.0*(outlPos[4]+outlPos[3])-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  outl[5] += 0.3535533905932737*(3.0*outlPos[7]-3.0*outlPos[6]+3.0*outlPos[5]-3.0*(outlPos[4]+outlPos[3])+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[6] += 0.3535533905932737*(3.0*(outlPos[7]+outlPos[6])-3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[7] += 0.3535533905932737*(5.196152422706631*outlPos[7]-5.196152422706631*(outlPos[6]+outlPos[5])+5.196152422706631*outlPos[4]-5.196152422706631*outlPos[3]+5.196152422706631*(outlPos[2]+outlPos[1])-5.196152422706631*outlPos[0]); 
  return std::abs(alpha0); 
} 
