#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionSurfPositivity3xSer_X1_P1(const double *cflFreqCtrlL, const double *cflFreqCtrlR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  double dfac3 = 2.0/dxv[2]; 
  double w3 = w[2]; 
  const double *v1 = &fr[8]; 
  const double *v2 = &fr[16]; 
  const double *v3 = &fr[24]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.1767766952966368*(1.732050807568877*v1[1]-1.0*v1[0]); 

  double alpha[4]; 
  alpha[0] = -0.5*(2.449489742783178*v1[1]-1.414213562373095*v1[0]); 
  alpha[1] = -0.5*(2.449489742783178*v1[4]-1.414213562373095*v1[2]); 
  alpha[2] = -0.5*(2.449489742783178*v1[5]-1.414213562373095*v1[3]); 
  alpha[3] = -0.5*(2.449489742783178*v1[7]-1.414213562373095*v1[6]); 
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
  // control node [y,z] = [-1/3,-1/3] 
  fCtrlL[0] = 0.04811252243246882*(2.449489742783178*fl[6]-4.242640687119286*(fl[3]+fl[2])+7.348469228349534*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = 0.04811252243246882*(2.449489742783178*fr[6]-4.242640687119286*(fr[3]+fr[2])+7.348469228349534*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [y,z] = [1/3,-1/3] 
  fCtrlL[1] = -0.04811252243246882*(2.449489742783178*fl[6]+4.242640687119286*fl[3]-4.242640687119286*fl[2]-7.348469228349534*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = -0.04811252243246882*(2.449489742783178*fr[6]+4.242640687119286*fr[3]-4.242640687119286*fr[2]-7.348469228349534*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [y,z] = [-1/3,1/3] 
  fCtrlL[2] = -0.04811252243246882*(2.449489742783178*fl[6]-4.242640687119286*fl[3]+4.242640687119286*fl[2]-7.348469228349534*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = -0.04811252243246882*(2.449489742783178*fr[6]-4.242640687119286*fr[3]+4.242640687119286*fr[2]-7.348469228349534*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [y,z] = [1/3,1/3] 
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
  alphaQuad = 0.5*alpha[3]-0.5*(alpha[2]+alpha[1])+0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = std::max(0., 0.5*(fL_AL[3]-1.0*(fL_AL[2]+fL_AL[1])+fL_AL[0])); 
  } else {
  fhatALQuad[0] = std::max(0., 0.5*(fR_AL[3]-1.0*(fR_AL[2]+fR_AL[1])+fR_AL[0])); 
  } 
  alphaQuad = 0.5*(alpha[1]+alpha[0])-0.5*(alpha[3]+alpha[2]); 
  if(alphaQuad > 0) {
  fhatALQuad[1] = std::max(0., -0.5*(fL_AL[3]+fL_AL[2]-1.0*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[1] = std::max(0., -0.5*(fR_AL[3]+fR_AL[2]-1.0*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = (-0.5*alpha[3])+0.5*alpha[2]-0.5*alpha[1]+0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = std::max(0., -0.5*(fL_AL[3]-1.0*fL_AL[2]+fL_AL[1]-1.0*fL_AL[0])); 
  } else {
  fhatALQuad[2] = std::max(0., -0.5*(fR_AL[3]-1.0*fR_AL[2]+fR_AL[1]-1.0*fR_AL[0])); 
  } 
  alphaQuad = 0.5*(alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = std::max(0., 0.5*(fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[3] = std::max(0., 0.5*(fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  double GhatCtrl[4];
  // control node [y,z] = [-1/3,-1/3] 
  GhatCtrl[0] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.1443375672974065*((alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac1; 
  if(std::abs(GhatCtrl[0]) < EPSILON) GhatCtrl[0] = 0.; 
  // control node [y,z] = [1/3,-1/3] 
  GhatCtrl[1] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.1443375672974065*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])-0.1443375672974065*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac1; 
  if(std::abs(GhatCtrl[1]) < EPSILON) GhatCtrl[1] = 0.; 
  // control node [y,z] = [-1/3,1/3] 
  GhatCtrl[2] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.1443375672974065*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.1443375672974065*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac1; 
  if(std::abs(GhatCtrl[2]) < EPSILON) GhatCtrl[2] = 0.; 
  // control node [y,z] = [1/3,1/3] 
  GhatCtrl[3] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.1443375672974065*((alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac1; 
  if(std::abs(GhatCtrl[3]) < EPSILON) GhatCtrl[3] = 0.; 
  double uFrac = 0., fCtrl = 0., alphaCtrl = 0.;
  if(GhatCtrl[0]<-EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]-0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac1)/cflFreqCtrlR[0]; 
    fCtrl = (-0.06804138174397717*fr[7])+0.1178511301977579*(fr[6]+fr[5]+fr[4])-0.2041241452319315*(fr[3]+fr[2]+fr[1])+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = -std::min(std::abs(GhatCtrl[0]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[0]>EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]-0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac1)/cflFreqCtrlL[1]; 
    fCtrl = 0.06804138174397717*fl[7]+0.1178511301977579*fl[6]-0.1178511301977579*(fl[5]+fl[4])-0.2041241452319315*(fl[3]+fl[2])+0.2041241452319315*fl[1]+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = std::min(GhatCtrl[0], std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[1]<-EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])-0.2886751345948129*alpha[2]+0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac1)/cflFreqCtrlR[2]; 
    fCtrl = 0.06804138174397717*fr[7]-0.1178511301977579*fr[6]+0.1178511301977579*fr[5]-0.1178511301977579*fr[4]-0.2041241452319315*fr[3]+0.2041241452319315*fr[2]-0.2041241452319315*fr[1]+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = -std::min(std::abs(GhatCtrl[1]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[1]>EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])-0.2886751345948129*alpha[2]+0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac1)/cflFreqCtrlL[3]; 
    fCtrl = (-0.06804138174397717*fl[7])-0.1178511301977579*(fl[6]+fl[5])+0.1178511301977579*fl[4]-0.2041241452319315*fl[3]+0.2041241452319315*(fl[2]+fl[1])+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = std::min(GhatCtrl[1], std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[2]<-EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])+0.2886751345948129*alpha[2]-0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac1)/cflFreqCtrlR[4]; 
    fCtrl = 0.06804138174397717*fr[7]-0.1178511301977579*(fr[6]+fr[5])+0.1178511301977579*fr[4]+0.2041241452319315*fr[3]-0.2041241452319315*(fr[2]+fr[1])+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = -std::min(std::abs(GhatCtrl[2]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[2]>EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])+0.2886751345948129*alpha[2]-0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac1)/cflFreqCtrlL[5]; 
    fCtrl = (-0.06804138174397717*fl[7])-0.1178511301977579*fl[6]+0.1178511301977579*fl[5]-0.1178511301977579*fl[4]+0.2041241452319315*fl[3]-0.2041241452319315*fl[2]+0.2041241452319315*fl[1]+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = std::min(GhatCtrl[2], std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[3]<-EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]+0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac1)/cflFreqCtrlR[6]; 
    fCtrl = (-0.06804138174397717*fr[7])+0.1178511301977579*fr[6]-0.1178511301977579*(fr[5]+fr[4])+0.2041241452319315*(fr[3]+fr[2])-0.2041241452319315*fr[1]+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = -std::min(std::abs(GhatCtrl[3]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[3]>EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]+0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac1)/cflFreqCtrlL[7]; 
    fCtrl = 0.06804138174397717*fl[7]+0.1178511301977579*(fl[6]+fl[5]+fl[4])+0.2041241452319315*(fl[3]+fl[2]+fl[1])+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = std::min(GhatCtrl[3], std::abs(uFrac)*fCtrl/dt);
  }

  incr[0] = 0.3535533905932737*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[1] = -0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[2] = 0.6123724356957944*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[3] = 0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[4] = -1.060660171779821*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[5] = -1.060660171779821*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[6] = 1.060660171779821*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[7] = -1.837117307087383*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 

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
double PassiveAdvectionSurfPositivity3xSer_X2_P1(const double *cflFreqCtrlL, const double *cflFreqCtrlR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  double dfac3 = 2.0/dxv[2]; 
  double w3 = w[2]; 
  const double *v1 = &fr[8]; 
  const double *v2 = &fr[16]; 
  const double *v3 = &fr[24]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.1767766952966368*(1.732050807568877*v2[2]-1.0*v2[0]); 

  double alpha[4]; 
  alpha[0] = -0.5*(2.449489742783178*v2[2]-1.414213562373095*v2[0]); 
  alpha[1] = -0.5*(2.449489742783178*v2[4]-1.414213562373095*v2[1]); 
  alpha[2] = -0.5*(2.449489742783178*v2[6]-1.414213562373095*v2[3]); 
  alpha[3] = -0.5*(2.449489742783178*v2[7]-1.414213562373095*v2[5]); 
  double rCtrlL[4], rCtrlR[4];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrlL[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rCtrlL[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlR[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rCtrlR[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  double fCtrlL[4], fCtrlR[4];  // fCtrl = anti-limited f evaluated at each control node on y surface 
  // control node [x,z] = [-1/3,-1/3] 
  fCtrlL[0] = 0.04811252243246882*(2.449489742783178*fl[5]-4.242640687119286*(fl[3]+fl[1])+7.348469228349534*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = 0.04811252243246882*(2.449489742783178*fr[5]-4.242640687119286*(fr[3]+fr[1])+7.348469228349534*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,z] = [1/3,-1/3] 
  fCtrlL[1] = -0.04811252243246882*(2.449489742783178*fl[5]+4.242640687119286*fl[3]-4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = -0.04811252243246882*(2.449489742783178*fr[5]+4.242640687119286*fr[3]-4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,z] = [-1/3,1/3] 
  fCtrlL[2] = -0.04811252243246882*(2.449489742783178*fl[5]-4.242640687119286*fl[3]+4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = -0.04811252243246882*(2.449489742783178*fr[5]-4.242640687119286*fr[3]+4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,z] = [1/3,1/3] 
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
  fhatALQuad[0] = std::max(0., 0.5*(fL_AL[3]-1.0*(fL_AL[2]+fL_AL[1])+fL_AL[0])); 
  } else {
  fhatALQuad[0] = std::max(0., 0.5*(fR_AL[3]-1.0*(fR_AL[2]+fR_AL[1])+fR_AL[0])); 
  } 
  alphaQuad = 0.5*(alpha[1]+alpha[0])-0.5*(alpha[3]+alpha[2]); 
  if(alphaQuad > 0) {
  fhatALQuad[1] = std::max(0., -0.5*(fL_AL[3]+fL_AL[2]-1.0*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[1] = std::max(0., -0.5*(fR_AL[3]+fR_AL[2]-1.0*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = (-0.5*alpha[3])+0.5*alpha[2]-0.5*alpha[1]+0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = std::max(0., -0.5*(fL_AL[3]-1.0*fL_AL[2]+fL_AL[1]-1.0*fL_AL[0])); 
  } else {
  fhatALQuad[2] = std::max(0., -0.5*(fR_AL[3]-1.0*fR_AL[2]+fR_AL[1]-1.0*fR_AL[0])); 
  } 
  alphaQuad = 0.5*(alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = std::max(0., 0.5*(fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[3] = std::max(0., 0.5*(fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  double GhatCtrl[4];
  // control node [x,z] = [-1/3,-1/3] 
  GhatCtrl[0] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.1443375672974065*((alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac2; 
  if(std::abs(GhatCtrl[0]) < EPSILON) GhatCtrl[0] = 0.; 
  // control node [x,z] = [1/3,-1/3] 
  GhatCtrl[1] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.1443375672974065*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])-0.1443375672974065*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac2; 
  if(std::abs(GhatCtrl[1]) < EPSILON) GhatCtrl[1] = 0.; 
  // control node [x,z] = [-1/3,1/3] 
  GhatCtrl[2] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.1443375672974065*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.1443375672974065*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac2; 
  if(std::abs(GhatCtrl[2]) < EPSILON) GhatCtrl[2] = 0.; 
  // control node [x,z] = [1/3,1/3] 
  GhatCtrl[3] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.1443375672974065*((alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac2; 
  if(std::abs(GhatCtrl[3]) < EPSILON) GhatCtrl[3] = 0.; 
  double uFrac = 0., fCtrl = 0., alphaCtrl = 0.;
  if(GhatCtrl[0]<-EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]-0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac2)/cflFreqCtrlR[0]; 
    fCtrl = (-0.06804138174397717*fr[7])+0.1178511301977579*(fr[6]+fr[5]+fr[4])-0.2041241452319315*(fr[3]+fr[2]+fr[1])+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = -std::min(std::abs(GhatCtrl[0]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[1]<-EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])-0.2886751345948129*alpha[2]+0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac2)/cflFreqCtrlR[1]; 
    fCtrl = 0.06804138174397717*fr[7]+0.1178511301977579*fr[6]-0.1178511301977579*(fr[5]+fr[4])-0.2041241452319315*(fr[3]+fr[2])+0.2041241452319315*fr[1]+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = -std::min(std::abs(GhatCtrl[1]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[0]>EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]-0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac2)/cflFreqCtrlL[2]; 
    fCtrl = 0.06804138174397717*fl[7]-0.1178511301977579*fl[6]+0.1178511301977579*fl[5]-0.1178511301977579*fl[4]-0.2041241452319315*fl[3]+0.2041241452319315*fl[2]-0.2041241452319315*fl[1]+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = std::min(GhatCtrl[0], std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[1]>EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])-0.2886751345948129*alpha[2]+0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac2)/cflFreqCtrlL[3]; 
    fCtrl = (-0.06804138174397717*fl[7])-0.1178511301977579*(fl[6]+fl[5])+0.1178511301977579*fl[4]-0.2041241452319315*fl[3]+0.2041241452319315*(fl[2]+fl[1])+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = std::min(GhatCtrl[1], std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[2]<-EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])+0.2886751345948129*alpha[2]-0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac2)/cflFreqCtrlR[4]; 
    fCtrl = 0.06804138174397717*fr[7]-0.1178511301977579*(fr[6]+fr[5])+0.1178511301977579*fr[4]+0.2041241452319315*fr[3]-0.2041241452319315*(fr[2]+fr[1])+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = -std::min(std::abs(GhatCtrl[2]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[3]<-EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]+0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac2)/cflFreqCtrlR[5]; 
    fCtrl = (-0.06804138174397717*fr[7])-0.1178511301977579*fr[6]+0.1178511301977579*fr[5]-0.1178511301977579*fr[4]+0.2041241452319315*fr[3]-0.2041241452319315*fr[2]+0.2041241452319315*fr[1]+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = -std::min(std::abs(GhatCtrl[3]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[2]>EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])+0.2886751345948129*alpha[2]-0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac2)/cflFreqCtrlL[6]; 
    fCtrl = (-0.06804138174397717*fl[7])+0.1178511301977579*fl[6]-0.1178511301977579*(fl[5]+fl[4])+0.2041241452319315*(fl[3]+fl[2])-0.2041241452319315*fl[1]+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = std::min(GhatCtrl[2], std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[3]>EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]+0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac2)/cflFreqCtrlL[7]; 
    fCtrl = 0.06804138174397717*fl[7]+0.1178511301977579*(fl[6]+fl[5]+fl[4])+0.2041241452319315*(fl[3]+fl[2]+fl[1])+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = std::min(GhatCtrl[3], std::abs(uFrac)*fCtrl/dt);
  }

  incr[0] = 0.3535533905932737*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[1] = 0.6123724356957944*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[2] = -0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[3] = 0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[4] = -1.060660171779821*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[5] = 1.060660171779821*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[6] = -1.060660171779821*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[7] = -1.837117307087383*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 

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
double PassiveAdvectionSurfPositivity3xSer_X3_P1(const double *cflFreqCtrlL, const double *cflFreqCtrlR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  double dfac3 = 2.0/dxv[2]; 
  double w3 = w[2]; 
  const double *v1 = &fr[8]; 
  const double *v2 = &fr[16]; 
  const double *v3 = &fr[24]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.1767766952966368*(1.732050807568877*v3[3]-1.0*v3[0]); 

  double alpha[4]; 
  alpha[0] = -0.5*(2.449489742783178*v3[3]-1.414213562373095*v3[0]); 
  alpha[1] = -0.5*(2.449489742783178*v3[5]-1.414213562373095*v3[1]); 
  alpha[2] = -0.5*(2.449489742783178*v3[6]-1.414213562373095*v3[2]); 
  alpha[3] = -0.5*(2.449489742783178*v3[7]-1.414213562373095*v3[4]); 
  double rCtrlL[4], rCtrlR[4];  // rCtrl=f1/f0 at each control node in dimensions other than z 
  rCtrlL[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[5])+9.0*fl[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[4]-3.0*(fl[2]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[5]+9.0*fl[3])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[4])+3.0*(fl[1]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rCtrlL[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]-1.0*fl[6])-9.0*fl[3]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[4])+3.0*(fl[2]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[5])+9.0*fl[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[4]+3.0*(fl[2]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlR[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[5])+9.0*fr[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[4]-3.0*(fr[2]+fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[5]+9.0*fr[3])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[4])+3.0*(fr[1]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rCtrlR[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]-1.0*fr[6])-9.0*fr[3]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[4])+3.0*(fr[2]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[5])+9.0*fr[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[4]+3.0*(fr[2]+fr[1])+5.196152422706631*fr[0])); 
  double fCtrlL[4], fCtrlR[4];  // fCtrl = anti-limited f evaluated at each control node on z surface 
  // control node [x,y] = [-1/3,-1/3] 
  fCtrlL[0] = 0.04811252243246882*(2.449489742783178*fl[4]-4.242640687119286*(fl[2]+fl[1])+7.348469228349534*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = 0.04811252243246882*(2.449489742783178*fr[4]-4.242640687119286*(fr[2]+fr[1])+7.348469228349534*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,y] = [1/3,-1/3] 
  fCtrlL[1] = -0.04811252243246882*(2.449489742783178*fl[4]+4.242640687119286*fl[2]-4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = -0.04811252243246882*(2.449489742783178*fr[4]+4.242640687119286*fr[2]-4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,y] = [-1/3,1/3] 
  fCtrlL[2] = -0.04811252243246882*(2.449489742783178*fl[4]-4.242640687119286*fl[2]+4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = -0.04811252243246882*(2.449489742783178*fr[4]-4.242640687119286*fr[2]+4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,y] = [1/3,1/3] 
  fCtrlL[3] = 0.04811252243246882*(2.449489742783178*fl[4]+4.242640687119286*(fl[2]+fl[1])+7.348469228349534*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = 0.04811252243246882*(2.449489742783178*fr[4]+4.242640687119286*(fr[2]+fr[1])+7.348469228349534*fr[0])*limTheta(rCtrlR[3],-1.0); 
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
  fhatALQuad[0] = std::max(0., 0.5*(fL_AL[3]-1.0*(fL_AL[2]+fL_AL[1])+fL_AL[0])); 
  } else {
  fhatALQuad[0] = std::max(0., 0.5*(fR_AL[3]-1.0*(fR_AL[2]+fR_AL[1])+fR_AL[0])); 
  } 
  alphaQuad = 0.5*(alpha[1]+alpha[0])-0.5*(alpha[3]+alpha[2]); 
  if(alphaQuad > 0) {
  fhatALQuad[1] = std::max(0., -0.5*(fL_AL[3]+fL_AL[2]-1.0*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[1] = std::max(0., -0.5*(fR_AL[3]+fR_AL[2]-1.0*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = (-0.5*alpha[3])+0.5*alpha[2]-0.5*alpha[1]+0.5*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = std::max(0., -0.5*(fL_AL[3]-1.0*fL_AL[2]+fL_AL[1]-1.0*fL_AL[0])); 
  } else {
  fhatALQuad[2] = std::max(0., -0.5*(fR_AL[3]-1.0*fR_AL[2]+fR_AL[1]-1.0*fR_AL[0])); 
  } 
  alphaQuad = 0.5*(alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = std::max(0., 0.5*(fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[3] = std::max(0., 0.5*(fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  double GhatCtrl[4];
  // control node [x,y] = [-1/3,-1/3] 
  GhatCtrl[0] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.1443375672974065*((alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac3; 
  if(std::abs(GhatCtrl[0]) < EPSILON) GhatCtrl[0] = 0.; 
  // control node [x,y] = [1/3,-1/3] 
  GhatCtrl[1] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.1443375672974065*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])-0.1443375672974065*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac3; 
  if(std::abs(GhatCtrl[1]) < EPSILON) GhatCtrl[1] = 0.; 
  // control node [x,y] = [-1/3,1/3] 
  GhatCtrl[2] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.1443375672974065*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.1443375672974065*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac3; 
  if(std::abs(GhatCtrl[2]) < EPSILON) GhatCtrl[2] = 0.; 
  // control node [x,y] = [1/3,1/3] 
  GhatCtrl[3] = (0.25*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.1443375672974065*((alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.08333333333333333*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2]))*dfac3; 
  if(std::abs(GhatCtrl[3]) < EPSILON) GhatCtrl[3] = 0.; 
  double uFrac = 0., fCtrl = 0., alphaCtrl = 0.;
  if(GhatCtrl[0]<-EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]-0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac3)/cflFreqCtrlR[0]; 
    fCtrl = (-0.06804138174397717*fr[7])+0.1178511301977579*(fr[6]+fr[5]+fr[4])-0.2041241452319315*(fr[3]+fr[2]+fr[1])+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = -std::min(std::abs(GhatCtrl[0]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[1]<-EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])-0.2886751345948129*alpha[2]+0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac3)/cflFreqCtrlR[1]; 
    fCtrl = 0.06804138174397717*fr[7]+0.1178511301977579*fr[6]-0.1178511301977579*(fr[5]+fr[4])-0.2041241452319315*(fr[3]+fr[2])+0.2041241452319315*fr[1]+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = -std::min(std::abs(GhatCtrl[1]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[2]<-EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])+0.2886751345948129*alpha[2]-0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac3)/cflFreqCtrlR[2]; 
    fCtrl = 0.06804138174397717*fr[7]-0.1178511301977579*fr[6]+0.1178511301977579*fr[5]-0.1178511301977579*fr[4]-0.2041241452319315*fr[3]+0.2041241452319315*fr[2]-0.2041241452319315*fr[1]+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = -std::min(std::abs(GhatCtrl[2]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[3]<-EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]+0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac3)/cflFreqCtrlR[3]; 
    fCtrl = (-0.06804138174397717*fr[7])-0.1178511301977579*(fr[6]+fr[5])+0.1178511301977579*fr[4]-0.2041241452319315*fr[3]+0.2041241452319315*(fr[2]+fr[1])+0.3535533905932737*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = -std::min(std::abs(GhatCtrl[3]), std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[0]>EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]-0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac3)/cflFreqCtrlL[4]; 
    fCtrl = 0.06804138174397717*fl[7]-0.1178511301977579*(fl[6]+fl[5])+0.1178511301977579*fl[4]+0.2041241452319315*fl[3]-0.2041241452319315*(fl[2]+fl[1])+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = std::min(GhatCtrl[0], std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[1]>EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])-0.2886751345948129*alpha[2]+0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac3)/cflFreqCtrlL[5]; 
    fCtrl = (-0.06804138174397717*fl[7])-0.1178511301977579*fl[6]+0.1178511301977579*fl[5]-0.1178511301977579*fl[4]+0.2041241452319315*fl[3]-0.2041241452319315*fl[2]+0.2041241452319315*fl[1]+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = std::min(GhatCtrl[1], std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[2]>EPSILON) {
    alphaCtrl = (-0.1666666666666667*alpha[3])+0.2886751345948129*alpha[2]-0.2886751345948129*alpha[1]+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac3)/cflFreqCtrlL[6]; 
    fCtrl = (-0.06804138174397717*fl[7])+0.1178511301977579*fl[6]-0.1178511301977579*(fl[5]+fl[4])+0.2041241452319315*(fl[3]+fl[2])-0.2041241452319315*fl[1]+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = std::min(GhatCtrl[2], std::abs(uFrac)*fCtrl/dt);
  }
  if(GhatCtrl[3]>EPSILON) {
    alphaCtrl = 0.1666666666666667*alpha[3]+0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0]; 
    uFrac = (alphaCtrl*dfac3)/cflFreqCtrlL[7]; 
    fCtrl = 0.06804138174397717*fl[7]+0.1178511301977579*(fl[6]+fl[5]+fl[4])+0.2041241452319315*(fl[3]+fl[2]+fl[1])+0.3535533905932737*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = std::min(GhatCtrl[3], std::abs(uFrac)*fCtrl/dt);
  }

  incr[0] = 0.3535533905932737*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[1] = 0.6123724356957944*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[2] = 0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[3] = -0.6123724356957944*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[4] = 1.060660171779821*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[5] = -1.060660171779821*(GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[6] = -1.060660171779821*(GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[7] = -1.837117307087383*(GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 

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
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
