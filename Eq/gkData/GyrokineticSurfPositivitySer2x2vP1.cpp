#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity2x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double *cflRateCtrlL, const double *cflRateCtrlR, const double *w, const double *dxv, const double dtApprox, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*(2.0*BdriftX[0]*m_*wv2+BmagInv[0]*(3.0*Phi[3]-1.732050807568877*Phi[2])*dfac_y*q_))/q_; 

  double alpha[8]; 
  alpha[0] = (0.5*(2.828427124746191*BdriftX[0]*m_*wv2+BmagInv[0]*(4.242640687119286*Phi[3]-2.449489742783178*Phi[2])*dfac_y*q_))/q_; 
  alpha[2] = (0.8164965809277261*BdriftX[0]*m_*wv)/(dfac_v*q_); 
  double rCtrlL[8], rCtrlR[8];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrlL[0] = -(1.0*(1.732050807568877*fl[15]-3.0*(fl[13]+fl[12]+fl[11])+5.196152422706631*(fl[8]+fl[6]+fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON-1.0*fl[14]+1.732050807568877*(fl[10]+fl[9]+fl[7])-3.0*(fl[4]+fl[3]+fl[2])+5.196152422706631*fl[0]); 
  rCtrlL[1] = (1.732050807568877*fl[15]+3.0*(fl[13]-1.0*(fl[12]+fl[11]))+5.196152422706631*(fl[5]-1.0*(fl[8]+fl[6]))+9.0*fl[1])/(20.78460969082652*EPSILON+fl[14]+1.732050807568877*(fl[10]-1.0*(fl[9]+fl[7]))+3.0*(fl[2]-1.0*(fl[4]+fl[3]))+5.196152422706631*fl[0]); 
  rCtrlL[2] = (1.732050807568877*fl[15]+3.0*((-1.0*fl[13])+fl[12]-1.0*fl[11])+5.196152422706631*((-1.0*fl[8])+fl[6]-1.0*fl[5])+9.0*fl[1])/(20.78460969082652*EPSILON+fl[14]+1.732050807568877*((-1.0*fl[10])+fl[9]-1.0*fl[7])+3.0*((-1.0*fl[4])+fl[3]-1.0*fl[2])+5.196152422706631*fl[0]); 
  rCtrlL[3] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[13]+fl[12]-1.0*fl[11])+5.196152422706631*fl[8]-1.0*(5.196152422706631*(fl[6]+fl[5])+9.0*fl[1])))/(20.78460969082652*EPSILON-1.0*fl[14]+1.732050807568877*(fl[7]-1.0*(fl[10]+fl[9]))+3.0*((-1.0*fl[4])+fl[3]+fl[2])+5.196152422706631*fl[0]); 
  rCtrlL[4] = (1.732050807568877*fl[15]+3.0*(fl[11]-1.0*(fl[13]+fl[12]))+5.196152422706631*(fl[8]-1.0*(fl[6]+fl[5]))+9.0*fl[1])/(20.78460969082652*EPSILON+fl[14]+1.732050807568877*(fl[7]-1.0*(fl[10]+fl[9]))+3.0*(fl[4]-1.0*(fl[3]+fl[2]))+5.196152422706631*fl[0]); 
  rCtrlL[5] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[13]-1.0*fl[12]+fl[11])+5.196152422706631*(fl[6]-1.0*fl[8])-1.0*(5.196152422706631*fl[5]+9.0*fl[1])))/(20.78460969082652*EPSILON-1.0*fl[14]+1.732050807568877*((-1.0*fl[10])+fl[9]-1.0*fl[7])+3.0*(fl[4]-1.0*fl[3]+fl[2])+5.196152422706631*fl[0]); 
  rCtrlL[6] = -(1.0*(1.732050807568877*fl[15]+3.0*((-1.0*fl[13])+fl[12]+fl[11])+5.196152422706631*(fl[5]-1.0*(fl[8]+fl[6]))-9.0*fl[1]))/(20.78460969082652*EPSILON-1.0*fl[14]+1.732050807568877*(fl[10]-1.0*(fl[9]+fl[7]))+3.0*(fl[4]+fl[3]-1.0*fl[2])+5.196152422706631*fl[0]); 
  rCtrlL[7] = (1.732050807568877*fl[15]+3.0*(fl[13]+fl[12]+fl[11])+5.196152422706631*(fl[8]+fl[6]+fl[5])+9.0*fl[1])/(20.78460969082652*EPSILON+fl[14]+1.732050807568877*(fl[10]+fl[9]+fl[7])+3.0*(fl[4]+fl[3]+fl[2])+5.196152422706631*fl[0]); 
  rCtrlR[0] = -(1.0*(1.732050807568877*fr[15]-3.0*(fr[13]+fr[12]+fr[11])+5.196152422706631*(fr[8]+fr[6]+fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON-1.0*fr[14]+1.732050807568877*(fr[10]+fr[9]+fr[7])-3.0*(fr[4]+fr[3]+fr[2])+5.196152422706631*fr[0]); 
  rCtrlR[1] = (1.732050807568877*fr[15]+3.0*(fr[13]-1.0*(fr[12]+fr[11]))+5.196152422706631*(fr[5]-1.0*(fr[8]+fr[6]))+9.0*fr[1])/(20.78460969082652*EPSILON+fr[14]+1.732050807568877*(fr[10]-1.0*(fr[9]+fr[7]))+3.0*(fr[2]-1.0*(fr[4]+fr[3]))+5.196152422706631*fr[0]); 
  rCtrlR[2] = (1.732050807568877*fr[15]+3.0*((-1.0*fr[13])+fr[12]-1.0*fr[11])+5.196152422706631*((-1.0*fr[8])+fr[6]-1.0*fr[5])+9.0*fr[1])/(20.78460969082652*EPSILON+fr[14]+1.732050807568877*((-1.0*fr[10])+fr[9]-1.0*fr[7])+3.0*((-1.0*fr[4])+fr[3]-1.0*fr[2])+5.196152422706631*fr[0]); 
  rCtrlR[3] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[13]+fr[12]-1.0*fr[11])+5.196152422706631*fr[8]-1.0*(5.196152422706631*(fr[6]+fr[5])+9.0*fr[1])))/(20.78460969082652*EPSILON-1.0*fr[14]+1.732050807568877*(fr[7]-1.0*(fr[10]+fr[9]))+3.0*((-1.0*fr[4])+fr[3]+fr[2])+5.196152422706631*fr[0]); 
  rCtrlR[4] = (1.732050807568877*fr[15]+3.0*(fr[11]-1.0*(fr[13]+fr[12]))+5.196152422706631*(fr[8]-1.0*(fr[6]+fr[5]))+9.0*fr[1])/(20.78460969082652*EPSILON+fr[14]+1.732050807568877*(fr[7]-1.0*(fr[10]+fr[9]))+3.0*(fr[4]-1.0*(fr[3]+fr[2]))+5.196152422706631*fr[0]); 
  rCtrlR[5] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[13]-1.0*fr[12]+fr[11])+5.196152422706631*(fr[6]-1.0*fr[8])-1.0*(5.196152422706631*fr[5]+9.0*fr[1])))/(20.78460969082652*EPSILON-1.0*fr[14]+1.732050807568877*((-1.0*fr[10])+fr[9]-1.0*fr[7])+3.0*(fr[4]-1.0*fr[3]+fr[2])+5.196152422706631*fr[0]); 
  rCtrlR[6] = -(1.0*(1.732050807568877*fr[15]+3.0*((-1.0*fr[13])+fr[12]+fr[11])+5.196152422706631*(fr[5]-1.0*(fr[8]+fr[6]))-9.0*fr[1]))/(20.78460969082652*EPSILON-1.0*fr[14]+1.732050807568877*(fr[10]-1.0*(fr[9]+fr[7]))+3.0*(fr[4]+fr[3]-1.0*fr[2])+5.196152422706631*fr[0]); 
  rCtrlR[7] = (1.732050807568877*fr[15]+3.0*(fr[13]+fr[12]+fr[11])+5.196152422706631*(fr[8]+fr[6]+fr[5])+9.0*fr[1])/(20.78460969082652*EPSILON+fr[14]+1.732050807568877*(fr[10]+fr[9]+fr[7])+3.0*(fr[4]+fr[3]+fr[2])+5.196152422706631*fr[0]); 
  double fCtrlL[8], fCtrlR[8];  // fCtrl = anti-limited f evaluated at each control node on x surface 
  // control node [y,vx,vy] = [-1/3,-1/3,-1/3] 
  fCtrlL[0] = -0.04811252243246882*(fl[14]-1.732050807568877*(fl[10]+fl[9]+fl[7])+3.0*(fl[4]+fl[3]+fl[2])-5.196152422706631*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.04811252243246882*(fr[14]-1.732050807568877*(fr[10]+fr[9]+fr[7])+3.0*(fr[4]+fr[3]+fr[2])-5.196152422706631*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [y,vx,vy] = [1/3,-1/3,-1/3] 
  fCtrlL[1] = 0.04811252243246882*(fl[14]+1.732050807568877*fl[10]-1.732050807568877*(fl[9]+fl[7])-3.0*(fl[4]+fl[3])+3.0*fl[2]+5.196152422706631*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.04811252243246882*(fr[14]+1.732050807568877*fr[10]-1.732050807568877*(fr[9]+fr[7])-3.0*(fr[4]+fr[3])+3.0*fr[2]+5.196152422706631*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [y,vx,vy] = [-1/3,1/3,-1/3] 
  fCtrlL[2] = 0.04811252243246882*(fl[14]-1.732050807568877*fl[10]+1.732050807568877*fl[9]-1.732050807568877*fl[7]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2]+5.196152422706631*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = 0.04811252243246882*(fr[14]-1.732050807568877*fr[10]+1.732050807568877*fr[9]-1.732050807568877*fr[7]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2]+5.196152422706631*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [y,vx,vy] = [1/3,1/3,-1/3] 
  fCtrlL[3] = -0.04811252243246882*(fl[14]+1.732050807568877*(fl[10]+fl[9])-1.732050807568877*fl[7]+3.0*fl[4]-3.0*(fl[3]+fl[2])-5.196152422706631*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = -0.04811252243246882*(fr[14]+1.732050807568877*(fr[10]+fr[9])-1.732050807568877*fr[7]+3.0*fr[4]-3.0*(fr[3]+fr[2])-5.196152422706631*fr[0])*limTheta(rCtrlR[3],-1.0); 
  // control node [y,vx,vy] = [-1/3,-1/3,1/3] 
  fCtrlL[4] = 0.04811252243246882*(fl[14]-1.732050807568877*(fl[10]+fl[9])+1.732050807568877*fl[7]+3.0*fl[4]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rCtrlL[4],1.0); 
  fCtrlR[4] = 0.04811252243246882*(fr[14]-1.732050807568877*(fr[10]+fr[9])+1.732050807568877*fr[7]+3.0*fr[4]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rCtrlR[4],-1.0); 
  // control node [y,vx,vy] = [1/3,-1/3,1/3] 
  fCtrlL[5] = -0.04811252243246882*(fl[14]+1.732050807568877*fl[10]-1.732050807568877*fl[9]+1.732050807568877*fl[7]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2]-5.196152422706631*fl[0])*limTheta(rCtrlL[5],1.0); 
  fCtrlR[5] = -0.04811252243246882*(fr[14]+1.732050807568877*fr[10]-1.732050807568877*fr[9]+1.732050807568877*fr[7]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2]-5.196152422706631*fr[0])*limTheta(rCtrlR[5],-1.0); 
  // control node [y,vx,vy] = [-1/3,1/3,1/3] 
  fCtrlL[6] = -0.04811252243246882*(fl[14]-1.732050807568877*fl[10]+1.732050807568877*(fl[9]+fl[7])-3.0*(fl[4]+fl[3])+3.0*fl[2]-5.196152422706631*fl[0])*limTheta(rCtrlL[6],1.0); 
  fCtrlR[6] = -0.04811252243246882*(fr[14]-1.732050807568877*fr[10]+1.732050807568877*(fr[9]+fr[7])-3.0*(fr[4]+fr[3])+3.0*fr[2]-5.196152422706631*fr[0])*limTheta(rCtrlR[6],-1.0); 
  // control node [y,vx,vy] = [1/3,1/3,1/3] 
  fCtrlL[7] = 0.04811252243246882*(fl[14]+1.732050807568877*(fl[10]+fl[9]+fl[7])+3.0*(fl[4]+fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rCtrlL[7],1.0); 
  fCtrlR[7] = 0.04811252243246882*(fr[14]+1.732050807568877*(fr[10]+fr[9]+fr[7])+3.0*(fr[4]+fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rCtrlR[7],-1.0); 
  double fL_AL[8], fR_AL[8];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.3535533905932737*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.6123724356957944*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*fCtrlL[4]+fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4])+fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]-1.0*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0])); 
  fL_AL[4] = 1.060660171779821*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]+fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fL_AL[5] = 1.060660171779821*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*(fCtrlL[4]+fCtrlL[3])+fCtrlL[2]-1.0*fCtrlL[1]+fCtrlL[0]); 
  fL_AL[6] = 1.060660171779821*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2])+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[7] = 1.837117307087383*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]-1.0*fCtrlL[3]+fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fR_AL[0] = 0.3535533905932737*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.6123724356957944*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*fCtrlR[4]+fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4])+fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]-1.0*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0])); 
  fR_AL[4] = 1.060660171779821*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]+fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  fR_AL[5] = 1.060660171779821*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*(fCtrlR[4]+fCtrlR[3])+fCtrlR[2]-1.0*fCtrlR[1]+fCtrlR[0]); 
  fR_AL[6] = 1.060660171779821*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2])+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[7] = 1.837117307087383*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]-1.0*fCtrlR[3]+fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[8], fhatAL[8]; 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double GhatCtrl[8];
  // control node [y,vpar,mu] = [-1/3,-1/3,-1/3] 
  GhatCtrl[0] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[0]) < EPSILON) GhatCtrl[0] = 0.; 
  // control node [y,vpar,mu] = [1/3,-1/3,-1/3] 
  GhatCtrl[1] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[1]) < EPSILON) GhatCtrl[1] = 0.; 
  // control node [y,vpar,mu] = [-1/3,1/3,-1/3] 
  GhatCtrl[2] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[2]) < EPSILON) GhatCtrl[2] = 0.; 
  // control node [y,vpar,mu] = [1/3,1/3,-1/3] 
  GhatCtrl[3] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[3]) < EPSILON) GhatCtrl[3] = 0.; 
  // control node [y,vpar,mu] = [-1/3,-1/3,1/3] 
  GhatCtrl[4] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[4]) < EPSILON) GhatCtrl[4] = 0.; 
  // control node [y,vpar,mu] = [1/3,-1/3,1/3] 
  GhatCtrl[5] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[5]) < EPSILON) GhatCtrl[5] = 0.; 
  // control node [y,vpar,mu] = [-1/3,1/3,1/3] 
  GhatCtrl[6] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[6]) < EPSILON) GhatCtrl[6] = 0.; 
  // control node [y,vpar,mu] = [1/3,1/3,1/3] 
  GhatCtrl[7] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[7]) < EPSILON) GhatCtrl[7] = 0.; 
  double uFrac = 0., fCtrl = 0., alphaCtrl = 0.;
  if(GhatCtrl[0]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[0]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*(fr[14]+fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.1443375672974065*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = -std::min(std::abs(GhatCtrl[0]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[0]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[1]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*fl[14]+0.04811252243246882*(fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9])-0.08333333333333333*fl[8]+0.08333333333333333*fl[7]-0.08333333333333333*(fl[6]+fl[5])-0.1443375672974065*(fl[4]+fl[3]+fl[2])+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = std::min(GhatCtrl[0], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[2]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*fr[14]-0.04811252243246882*fr[13]+0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*fr[9]+0.08333333333333333*fr[8]-0.08333333333333333*fr[7]+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]-0.1443375672974065*(fr[4]+fr[3])+0.1443375672974065*fr[2]-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = -std::min(std::abs(GhatCtrl[1]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[3]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*(fl[14]+fl[13])-0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*(fl[9]+fl[8]+fl[7]+fl[6])+0.08333333333333333*fl[5]-0.1443375672974065*(fl[4]+fl[3])+0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = std::min(GhatCtrl[1], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[4]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*(fr[14]+fr[13])-0.04811252243246882*fr[12]+0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*(fr[9]+fr[8])-0.08333333333333333*(fr[7]+fr[6])+0.08333333333333333*fr[5]-0.1443375672974065*fr[4]+0.1443375672974065*fr[3]-0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = -std::min(std::abs(GhatCtrl[2]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[5]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*fl[14]-0.04811252243246882*fl[13]+0.04811252243246882*fl[12]-0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*fl[9]-0.08333333333333333*(fl[8]+fl[7])+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]-0.1443375672974065*fl[4]+0.1443375672974065*fl[3]-0.1443375672974065*fl[2]+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = std::min(GhatCtrl[2], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[6]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*fr[14]+0.04811252243246882*(fr[13]+fr[12])-0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9])+0.08333333333333333*(fr[8]+fr[7])-0.08333333333333333*(fr[6]+fr[5])-0.1443375672974065*fr[4]+0.1443375672974065*(fr[3]+fr[2])-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = -std::min(std::abs(GhatCtrl[3]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[7]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*(fl[14]+fl[13]+fl[12])+0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9]+fl[8])+0.08333333333333333*(fl[7]+fl[6]+fl[5])-0.1443375672974065*fl[4]+0.1443375672974065*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = std::min(GhatCtrl[3], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[8]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*(fr[14]+fr[13]+fr[12])-0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9]+fr[8])+0.08333333333333333*(fr[7]+fr[6]+fr[5])+0.1443375672974065*fr[4]-0.1443375672974065*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = -std::min(std::abs(GhatCtrl[4]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[9]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*fl[14]-0.04811252243246882*(fl[13]+fl[12])+0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9])+0.08333333333333333*(fl[8]+fl[7])-0.08333333333333333*(fl[6]+fl[5])+0.1443375672974065*fl[4]-0.1443375672974065*(fl[3]+fl[2])+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = std::min(GhatCtrl[4], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[10]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*fr[14]+0.04811252243246882*fr[13]-0.04811252243246882*fr[12]+0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*fr[9]-0.08333333333333333*(fr[8]+fr[7])+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]+0.1443375672974065*fr[4]-0.1443375672974065*fr[3]+0.1443375672974065*fr[2]-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = -std::min(std::abs(GhatCtrl[5]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[11]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*(fl[14]+fl[13])+0.04811252243246882*fl[12]-0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*(fl[9]+fl[8])-0.08333333333333333*(fl[7]+fl[6])+0.08333333333333333*fl[5]+0.1443375672974065*fl[4]-0.1443375672974065*fl[3]+0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = std::min(GhatCtrl[5], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[12]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*(fr[14]+fr[13])+0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*(fr[9]+fr[8]+fr[7]+fr[6])+0.08333333333333333*fr[5]+0.1443375672974065*(fr[4]+fr[3])-0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = -std::min(std::abs(GhatCtrl[6]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[13]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*fl[14]+0.04811252243246882*fl[13]-0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*fl[9]+0.08333333333333333*fl[8]-0.08333333333333333*fl[7]+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]+0.1443375672974065*(fl[4]+fl[3])-0.1443375672974065*fl[2]+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = std::min(GhatCtrl[6], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[14]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*fr[14]-0.04811252243246882*(fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9])-0.08333333333333333*fr[8]+0.08333333333333333*fr[7]-0.08333333333333333*(fr[6]+fr[5])+0.1443375672974065*(fr[4]+fr[3]+fr[2])-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = -std::min(std::abs(GhatCtrl[7]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[15]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*(fl[14]+fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+0.1443375672974065*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = std::min(GhatCtrl[7], std::abs(uFrac)*fCtrl/dtApprox);
  }

  incr[0] = 0.25*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[1] = -0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[2] = 0.4330127018922193*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[3] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[4] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[5] = -0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[6] = -0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[7] = 0.75*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[8] = -0.75*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[9] = 0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[10] = 0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[11] = -1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[12] = -1.299038105676658*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[13] = -1.299038105676658*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[14] = 1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[15] = -2.25*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 

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
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Y_P1_Bvars_0(const double q_, const double m_, const double *cflRateCtrlL, const double *cflRateCtrlR, const double *w, const double *dxv, const double dtApprox, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*(2.0*BdriftY[0]*m_*wv2+BmagInv[0]*(1.732050807568877*Phi[1]-3.0*Phi[3])*dfac_x*q_))/q_; 

  double alpha[8]; 
  alpha[0] = (0.5*(2.828427124746191*BdriftY[0]*m_*wv2+BmagInv[0]*(2.449489742783178*Phi[1]-4.242640687119286*Phi[3])*dfac_x*q_))/q_; 
  alpha[2] = (0.8164965809277261*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  double rCtrlL[8], rCtrlR[8];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrlL[0] = -(1.0*(1.732050807568877*fl[15]-3.0*(fl[14]+fl[12]+fl[11])+5.196152422706631*(fl[9]+fl[7]+fl[5])-9.0*fl[2]))/(20.78460969082652*EPSILON-1.0*fl[13]+1.732050807568877*(fl[10]+fl[8]+fl[6])-3.0*(fl[4]+fl[3]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[1] = (1.732050807568877*fl[15]+3.0*(fl[14]-1.0*(fl[12]+fl[11]))+5.196152422706631*(fl[5]-1.0*(fl[9]+fl[7]))+9.0*fl[2])/(20.78460969082652*EPSILON+fl[13]+1.732050807568877*(fl[10]-1.0*(fl[8]+fl[6]))+3.0*(fl[1]-1.0*(fl[4]+fl[3]))+5.196152422706631*fl[0]); 
  rCtrlL[2] = (1.732050807568877*fl[15]+3.0*((-1.0*fl[14])+fl[12]-1.0*fl[11])+5.196152422706631*((-1.0*fl[9])+fl[7]-1.0*fl[5])+9.0*fl[2])/(20.78460969082652*EPSILON+fl[13]+1.732050807568877*((-1.0*fl[10])+fl[8]-1.0*fl[6])+3.0*((-1.0*fl[4])+fl[3]-1.0*fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[3] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[14]+fl[12]-1.0*fl[11])+5.196152422706631*fl[9]-1.0*(5.196152422706631*(fl[7]+fl[5])+9.0*fl[2])))/(20.78460969082652*EPSILON-1.0*fl[13]+1.732050807568877*(fl[6]-1.0*(fl[10]+fl[8]))+3.0*((-1.0*fl[4])+fl[3]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[4] = (1.732050807568877*fl[15]+3.0*(fl[11]-1.0*(fl[14]+fl[12]))+5.196152422706631*(fl[9]-1.0*(fl[7]+fl[5]))+9.0*fl[2])/(20.78460969082652*EPSILON+fl[13]+1.732050807568877*(fl[6]-1.0*(fl[10]+fl[8]))+3.0*(fl[4]-1.0*(fl[3]+fl[1]))+5.196152422706631*fl[0]); 
  rCtrlL[5] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[14]-1.0*fl[12]+fl[11])+5.196152422706631*(fl[7]-1.0*fl[9])-1.0*(5.196152422706631*fl[5]+9.0*fl[2])))/(20.78460969082652*EPSILON-1.0*fl[13]+1.732050807568877*((-1.0*fl[10])+fl[8]-1.0*fl[6])+3.0*(fl[4]-1.0*fl[3]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[6] = -(1.0*(1.732050807568877*fl[15]+3.0*((-1.0*fl[14])+fl[12]+fl[11])+5.196152422706631*(fl[5]-1.0*(fl[9]+fl[7]))-9.0*fl[2]))/(20.78460969082652*EPSILON-1.0*fl[13]+1.732050807568877*(fl[10]-1.0*(fl[8]+fl[6]))+3.0*(fl[4]+fl[3]-1.0*fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[7] = (1.732050807568877*fl[15]+3.0*(fl[14]+fl[12]+fl[11])+5.196152422706631*(fl[9]+fl[7]+fl[5])+9.0*fl[2])/(20.78460969082652*EPSILON+fl[13]+1.732050807568877*(fl[10]+fl[8]+fl[6])+3.0*(fl[4]+fl[3]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlR[0] = -(1.0*(1.732050807568877*fr[15]-3.0*(fr[14]+fr[12]+fr[11])+5.196152422706631*(fr[9]+fr[7]+fr[5])-9.0*fr[2]))/(20.78460969082652*EPSILON-1.0*fr[13]+1.732050807568877*(fr[10]+fr[8]+fr[6])-3.0*(fr[4]+fr[3]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[1] = (1.732050807568877*fr[15]+3.0*(fr[14]-1.0*(fr[12]+fr[11]))+5.196152422706631*(fr[5]-1.0*(fr[9]+fr[7]))+9.0*fr[2])/(20.78460969082652*EPSILON+fr[13]+1.732050807568877*(fr[10]-1.0*(fr[8]+fr[6]))+3.0*(fr[1]-1.0*(fr[4]+fr[3]))+5.196152422706631*fr[0]); 
  rCtrlR[2] = (1.732050807568877*fr[15]+3.0*((-1.0*fr[14])+fr[12]-1.0*fr[11])+5.196152422706631*((-1.0*fr[9])+fr[7]-1.0*fr[5])+9.0*fr[2])/(20.78460969082652*EPSILON+fr[13]+1.732050807568877*((-1.0*fr[10])+fr[8]-1.0*fr[6])+3.0*((-1.0*fr[4])+fr[3]-1.0*fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[3] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[14]+fr[12]-1.0*fr[11])+5.196152422706631*fr[9]-1.0*(5.196152422706631*(fr[7]+fr[5])+9.0*fr[2])))/(20.78460969082652*EPSILON-1.0*fr[13]+1.732050807568877*(fr[6]-1.0*(fr[10]+fr[8]))+3.0*((-1.0*fr[4])+fr[3]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[4] = (1.732050807568877*fr[15]+3.0*(fr[11]-1.0*(fr[14]+fr[12]))+5.196152422706631*(fr[9]-1.0*(fr[7]+fr[5]))+9.0*fr[2])/(20.78460969082652*EPSILON+fr[13]+1.732050807568877*(fr[6]-1.0*(fr[10]+fr[8]))+3.0*(fr[4]-1.0*(fr[3]+fr[1]))+5.196152422706631*fr[0]); 
  rCtrlR[5] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[14]-1.0*fr[12]+fr[11])+5.196152422706631*(fr[7]-1.0*fr[9])-1.0*(5.196152422706631*fr[5]+9.0*fr[2])))/(20.78460969082652*EPSILON-1.0*fr[13]+1.732050807568877*((-1.0*fr[10])+fr[8]-1.0*fr[6])+3.0*(fr[4]-1.0*fr[3]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[6] = -(1.0*(1.732050807568877*fr[15]+3.0*((-1.0*fr[14])+fr[12]+fr[11])+5.196152422706631*(fr[5]-1.0*(fr[9]+fr[7]))-9.0*fr[2]))/(20.78460969082652*EPSILON-1.0*fr[13]+1.732050807568877*(fr[10]-1.0*(fr[8]+fr[6]))+3.0*(fr[4]+fr[3]-1.0*fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[7] = (1.732050807568877*fr[15]+3.0*(fr[14]+fr[12]+fr[11])+5.196152422706631*(fr[9]+fr[7]+fr[5])+9.0*fr[2])/(20.78460969082652*EPSILON+fr[13]+1.732050807568877*(fr[10]+fr[8]+fr[6])+3.0*(fr[4]+fr[3]+fr[1])+5.196152422706631*fr[0]); 
  double fCtrlL[8], fCtrlR[8];  // fCtrl = anti-limited f evaluated at each control node on y surface 
  // control node [x,vx,vy] = [-1/3,-1/3,-1/3] 
  fCtrlL[0] = -0.04811252243246882*(fl[13]-1.732050807568877*(fl[10]+fl[8]+fl[6])+3.0*(fl[4]+fl[3]+fl[1])-5.196152422706631*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.04811252243246882*(fr[13]-1.732050807568877*(fr[10]+fr[8]+fr[6])+3.0*(fr[4]+fr[3]+fr[1])-5.196152422706631*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,vx,vy] = [1/3,-1/3,-1/3] 
  fCtrlL[1] = 0.04811252243246882*(fl[13]+1.732050807568877*fl[10]-1.732050807568877*(fl[8]+fl[6])-3.0*(fl[4]+fl[3])+3.0*fl[1]+5.196152422706631*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.04811252243246882*(fr[13]+1.732050807568877*fr[10]-1.732050807568877*(fr[8]+fr[6])-3.0*(fr[4]+fr[3])+3.0*fr[1]+5.196152422706631*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,vx,vy] = [-1/3,1/3,-1/3] 
  fCtrlL[2] = 0.04811252243246882*(fl[13]-1.732050807568877*fl[10]+1.732050807568877*fl[8]-1.732050807568877*fl[6]-3.0*fl[4]+3.0*fl[3]-3.0*fl[1]+5.196152422706631*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = 0.04811252243246882*(fr[13]-1.732050807568877*fr[10]+1.732050807568877*fr[8]-1.732050807568877*fr[6]-3.0*fr[4]+3.0*fr[3]-3.0*fr[1]+5.196152422706631*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,vx,vy] = [1/3,1/3,-1/3] 
  fCtrlL[3] = -0.04811252243246882*(fl[13]+1.732050807568877*(fl[10]+fl[8])-1.732050807568877*fl[6]+3.0*fl[4]-3.0*(fl[3]+fl[1])-5.196152422706631*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = -0.04811252243246882*(fr[13]+1.732050807568877*(fr[10]+fr[8])-1.732050807568877*fr[6]+3.0*fr[4]-3.0*(fr[3]+fr[1])-5.196152422706631*fr[0])*limTheta(rCtrlR[3],-1.0); 
  // control node [x,vx,vy] = [-1/3,-1/3,1/3] 
  fCtrlL[4] = 0.04811252243246882*(fl[13]-1.732050807568877*(fl[10]+fl[8])+1.732050807568877*fl[6]+3.0*fl[4]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rCtrlL[4],1.0); 
  fCtrlR[4] = 0.04811252243246882*(fr[13]-1.732050807568877*(fr[10]+fr[8])+1.732050807568877*fr[6]+3.0*fr[4]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rCtrlR[4],-1.0); 
  // control node [x,vx,vy] = [1/3,-1/3,1/3] 
  fCtrlL[5] = -0.04811252243246882*(fl[13]+1.732050807568877*fl[10]-1.732050807568877*fl[8]+1.732050807568877*fl[6]-3.0*fl[4]+3.0*fl[3]-3.0*fl[1]-5.196152422706631*fl[0])*limTheta(rCtrlL[5],1.0); 
  fCtrlR[5] = -0.04811252243246882*(fr[13]+1.732050807568877*fr[10]-1.732050807568877*fr[8]+1.732050807568877*fr[6]-3.0*fr[4]+3.0*fr[3]-3.0*fr[1]-5.196152422706631*fr[0])*limTheta(rCtrlR[5],-1.0); 
  // control node [x,vx,vy] = [-1/3,1/3,1/3] 
  fCtrlL[6] = -0.04811252243246882*(fl[13]-1.732050807568877*fl[10]+1.732050807568877*(fl[8]+fl[6])-3.0*(fl[4]+fl[3])+3.0*fl[1]-5.196152422706631*fl[0])*limTheta(rCtrlL[6],1.0); 
  fCtrlR[6] = -0.04811252243246882*(fr[13]-1.732050807568877*fr[10]+1.732050807568877*(fr[8]+fr[6])-3.0*(fr[4]+fr[3])+3.0*fr[1]-5.196152422706631*fr[0])*limTheta(rCtrlR[6],-1.0); 
  // control node [x,vx,vy] = [1/3,1/3,1/3] 
  fCtrlL[7] = 0.04811252243246882*(fl[13]+1.732050807568877*(fl[10]+fl[8]+fl[6])+3.0*(fl[4]+fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rCtrlL[7],1.0); 
  fCtrlR[7] = 0.04811252243246882*(fr[13]+1.732050807568877*(fr[10]+fr[8]+fr[6])+3.0*(fr[4]+fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rCtrlR[7],-1.0); 
  double fL_AL[8], fR_AL[8];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.3535533905932737*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.6123724356957944*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*fCtrlL[4]+fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4])+fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]-1.0*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0])); 
  fL_AL[4] = 1.060660171779821*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]+fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fL_AL[5] = 1.060660171779821*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*(fCtrlL[4]+fCtrlL[3])+fCtrlL[2]-1.0*fCtrlL[1]+fCtrlL[0]); 
  fL_AL[6] = 1.060660171779821*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2])+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[7] = 1.837117307087383*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]-1.0*fCtrlL[3]+fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fR_AL[0] = 0.3535533905932737*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.6123724356957944*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*fCtrlR[4]+fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4])+fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]-1.0*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0])); 
  fR_AL[4] = 1.060660171779821*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]+fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  fR_AL[5] = 1.060660171779821*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*(fCtrlR[4]+fCtrlR[3])+fCtrlR[2]-1.0*fCtrlR[1]+fCtrlR[0]); 
  fR_AL[6] = 1.060660171779821*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2])+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[7] = 1.837117307087383*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]-1.0*fCtrlR[3]+fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[8], fhatAL[8]; 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double GhatCtrl[8];
  // control node [x,vpar,mu] = [-1/3,-1/3,-1/3] 
  GhatCtrl[0] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_y; 
  if(std::abs(GhatCtrl[0]) < EPSILON) GhatCtrl[0] = 0.; 
  // control node [x,vpar,mu] = [1/3,-1/3,-1/3] 
  GhatCtrl[1] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_y; 
  if(std::abs(GhatCtrl[1]) < EPSILON) GhatCtrl[1] = 0.; 
  // control node [x,vpar,mu] = [-1/3,1/3,-1/3] 
  GhatCtrl[2] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_y; 
  if(std::abs(GhatCtrl[2]) < EPSILON) GhatCtrl[2] = 0.; 
  // control node [x,vpar,mu] = [1/3,1/3,-1/3] 
  GhatCtrl[3] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_y; 
  if(std::abs(GhatCtrl[3]) < EPSILON) GhatCtrl[3] = 0.; 
  // control node [x,vpar,mu] = [-1/3,-1/3,1/3] 
  GhatCtrl[4] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_y; 
  if(std::abs(GhatCtrl[4]) < EPSILON) GhatCtrl[4] = 0.; 
  // control node [x,vpar,mu] = [1/3,-1/3,1/3] 
  GhatCtrl[5] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_y; 
  if(std::abs(GhatCtrl[5]) < EPSILON) GhatCtrl[5] = 0.; 
  // control node [x,vpar,mu] = [-1/3,1/3,1/3] 
  GhatCtrl[6] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_y; 
  if(std::abs(GhatCtrl[6]) < EPSILON) GhatCtrl[6] = 0.; 
  // control node [x,vpar,mu] = [1/3,1/3,1/3] 
  GhatCtrl[7] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_y; 
  if(std::abs(GhatCtrl[7]) < EPSILON) GhatCtrl[7] = 0.; 
  double uFrac = 0., fCtrl = 0., alphaCtrl = 0.;
  if(GhatCtrl[0]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[0]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*(fr[14]+fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.1443375672974065*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = -std::min(std::abs(GhatCtrl[0]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[1]; 
    fCtrl = (-0.02777777777777778*fr[15])-0.04811252243246882*fr[14]+0.04811252243246882*(fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9])-0.08333333333333333*fr[8]+0.08333333333333333*fr[7]-0.08333333333333333*(fr[6]+fr[5])-0.1443375672974065*(fr[4]+fr[3]+fr[2])+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = -std::min(std::abs(GhatCtrl[1]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[0]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[2]; 
    fCtrl = (-0.02777777777777778*fl[15])+0.04811252243246882*fl[14]-0.04811252243246882*fl[13]+0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*fl[9]+0.08333333333333333*fl[8]-0.08333333333333333*fl[7]+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]-0.1443375672974065*(fl[4]+fl[3])+0.1443375672974065*fl[2]-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = std::min(GhatCtrl[0], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[3]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*(fl[14]+fl[13])-0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*(fl[9]+fl[8]+fl[7]+fl[6])+0.08333333333333333*fl[5]-0.1443375672974065*(fl[4]+fl[3])+0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = std::min(GhatCtrl[1], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[4]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*(fr[14]+fr[13])-0.04811252243246882*fr[12]+0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*(fr[9]+fr[8])-0.08333333333333333*(fr[7]+fr[6])+0.08333333333333333*fr[5]-0.1443375672974065*fr[4]+0.1443375672974065*fr[3]-0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = -std::min(std::abs(GhatCtrl[2]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[5]; 
    fCtrl = 0.02777777777777778*fr[15]+0.04811252243246882*fr[14]-0.04811252243246882*fr[13]+0.04811252243246882*fr[12]-0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*fr[9]-0.08333333333333333*(fr[8]+fr[7])+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]-0.1443375672974065*fr[4]+0.1443375672974065*fr[3]-0.1443375672974065*fr[2]+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = -std::min(std::abs(GhatCtrl[3]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[6]; 
    fCtrl = 0.02777777777777778*fl[15]-0.04811252243246882*fl[14]+0.04811252243246882*(fl[13]+fl[12])-0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9])+0.08333333333333333*(fl[8]+fl[7])-0.08333333333333333*(fl[6]+fl[5])-0.1443375672974065*fl[4]+0.1443375672974065*(fl[3]+fl[2])-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = std::min(GhatCtrl[2], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[7]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*(fl[14]+fl[13]+fl[12])+0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9]+fl[8])+0.08333333333333333*(fl[7]+fl[6]+fl[5])-0.1443375672974065*fl[4]+0.1443375672974065*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = std::min(GhatCtrl[3], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[8]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*(fr[14]+fr[13]+fr[12])-0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9]+fr[8])+0.08333333333333333*(fr[7]+fr[6]+fr[5])+0.1443375672974065*fr[4]-0.1443375672974065*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = -std::min(std::abs(GhatCtrl[4]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[9]; 
    fCtrl = 0.02777777777777778*fr[15]+0.04811252243246882*fr[14]-0.04811252243246882*(fr[13]+fr[12])+0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9])+0.08333333333333333*(fr[8]+fr[7])-0.08333333333333333*(fr[6]+fr[5])+0.1443375672974065*fr[4]-0.1443375672974065*(fr[3]+fr[2])+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = -std::min(std::abs(GhatCtrl[5]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[10]; 
    fCtrl = 0.02777777777777778*fl[15]-0.04811252243246882*fl[14]+0.04811252243246882*fl[13]-0.04811252243246882*fl[12]+0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*fl[9]-0.08333333333333333*(fl[8]+fl[7])+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]+0.1443375672974065*fl[4]-0.1443375672974065*fl[3]+0.1443375672974065*fl[2]-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = std::min(GhatCtrl[4], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[11]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*(fl[14]+fl[13])+0.04811252243246882*fl[12]-0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*(fl[9]+fl[8])-0.08333333333333333*(fl[7]+fl[6])+0.08333333333333333*fl[5]+0.1443375672974065*fl[4]-0.1443375672974065*fl[3]+0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = std::min(GhatCtrl[5], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[12]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*(fr[14]+fr[13])+0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*(fr[9]+fr[8]+fr[7]+fr[6])+0.08333333333333333*fr[5]+0.1443375672974065*(fr[4]+fr[3])-0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = -std::min(std::abs(GhatCtrl[6]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[13]; 
    fCtrl = (-0.02777777777777778*fr[15])-0.04811252243246882*fr[14]+0.04811252243246882*fr[13]-0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*fr[9]+0.08333333333333333*fr[8]-0.08333333333333333*fr[7]+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]+0.1443375672974065*(fr[4]+fr[3])-0.1443375672974065*fr[2]+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = -std::min(std::abs(GhatCtrl[7]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[14]; 
    fCtrl = (-0.02777777777777778*fl[15])+0.04811252243246882*fl[14]-0.04811252243246882*(fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9])-0.08333333333333333*fl[8]+0.08333333333333333*fl[7]-0.08333333333333333*(fl[6]+fl[5])+0.1443375672974065*(fl[4]+fl[3]+fl[2])-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = std::min(GhatCtrl[6], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[15]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*(fl[14]+fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+0.1443375672974065*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = std::min(GhatCtrl[7], std::abs(uFrac)*fCtrl/dtApprox);
  }

  incr[0] = 0.25*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[1] = 0.4330127018922193*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[2] = -0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[3] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[4] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[5] = -0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[6] = 0.75*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[7] = -0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[8] = 0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[9] = -0.75*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[10] = 0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[11] = -1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[12] = -1.299038105676658*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[13] = 1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[14] = -1.299038105676658*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[15] = -2.25*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 

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
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double *cflRateCtrlL, const double *cflRateCtrlR, const double *w, const double *dxv, const double dtApprox, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.2165063509461096*(dfac_v*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*wv-1.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)))/dfac_v; 

  double alpha[8]; 
  alpha[0] = -(0.7071067811865475*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*(1.732050807568877*dfac_v*wv-1.732050807568877))/dfac_v; 
  alpha[1] = -(0.7071067811865475*BdriftY[0]*Phi[3]*dfac_y*(1.732050807568877*dfac_v*wv-1.732050807568877))/dfac_v; 
  alpha[2] = -(0.7071067811865475*BdriftX[0]*Phi[3]*dfac_x*(1.732050807568877*dfac_v*wv-1.732050807568877))/dfac_v; 
  double rCtrlL[8], rCtrlR[8];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrlL[0] = -(1.0*(1.732050807568877*fl[15]-3.0*(fl[14]+fl[13]+fl[11])+5.196152422706631*(fl[10]+fl[7]+fl[6])-9.0*fl[3]))/(20.78460969082652*EPSILON-1.0*fl[12]+1.732050807568877*(fl[9]+fl[8]+fl[5])-3.0*(fl[4]+fl[2]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[1] = (1.732050807568877*fl[15]+3.0*(fl[14]-1.0*(fl[13]+fl[11]))+5.196152422706631*(fl[6]-1.0*(fl[10]+fl[7]))+9.0*fl[3])/(20.78460969082652*EPSILON+fl[12]+1.732050807568877*(fl[9]-1.0*(fl[8]+fl[5]))+3.0*(fl[1]-1.0*(fl[4]+fl[2]))+5.196152422706631*fl[0]); 
  rCtrlL[2] = (1.732050807568877*fl[15]+3.0*((-1.0*fl[14])+fl[13]-1.0*fl[11])+5.196152422706631*((-1.0*fl[10])+fl[7]-1.0*fl[6])+9.0*fl[3])/(20.78460969082652*EPSILON+fl[12]+1.732050807568877*((-1.0*fl[9])+fl[8]-1.0*fl[5])+3.0*((-1.0*fl[4])+fl[2]-1.0*fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[3] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[14]+fl[13]-1.0*fl[11])+5.196152422706631*fl[10]-1.0*(5.196152422706631*(fl[7]+fl[6])+9.0*fl[3])))/(20.78460969082652*EPSILON-1.0*fl[12]+1.732050807568877*(fl[5]-1.0*(fl[9]+fl[8]))+3.0*((-1.0*fl[4])+fl[2]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[4] = (1.732050807568877*fl[15]+3.0*(fl[11]-1.0*(fl[14]+fl[13]))+5.196152422706631*(fl[10]-1.0*(fl[7]+fl[6]))+9.0*fl[3])/(20.78460969082652*EPSILON+fl[12]+1.732050807568877*(fl[5]-1.0*(fl[9]+fl[8]))+3.0*(fl[4]-1.0*(fl[2]+fl[1]))+5.196152422706631*fl[0]); 
  rCtrlL[5] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[14]-1.0*fl[13]+fl[11])+5.196152422706631*(fl[7]-1.0*fl[10])-1.0*(5.196152422706631*fl[6]+9.0*fl[3])))/(20.78460969082652*EPSILON-1.0*fl[12]+1.732050807568877*((-1.0*fl[9])+fl[8]-1.0*fl[5])+3.0*(fl[4]-1.0*fl[2]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[6] = -(1.0*(1.732050807568877*fl[15]+3.0*((-1.0*fl[14])+fl[13]+fl[11])+5.196152422706631*(fl[6]-1.0*(fl[10]+fl[7]))-9.0*fl[3]))/(20.78460969082652*EPSILON-1.0*fl[12]+1.732050807568877*(fl[9]-1.0*(fl[8]+fl[5]))+3.0*(fl[4]+fl[2]-1.0*fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[7] = (1.732050807568877*fl[15]+3.0*(fl[14]+fl[13]+fl[11])+5.196152422706631*(fl[10]+fl[7]+fl[6])+9.0*fl[3])/(20.78460969082652*EPSILON+fl[12]+1.732050807568877*(fl[9]+fl[8]+fl[5])+3.0*(fl[4]+fl[2]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlR[0] = -(1.0*(1.732050807568877*fr[15]-3.0*(fr[14]+fr[13]+fr[11])+5.196152422706631*(fr[10]+fr[7]+fr[6])-9.0*fr[3]))/(20.78460969082652*EPSILON-1.0*fr[12]+1.732050807568877*(fr[9]+fr[8]+fr[5])-3.0*(fr[4]+fr[2]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[1] = (1.732050807568877*fr[15]+3.0*(fr[14]-1.0*(fr[13]+fr[11]))+5.196152422706631*(fr[6]-1.0*(fr[10]+fr[7]))+9.0*fr[3])/(20.78460969082652*EPSILON+fr[12]+1.732050807568877*(fr[9]-1.0*(fr[8]+fr[5]))+3.0*(fr[1]-1.0*(fr[4]+fr[2]))+5.196152422706631*fr[0]); 
  rCtrlR[2] = (1.732050807568877*fr[15]+3.0*((-1.0*fr[14])+fr[13]-1.0*fr[11])+5.196152422706631*((-1.0*fr[10])+fr[7]-1.0*fr[6])+9.0*fr[3])/(20.78460969082652*EPSILON+fr[12]+1.732050807568877*((-1.0*fr[9])+fr[8]-1.0*fr[5])+3.0*((-1.0*fr[4])+fr[2]-1.0*fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[3] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[14]+fr[13]-1.0*fr[11])+5.196152422706631*fr[10]-1.0*(5.196152422706631*(fr[7]+fr[6])+9.0*fr[3])))/(20.78460969082652*EPSILON-1.0*fr[12]+1.732050807568877*(fr[5]-1.0*(fr[9]+fr[8]))+3.0*((-1.0*fr[4])+fr[2]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[4] = (1.732050807568877*fr[15]+3.0*(fr[11]-1.0*(fr[14]+fr[13]))+5.196152422706631*(fr[10]-1.0*(fr[7]+fr[6]))+9.0*fr[3])/(20.78460969082652*EPSILON+fr[12]+1.732050807568877*(fr[5]-1.0*(fr[9]+fr[8]))+3.0*(fr[4]-1.0*(fr[2]+fr[1]))+5.196152422706631*fr[0]); 
  rCtrlR[5] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[14]-1.0*fr[13]+fr[11])+5.196152422706631*(fr[7]-1.0*fr[10])-1.0*(5.196152422706631*fr[6]+9.0*fr[3])))/(20.78460969082652*EPSILON-1.0*fr[12]+1.732050807568877*((-1.0*fr[9])+fr[8]-1.0*fr[5])+3.0*(fr[4]-1.0*fr[2]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[6] = -(1.0*(1.732050807568877*fr[15]+3.0*((-1.0*fr[14])+fr[13]+fr[11])+5.196152422706631*(fr[6]-1.0*(fr[10]+fr[7]))-9.0*fr[3]))/(20.78460969082652*EPSILON-1.0*fr[12]+1.732050807568877*(fr[9]-1.0*(fr[8]+fr[5]))+3.0*(fr[4]+fr[2]-1.0*fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[7] = (1.732050807568877*fr[15]+3.0*(fr[14]+fr[13]+fr[11])+5.196152422706631*(fr[10]+fr[7]+fr[6])+9.0*fr[3])/(20.78460969082652*EPSILON+fr[12]+1.732050807568877*(fr[9]+fr[8]+fr[5])+3.0*(fr[4]+fr[2]+fr[1])+5.196152422706631*fr[0]); 
  double fCtrlL[8], fCtrlR[8];  // fCtrl = anti-limited f evaluated at each control node on vx surface 
  // control node [x,y,vy] = [-1/3,-1/3,-1/3] 
  fCtrlL[0] = -0.04811252243246882*(fl[12]-1.732050807568877*(fl[9]+fl[8]+fl[5])+3.0*(fl[4]+fl[2]+fl[1])-5.196152422706631*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.04811252243246882*(fr[12]-1.732050807568877*(fr[9]+fr[8]+fr[5])+3.0*(fr[4]+fr[2]+fr[1])-5.196152422706631*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,y,vy] = [1/3,-1/3,-1/3] 
  fCtrlL[1] = 0.04811252243246882*(fl[12]+1.732050807568877*fl[9]-1.732050807568877*(fl[8]+fl[5])-3.0*(fl[4]+fl[2])+3.0*fl[1]+5.196152422706631*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.04811252243246882*(fr[12]+1.732050807568877*fr[9]-1.732050807568877*(fr[8]+fr[5])-3.0*(fr[4]+fr[2])+3.0*fr[1]+5.196152422706631*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,y,vy] = [-1/3,1/3,-1/3] 
  fCtrlL[2] = 0.04811252243246882*(fl[12]-1.732050807568877*fl[9]+1.732050807568877*fl[8]-1.732050807568877*fl[5]-3.0*fl[4]+3.0*fl[2]-3.0*fl[1]+5.196152422706631*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = 0.04811252243246882*(fr[12]-1.732050807568877*fr[9]+1.732050807568877*fr[8]-1.732050807568877*fr[5]-3.0*fr[4]+3.0*fr[2]-3.0*fr[1]+5.196152422706631*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,y,vy] = [1/3,1/3,-1/3] 
  fCtrlL[3] = -0.04811252243246882*(fl[12]+1.732050807568877*(fl[9]+fl[8])-1.732050807568877*fl[5]+3.0*fl[4]-3.0*(fl[2]+fl[1])-5.196152422706631*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = -0.04811252243246882*(fr[12]+1.732050807568877*(fr[9]+fr[8])-1.732050807568877*fr[5]+3.0*fr[4]-3.0*(fr[2]+fr[1])-5.196152422706631*fr[0])*limTheta(rCtrlR[3],-1.0); 
  // control node [x,y,vy] = [-1/3,-1/3,1/3] 
  fCtrlL[4] = 0.04811252243246882*(fl[12]-1.732050807568877*(fl[9]+fl[8])+1.732050807568877*fl[5]+3.0*fl[4]-3.0*(fl[2]+fl[1])+5.196152422706631*fl[0])*limTheta(rCtrlL[4],1.0); 
  fCtrlR[4] = 0.04811252243246882*(fr[12]-1.732050807568877*(fr[9]+fr[8])+1.732050807568877*fr[5]+3.0*fr[4]-3.0*(fr[2]+fr[1])+5.196152422706631*fr[0])*limTheta(rCtrlR[4],-1.0); 
  // control node [x,y,vy] = [1/3,-1/3,1/3] 
  fCtrlL[5] = -0.04811252243246882*(fl[12]+1.732050807568877*fl[9]-1.732050807568877*fl[8]+1.732050807568877*fl[5]-3.0*fl[4]+3.0*fl[2]-3.0*fl[1]-5.196152422706631*fl[0])*limTheta(rCtrlL[5],1.0); 
  fCtrlR[5] = -0.04811252243246882*(fr[12]+1.732050807568877*fr[9]-1.732050807568877*fr[8]+1.732050807568877*fr[5]-3.0*fr[4]+3.0*fr[2]-3.0*fr[1]-5.196152422706631*fr[0])*limTheta(rCtrlR[5],-1.0); 
  // control node [x,y,vy] = [-1/3,1/3,1/3] 
  fCtrlL[6] = -0.04811252243246882*(fl[12]-1.732050807568877*fl[9]+1.732050807568877*(fl[8]+fl[5])-3.0*(fl[4]+fl[2])+3.0*fl[1]-5.196152422706631*fl[0])*limTheta(rCtrlL[6],1.0); 
  fCtrlR[6] = -0.04811252243246882*(fr[12]-1.732050807568877*fr[9]+1.732050807568877*(fr[8]+fr[5])-3.0*(fr[4]+fr[2])+3.0*fr[1]-5.196152422706631*fr[0])*limTheta(rCtrlR[6],-1.0); 
  // control node [x,y,vy] = [1/3,1/3,1/3] 
  fCtrlL[7] = 0.04811252243246882*(fl[12]+1.732050807568877*(fl[9]+fl[8]+fl[5])+3.0*(fl[4]+fl[2]+fl[1])+5.196152422706631*fl[0])*limTheta(rCtrlL[7],1.0); 
  fCtrlR[7] = 0.04811252243246882*(fr[12]+1.732050807568877*(fr[9]+fr[8]+fr[5])+3.0*(fr[4]+fr[2]+fr[1])+5.196152422706631*fr[0])*limTheta(rCtrlR[7],-1.0); 
  double fL_AL[8], fR_AL[8];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.3535533905932737*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.6123724356957944*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*fCtrlL[4]+fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4])+fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]-1.0*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0])); 
  fL_AL[4] = 1.060660171779821*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]+fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fL_AL[5] = 1.060660171779821*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*(fCtrlL[4]+fCtrlL[3])+fCtrlL[2]-1.0*fCtrlL[1]+fCtrlL[0]); 
  fL_AL[6] = 1.060660171779821*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2])+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[7] = 1.837117307087383*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]-1.0*fCtrlL[3]+fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fR_AL[0] = 0.3535533905932737*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.6123724356957944*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*fCtrlR[4]+fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4])+fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]-1.0*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0])); 
  fR_AL[4] = 1.060660171779821*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]+fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  fR_AL[5] = 1.060660171779821*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*(fCtrlR[4]+fCtrlR[3])+fCtrlR[2]-1.0*fCtrlR[1]+fCtrlR[0]); 
  fR_AL[6] = 1.060660171779821*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2])+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[7] = 1.837117307087383*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]-1.0*fCtrlR[3]+fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[8], fhatAL[8]; 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*(alpha[2]+alpha[1]); 
  if(alphaQuad > 0) {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0])-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*(alpha[2]+alpha[1]); 
  if(alphaQuad > 0) {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0])-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double GhatCtrl[8];
  // control node [x,y,mu] = [-1/3,-1/3,-1/3] 
  GhatCtrl[0] = (0.04166666666666666*((alpha[2]+alpha[1])*fhatAL[7]+alpha[0]*(fhatAL[6]+fhatAL[5])+(alpha[2]+alpha[1])*fhatAL[3])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+(alpha[2]+alpha[1])*fhatAL[4]+alpha[0]*(fhatAL[3]+fhatAL[2])+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0]))*dfac_v; 
  if(std::abs(GhatCtrl[0]) < EPSILON) GhatCtrl[0] = 0.; 
  // control node [x,y,mu] = [1/3,-1/3,-1/3] 
  GhatCtrl[1] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]+alpha[1]*fhatAL[3]))+0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+alpha[0]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])-0.07216878364870323*(alpha[1]*fhatAL[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.04166666666666666*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0]))*dfac_v; 
  if(std::abs(GhatCtrl[1]) < EPSILON) GhatCtrl[1] = 0.; 
  // control node [x,y,mu] = [-1/3,1/3,-1/3] 
  GhatCtrl[2] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]+alpha[1]*fhatAL[3])-0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+alpha[2]*fhatAL[4]+alpha[0]*(fhatAL[3]+fhatAL[1])+fhatAL[0]*alpha[1])+0.07216878364870323*(alpha[1]*fhatAL[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.04166666666666666*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0]))*dfac_v; 
  if(std::abs(GhatCtrl[2]) < EPSILON) GhatCtrl[2] = 0.; 
  // control node [x,y,mu] = [1/3,1/3,-1/3] 
  GhatCtrl[3] = ((-0.04166666666666666*((alpha[2]+alpha[1])*fhatAL[7]+alpha[0]*(fhatAL[6]+fhatAL[5])+(alpha[2]+alpha[1])*fhatAL[3]))-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+alpha[0]*fhatAL[3])+0.07216878364870323*((alpha[2]+alpha[1])*fhatAL[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0]))*dfac_v; 
  if(std::abs(GhatCtrl[3]) < EPSILON) GhatCtrl[3] = 0.; 
  // control node [x,y,mu] = [-1/3,-1/3,1/3] 
  GhatCtrl[4] = ((-0.04166666666666666*((alpha[2]+alpha[1])*fhatAL[7]+alpha[0]*(fhatAL[6]+fhatAL[5])+(alpha[2]+alpha[1])*fhatAL[3]))+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+alpha[0]*fhatAL[3])-0.07216878364870323*((alpha[2]+alpha[1])*fhatAL[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0]))*dfac_v; 
  if(std::abs(GhatCtrl[4]) < EPSILON) GhatCtrl[4] = 0.; 
  // control node [x,y,mu] = [1/3,-1/3,1/3] 
  GhatCtrl[5] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]+alpha[1]*fhatAL[3])-0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+alpha[2]*fhatAL[4]+alpha[0]*(fhatAL[3]+fhatAL[1])+fhatAL[0]*alpha[1])-0.07216878364870323*(alpha[1]*fhatAL[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.04166666666666666*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0]))*dfac_v; 
  if(std::abs(GhatCtrl[5]) < EPSILON) GhatCtrl[5] = 0.; 
  // control node [x,y,mu] = [-1/3,1/3,1/3] 
  GhatCtrl[6] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]+alpha[1]*fhatAL[3]))+0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+alpha[0]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.07216878364870323*(alpha[1]*fhatAL[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.04166666666666666*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0]))*dfac_v; 
  if(std::abs(GhatCtrl[6]) < EPSILON) GhatCtrl[6] = 0.; 
  // control node [x,y,mu] = [1/3,1/3,1/3] 
  GhatCtrl[7] = (0.04166666666666666*((alpha[2]+alpha[1])*fhatAL[7]+alpha[0]*(fhatAL[6]+fhatAL[5])+(alpha[2]+alpha[1])*fhatAL[3])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+(alpha[2]+alpha[1])*fhatAL[4]+alpha[0]*(fhatAL[3]+fhatAL[2])+fhatAL[0]*alpha[2]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0]))*dfac_v; 
  if(std::abs(GhatCtrl[7]) < EPSILON) GhatCtrl[7] = 0.; 
  double uFrac = 0., fCtrl = 0., alphaCtrl = 0.;
  if(GhatCtrl[0]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*(alpha[2]+alpha[1]); 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[0]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*(fr[14]+fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.1443375672974065*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = -std::min(std::abs(GhatCtrl[0]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]<-EPSILON) {
    alphaCtrl = (-0.2041241452319315*alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[1]; 
    fCtrl = (-0.02777777777777778*fr[15])-0.04811252243246882*fr[14]+0.04811252243246882*(fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9])-0.08333333333333333*fr[8]+0.08333333333333333*fr[7]-0.08333333333333333*(fr[6]+fr[5])-0.1443375672974065*(fr[4]+fr[3]+fr[2])+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = -std::min(std::abs(GhatCtrl[1]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[2]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*fr[14]-0.04811252243246882*fr[13]+0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*fr[9]+0.08333333333333333*fr[8]-0.08333333333333333*fr[7]+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]-0.1443375672974065*(fr[4]+fr[3])+0.1443375672974065*fr[2]-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = -std::min(std::abs(GhatCtrl[2]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[3]; 
    fCtrl = 0.02777777777777778*fr[15]+0.04811252243246882*(fr[14]+fr[13])-0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*(fr[9]+fr[8]+fr[7]+fr[6])+0.08333333333333333*fr[5]-0.1443375672974065*(fr[4]+fr[3])+0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = -std::min(std::abs(GhatCtrl[3]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[0]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*(alpha[2]+alpha[1]); 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[4]; 
    fCtrl = (-0.02777777777777778*fl[15])+0.04811252243246882*(fl[14]+fl[13])-0.04811252243246882*fl[12]+0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*(fl[9]+fl[8])-0.08333333333333333*(fl[7]+fl[6])+0.08333333333333333*fl[5]-0.1443375672974065*fl[4]+0.1443375672974065*fl[3]-0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = std::min(GhatCtrl[0], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]>EPSILON) {
    alphaCtrl = (-0.2041241452319315*alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[5]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*fl[14]-0.04811252243246882*fl[13]+0.04811252243246882*fl[12]-0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*fl[9]-0.08333333333333333*(fl[8]+fl[7])+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]-0.1443375672974065*fl[4]+0.1443375672974065*fl[3]-0.1443375672974065*fl[2]+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = std::min(GhatCtrl[1], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[6]; 
    fCtrl = 0.02777777777777778*fl[15]-0.04811252243246882*fl[14]+0.04811252243246882*(fl[13]+fl[12])-0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9])+0.08333333333333333*(fl[8]+fl[7])-0.08333333333333333*(fl[6]+fl[5])-0.1443375672974065*fl[4]+0.1443375672974065*(fl[3]+fl[2])-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = std::min(GhatCtrl[2], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]>EPSILON) {
    alphaCtrl = 0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[7]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*(fl[14]+fl[13]+fl[12])+0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9]+fl[8])+0.08333333333333333*(fl[7]+fl[6]+fl[5])-0.1443375672974065*fl[4]+0.1443375672974065*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = std::min(GhatCtrl[3], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*(alpha[2]+alpha[1]); 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[8]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*(fr[14]+fr[13]+fr[12])-0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9]+fr[8])+0.08333333333333333*(fr[7]+fr[6]+fr[5])+0.1443375672974065*fr[4]-0.1443375672974065*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = -std::min(std::abs(GhatCtrl[4]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]<-EPSILON) {
    alphaCtrl = (-0.2041241452319315*alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[9]; 
    fCtrl = 0.02777777777777778*fr[15]+0.04811252243246882*fr[14]-0.04811252243246882*(fr[13]+fr[12])+0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9])+0.08333333333333333*(fr[8]+fr[7])-0.08333333333333333*(fr[6]+fr[5])+0.1443375672974065*fr[4]-0.1443375672974065*(fr[3]+fr[2])+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = -std::min(std::abs(GhatCtrl[5]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[10]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*fr[14]+0.04811252243246882*fr[13]-0.04811252243246882*fr[12]+0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*fr[9]-0.08333333333333333*(fr[8]+fr[7])+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]+0.1443375672974065*fr[4]-0.1443375672974065*fr[3]+0.1443375672974065*fr[2]-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = -std::min(std::abs(GhatCtrl[6]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[11]; 
    fCtrl = (-0.02777777777777778*fr[15])-0.04811252243246882*(fr[14]+fr[13])+0.04811252243246882*fr[12]-0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*(fr[9]+fr[8])-0.08333333333333333*(fr[7]+fr[6])+0.08333333333333333*fr[5]+0.1443375672974065*fr[4]-0.1443375672974065*fr[3]+0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = -std::min(std::abs(GhatCtrl[7]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*(alpha[2]+alpha[1]); 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[12]; 
    fCtrl = 0.02777777777777778*fl[15]-0.04811252243246882*(fl[14]+fl[13])+0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*(fl[9]+fl[8]+fl[7]+fl[6])+0.08333333333333333*fl[5]+0.1443375672974065*(fl[4]+fl[3])-0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = std::min(GhatCtrl[4], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]>EPSILON) {
    alphaCtrl = (-0.2041241452319315*alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[13]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*fl[14]+0.04811252243246882*fl[13]-0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*fl[9]+0.08333333333333333*fl[8]-0.08333333333333333*fl[7]+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]+0.1443375672974065*(fl[4]+fl[3])-0.1443375672974065*fl[2]+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = std::min(GhatCtrl[5], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[14]; 
    fCtrl = (-0.02777777777777778*fl[15])+0.04811252243246882*fl[14]-0.04811252243246882*(fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9])-0.08333333333333333*fl[8]+0.08333333333333333*fl[7]-0.08333333333333333*(fl[6]+fl[5])+0.1443375672974065*(fl[4]+fl[3]+fl[2])-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = std::min(GhatCtrl[6], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]>EPSILON) {
    alphaCtrl = 0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[15]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*(fl[14]+fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+0.1443375672974065*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = std::min(GhatCtrl[7], std::abs(uFrac)*fCtrl/dtApprox);
  }

  incr[0] = 0.25*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[1] = 0.4330127018922193*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[2] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[3] = -0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[4] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[5] = 0.75*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[6] = -0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[7] = -0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[8] = 0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[9] = 0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[10] = -0.75*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[11] = -1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[12] = 1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[13] = -1.299038105676658*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[14] = -1.299038105676658*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[15] = -2.25*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 

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
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double *cflRateCtrlL, const double *cflRateCtrlR, const double *w, const double *dxv, const double dtApprox, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.125*(2.0*(1.732050807568877*BdriftX[1]-1.0*BdriftX[0])*m_*wv2+((5.196152422706631*BmagInv[1]-3.0*BmagInv[0])*Phi[3]+(1.732050807568877*BmagInv[0]-3.0*BmagInv[1])*Phi[2])*dfac_y*q_))/q_; 

  double alpha[8]; 
  alpha[0] = -(0.5*((4.898979485566357*BdriftX[1]-2.828427124746191*BdriftX[0])*m_*wv2+((7.348469228349534*BmagInv[1]-4.242640687119286*BmagInv[0])*Phi[3]+(2.449489742783178*BmagInv[0]-4.242640687119286*BmagInv[1])*Phi[2])*dfac_y*q_))/q_; 
  alpha[2] = -(0.3333333333333333*(4.242640687119286*BdriftX[1]-2.449489742783178*BdriftX[0])*m_*wv)/(dfac_v*q_); 
  double rCtrlL[8], rCtrlR[8];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrlL[0] = -(1.0*(1.732050807568877*fl[15]-3.0*(fl[13]+fl[12]+fl[11])+5.196152422706631*(fl[8]+fl[6]+fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON-1.0*fl[14]+1.732050807568877*(fl[10]+fl[9]+fl[7])-3.0*(fl[4]+fl[3]+fl[2])+5.196152422706631*fl[0]); 
  rCtrlL[1] = (1.732050807568877*fl[15]+3.0*(fl[13]-1.0*(fl[12]+fl[11]))+5.196152422706631*(fl[5]-1.0*(fl[8]+fl[6]))+9.0*fl[1])/(20.78460969082652*EPSILON+fl[14]+1.732050807568877*(fl[10]-1.0*(fl[9]+fl[7]))+3.0*(fl[2]-1.0*(fl[4]+fl[3]))+5.196152422706631*fl[0]); 
  rCtrlL[2] = (1.732050807568877*fl[15]+3.0*((-1.0*fl[13])+fl[12]-1.0*fl[11])+5.196152422706631*((-1.0*fl[8])+fl[6]-1.0*fl[5])+9.0*fl[1])/(20.78460969082652*EPSILON+fl[14]+1.732050807568877*((-1.0*fl[10])+fl[9]-1.0*fl[7])+3.0*((-1.0*fl[4])+fl[3]-1.0*fl[2])+5.196152422706631*fl[0]); 
  rCtrlL[3] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[13]+fl[12]-1.0*fl[11])+5.196152422706631*fl[8]-1.0*(5.196152422706631*(fl[6]+fl[5])+9.0*fl[1])))/(20.78460969082652*EPSILON-1.0*fl[14]+1.732050807568877*(fl[7]-1.0*(fl[10]+fl[9]))+3.0*((-1.0*fl[4])+fl[3]+fl[2])+5.196152422706631*fl[0]); 
  rCtrlL[4] = (1.732050807568877*fl[15]+3.0*(fl[11]-1.0*(fl[13]+fl[12]))+5.196152422706631*(fl[8]-1.0*(fl[6]+fl[5]))+9.0*fl[1])/(20.78460969082652*EPSILON+fl[14]+1.732050807568877*(fl[7]-1.0*(fl[10]+fl[9]))+3.0*(fl[4]-1.0*(fl[3]+fl[2]))+5.196152422706631*fl[0]); 
  rCtrlL[5] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[13]-1.0*fl[12]+fl[11])+5.196152422706631*(fl[6]-1.0*fl[8])-1.0*(5.196152422706631*fl[5]+9.0*fl[1])))/(20.78460969082652*EPSILON-1.0*fl[14]+1.732050807568877*((-1.0*fl[10])+fl[9]-1.0*fl[7])+3.0*(fl[4]-1.0*fl[3]+fl[2])+5.196152422706631*fl[0]); 
  rCtrlL[6] = -(1.0*(1.732050807568877*fl[15]+3.0*((-1.0*fl[13])+fl[12]+fl[11])+5.196152422706631*(fl[5]-1.0*(fl[8]+fl[6]))-9.0*fl[1]))/(20.78460969082652*EPSILON-1.0*fl[14]+1.732050807568877*(fl[10]-1.0*(fl[9]+fl[7]))+3.0*(fl[4]+fl[3]-1.0*fl[2])+5.196152422706631*fl[0]); 
  rCtrlL[7] = (1.732050807568877*fl[15]+3.0*(fl[13]+fl[12]+fl[11])+5.196152422706631*(fl[8]+fl[6]+fl[5])+9.0*fl[1])/(20.78460969082652*EPSILON+fl[14]+1.732050807568877*(fl[10]+fl[9]+fl[7])+3.0*(fl[4]+fl[3]+fl[2])+5.196152422706631*fl[0]); 
  rCtrlR[0] = -(1.0*(1.732050807568877*fr[15]-3.0*(fr[13]+fr[12]+fr[11])+5.196152422706631*(fr[8]+fr[6]+fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON-1.0*fr[14]+1.732050807568877*(fr[10]+fr[9]+fr[7])-3.0*(fr[4]+fr[3]+fr[2])+5.196152422706631*fr[0]); 
  rCtrlR[1] = (1.732050807568877*fr[15]+3.0*(fr[13]-1.0*(fr[12]+fr[11]))+5.196152422706631*(fr[5]-1.0*(fr[8]+fr[6]))+9.0*fr[1])/(20.78460969082652*EPSILON+fr[14]+1.732050807568877*(fr[10]-1.0*(fr[9]+fr[7]))+3.0*(fr[2]-1.0*(fr[4]+fr[3]))+5.196152422706631*fr[0]); 
  rCtrlR[2] = (1.732050807568877*fr[15]+3.0*((-1.0*fr[13])+fr[12]-1.0*fr[11])+5.196152422706631*((-1.0*fr[8])+fr[6]-1.0*fr[5])+9.0*fr[1])/(20.78460969082652*EPSILON+fr[14]+1.732050807568877*((-1.0*fr[10])+fr[9]-1.0*fr[7])+3.0*((-1.0*fr[4])+fr[3]-1.0*fr[2])+5.196152422706631*fr[0]); 
  rCtrlR[3] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[13]+fr[12]-1.0*fr[11])+5.196152422706631*fr[8]-1.0*(5.196152422706631*(fr[6]+fr[5])+9.0*fr[1])))/(20.78460969082652*EPSILON-1.0*fr[14]+1.732050807568877*(fr[7]-1.0*(fr[10]+fr[9]))+3.0*((-1.0*fr[4])+fr[3]+fr[2])+5.196152422706631*fr[0]); 
  rCtrlR[4] = (1.732050807568877*fr[15]+3.0*(fr[11]-1.0*(fr[13]+fr[12]))+5.196152422706631*(fr[8]-1.0*(fr[6]+fr[5]))+9.0*fr[1])/(20.78460969082652*EPSILON+fr[14]+1.732050807568877*(fr[7]-1.0*(fr[10]+fr[9]))+3.0*(fr[4]-1.0*(fr[3]+fr[2]))+5.196152422706631*fr[0]); 
  rCtrlR[5] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[13]-1.0*fr[12]+fr[11])+5.196152422706631*(fr[6]-1.0*fr[8])-1.0*(5.196152422706631*fr[5]+9.0*fr[1])))/(20.78460969082652*EPSILON-1.0*fr[14]+1.732050807568877*((-1.0*fr[10])+fr[9]-1.0*fr[7])+3.0*(fr[4]-1.0*fr[3]+fr[2])+5.196152422706631*fr[0]); 
  rCtrlR[6] = -(1.0*(1.732050807568877*fr[15]+3.0*((-1.0*fr[13])+fr[12]+fr[11])+5.196152422706631*(fr[5]-1.0*(fr[8]+fr[6]))-9.0*fr[1]))/(20.78460969082652*EPSILON-1.0*fr[14]+1.732050807568877*(fr[10]-1.0*(fr[9]+fr[7]))+3.0*(fr[4]+fr[3]-1.0*fr[2])+5.196152422706631*fr[0]); 
  rCtrlR[7] = (1.732050807568877*fr[15]+3.0*(fr[13]+fr[12]+fr[11])+5.196152422706631*(fr[8]+fr[6]+fr[5])+9.0*fr[1])/(20.78460969082652*EPSILON+fr[14]+1.732050807568877*(fr[10]+fr[9]+fr[7])+3.0*(fr[4]+fr[3]+fr[2])+5.196152422706631*fr[0]); 
  double fCtrlL[8], fCtrlR[8];  // fCtrl = anti-limited f evaluated at each control node on x surface 
  // control node [y,vx,vy] = [-1/3,-1/3,-1/3] 
  fCtrlL[0] = -0.04811252243246882*(fl[14]-1.732050807568877*(fl[10]+fl[9]+fl[7])+3.0*(fl[4]+fl[3]+fl[2])-5.196152422706631*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.04811252243246882*(fr[14]-1.732050807568877*(fr[10]+fr[9]+fr[7])+3.0*(fr[4]+fr[3]+fr[2])-5.196152422706631*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [y,vx,vy] = [1/3,-1/3,-1/3] 
  fCtrlL[1] = 0.04811252243246882*(fl[14]+1.732050807568877*fl[10]-1.732050807568877*(fl[9]+fl[7])-3.0*(fl[4]+fl[3])+3.0*fl[2]+5.196152422706631*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.04811252243246882*(fr[14]+1.732050807568877*fr[10]-1.732050807568877*(fr[9]+fr[7])-3.0*(fr[4]+fr[3])+3.0*fr[2]+5.196152422706631*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [y,vx,vy] = [-1/3,1/3,-1/3] 
  fCtrlL[2] = 0.04811252243246882*(fl[14]-1.732050807568877*fl[10]+1.732050807568877*fl[9]-1.732050807568877*fl[7]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2]+5.196152422706631*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = 0.04811252243246882*(fr[14]-1.732050807568877*fr[10]+1.732050807568877*fr[9]-1.732050807568877*fr[7]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2]+5.196152422706631*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [y,vx,vy] = [1/3,1/3,-1/3] 
  fCtrlL[3] = -0.04811252243246882*(fl[14]+1.732050807568877*(fl[10]+fl[9])-1.732050807568877*fl[7]+3.0*fl[4]-3.0*(fl[3]+fl[2])-5.196152422706631*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = -0.04811252243246882*(fr[14]+1.732050807568877*(fr[10]+fr[9])-1.732050807568877*fr[7]+3.0*fr[4]-3.0*(fr[3]+fr[2])-5.196152422706631*fr[0])*limTheta(rCtrlR[3],-1.0); 
  // control node [y,vx,vy] = [-1/3,-1/3,1/3] 
  fCtrlL[4] = 0.04811252243246882*(fl[14]-1.732050807568877*(fl[10]+fl[9])+1.732050807568877*fl[7]+3.0*fl[4]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rCtrlL[4],1.0); 
  fCtrlR[4] = 0.04811252243246882*(fr[14]-1.732050807568877*(fr[10]+fr[9])+1.732050807568877*fr[7]+3.0*fr[4]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rCtrlR[4],-1.0); 
  // control node [y,vx,vy] = [1/3,-1/3,1/3] 
  fCtrlL[5] = -0.04811252243246882*(fl[14]+1.732050807568877*fl[10]-1.732050807568877*fl[9]+1.732050807568877*fl[7]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2]-5.196152422706631*fl[0])*limTheta(rCtrlL[5],1.0); 
  fCtrlR[5] = -0.04811252243246882*(fr[14]+1.732050807568877*fr[10]-1.732050807568877*fr[9]+1.732050807568877*fr[7]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2]-5.196152422706631*fr[0])*limTheta(rCtrlR[5],-1.0); 
  // control node [y,vx,vy] = [-1/3,1/3,1/3] 
  fCtrlL[6] = -0.04811252243246882*(fl[14]-1.732050807568877*fl[10]+1.732050807568877*(fl[9]+fl[7])-3.0*(fl[4]+fl[3])+3.0*fl[2]-5.196152422706631*fl[0])*limTheta(rCtrlL[6],1.0); 
  fCtrlR[6] = -0.04811252243246882*(fr[14]-1.732050807568877*fr[10]+1.732050807568877*(fr[9]+fr[7])-3.0*(fr[4]+fr[3])+3.0*fr[2]-5.196152422706631*fr[0])*limTheta(rCtrlR[6],-1.0); 
  // control node [y,vx,vy] = [1/3,1/3,1/3] 
  fCtrlL[7] = 0.04811252243246882*(fl[14]+1.732050807568877*(fl[10]+fl[9]+fl[7])+3.0*(fl[4]+fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rCtrlL[7],1.0); 
  fCtrlR[7] = 0.04811252243246882*(fr[14]+1.732050807568877*(fr[10]+fr[9]+fr[7])+3.0*(fr[4]+fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rCtrlR[7],-1.0); 
  double fL_AL[8], fR_AL[8];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.3535533905932737*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.6123724356957944*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*fCtrlL[4]+fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4])+fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]-1.0*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0])); 
  fL_AL[4] = 1.060660171779821*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]+fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fL_AL[5] = 1.060660171779821*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*(fCtrlL[4]+fCtrlL[3])+fCtrlL[2]-1.0*fCtrlL[1]+fCtrlL[0]); 
  fL_AL[6] = 1.060660171779821*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2])+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[7] = 1.837117307087383*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]-1.0*fCtrlL[3]+fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fR_AL[0] = 0.3535533905932737*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.6123724356957944*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*fCtrlR[4]+fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4])+fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]-1.0*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0])); 
  fR_AL[4] = 1.060660171779821*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]+fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  fR_AL[5] = 1.060660171779821*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*(fCtrlR[4]+fCtrlR[3])+fCtrlR[2]-1.0*fCtrlR[1]+fCtrlR[0]); 
  fR_AL[6] = 1.060660171779821*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2])+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[7] = 1.837117307087383*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]-1.0*fCtrlR[3]+fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[8], fhatAL[8]; 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double GhatCtrl[8];
  // control node [y,vpar,mu] = [-1/3,-1/3,-1/3] 
  GhatCtrl[0] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[0]) < EPSILON) GhatCtrl[0] = 0.; 
  // control node [y,vpar,mu] = [1/3,-1/3,-1/3] 
  GhatCtrl[1] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[1]) < EPSILON) GhatCtrl[1] = 0.; 
  // control node [y,vpar,mu] = [-1/3,1/3,-1/3] 
  GhatCtrl[2] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[2]) < EPSILON) GhatCtrl[2] = 0.; 
  // control node [y,vpar,mu] = [1/3,1/3,-1/3] 
  GhatCtrl[3] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])-0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[3]) < EPSILON) GhatCtrl[3] = 0.; 
  // control node [y,vpar,mu] = [-1/3,-1/3,1/3] 
  GhatCtrl[4] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[4]) < EPSILON) GhatCtrl[4] = 0.; 
  // control node [y,vpar,mu] = [1/3,-1/3,1/3] 
  GhatCtrl[5] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])-0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[5]) < EPSILON) GhatCtrl[5] = 0.; 
  // control node [y,vpar,mu] = [-1/3,1/3,1/3] 
  GhatCtrl[6] = ((-0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]))-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])-0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])-0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[6]) < EPSILON) GhatCtrl[6] = 0.; 
  // control node [y,vpar,mu] = [1/3,1/3,1/3] 
  GhatCtrl[7] = (0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])+0.07216878364870323*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])+0.04166666666666666*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])+0.07216878364870323*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])+0.04166666666666666*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])+0.125*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))*dfac_x; 
  if(std::abs(GhatCtrl[7]) < EPSILON) GhatCtrl[7] = 0.; 
  double uFrac = 0., fCtrl = 0., alphaCtrl = 0.;
  if(GhatCtrl[0]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[0]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*(fr[14]+fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.1443375672974065*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = -std::min(std::abs(GhatCtrl[0]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[0]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[1]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*fl[14]+0.04811252243246882*(fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9])-0.08333333333333333*fl[8]+0.08333333333333333*fl[7]-0.08333333333333333*(fl[6]+fl[5])-0.1443375672974065*(fl[4]+fl[3]+fl[2])+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = std::min(GhatCtrl[0], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[2]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*fr[14]-0.04811252243246882*fr[13]+0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*fr[9]+0.08333333333333333*fr[8]-0.08333333333333333*fr[7]+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]-0.1443375672974065*(fr[4]+fr[3])+0.1443375672974065*fr[2]-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = -std::min(std::abs(GhatCtrl[1]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[3]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*(fl[14]+fl[13])-0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*(fl[9]+fl[8]+fl[7]+fl[6])+0.08333333333333333*fl[5]-0.1443375672974065*(fl[4]+fl[3])+0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = std::min(GhatCtrl[1], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[4]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*(fr[14]+fr[13])-0.04811252243246882*fr[12]+0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*(fr[9]+fr[8])-0.08333333333333333*(fr[7]+fr[6])+0.08333333333333333*fr[5]-0.1443375672974065*fr[4]+0.1443375672974065*fr[3]-0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = -std::min(std::abs(GhatCtrl[2]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[5]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*fl[14]-0.04811252243246882*fl[13]+0.04811252243246882*fl[12]-0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*fl[9]-0.08333333333333333*(fl[8]+fl[7])+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]-0.1443375672974065*fl[4]+0.1443375672974065*fl[3]-0.1443375672974065*fl[2]+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = std::min(GhatCtrl[2], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[6]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*fr[14]+0.04811252243246882*(fr[13]+fr[12])-0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9])+0.08333333333333333*(fr[8]+fr[7])-0.08333333333333333*(fr[6]+fr[5])-0.1443375672974065*fr[4]+0.1443375672974065*(fr[3]+fr[2])-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = -std::min(std::abs(GhatCtrl[3]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[7]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*(fl[14]+fl[13]+fl[12])+0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9]+fl[8])+0.08333333333333333*(fl[7]+fl[6]+fl[5])-0.1443375672974065*fl[4]+0.1443375672974065*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = std::min(GhatCtrl[3], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[8]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*(fr[14]+fr[13]+fr[12])-0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9]+fr[8])+0.08333333333333333*(fr[7]+fr[6]+fr[5])+0.1443375672974065*fr[4]-0.1443375672974065*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = -std::min(std::abs(GhatCtrl[4]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[9]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*fl[14]-0.04811252243246882*(fl[13]+fl[12])+0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9])+0.08333333333333333*(fl[8]+fl[7])-0.08333333333333333*(fl[6]+fl[5])+0.1443375672974065*fl[4]-0.1443375672974065*(fl[3]+fl[2])+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = std::min(GhatCtrl[4], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]<-EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[10]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*fr[14]+0.04811252243246882*fr[13]-0.04811252243246882*fr[12]+0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*fr[9]-0.08333333333333333*(fr[8]+fr[7])+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]+0.1443375672974065*fr[4]-0.1443375672974065*fr[3]+0.1443375672974065*fr[2]-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = -std::min(std::abs(GhatCtrl[5]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]>EPSILON) {
    alphaCtrl = 0.3535533905932737*alpha[0]-0.2041241452319315*alpha[2]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[11]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*(fl[14]+fl[13])+0.04811252243246882*fl[12]-0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*(fl[9]+fl[8])-0.08333333333333333*(fl[7]+fl[6])+0.08333333333333333*fl[5]+0.1443375672974065*fl[4]-0.1443375672974065*fl[3]+0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = std::min(GhatCtrl[5], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[12]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*(fr[14]+fr[13])+0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*(fr[9]+fr[8]+fr[7]+fr[6])+0.08333333333333333*fr[5]+0.1443375672974065*(fr[4]+fr[3])-0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = -std::min(std::abs(GhatCtrl[6]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[13]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*fl[14]+0.04811252243246882*fl[13]-0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*fl[9]+0.08333333333333333*fl[8]-0.08333333333333333*fl[7]+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]+0.1443375672974065*(fl[4]+fl[3])-0.1443375672974065*fl[2]+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = std::min(GhatCtrl[6], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]<-EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlR[14]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*fr[14]-0.04811252243246882*(fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9])-0.08333333333333333*fr[8]+0.08333333333333333*fr[7]-0.08333333333333333*(fr[6]+fr[5])+0.1443375672974065*(fr[4]+fr[3]+fr[2])-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = -std::min(std::abs(GhatCtrl[7]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]>EPSILON) {
    alphaCtrl = 0.2041241452319315*alpha[2]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_x)/cflRateCtrlL[15]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*(fl[14]+fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+0.1443375672974065*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = std::min(GhatCtrl[7], std::abs(uFrac)*fCtrl/dtApprox);
  }

  incr[0] = 0.25*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[1] = -0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[2] = 0.4330127018922193*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[3] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[4] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[5] = -0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[6] = -0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[7] = 0.75*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[8] = -0.75*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[9] = 0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[10] = 0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[11] = -1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[12] = -1.299038105676658*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[13] = -1.299038105676658*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[14] = 1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[15] = -2.25*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 

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
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Y_P1_Bvars_1(const double q_, const double m_, const double *cflRateCtrlL, const double *cflRateCtrlR, const double *w, const double *dxv, const double dtApprox, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*(2.0*BdriftY[0]*m_*wv2+BmagInv[0]*dfac_x*(1.732050807568877*Bmag[1]*wm+(1.732050807568877*Phi[1]-3.0*Phi[3])*q_)))/q_; 

  double alpha[8]; 
  alpha[0] = (0.7071067811865475*(2.0*BdriftY[0]*m_*wv2+BmagInv[0]*dfac_x*(1.732050807568877*Bmag[1]*wm+(1.732050807568877*Phi[1]-3.0*Phi[3])*q_)))/q_; 
  alpha[1] = (0.7071067811865475*(2.0*BdriftY[1]*m_*wv2+BmagInv[1]*dfac_x*(1.732050807568877*Bmag[1]*wm+(1.732050807568877*Phi[1]-3.0*Phi[3])*q_)))/q_; 
  alpha[2] = (0.8164965809277261*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  alpha[3] = (0.7071067811865475*BmagInv[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alpha[4] = (0.8164965809277261*BdriftY[1]*m_*wv)/(dfac_v*q_); 
  alpha[5] = (0.7071067811865475*Bmag[1]*BmagInv[1]*dfac_x)/(dfac_m*q_); 
  double rCtrlL[8], rCtrlR[8];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrlL[0] = -(1.0*(1.732050807568877*fl[15]-3.0*(fl[14]+fl[12]+fl[11])+5.196152422706631*(fl[9]+fl[7]+fl[5])-9.0*fl[2]))/(20.78460969082652*EPSILON-1.0*fl[13]+1.732050807568877*(fl[10]+fl[8]+fl[6])-3.0*(fl[4]+fl[3]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[1] = (1.732050807568877*fl[15]+3.0*(fl[14]-1.0*(fl[12]+fl[11]))+5.196152422706631*(fl[5]-1.0*(fl[9]+fl[7]))+9.0*fl[2])/(20.78460969082652*EPSILON+fl[13]+1.732050807568877*(fl[10]-1.0*(fl[8]+fl[6]))+3.0*(fl[1]-1.0*(fl[4]+fl[3]))+5.196152422706631*fl[0]); 
  rCtrlL[2] = (1.732050807568877*fl[15]+3.0*((-1.0*fl[14])+fl[12]-1.0*fl[11])+5.196152422706631*((-1.0*fl[9])+fl[7]-1.0*fl[5])+9.0*fl[2])/(20.78460969082652*EPSILON+fl[13]+1.732050807568877*((-1.0*fl[10])+fl[8]-1.0*fl[6])+3.0*((-1.0*fl[4])+fl[3]-1.0*fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[3] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[14]+fl[12]-1.0*fl[11])+5.196152422706631*fl[9]-1.0*(5.196152422706631*(fl[7]+fl[5])+9.0*fl[2])))/(20.78460969082652*EPSILON-1.0*fl[13]+1.732050807568877*(fl[6]-1.0*(fl[10]+fl[8]))+3.0*((-1.0*fl[4])+fl[3]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[4] = (1.732050807568877*fl[15]+3.0*(fl[11]-1.0*(fl[14]+fl[12]))+5.196152422706631*(fl[9]-1.0*(fl[7]+fl[5]))+9.0*fl[2])/(20.78460969082652*EPSILON+fl[13]+1.732050807568877*(fl[6]-1.0*(fl[10]+fl[8]))+3.0*(fl[4]-1.0*(fl[3]+fl[1]))+5.196152422706631*fl[0]); 
  rCtrlL[5] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[14]-1.0*fl[12]+fl[11])+5.196152422706631*(fl[7]-1.0*fl[9])-1.0*(5.196152422706631*fl[5]+9.0*fl[2])))/(20.78460969082652*EPSILON-1.0*fl[13]+1.732050807568877*((-1.0*fl[10])+fl[8]-1.0*fl[6])+3.0*(fl[4]-1.0*fl[3]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[6] = -(1.0*(1.732050807568877*fl[15]+3.0*((-1.0*fl[14])+fl[12]+fl[11])+5.196152422706631*(fl[5]-1.0*(fl[9]+fl[7]))-9.0*fl[2]))/(20.78460969082652*EPSILON-1.0*fl[13]+1.732050807568877*(fl[10]-1.0*(fl[8]+fl[6]))+3.0*(fl[4]+fl[3]-1.0*fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[7] = (1.732050807568877*fl[15]+3.0*(fl[14]+fl[12]+fl[11])+5.196152422706631*(fl[9]+fl[7]+fl[5])+9.0*fl[2])/(20.78460969082652*EPSILON+fl[13]+1.732050807568877*(fl[10]+fl[8]+fl[6])+3.0*(fl[4]+fl[3]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlR[0] = -(1.0*(1.732050807568877*fr[15]-3.0*(fr[14]+fr[12]+fr[11])+5.196152422706631*(fr[9]+fr[7]+fr[5])-9.0*fr[2]))/(20.78460969082652*EPSILON-1.0*fr[13]+1.732050807568877*(fr[10]+fr[8]+fr[6])-3.0*(fr[4]+fr[3]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[1] = (1.732050807568877*fr[15]+3.0*(fr[14]-1.0*(fr[12]+fr[11]))+5.196152422706631*(fr[5]-1.0*(fr[9]+fr[7]))+9.0*fr[2])/(20.78460969082652*EPSILON+fr[13]+1.732050807568877*(fr[10]-1.0*(fr[8]+fr[6]))+3.0*(fr[1]-1.0*(fr[4]+fr[3]))+5.196152422706631*fr[0]); 
  rCtrlR[2] = (1.732050807568877*fr[15]+3.0*((-1.0*fr[14])+fr[12]-1.0*fr[11])+5.196152422706631*((-1.0*fr[9])+fr[7]-1.0*fr[5])+9.0*fr[2])/(20.78460969082652*EPSILON+fr[13]+1.732050807568877*((-1.0*fr[10])+fr[8]-1.0*fr[6])+3.0*((-1.0*fr[4])+fr[3]-1.0*fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[3] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[14]+fr[12]-1.0*fr[11])+5.196152422706631*fr[9]-1.0*(5.196152422706631*(fr[7]+fr[5])+9.0*fr[2])))/(20.78460969082652*EPSILON-1.0*fr[13]+1.732050807568877*(fr[6]-1.0*(fr[10]+fr[8]))+3.0*((-1.0*fr[4])+fr[3]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[4] = (1.732050807568877*fr[15]+3.0*(fr[11]-1.0*(fr[14]+fr[12]))+5.196152422706631*(fr[9]-1.0*(fr[7]+fr[5]))+9.0*fr[2])/(20.78460969082652*EPSILON+fr[13]+1.732050807568877*(fr[6]-1.0*(fr[10]+fr[8]))+3.0*(fr[4]-1.0*(fr[3]+fr[1]))+5.196152422706631*fr[0]); 
  rCtrlR[5] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[14]-1.0*fr[12]+fr[11])+5.196152422706631*(fr[7]-1.0*fr[9])-1.0*(5.196152422706631*fr[5]+9.0*fr[2])))/(20.78460969082652*EPSILON-1.0*fr[13]+1.732050807568877*((-1.0*fr[10])+fr[8]-1.0*fr[6])+3.0*(fr[4]-1.0*fr[3]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[6] = -(1.0*(1.732050807568877*fr[15]+3.0*((-1.0*fr[14])+fr[12]+fr[11])+5.196152422706631*(fr[5]-1.0*(fr[9]+fr[7]))-9.0*fr[2]))/(20.78460969082652*EPSILON-1.0*fr[13]+1.732050807568877*(fr[10]-1.0*(fr[8]+fr[6]))+3.0*(fr[4]+fr[3]-1.0*fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[7] = (1.732050807568877*fr[15]+3.0*(fr[14]+fr[12]+fr[11])+5.196152422706631*(fr[9]+fr[7]+fr[5])+9.0*fr[2])/(20.78460969082652*EPSILON+fr[13]+1.732050807568877*(fr[10]+fr[8]+fr[6])+3.0*(fr[4]+fr[3]+fr[1])+5.196152422706631*fr[0]); 
  double fCtrlL[8], fCtrlR[8];  // fCtrl = anti-limited f evaluated at each control node on y surface 
  // control node [x,vx,vy] = [-1/3,-1/3,-1/3] 
  fCtrlL[0] = -0.04811252243246882*(fl[13]-1.732050807568877*(fl[10]+fl[8]+fl[6])+3.0*(fl[4]+fl[3]+fl[1])-5.196152422706631*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.04811252243246882*(fr[13]-1.732050807568877*(fr[10]+fr[8]+fr[6])+3.0*(fr[4]+fr[3]+fr[1])-5.196152422706631*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,vx,vy] = [1/3,-1/3,-1/3] 
  fCtrlL[1] = 0.04811252243246882*(fl[13]+1.732050807568877*fl[10]-1.732050807568877*(fl[8]+fl[6])-3.0*(fl[4]+fl[3])+3.0*fl[1]+5.196152422706631*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.04811252243246882*(fr[13]+1.732050807568877*fr[10]-1.732050807568877*(fr[8]+fr[6])-3.0*(fr[4]+fr[3])+3.0*fr[1]+5.196152422706631*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,vx,vy] = [-1/3,1/3,-1/3] 
  fCtrlL[2] = 0.04811252243246882*(fl[13]-1.732050807568877*fl[10]+1.732050807568877*fl[8]-1.732050807568877*fl[6]-3.0*fl[4]+3.0*fl[3]-3.0*fl[1]+5.196152422706631*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = 0.04811252243246882*(fr[13]-1.732050807568877*fr[10]+1.732050807568877*fr[8]-1.732050807568877*fr[6]-3.0*fr[4]+3.0*fr[3]-3.0*fr[1]+5.196152422706631*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,vx,vy] = [1/3,1/3,-1/3] 
  fCtrlL[3] = -0.04811252243246882*(fl[13]+1.732050807568877*(fl[10]+fl[8])-1.732050807568877*fl[6]+3.0*fl[4]-3.0*(fl[3]+fl[1])-5.196152422706631*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = -0.04811252243246882*(fr[13]+1.732050807568877*(fr[10]+fr[8])-1.732050807568877*fr[6]+3.0*fr[4]-3.0*(fr[3]+fr[1])-5.196152422706631*fr[0])*limTheta(rCtrlR[3],-1.0); 
  // control node [x,vx,vy] = [-1/3,-1/3,1/3] 
  fCtrlL[4] = 0.04811252243246882*(fl[13]-1.732050807568877*(fl[10]+fl[8])+1.732050807568877*fl[6]+3.0*fl[4]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rCtrlL[4],1.0); 
  fCtrlR[4] = 0.04811252243246882*(fr[13]-1.732050807568877*(fr[10]+fr[8])+1.732050807568877*fr[6]+3.0*fr[4]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rCtrlR[4],-1.0); 
  // control node [x,vx,vy] = [1/3,-1/3,1/3] 
  fCtrlL[5] = -0.04811252243246882*(fl[13]+1.732050807568877*fl[10]-1.732050807568877*fl[8]+1.732050807568877*fl[6]-3.0*fl[4]+3.0*fl[3]-3.0*fl[1]-5.196152422706631*fl[0])*limTheta(rCtrlL[5],1.0); 
  fCtrlR[5] = -0.04811252243246882*(fr[13]+1.732050807568877*fr[10]-1.732050807568877*fr[8]+1.732050807568877*fr[6]-3.0*fr[4]+3.0*fr[3]-3.0*fr[1]-5.196152422706631*fr[0])*limTheta(rCtrlR[5],-1.0); 
  // control node [x,vx,vy] = [-1/3,1/3,1/3] 
  fCtrlL[6] = -0.04811252243246882*(fl[13]-1.732050807568877*fl[10]+1.732050807568877*(fl[8]+fl[6])-3.0*(fl[4]+fl[3])+3.0*fl[1]-5.196152422706631*fl[0])*limTheta(rCtrlL[6],1.0); 
  fCtrlR[6] = -0.04811252243246882*(fr[13]-1.732050807568877*fr[10]+1.732050807568877*(fr[8]+fr[6])-3.0*(fr[4]+fr[3])+3.0*fr[1]-5.196152422706631*fr[0])*limTheta(rCtrlR[6],-1.0); 
  // control node [x,vx,vy] = [1/3,1/3,1/3] 
  fCtrlL[7] = 0.04811252243246882*(fl[13]+1.732050807568877*(fl[10]+fl[8]+fl[6])+3.0*(fl[4]+fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rCtrlL[7],1.0); 
  fCtrlR[7] = 0.04811252243246882*(fr[13]+1.732050807568877*(fr[10]+fr[8]+fr[6])+3.0*(fr[4]+fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rCtrlR[7],-1.0); 
  double fL_AL[8], fR_AL[8];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.3535533905932737*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.6123724356957944*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*fCtrlL[4]+fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4])+fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]-1.0*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0])); 
  fL_AL[4] = 1.060660171779821*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]+fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fL_AL[5] = 1.060660171779821*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*(fCtrlL[4]+fCtrlL[3])+fCtrlL[2]-1.0*fCtrlL[1]+fCtrlL[0]); 
  fL_AL[6] = 1.060660171779821*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2])+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[7] = 1.837117307087383*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]-1.0*fCtrlL[3]+fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fR_AL[0] = 0.3535533905932737*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.6123724356957944*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*fCtrlR[4]+fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4])+fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]-1.0*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0])); 
  fR_AL[4] = 1.060660171779821*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]+fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  fR_AL[5] = 1.060660171779821*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*(fCtrlR[4]+fCtrlR[3])+fCtrlR[2]-1.0*fCtrlR[1]+fCtrlR[0]); 
  fR_AL[6] = 1.060660171779821*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2])+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[7] = 1.837117307087383*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]-1.0*fCtrlR[3]+fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[8], fhatAL[8]; 
  alphaQuad = 0.3535533905932737*(alpha[5]+alpha[4])-0.3535533905932737*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0])-0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2]); 
  if(alphaQuad > 0) {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*alpha[5]-0.3535533905932737*(alpha[4]+alpha[3])+0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = (-0.3535533905932737*alpha[5])+0.3535533905932737*alpha[4]-0.3535533905932737*alpha[3]+0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = (-0.3535533905932737*alpha[5])+0.3535533905932737*(alpha[4]+alpha[3])-0.3535533905932737*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[5]-0.3535533905932737*alpha[4]+0.3535533905932737*alpha[3]-0.3535533905932737*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = (-0.3535533905932737*(alpha[5]+alpha[4]))+0.3535533905932737*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double GhatCtrl[8];
  // control node [x,vpar,mu] = [-1/3,-1/3,-1/3] 
  GhatCtrl[0] = ((-0.07216878364870323*((alpha[5]+alpha[4])*fhatAL[7]+(alpha[3]+alpha[2])*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))+0.04166666666666666*((alpha[3]+alpha[2]+alpha[1])*fhatAL[7]+alpha[5]*fhatAL[6]+(alpha[4]+alpha[0])*(fhatAL[6]+fhatAL[5])+(fhatAL[4]+fhatAL[0])*alpha[5]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+(alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_y; 
  if(std::abs(GhatCtrl[0]) < EPSILON) GhatCtrl[0] = 0.; 
  // control node [x,vpar,mu] = [1/3,-1/3,-1/3] 
  GhatCtrl[1] = ((-0.07216878364870323*((alpha[5]+alpha[4])*fhatAL[7]+(alpha[3]+alpha[2])*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))-0.04166666666666666*((alpha[3]+alpha[2])*fhatAL[7]+(alpha[5]+alpha[4])*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_y; 
  if(std::abs(GhatCtrl[1]) < EPSILON) GhatCtrl[1] = 0.; 
  // control node [x,vpar,mu] = [-1/3,1/3,-1/3] 
  GhatCtrl[2] = (0.07216878364870323*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.07216878364870323*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])-0.04166666666666666*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[4]*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])-0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_y; 
  if(std::abs(GhatCtrl[2]) < EPSILON) GhatCtrl[2] = 0.; 
  // control node [x,vpar,mu] = [1/3,1/3,-1/3] 
  GhatCtrl[3] = (0.07216878364870323*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.07216878364870323*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])+0.04166666666666666*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])-0.04166666666666666*((alpha[2]+alpha[1])*fhatAL[7]+(alpha[4]+alpha[0])*fhatAL[6]+(alpha[4]+alpha[0])*fhatAL[5]+(fhatAL[4]+fhatAL[0])*alpha[5]+(alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_y; 
  if(std::abs(GhatCtrl[3]) < EPSILON) GhatCtrl[3] = 0.; 
  // control node [x,vpar,mu] = [-1/3,-1/3,1/3] 
  GhatCtrl[4] = ((-0.07216878364870323*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))+0.07216878364870323*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])+0.04166666666666666*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])-0.04166666666666666*((alpha[2]+alpha[1])*fhatAL[7]+(alpha[4]+alpha[0])*fhatAL[6]+(alpha[4]+alpha[0])*fhatAL[5]+(fhatAL[4]+fhatAL[0])*alpha[5]+(alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_y; 
  if(std::abs(GhatCtrl[4]) < EPSILON) GhatCtrl[4] = 0.; 
  // control node [x,vpar,mu] = [1/3,-1/3,1/3] 
  GhatCtrl[5] = ((-0.07216878364870323*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))+0.07216878364870323*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])-0.04166666666666666*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[4]*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])-0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_y; 
  if(std::abs(GhatCtrl[5]) < EPSILON) GhatCtrl[5] = 0.; 
  // control node [x,vpar,mu] = [-1/3,1/3,1/3] 
  GhatCtrl[6] = (0.07216878364870323*((alpha[5]+alpha[4])*fhatAL[7]+(alpha[3]+alpha[2])*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.04166666666666666*((alpha[3]+alpha[2])*fhatAL[7]+(alpha[5]+alpha[4])*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_y; 
  if(std::abs(GhatCtrl[6]) < EPSILON) GhatCtrl[6] = 0.; 
  // control node [x,vpar,mu] = [1/3,1/3,1/3] 
  GhatCtrl[7] = (0.07216878364870323*((alpha[5]+alpha[4])*fhatAL[7]+(alpha[3]+alpha[2])*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])+0.04166666666666666*((alpha[3]+alpha[2]+alpha[1])*fhatAL[7]+alpha[5]*fhatAL[6]+(alpha[4]+alpha[0])*(fhatAL[6]+fhatAL[5])+(fhatAL[4]+fhatAL[0])*alpha[5]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+(alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_y; 
  if(std::abs(GhatCtrl[7]) < EPSILON) GhatCtrl[7] = 0.; 
  double uFrac = 0., fCtrl = 0., alphaCtrl = 0.;
  if(GhatCtrl[0]<-EPSILON) {
    alphaCtrl = 0.1178511301977579*(alpha[5]+alpha[4])-0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[0]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*(fr[14]+fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.1443375672974065*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = -std::min(std::abs(GhatCtrl[0]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]<-EPSILON) {
    alphaCtrl = (-0.1178511301977579*(alpha[5]+alpha[4]))-0.2041241452319315*(alpha[3]+alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[1]; 
    fCtrl = (-0.02777777777777778*fr[15])-0.04811252243246882*fr[14]+0.04811252243246882*(fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9])-0.08333333333333333*fr[8]+0.08333333333333333*fr[7]-0.08333333333333333*(fr[6]+fr[5])-0.1443375672974065*(fr[4]+fr[3]+fr[2])+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = -std::min(std::abs(GhatCtrl[1]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[0]>EPSILON) {
    alphaCtrl = 0.1178511301977579*(alpha[5]+alpha[4])-0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[2]; 
    fCtrl = (-0.02777777777777778*fl[15])+0.04811252243246882*fl[14]-0.04811252243246882*fl[13]+0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*fl[9]+0.08333333333333333*fl[8]-0.08333333333333333*fl[7]+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]-0.1443375672974065*(fl[4]+fl[3])+0.1443375672974065*fl[2]-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = std::min(GhatCtrl[0], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]>EPSILON) {
    alphaCtrl = (-0.1178511301977579*(alpha[5]+alpha[4]))-0.2041241452319315*(alpha[3]+alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[3]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*(fl[14]+fl[13])-0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*(fl[9]+fl[8]+fl[7]+fl[6])+0.08333333333333333*fl[5]-0.1443375672974065*(fl[4]+fl[3])+0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = std::min(GhatCtrl[1], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]<-EPSILON) {
    alphaCtrl = 0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[4]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*(fr[14]+fr[13])-0.04811252243246882*fr[12]+0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*(fr[9]+fr[8])-0.08333333333333333*(fr[7]+fr[6])+0.08333333333333333*fr[5]-0.1443375672974065*fr[4]+0.1443375672974065*fr[3]-0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = -std::min(std::abs(GhatCtrl[2]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]<-EPSILON) {
    alphaCtrl = (-0.1178511301977579*alpha[5])+0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[5]; 
    fCtrl = 0.02777777777777778*fr[15]+0.04811252243246882*fr[14]-0.04811252243246882*fr[13]+0.04811252243246882*fr[12]-0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*fr[9]-0.08333333333333333*(fr[8]+fr[7])+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]-0.1443375672974065*fr[4]+0.1443375672974065*fr[3]-0.1443375672974065*fr[2]+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = -std::min(std::abs(GhatCtrl[3]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]>EPSILON) {
    alphaCtrl = 0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[6]; 
    fCtrl = 0.02777777777777778*fl[15]-0.04811252243246882*fl[14]+0.04811252243246882*(fl[13]+fl[12])-0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9])+0.08333333333333333*(fl[8]+fl[7])-0.08333333333333333*(fl[6]+fl[5])-0.1443375672974065*fl[4]+0.1443375672974065*(fl[3]+fl[2])-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = std::min(GhatCtrl[2], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]>EPSILON) {
    alphaCtrl = (-0.1178511301977579*alpha[5])+0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[7]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*(fl[14]+fl[13]+fl[12])+0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9]+fl[8])+0.08333333333333333*(fl[7]+fl[6]+fl[5])-0.1443375672974065*fl[4]+0.1443375672974065*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = std::min(GhatCtrl[3], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]<-EPSILON) {
    alphaCtrl = (-0.1178511301977579*alpha[5])+0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[8]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*(fr[14]+fr[13]+fr[12])-0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9]+fr[8])+0.08333333333333333*(fr[7]+fr[6]+fr[5])+0.1443375672974065*fr[4]-0.1443375672974065*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = -std::min(std::abs(GhatCtrl[4]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]<-EPSILON) {
    alphaCtrl = 0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*alpha[2]+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[9]; 
    fCtrl = 0.02777777777777778*fr[15]+0.04811252243246882*fr[14]-0.04811252243246882*(fr[13]+fr[12])+0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9])+0.08333333333333333*(fr[8]+fr[7])-0.08333333333333333*(fr[6]+fr[5])+0.1443375672974065*fr[4]-0.1443375672974065*(fr[3]+fr[2])+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = -std::min(std::abs(GhatCtrl[5]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]>EPSILON) {
    alphaCtrl = (-0.1178511301977579*alpha[5])+0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[10]; 
    fCtrl = 0.02777777777777778*fl[15]-0.04811252243246882*fl[14]+0.04811252243246882*fl[13]-0.04811252243246882*fl[12]+0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*fl[9]-0.08333333333333333*(fl[8]+fl[7])+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]+0.1443375672974065*fl[4]-0.1443375672974065*fl[3]+0.1443375672974065*fl[2]-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = std::min(GhatCtrl[4], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]>EPSILON) {
    alphaCtrl = 0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*alpha[2]+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[11]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*(fl[14]+fl[13])+0.04811252243246882*fl[12]-0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*(fl[9]+fl[8])-0.08333333333333333*(fl[7]+fl[6])+0.08333333333333333*fl[5]+0.1443375672974065*fl[4]-0.1443375672974065*fl[3]+0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = std::min(GhatCtrl[5], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]<-EPSILON) {
    alphaCtrl = (-0.1178511301977579*(alpha[5]+alpha[4]))+0.2041241452319315*(alpha[3]+alpha[2])-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[12]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*(fr[14]+fr[13])+0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*(fr[9]+fr[8]+fr[7]+fr[6])+0.08333333333333333*fr[5]+0.1443375672974065*(fr[4]+fr[3])-0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = -std::min(std::abs(GhatCtrl[6]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]<-EPSILON) {
    alphaCtrl = 0.1178511301977579*(alpha[5]+alpha[4])+0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlR[13]; 
    fCtrl = (-0.02777777777777778*fr[15])-0.04811252243246882*fr[14]+0.04811252243246882*fr[13]-0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*fr[9]+0.08333333333333333*fr[8]-0.08333333333333333*fr[7]+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]+0.1443375672974065*(fr[4]+fr[3])-0.1443375672974065*fr[2]+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = -std::min(std::abs(GhatCtrl[7]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]>EPSILON) {
    alphaCtrl = (-0.1178511301977579*(alpha[5]+alpha[4]))+0.2041241452319315*(alpha[3]+alpha[2])-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[14]; 
    fCtrl = (-0.02777777777777778*fl[15])+0.04811252243246882*fl[14]-0.04811252243246882*(fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9])-0.08333333333333333*fl[8]+0.08333333333333333*fl[7]-0.08333333333333333*(fl[6]+fl[5])+0.1443375672974065*(fl[4]+fl[3]+fl[2])-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = std::min(GhatCtrl[6], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]>EPSILON) {
    alphaCtrl = 0.1178511301977579*(alpha[5]+alpha[4])+0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_y)/cflRateCtrlL[15]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*(fl[14]+fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+0.1443375672974065*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = std::min(GhatCtrl[7], std::abs(uFrac)*fCtrl/dtApprox);
  }

  incr[0] = 0.25*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[1] = 0.4330127018922193*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[2] = -0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[3] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[4] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[5] = -0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[6] = 0.75*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[7] = -0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[8] = 0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[9] = -0.75*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[10] = 0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[11] = -1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[12] = -1.299038105676658*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[13] = 1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[14] = -1.299038105676658*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[15] = -2.25*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 

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
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double *cflRateCtrlL, const double *cflRateCtrlR, const double *w, const double *dxv, const double dtApprox, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.2165063509461096*(dfac_v*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*wv-1.0*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)))/(dfac_v*q_); 

  double alpha[8]; 
  alpha[0] = -(0.7071067811865475*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*(1.732050807568877*dfac_v*wv-1.732050807568877))/(dfac_v*q_); 
  alpha[1] = -(0.7071067811865475*(BdriftX[1]*Bmag[1]*dfac_x*wm+((BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*q_)*(1.732050807568877*dfac_v*wv-1.732050807568877))/(dfac_v*q_); 
  alpha[2] = -(0.7071067811865475*BdriftX[0]*Phi[3]*dfac_x*(1.732050807568877*dfac_v*wv-1.732050807568877))/dfac_v; 
  alpha[3] = -(0.7071067811865475*BdriftX[0]*Bmag[1]*dfac_x*(dfac_v*wv-1.0))/(dfac_m*dfac_v*q_); 
  alpha[4] = -(0.7071067811865475*BdriftX[1]*Phi[3]*dfac_x*(1.732050807568877*dfac_v*wv-1.732050807568877))/dfac_v; 
  alpha[5] = -(0.7071067811865475*BdriftX[1]*Bmag[1]*dfac_x*(dfac_v*wv-1.0))/(dfac_m*dfac_v*q_); 
  double rCtrlL[8], rCtrlR[8];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrlL[0] = -(1.0*(1.732050807568877*fl[15]-3.0*(fl[14]+fl[13]+fl[11])+5.196152422706631*(fl[10]+fl[7]+fl[6])-9.0*fl[3]))/(20.78460969082652*EPSILON-1.0*fl[12]+1.732050807568877*(fl[9]+fl[8]+fl[5])-3.0*(fl[4]+fl[2]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[1] = (1.732050807568877*fl[15]+3.0*(fl[14]-1.0*(fl[13]+fl[11]))+5.196152422706631*(fl[6]-1.0*(fl[10]+fl[7]))+9.0*fl[3])/(20.78460969082652*EPSILON+fl[12]+1.732050807568877*(fl[9]-1.0*(fl[8]+fl[5]))+3.0*(fl[1]-1.0*(fl[4]+fl[2]))+5.196152422706631*fl[0]); 
  rCtrlL[2] = (1.732050807568877*fl[15]+3.0*((-1.0*fl[14])+fl[13]-1.0*fl[11])+5.196152422706631*((-1.0*fl[10])+fl[7]-1.0*fl[6])+9.0*fl[3])/(20.78460969082652*EPSILON+fl[12]+1.732050807568877*((-1.0*fl[9])+fl[8]-1.0*fl[5])+3.0*((-1.0*fl[4])+fl[2]-1.0*fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[3] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[14]+fl[13]-1.0*fl[11])+5.196152422706631*fl[10]-1.0*(5.196152422706631*(fl[7]+fl[6])+9.0*fl[3])))/(20.78460969082652*EPSILON-1.0*fl[12]+1.732050807568877*(fl[5]-1.0*(fl[9]+fl[8]))+3.0*((-1.0*fl[4])+fl[2]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[4] = (1.732050807568877*fl[15]+3.0*(fl[11]-1.0*(fl[14]+fl[13]))+5.196152422706631*(fl[10]-1.0*(fl[7]+fl[6]))+9.0*fl[3])/(20.78460969082652*EPSILON+fl[12]+1.732050807568877*(fl[5]-1.0*(fl[9]+fl[8]))+3.0*(fl[4]-1.0*(fl[2]+fl[1]))+5.196152422706631*fl[0]); 
  rCtrlL[5] = -(1.0*(1.732050807568877*fl[15]+3.0*(fl[14]-1.0*fl[13]+fl[11])+5.196152422706631*(fl[7]-1.0*fl[10])-1.0*(5.196152422706631*fl[6]+9.0*fl[3])))/(20.78460969082652*EPSILON-1.0*fl[12]+1.732050807568877*((-1.0*fl[9])+fl[8]-1.0*fl[5])+3.0*(fl[4]-1.0*fl[2]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[6] = -(1.0*(1.732050807568877*fl[15]+3.0*((-1.0*fl[14])+fl[13]+fl[11])+5.196152422706631*(fl[6]-1.0*(fl[10]+fl[7]))-9.0*fl[3]))/(20.78460969082652*EPSILON-1.0*fl[12]+1.732050807568877*(fl[9]-1.0*(fl[8]+fl[5]))+3.0*(fl[4]+fl[2]-1.0*fl[1])+5.196152422706631*fl[0]); 
  rCtrlL[7] = (1.732050807568877*fl[15]+3.0*(fl[14]+fl[13]+fl[11])+5.196152422706631*(fl[10]+fl[7]+fl[6])+9.0*fl[3])/(20.78460969082652*EPSILON+fl[12]+1.732050807568877*(fl[9]+fl[8]+fl[5])+3.0*(fl[4]+fl[2]+fl[1])+5.196152422706631*fl[0]); 
  rCtrlR[0] = -(1.0*(1.732050807568877*fr[15]-3.0*(fr[14]+fr[13]+fr[11])+5.196152422706631*(fr[10]+fr[7]+fr[6])-9.0*fr[3]))/(20.78460969082652*EPSILON-1.0*fr[12]+1.732050807568877*(fr[9]+fr[8]+fr[5])-3.0*(fr[4]+fr[2]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[1] = (1.732050807568877*fr[15]+3.0*(fr[14]-1.0*(fr[13]+fr[11]))+5.196152422706631*(fr[6]-1.0*(fr[10]+fr[7]))+9.0*fr[3])/(20.78460969082652*EPSILON+fr[12]+1.732050807568877*(fr[9]-1.0*(fr[8]+fr[5]))+3.0*(fr[1]-1.0*(fr[4]+fr[2]))+5.196152422706631*fr[0]); 
  rCtrlR[2] = (1.732050807568877*fr[15]+3.0*((-1.0*fr[14])+fr[13]-1.0*fr[11])+5.196152422706631*((-1.0*fr[10])+fr[7]-1.0*fr[6])+9.0*fr[3])/(20.78460969082652*EPSILON+fr[12]+1.732050807568877*((-1.0*fr[9])+fr[8]-1.0*fr[5])+3.0*((-1.0*fr[4])+fr[2]-1.0*fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[3] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[14]+fr[13]-1.0*fr[11])+5.196152422706631*fr[10]-1.0*(5.196152422706631*(fr[7]+fr[6])+9.0*fr[3])))/(20.78460969082652*EPSILON-1.0*fr[12]+1.732050807568877*(fr[5]-1.0*(fr[9]+fr[8]))+3.0*((-1.0*fr[4])+fr[2]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[4] = (1.732050807568877*fr[15]+3.0*(fr[11]-1.0*(fr[14]+fr[13]))+5.196152422706631*(fr[10]-1.0*(fr[7]+fr[6]))+9.0*fr[3])/(20.78460969082652*EPSILON+fr[12]+1.732050807568877*(fr[5]-1.0*(fr[9]+fr[8]))+3.0*(fr[4]-1.0*(fr[2]+fr[1]))+5.196152422706631*fr[0]); 
  rCtrlR[5] = -(1.0*(1.732050807568877*fr[15]+3.0*(fr[14]-1.0*fr[13]+fr[11])+5.196152422706631*(fr[7]-1.0*fr[10])-1.0*(5.196152422706631*fr[6]+9.0*fr[3])))/(20.78460969082652*EPSILON-1.0*fr[12]+1.732050807568877*((-1.0*fr[9])+fr[8]-1.0*fr[5])+3.0*(fr[4]-1.0*fr[2]+fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[6] = -(1.0*(1.732050807568877*fr[15]+3.0*((-1.0*fr[14])+fr[13]+fr[11])+5.196152422706631*(fr[6]-1.0*(fr[10]+fr[7]))-9.0*fr[3]))/(20.78460969082652*EPSILON-1.0*fr[12]+1.732050807568877*(fr[9]-1.0*(fr[8]+fr[5]))+3.0*(fr[4]+fr[2]-1.0*fr[1])+5.196152422706631*fr[0]); 
  rCtrlR[7] = (1.732050807568877*fr[15]+3.0*(fr[14]+fr[13]+fr[11])+5.196152422706631*(fr[10]+fr[7]+fr[6])+9.0*fr[3])/(20.78460969082652*EPSILON+fr[12]+1.732050807568877*(fr[9]+fr[8]+fr[5])+3.0*(fr[4]+fr[2]+fr[1])+5.196152422706631*fr[0]); 
  double fCtrlL[8], fCtrlR[8];  // fCtrl = anti-limited f evaluated at each control node on vx surface 
  // control node [x,y,vy] = [-1/3,-1/3,-1/3] 
  fCtrlL[0] = -0.04811252243246882*(fl[12]-1.732050807568877*(fl[9]+fl[8]+fl[5])+3.0*(fl[4]+fl[2]+fl[1])-5.196152422706631*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.04811252243246882*(fr[12]-1.732050807568877*(fr[9]+fr[8]+fr[5])+3.0*(fr[4]+fr[2]+fr[1])-5.196152422706631*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,y,vy] = [1/3,-1/3,-1/3] 
  fCtrlL[1] = 0.04811252243246882*(fl[12]+1.732050807568877*fl[9]-1.732050807568877*(fl[8]+fl[5])-3.0*(fl[4]+fl[2])+3.0*fl[1]+5.196152422706631*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.04811252243246882*(fr[12]+1.732050807568877*fr[9]-1.732050807568877*(fr[8]+fr[5])-3.0*(fr[4]+fr[2])+3.0*fr[1]+5.196152422706631*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,y,vy] = [-1/3,1/3,-1/3] 
  fCtrlL[2] = 0.04811252243246882*(fl[12]-1.732050807568877*fl[9]+1.732050807568877*fl[8]-1.732050807568877*fl[5]-3.0*fl[4]+3.0*fl[2]-3.0*fl[1]+5.196152422706631*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = 0.04811252243246882*(fr[12]-1.732050807568877*fr[9]+1.732050807568877*fr[8]-1.732050807568877*fr[5]-3.0*fr[4]+3.0*fr[2]-3.0*fr[1]+5.196152422706631*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,y,vy] = [1/3,1/3,-1/3] 
  fCtrlL[3] = -0.04811252243246882*(fl[12]+1.732050807568877*(fl[9]+fl[8])-1.732050807568877*fl[5]+3.0*fl[4]-3.0*(fl[2]+fl[1])-5.196152422706631*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = -0.04811252243246882*(fr[12]+1.732050807568877*(fr[9]+fr[8])-1.732050807568877*fr[5]+3.0*fr[4]-3.0*(fr[2]+fr[1])-5.196152422706631*fr[0])*limTheta(rCtrlR[3],-1.0); 
  // control node [x,y,vy] = [-1/3,-1/3,1/3] 
  fCtrlL[4] = 0.04811252243246882*(fl[12]-1.732050807568877*(fl[9]+fl[8])+1.732050807568877*fl[5]+3.0*fl[4]-3.0*(fl[2]+fl[1])+5.196152422706631*fl[0])*limTheta(rCtrlL[4],1.0); 
  fCtrlR[4] = 0.04811252243246882*(fr[12]-1.732050807568877*(fr[9]+fr[8])+1.732050807568877*fr[5]+3.0*fr[4]-3.0*(fr[2]+fr[1])+5.196152422706631*fr[0])*limTheta(rCtrlR[4],-1.0); 
  // control node [x,y,vy] = [1/3,-1/3,1/3] 
  fCtrlL[5] = -0.04811252243246882*(fl[12]+1.732050807568877*fl[9]-1.732050807568877*fl[8]+1.732050807568877*fl[5]-3.0*fl[4]+3.0*fl[2]-3.0*fl[1]-5.196152422706631*fl[0])*limTheta(rCtrlL[5],1.0); 
  fCtrlR[5] = -0.04811252243246882*(fr[12]+1.732050807568877*fr[9]-1.732050807568877*fr[8]+1.732050807568877*fr[5]-3.0*fr[4]+3.0*fr[2]-3.0*fr[1]-5.196152422706631*fr[0])*limTheta(rCtrlR[5],-1.0); 
  // control node [x,y,vy] = [-1/3,1/3,1/3] 
  fCtrlL[6] = -0.04811252243246882*(fl[12]-1.732050807568877*fl[9]+1.732050807568877*(fl[8]+fl[5])-3.0*(fl[4]+fl[2])+3.0*fl[1]-5.196152422706631*fl[0])*limTheta(rCtrlL[6],1.0); 
  fCtrlR[6] = -0.04811252243246882*(fr[12]-1.732050807568877*fr[9]+1.732050807568877*(fr[8]+fr[5])-3.0*(fr[4]+fr[2])+3.0*fr[1]-5.196152422706631*fr[0])*limTheta(rCtrlR[6],-1.0); 
  // control node [x,y,vy] = [1/3,1/3,1/3] 
  fCtrlL[7] = 0.04811252243246882*(fl[12]+1.732050807568877*(fl[9]+fl[8]+fl[5])+3.0*(fl[4]+fl[2]+fl[1])+5.196152422706631*fl[0])*limTheta(rCtrlL[7],1.0); 
  fCtrlR[7] = 0.04811252243246882*(fr[12]+1.732050807568877*(fr[9]+fr[8]+fr[5])+3.0*(fr[4]+fr[2]+fr[1])+5.196152422706631*fr[0])*limTheta(rCtrlR[7],-1.0); 
  double fL_AL[8], fR_AL[8];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.3535533905932737*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.6123724356957944*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*fCtrlL[4]+fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4])+fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 0.6123724356957944*(fCtrlL[7]+fCtrlL[6]+fCtrlL[5]+fCtrlL[4]-1.0*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0])); 
  fL_AL[4] = 1.060660171779821*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]+fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fL_AL[5] = 1.060660171779821*(fCtrlL[7]-1.0*fCtrlL[6]+fCtrlL[5]-1.0*(fCtrlL[4]+fCtrlL[3])+fCtrlL[2]-1.0*fCtrlL[1]+fCtrlL[0]); 
  fL_AL[6] = 1.060660171779821*(fCtrlL[7]+fCtrlL[6]-1.0*(fCtrlL[5]+fCtrlL[4]+fCtrlL[3]+fCtrlL[2])+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[7] = 1.837117307087383*(fCtrlL[7]-1.0*(fCtrlL[6]+fCtrlL[5])+fCtrlL[4]-1.0*fCtrlL[3]+fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fR_AL[0] = 0.3535533905932737*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.6123724356957944*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*fCtrlR[4]+fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4])+fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 0.6123724356957944*(fCtrlR[7]+fCtrlR[6]+fCtrlR[5]+fCtrlR[4]-1.0*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0])); 
  fR_AL[4] = 1.060660171779821*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]+fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  fR_AL[5] = 1.060660171779821*(fCtrlR[7]-1.0*fCtrlR[6]+fCtrlR[5]-1.0*(fCtrlR[4]+fCtrlR[3])+fCtrlR[2]-1.0*fCtrlR[1]+fCtrlR[0]); 
  fR_AL[6] = 1.060660171779821*(fCtrlR[7]+fCtrlR[6]-1.0*(fCtrlR[5]+fCtrlR[4]+fCtrlR[3]+fCtrlR[2])+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[7] = 1.837117307087383*(fCtrlR[7]-1.0*(fCtrlR[6]+fCtrlR[5])+fCtrlR[4]-1.0*fCtrlR[3]+fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[8], fhatAL[8]; 
  alphaQuad = 0.3535533905932737*(alpha[5]+alpha[4])-0.3535533905932737*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[0] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0])-0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2]); 
  if(alphaQuad > 0) {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[1] = std::max(0., 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = 0.3535533905932737*alpha[5]-0.3535533905932737*(alpha[4]+alpha[3])+0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[2] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = (-0.3535533905932737*alpha[5])+0.3535533905932737*alpha[4]-0.3535533905932737*alpha[3]+0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[3] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = (-0.3535533905932737*alpha[5])+0.3535533905932737*(alpha[4]+alpha[3])-0.3535533905932737*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[4] = std::max(0., 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[5]-0.3535533905932737*alpha[4]+0.3535533905932737*alpha[3]-0.3535533905932737*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0]))); 
  } else {
  fhatALQuad[5] = std::max(0., -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0]))); 
  } 
  alphaQuad = (-0.3535533905932737*(alpha[5]+alpha[4]))+0.3535533905932737*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0])); 
  } else {
  fhatALQuad[6] = std::max(0., -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[7] = std::max(0., 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double GhatCtrl[8];
  // control node [x,y,mu] = [-1/3,-1/3,-1/3] 
  GhatCtrl[0] = ((-0.07216878364870323*((alpha[5]+alpha[4])*fhatAL[7]+(alpha[3]+alpha[2])*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))+0.04166666666666666*((alpha[3]+alpha[2]+alpha[1])*fhatAL[7]+alpha[5]*fhatAL[6]+(alpha[4]+alpha[0])*(fhatAL[6]+fhatAL[5])+(fhatAL[4]+fhatAL[0])*alpha[5]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+(alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_v; 
  if(std::abs(GhatCtrl[0]) < EPSILON) GhatCtrl[0] = 0.; 
  // control node [x,y,mu] = [1/3,-1/3,-1/3] 
  GhatCtrl[1] = ((-0.07216878364870323*((alpha[5]+alpha[4])*fhatAL[7]+(alpha[3]+alpha[2])*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))-0.04166666666666666*((alpha[3]+alpha[2])*fhatAL[7]+(alpha[5]+alpha[4])*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_v; 
  if(std::abs(GhatCtrl[1]) < EPSILON) GhatCtrl[1] = 0.; 
  // control node [x,y,mu] = [-1/3,1/3,-1/3] 
  GhatCtrl[2] = (0.07216878364870323*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.07216878364870323*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])-0.04166666666666666*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[4]*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])-0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_v; 
  if(std::abs(GhatCtrl[2]) < EPSILON) GhatCtrl[2] = 0.; 
  // control node [x,y,mu] = [1/3,1/3,-1/3] 
  GhatCtrl[3] = (0.07216878364870323*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.07216878364870323*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])+0.04166666666666666*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])-0.04166666666666666*((alpha[2]+alpha[1])*fhatAL[7]+(alpha[4]+alpha[0])*fhatAL[6]+(alpha[4]+alpha[0])*fhatAL[5]+(fhatAL[4]+fhatAL[0])*alpha[5]+(alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_v; 
  if(std::abs(GhatCtrl[3]) < EPSILON) GhatCtrl[3] = 0.; 
  // control node [x,y,mu] = [-1/3,-1/3,1/3] 
  GhatCtrl[4] = ((-0.07216878364870323*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))+0.07216878364870323*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])+0.04166666666666666*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])-0.04166666666666666*((alpha[2]+alpha[1])*fhatAL[7]+(alpha[4]+alpha[0])*fhatAL[6]+(alpha[4]+alpha[0])*fhatAL[5]+(fhatAL[4]+fhatAL[0])*alpha[5]+(alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_v; 
  if(std::abs(GhatCtrl[4]) < EPSILON) GhatCtrl[4] = 0.; 
  // control node [x,y,mu] = [1/3,-1/3,1/3] 
  GhatCtrl[5] = ((-0.07216878364870323*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2]))+0.07216878364870323*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])-0.04166666666666666*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.04166666666666666*(alpha[2]*fhatAL[7]+alpha[4]*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])-0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_v; 
  if(std::abs(GhatCtrl[5]) < EPSILON) GhatCtrl[5] = 0.; 
  // control node [x,y,mu] = [-1/3,1/3,1/3] 
  GhatCtrl[6] = (0.07216878364870323*((alpha[5]+alpha[4])*fhatAL[7]+(alpha[3]+alpha[2])*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])-0.04166666666666666*((alpha[3]+alpha[2])*fhatAL[7]+(alpha[5]+alpha[4])*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.04166666666666666*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])-0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])-0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_v; 
  if(std::abs(GhatCtrl[6]) < EPSILON) GhatCtrl[6] = 0.; 
  // control node [x,y,mu] = [1/3,1/3,1/3] 
  GhatCtrl[7] = (0.07216878364870323*((alpha[5]+alpha[4])*fhatAL[7]+(alpha[3]+alpha[2])*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])+0.04166666666666666*((alpha[3]+alpha[2]+alpha[1])*fhatAL[7]+alpha[5]*fhatAL[6]+(alpha[4]+alpha[0])*(fhatAL[6]+fhatAL[5])+(fhatAL[4]+fhatAL[0])*alpha[5]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+(alpha[2]+alpha[1])*fhatAL[3]+(fhatAL[2]+fhatAL[1])*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])+0.02405626121623441*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])+0.125*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])+0.07216878364870323*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1]))*dfac_v; 
  if(std::abs(GhatCtrl[7]) < EPSILON) GhatCtrl[7] = 0.; 
  double uFrac = 0., fCtrl = 0., alphaCtrl = 0.;
  if(GhatCtrl[0]<-EPSILON) {
    alphaCtrl = 0.1178511301977579*(alpha[5]+alpha[4])-0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[0]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*(fr[14]+fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-0.1443375672974065*(fr[4]+fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = -std::min(std::abs(GhatCtrl[0]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]<-EPSILON) {
    alphaCtrl = (-0.1178511301977579*(alpha[5]+alpha[4]))-0.2041241452319315*(alpha[3]+alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[1]; 
    fCtrl = (-0.02777777777777778*fr[15])-0.04811252243246882*fr[14]+0.04811252243246882*(fr[13]+fr[12]+fr[11])+0.08333333333333333*(fr[10]+fr[9])-0.08333333333333333*fr[8]+0.08333333333333333*fr[7]-0.08333333333333333*(fr[6]+fr[5])-0.1443375672974065*(fr[4]+fr[3]+fr[2])+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = -std::min(std::abs(GhatCtrl[1]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]<-EPSILON) {
    alphaCtrl = 0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[2]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*fr[14]-0.04811252243246882*fr[13]+0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*fr[9]+0.08333333333333333*fr[8]-0.08333333333333333*fr[7]+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]-0.1443375672974065*(fr[4]+fr[3])+0.1443375672974065*fr[2]-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = -std::min(std::abs(GhatCtrl[2]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]<-EPSILON) {
    alphaCtrl = (-0.1178511301977579*alpha[5])+0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[3]; 
    fCtrl = 0.02777777777777778*fr[15]+0.04811252243246882*(fr[14]+fr[13])-0.04811252243246882*(fr[12]+fr[11])+0.08333333333333333*fr[10]-0.08333333333333333*(fr[9]+fr[8]+fr[7]+fr[6])+0.08333333333333333*fr[5]-0.1443375672974065*(fr[4]+fr[3])+0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = -std::min(std::abs(GhatCtrl[3]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[0]>EPSILON) {
    alphaCtrl = 0.1178511301977579*(alpha[5]+alpha[4])-0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[4]; 
    fCtrl = (-0.02777777777777778*fl[15])+0.04811252243246882*(fl[14]+fl[13])-0.04811252243246882*fl[12]+0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*(fl[9]+fl[8])-0.08333333333333333*(fl[7]+fl[6])+0.08333333333333333*fl[5]-0.1443375672974065*fl[4]+0.1443375672974065*fl[3]-0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[0] = 0.;
    else GhatCtrl[0] = std::min(GhatCtrl[0], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[1]>EPSILON) {
    alphaCtrl = (-0.1178511301977579*(alpha[5]+alpha[4]))-0.2041241452319315*(alpha[3]+alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[5]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*fl[14]-0.04811252243246882*fl[13]+0.04811252243246882*fl[12]-0.04811252243246882*fl[11]-0.08333333333333333*fl[10]+0.08333333333333333*fl[9]-0.08333333333333333*(fl[8]+fl[7])+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]-0.1443375672974065*fl[4]+0.1443375672974065*fl[3]-0.1443375672974065*fl[2]+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[1] = 0.;
    else GhatCtrl[1] = std::min(GhatCtrl[1], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[2]>EPSILON) {
    alphaCtrl = 0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[6]; 
    fCtrl = 0.02777777777777778*fl[15]-0.04811252243246882*fl[14]+0.04811252243246882*(fl[13]+fl[12])-0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9])+0.08333333333333333*(fl[8]+fl[7])-0.08333333333333333*(fl[6]+fl[5])-0.1443375672974065*fl[4]+0.1443375672974065*(fl[3]+fl[2])-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[2] = 0.;
    else GhatCtrl[2] = std::min(GhatCtrl[2], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[3]>EPSILON) {
    alphaCtrl = (-0.1178511301977579*alpha[5])+0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[7]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*(fl[14]+fl[13]+fl[12])+0.04811252243246882*fl[11]-0.08333333333333333*(fl[10]+fl[9]+fl[8])+0.08333333333333333*(fl[7]+fl[6]+fl[5])-0.1443375672974065*fl[4]+0.1443375672974065*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[3] = 0.;
    else GhatCtrl[3] = std::min(GhatCtrl[3], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]<-EPSILON) {
    alphaCtrl = (-0.1178511301977579*alpha[5])+0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[8]; 
    fCtrl = (-0.02777777777777778*fr[15])+0.04811252243246882*(fr[14]+fr[13]+fr[12])-0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9]+fr[8])+0.08333333333333333*(fr[7]+fr[6]+fr[5])+0.1443375672974065*fr[4]-0.1443375672974065*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = -std::min(std::abs(GhatCtrl[4]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]<-EPSILON) {
    alphaCtrl = 0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*alpha[2]+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[9]; 
    fCtrl = 0.02777777777777778*fr[15]+0.04811252243246882*fr[14]-0.04811252243246882*(fr[13]+fr[12])+0.04811252243246882*fr[11]-0.08333333333333333*(fr[10]+fr[9])+0.08333333333333333*(fr[8]+fr[7])-0.08333333333333333*(fr[6]+fr[5])+0.1443375672974065*fr[4]-0.1443375672974065*(fr[3]+fr[2])+0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = -std::min(std::abs(GhatCtrl[5]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]<-EPSILON) {
    alphaCtrl = (-0.1178511301977579*(alpha[5]+alpha[4]))+0.2041241452319315*(alpha[3]+alpha[2])-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[10]; 
    fCtrl = 0.02777777777777778*fr[15]-0.04811252243246882*fr[14]+0.04811252243246882*fr[13]-0.04811252243246882*fr[12]+0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*fr[9]-0.08333333333333333*(fr[8]+fr[7])+0.08333333333333333*fr[6]-0.08333333333333333*fr[5]+0.1443375672974065*fr[4]-0.1443375672974065*fr[3]+0.1443375672974065*fr[2]-0.1443375672974065*fr[1]+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = -std::min(std::abs(GhatCtrl[6]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]<-EPSILON) {
    alphaCtrl = 0.1178511301977579*(alpha[5]+alpha[4])+0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlR[11]; 
    fCtrl = (-0.02777777777777778*fr[15])-0.04811252243246882*(fr[14]+fr[13])+0.04811252243246882*fr[12]-0.04811252243246882*fr[11]-0.08333333333333333*fr[10]+0.08333333333333333*(fr[9]+fr[8])-0.08333333333333333*(fr[7]+fr[6])+0.08333333333333333*fr[5]+0.1443375672974065*fr[4]-0.1443375672974065*fr[3]+0.1443375672974065*(fr[2]+fr[1])+0.25*fr[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = -std::min(std::abs(GhatCtrl[7]), std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[4]>EPSILON) {
    alphaCtrl = (-0.1178511301977579*alpha[5])+0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[12]; 
    fCtrl = 0.02777777777777778*fl[15]-0.04811252243246882*(fl[14]+fl[13])+0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*(fl[9]+fl[8]+fl[7]+fl[6])+0.08333333333333333*fl[5]+0.1443375672974065*(fl[4]+fl[3])-0.1443375672974065*(fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[4] = 0.;
    else GhatCtrl[4] = std::min(GhatCtrl[4], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[5]>EPSILON) {
    alphaCtrl = 0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*alpha[2]+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[13]; 
    fCtrl = (-0.02777777777777778*fl[15])-0.04811252243246882*fl[14]+0.04811252243246882*fl[13]-0.04811252243246882*(fl[12]+fl[11])+0.08333333333333333*fl[10]-0.08333333333333333*fl[9]+0.08333333333333333*fl[8]-0.08333333333333333*fl[7]+0.08333333333333333*fl[6]-0.08333333333333333*fl[5]+0.1443375672974065*(fl[4]+fl[3])-0.1443375672974065*fl[2]+0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[5] = 0.;
    else GhatCtrl[5] = std::min(GhatCtrl[5], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[6]>EPSILON) {
    alphaCtrl = (-0.1178511301977579*(alpha[5]+alpha[4]))+0.2041241452319315*(alpha[3]+alpha[2])-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[14]; 
    fCtrl = (-0.02777777777777778*fl[15])+0.04811252243246882*fl[14]-0.04811252243246882*(fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9])-0.08333333333333333*fl[8]+0.08333333333333333*fl[7]-0.08333333333333333*(fl[6]+fl[5])+0.1443375672974065*(fl[4]+fl[3]+fl[2])-0.1443375672974065*fl[1]+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[6] = 0.;
    else GhatCtrl[6] = std::min(GhatCtrl[6], std::abs(uFrac)*fCtrl/dtApprox);
  }
  if(GhatCtrl[7]>EPSILON) {
    alphaCtrl = 0.1178511301977579*(alpha[5]+alpha[4])+0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
    uFrac = (alphaCtrl*dfac_v)/cflRateCtrlL[15]; 
    fCtrl = 0.02777777777777778*fl[15]+0.04811252243246882*(fl[14]+fl[13]+fl[12]+fl[11])+0.08333333333333333*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+0.1443375672974065*(fl[4]+fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
    if(std::abs(alphaCtrl)<EPSILON || std::abs(fCtrl)<EPSILON) GhatCtrl[7] = 0.;
    else GhatCtrl[7] = std::min(GhatCtrl[7], std::abs(uFrac)*fCtrl/dtApprox);
  }

  incr[0] = 0.25*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[1] = 0.4330127018922193*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[2] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[3] = -0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0]); 
  incr[4] = 0.4330127018922193*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[5] = 0.75*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[6] = -0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*GhatCtrl[4]+GhatCtrl[3]-1.0*GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[7] = -0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4])+GhatCtrl[3]+GhatCtrl[2]-1.0*(GhatCtrl[1]+GhatCtrl[0])); 
  incr[8] = 0.75*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[9] = 0.75*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[10] = -0.75*(GhatCtrl[7]+GhatCtrl[6]+GhatCtrl[5]+GhatCtrl[4]-1.0*(GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]+GhatCtrl[0])); 
  incr[11] = -1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]+GhatCtrl[3]-1.0*(GhatCtrl[2]+GhatCtrl[1])+GhatCtrl[0]); 
  incr[12] = 1.299038105676658*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 
  incr[13] = -1.299038105676658*(GhatCtrl[7]-1.0*GhatCtrl[6]+GhatCtrl[5]-1.0*(GhatCtrl[4]+GhatCtrl[3])+GhatCtrl[2]-1.0*GhatCtrl[1]+GhatCtrl[0]); 
  incr[14] = -1.299038105676658*(GhatCtrl[7]+GhatCtrl[6]-1.0*(GhatCtrl[5]+GhatCtrl[4]+GhatCtrl[3]+GhatCtrl[2])+GhatCtrl[1]+GhatCtrl[0]); 
  incr[15] = -2.25*(GhatCtrl[7]-1.0*(GhatCtrl[6]+GhatCtrl[5])+GhatCtrl[4]-1.0*GhatCtrl[3]+GhatCtrl[2]+GhatCtrl[1]-1.0*GhatCtrl[0]); 

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
  return std::abs(alpha0); 
} 
