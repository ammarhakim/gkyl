#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity2x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
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
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[6]+fl[3])+alpha[0]*(1.732050807568877*fl[1]+fl[0]))*dfac_x; 
  incr[1] = -0.1767766952966368*(alpha[2]*(3.0*fl[6]+1.732050807568877*fl[3])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[2] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[11]+fl[7])+alpha[0]*(1.732050807568877*fl[5]+fl[2]))*dfac_x; 
  incr[3] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[6]+fl[3])+(1.732050807568877*fl[1]+fl[0])*alpha[2])*dfac_x; 
  incr[4] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[13]+fl[10])+alpha[0]*(1.732050807568877*fl[8]+fl[4]))*dfac_x; 
  incr[5] = -0.1767766952966368*(alpha[2]*(3.0*fl[11]+1.732050807568877*fl[7])+alpha[0]*(3.0*fl[5]+1.732050807568877*fl[2]))*dfac_x; 
  incr[6] = -0.1767766952966368*(alpha[0]*(3.0*fl[6]+1.732050807568877*fl[3])+(3.0*fl[1]+1.732050807568877*fl[0])*alpha[2])*dfac_x; 
  incr[7] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[11]+fl[7])+alpha[2]*(1.732050807568877*fl[5]+fl[2]))*dfac_x; 
  incr[8] = -0.1767766952966368*(alpha[2]*(3.0*fl[13]+1.732050807568877*fl[10])+alpha[0]*(3.0*fl[8]+1.732050807568877*fl[4]))*dfac_x; 
  incr[9] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[15]+fl[14])+alpha[0]*(1.732050807568877*fl[12]+fl[9]))*dfac_x; 
  incr[10] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[13]+fl[10])+alpha[2]*(1.732050807568877*fl[8]+fl[4]))*dfac_x; 
  incr[11] = -0.1767766952966368*(alpha[0]*(3.0*fl[11]+1.732050807568877*fl[7])+alpha[2]*(3.0*fl[5]+1.732050807568877*fl[2]))*dfac_x; 
  incr[12] = -0.1767766952966368*(alpha[2]*(3.0*fl[15]+1.732050807568877*fl[14])+alpha[0]*(3.0*fl[12]+1.732050807568877*fl[9]))*dfac_x; 
  incr[13] = -0.1767766952966368*(alpha[0]*(3.0*fl[13]+1.732050807568877*fl[10])+alpha[2]*(3.0*fl[8]+1.732050807568877*fl[4]))*dfac_x; 
  incr[14] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[15]+fl[14])+alpha[2]*(1.732050807568877*fl[12]+fl[9]))*dfac_x; 
  incr[15] = -0.1767766952966368*(alpha[0]*(3.0*fl[15]+1.732050807568877*fl[14])+alpha[2]*(3.0*fl[12]+1.732050807568877*fl[9]))*dfac_x; 
  } else { 
  incr[0] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[6]-1.0*fr[3])+alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac_x; 
  incr[1] = 0.1767766952966368*(alpha[2]*(3.0*fr[6]-1.732050807568877*fr[3])+alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac_x; 
  incr[2] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[11]-1.0*fr[7])+alpha[0]*(1.732050807568877*fr[5]-1.0*fr[2]))*dfac_x; 
  incr[3] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[6]-1.0*fr[3])+(1.732050807568877*fr[1]-1.0*fr[0])*alpha[2])*dfac_x; 
  incr[4] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[13]-1.0*fr[10])+alpha[0]*(1.732050807568877*fr[8]-1.0*fr[4]))*dfac_x; 
  incr[5] = 0.1767766952966368*(alpha[2]*(3.0*fr[11]-1.732050807568877*fr[7])+alpha[0]*(3.0*fr[5]-1.732050807568877*fr[2]))*dfac_x; 
  incr[6] = 0.1767766952966368*(alpha[0]*(3.0*fr[6]-1.732050807568877*fr[3])+(3.0*fr[1]-1.732050807568877*fr[0])*alpha[2])*dfac_x; 
  incr[7] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[11]-1.0*fr[7])+alpha[2]*(1.732050807568877*fr[5]-1.0*fr[2]))*dfac_x; 
  incr[8] = 0.1767766952966368*(alpha[2]*(3.0*fr[13]-1.732050807568877*fr[10])+alpha[0]*(3.0*fr[8]-1.732050807568877*fr[4]))*dfac_x; 
  incr[9] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[15]-1.0*fr[14])+alpha[0]*(1.732050807568877*fr[12]-1.0*fr[9]))*dfac_x; 
  incr[10] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[13]-1.0*fr[10])+alpha[2]*(1.732050807568877*fr[8]-1.0*fr[4]))*dfac_x; 
  incr[11] = 0.1767766952966368*(alpha[0]*(3.0*fr[11]-1.732050807568877*fr[7])+alpha[2]*(3.0*fr[5]-1.732050807568877*fr[2]))*dfac_x; 
  incr[12] = 0.1767766952966368*(alpha[2]*(3.0*fr[15]-1.732050807568877*fr[14])+alpha[0]*(3.0*fr[12]-1.732050807568877*fr[9]))*dfac_x; 
  incr[13] = 0.1767766952966368*(alpha[0]*(3.0*fr[13]-1.732050807568877*fr[10])+alpha[2]*(3.0*fr[8]-1.732050807568877*fr[4]))*dfac_x; 
  incr[14] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[15]-1.0*fr[14])+alpha[2]*(1.732050807568877*fr[12]-1.0*fr[9]))*dfac_x; 
  incr[15] = 0.1767766952966368*(alpha[0]*(3.0*fr[15]-1.732050807568877*fr[14])+alpha[2]*(3.0*fr[12]-1.732050807568877*fr[9]))*dfac_x; 
  }
#elif upwindType == QUAD 
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
  fhatALQuad[0] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[0] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[2] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[2] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[4] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[4] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[6] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[6] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[7] = 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
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
 
  incr[0] = 0.25*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_x; 
  incr[1] = -0.4330127018922193*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_x; 
  incr[2] = 0.25*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])*dfac_x; 
  incr[3] = 0.25*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_x; 
  incr[4] = 0.25*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])*dfac_x; 
  incr[5] = -0.4330127018922193*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])*dfac_x; 
  incr[6] = -0.4330127018922193*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_x; 
  incr[7] = 0.25*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])*dfac_x; 
  incr[8] = -0.4330127018922193*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])*dfac_x; 
  incr[9] = 0.25*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])*dfac_x; 
  incr[10] = 0.25*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])*dfac_x; 
  incr[11] = -0.4330127018922193*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])*dfac_x; 
  incr[12] = -0.4330127018922193*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])*dfac_x; 
  incr[13] = -0.4330127018922193*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])*dfac_x; 
  incr[14] = 0.25*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])*dfac_x; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])*dfac_x; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[16], outrPos[16]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.25 : positivityWeightByDirL[1]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.25 : positivityWeightByDirR[1]/positivityWeightByDirR[0]; 
  outlPos[0] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[1] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[2] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[3] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[4] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[5] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[6] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[7] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[8] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[9] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[10] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outlPos[11] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outlPos[12] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outlPos[13] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outlPos[14] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[15] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[0] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[1] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[2] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[3] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[4] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[5] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[6] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[7] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[8] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[9] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[10] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[11] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[12] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[13] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[14] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[15] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  if(outlPos[1] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12]+fl[11])+3.0*(fl[4]+fl[3]+fl[2])-3.0*fl[1])-3.0*(fl[10]+fl[9])+3.0*fl[8]-3.0*fl[7]+3.0*(fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[1]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]+fr[11]+3.0*(fr[4]+fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[1] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[1] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[3] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]-1.0*(fl[12]+fl[11])-3.0*(fl[4]+fl[3])+3.0*(fl[2]+fl[1]))+3.0*fl[10]-3.0*(fl[9]+fl[8]+fl[7]+fl[6])+3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[2] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]+fr[11]-3.0*(fr[4]+fr[3])+3.0*fr[2]-3.0*fr[1])+3.0*fr[10]-3.0*fr[9]+3.0*fr[8]-3.0*fr[7]+3.0*fr[6]-3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[2]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[2] *= limFac; 
  outrPos[3] *= limFac; 
  outrPos[2] *= limFac; 
  if(outlPos[5] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]-1.0*fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2]+3.0*fl[1])-3.0*fl[10]+3.0*fl[9]-3.0*(fl[8]+fl[7])+3.0*fl[6]-3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[5]); 
  else limFac = 1.0; 
  if(outrPos[4] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]-1.0*fr[12]+fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*(fr[2]+fr[1]))-3.0*fr[10]+3.0*(fr[9]+fr[8])-3.0*(fr[7]+fr[6])+3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[4]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[5] *= limFac; 
  outlPos[4] *= limFac; 
  outrPos[5] *= limFac; 
  outrPos[4] *= limFac; 
  if(outlPos[7] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]-1.0*fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8])-3.0*(fl[7]+fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[6] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12])+fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2])+3.0*fr[1])+3.0*(fr[10]+fr[9])-3.0*(fr[8]+fr[7])+3.0*(fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[6]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[6] *= limFac; 
  outrPos[7] *= limFac; 
  outrPos[6] *= limFac; 
  if(outlPos[9] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12])+fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2])+3.0*fl[1])-3.0*(fl[10]+fl[9])+3.0*(fl[8]+fl[7])-3.0*(fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[9]); 
  else limFac = 1.0; 
  if(outrPos[8] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]-1.0*fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8])+3.0*(fr[7]+fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[8]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[9] *= limFac; 
  outlPos[8] *= limFac; 
  outrPos[9] *= limFac; 
  outrPos[8] *= limFac; 
  if(outlPos[11] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]-1.0*fl[12]+fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*(fl[2]+fl[1]))+3.0*fl[10]-3.0*(fl[9]+fl[8])+3.0*(fl[7]+fl[6])-3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[11]); 
  else limFac = 1.0; 
  if(outrPos[10] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]-1.0*fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2]+3.0*fr[1])+3.0*fr[10]-3.0*fr[9]+3.0*(fr[8]+fr[7])-3.0*fr[6]+3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[10]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[11] *= limFac; 
  outlPos[10] *= limFac; 
  outrPos[11] *= limFac; 
  outrPos[10] *= limFac; 
  if(outlPos[13] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]+fl[11]-3.0*(fl[4]+fl[3])+3.0*fl[2]-3.0*fl[1])-3.0*fl[10]+3.0*fl[9]-3.0*fl[8]+3.0*fl[7]-3.0*fl[6]+3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[13]); 
  else limFac = 1.0; 
  if(outrPos[12] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]-1.0*(fr[12]+fr[11])-3.0*(fr[4]+fr[3])+3.0*(fr[2]+fr[1]))-3.0*fr[10]+3.0*(fr[9]+fr[8]+fr[7]+fr[6])-3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[12]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[13] *= limFac; 
  outlPos[12] *= limFac; 
  outrPos[13] *= limFac; 
  outrPos[12] *= limFac; 
  if(outlPos[15] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]+fl[11]+3.0*(fl[4]+fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[15]); 
  else limFac = 1.0; 
  if(outrPos[14] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12]+fr[11])+3.0*(fr[4]+fr[3]+fr[2])-3.0*fr[1])+3.0*(fr[10]+fr[9])-3.0*fr[8]+3.0*fr[7]-3.0*(fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[14]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[15] *= limFac; 
  outlPos[14] *= limFac; 
  outrPos[15] *= limFac; 
  outrPos[14] *= limFac; 
  outr[0] += 0.25*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.4330127018922193*(outrPos[15]-1.0*outrPos[14]+outrPos[13]-1.0*outrPos[12]+outrPos[11]-1.0*outrPos[10]+outrPos[9]-1.0*outrPos[8]+outrPos[7]-1.0*outrPos[6]+outrPos[5]-1.0*outrPos[4]+outrPos[3]-1.0*outrPos[2]+outrPos[1]-1.0*outrPos[0]); 
  outr[2] += 0.4330127018922193*(outrPos[15]+outrPos[14]-1.0*(outrPos[13]+outrPos[12])+outrPos[11]+outrPos[10]-1.0*(outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]-1.0*(outrPos[5]+outrPos[4])+outrPos[3]+outrPos[2]-1.0*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]-1.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]-1.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[4] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]-1.0*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[5] += 0.25*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*(outrPos[8]+outrPos[7])-3.0*(outrPos[6]+outrPos[5])+3.0*(outrPos[4]+outrPos[3])-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 
  outr[6] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*(outrPos[8]+outrPos[7])-3.0*outrPos[6]+3.0*outrPos[5]-3.0*(outrPos[4]+outrPos[3])+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[7] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])-3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[8] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*outrPos[12]+3.0*outrPos[11]-3.0*outrPos[10]+3.0*outrPos[9]-3.0*(outrPos[8]+outrPos[7])+3.0*outrPos[6]-3.0*outrPos[5]+3.0*outrPos[4]-3.0*outrPos[3]+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[9] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12])+3.0*(outrPos[11]+outrPos[10])-3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4])-3.0*(outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[10] += 0.25*(3.0*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12])-3.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])+3.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[11] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*outrPos[12]-3.0*outrPos[11]+3.0*(outrPos[10]+outrPos[9])-3.0*outrPos[8]+3.0*outrPos[7]-3.0*(outrPos[6]+outrPos[5])+3.0*outrPos[4]-3.0*outrPos[3]+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[12] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*outrPos[8]-3.0*outrPos[7]+3.0*(outrPos[6]+outrPos[5])-3.0*(outrPos[4]+outrPos[3])+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[13] += 0.4330127018922193*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*outrPos[8]-3.0*outrPos[7]+3.0*outrPos[6]-3.0*outrPos[5]+3.0*(outrPos[4]+outrPos[3])-3.0*outrPos[2]+3.0*outrPos[1]-3.0*outrPos[0]); 
  outr[14] += 0.4330127018922193*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8])-3.0*(outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])-3.0*(outrPos[1]+outrPos[0])); 
  outr[15] += 0.25*(9.0*outrPos[15]-9.0*(outrPos[14]+outrPos[13])+9.0*outrPos[12]-9.0*outrPos[11]+9.0*(outrPos[10]+outrPos[9])-9.0*(outrPos[8]+outrPos[7])+9.0*(outrPos[6]+outrPos[5])-9.0*outrPos[4]+9.0*outrPos[3]-9.0*(outrPos[2]+outrPos[1])+9.0*outrPos[0]); 

  outl[0] += 0.25*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.4330127018922193*(outlPos[15]-1.0*outlPos[14]+outlPos[13]-1.0*outlPos[12]+outlPos[11]-1.0*outlPos[10]+outlPos[9]-1.0*outlPos[8]+outlPos[7]-1.0*outlPos[6]+outlPos[5]-1.0*outlPos[4]+outlPos[3]-1.0*outlPos[2]+outlPos[1]-1.0*outlPos[0]); 
  outl[2] += 0.4330127018922193*(outlPos[15]+outlPos[14]-1.0*(outlPos[13]+outlPos[12])+outlPos[11]+outlPos[10]-1.0*(outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]-1.0*(outlPos[5]+outlPos[4])+outlPos[3]+outlPos[2]-1.0*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]-1.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]-1.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[4] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]-1.0*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[5] += 0.25*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*(outlPos[8]+outlPos[7])-3.0*(outlPos[6]+outlPos[5])+3.0*(outlPos[4]+outlPos[3])-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  outl[6] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*(outlPos[8]+outlPos[7])-3.0*outlPos[6]+3.0*outlPos[5]-3.0*(outlPos[4]+outlPos[3])+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[7] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])-3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[8] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*outlPos[12]+3.0*outlPos[11]-3.0*outlPos[10]+3.0*outlPos[9]-3.0*(outlPos[8]+outlPos[7])+3.0*outlPos[6]-3.0*outlPos[5]+3.0*outlPos[4]-3.0*outlPos[3]+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[9] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12])+3.0*(outlPos[11]+outlPos[10])-3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4])-3.0*(outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[10] += 0.25*(3.0*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12])-3.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])+3.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[11] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*outlPos[12]-3.0*outlPos[11]+3.0*(outlPos[10]+outlPos[9])-3.0*outlPos[8]+3.0*outlPos[7]-3.0*(outlPos[6]+outlPos[5])+3.0*outlPos[4]-3.0*outlPos[3]+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[12] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*outlPos[8]-3.0*outlPos[7]+3.0*(outlPos[6]+outlPos[5])-3.0*(outlPos[4]+outlPos[3])+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[13] += 0.4330127018922193*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*outlPos[8]-3.0*outlPos[7]+3.0*outlPos[6]-3.0*outlPos[5]+3.0*(outlPos[4]+outlPos[3])-3.0*outlPos[2]+3.0*outlPos[1]-3.0*outlPos[0]); 
  outl[14] += 0.4330127018922193*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8])-3.0*(outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])-3.0*(outlPos[1]+outlPos[0])); 
  outl[15] += 0.25*(9.0*outlPos[15]-9.0*(outlPos[14]+outlPos[13])+9.0*outlPos[12]-9.0*outlPos[11]+9.0*(outlPos[10]+outlPos[9])-9.0*(outlPos[8]+outlPos[7])+9.0*(outlPos[6]+outlPos[5])-9.0*outlPos[4]+9.0*outlPos[3]-9.0*(outlPos[2]+outlPos[1])+9.0*outlPos[0]); 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Y_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
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
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[7]+fl[3])+alpha[0]*(1.732050807568877*fl[2]+fl[0]))*dfac_y; 
  incr[1] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[11]+fl[6])+alpha[0]*(1.732050807568877*fl[5]+fl[1]))*dfac_y; 
  incr[2] = -0.1767766952966368*(alpha[2]*(3.0*fl[7]+1.732050807568877*fl[3])+alpha[0]*(3.0*fl[2]+1.732050807568877*fl[0]))*dfac_y; 
  incr[3] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[7]+fl[3])+alpha[2]*(1.732050807568877*fl[2]+fl[0]))*dfac_y; 
  incr[4] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[14]+fl[10])+alpha[0]*(1.732050807568877*fl[9]+fl[4]))*dfac_y; 
  incr[5] = -0.1767766952966368*(alpha[2]*(3.0*fl[11]+1.732050807568877*fl[6])+alpha[0]*(3.0*fl[5]+1.732050807568877*fl[1]))*dfac_y; 
  incr[6] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[11]+fl[6])+alpha[2]*(1.732050807568877*fl[5]+fl[1]))*dfac_y; 
  incr[7] = -0.1767766952966368*(alpha[0]*(3.0*fl[7]+1.732050807568877*fl[3])+alpha[2]*(3.0*fl[2]+1.732050807568877*fl[0]))*dfac_y; 
  incr[8] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[15]+fl[13])+alpha[0]*(1.732050807568877*fl[12]+fl[8]))*dfac_y; 
  incr[9] = -0.1767766952966368*(alpha[2]*(3.0*fl[14]+1.732050807568877*fl[10])+alpha[0]*(3.0*fl[9]+1.732050807568877*fl[4]))*dfac_y; 
  incr[10] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[14]+fl[10])+alpha[2]*(1.732050807568877*fl[9]+fl[4]))*dfac_y; 
  incr[11] = -0.1767766952966368*(alpha[0]*(3.0*fl[11]+1.732050807568877*fl[6])+alpha[2]*(3.0*fl[5]+1.732050807568877*fl[1]))*dfac_y; 
  incr[12] = -0.1767766952966368*(alpha[2]*(3.0*fl[15]+1.732050807568877*fl[13])+alpha[0]*(3.0*fl[12]+1.732050807568877*fl[8]))*dfac_y; 
  incr[13] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[15]+fl[13])+alpha[2]*(1.732050807568877*fl[12]+fl[8]))*dfac_y; 
  incr[14] = -0.1767766952966368*(alpha[0]*(3.0*fl[14]+1.732050807568877*fl[10])+alpha[2]*(3.0*fl[9]+1.732050807568877*fl[4]))*dfac_y; 
  incr[15] = -0.1767766952966368*(alpha[0]*(3.0*fl[15]+1.732050807568877*fl[13])+alpha[2]*(3.0*fl[12]+1.732050807568877*fl[8]))*dfac_y; 
  } else { 
  incr[0] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[7]-1.0*fr[3])+alpha[0]*(1.732050807568877*fr[2]-1.0*fr[0]))*dfac_y; 
  incr[1] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[11]-1.0*fr[6])+alpha[0]*(1.732050807568877*fr[5]-1.0*fr[1]))*dfac_y; 
  incr[2] = 0.1767766952966368*(alpha[2]*(3.0*fr[7]-1.732050807568877*fr[3])+alpha[0]*(3.0*fr[2]-1.732050807568877*fr[0]))*dfac_y; 
  incr[3] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[7]-1.0*fr[3])+alpha[2]*(1.732050807568877*fr[2]-1.0*fr[0]))*dfac_y; 
  incr[4] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[14]-1.0*fr[10])+alpha[0]*(1.732050807568877*fr[9]-1.0*fr[4]))*dfac_y; 
  incr[5] = 0.1767766952966368*(alpha[2]*(3.0*fr[11]-1.732050807568877*fr[6])+alpha[0]*(3.0*fr[5]-1.732050807568877*fr[1]))*dfac_y; 
  incr[6] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[11]-1.0*fr[6])+alpha[2]*(1.732050807568877*fr[5]-1.0*fr[1]))*dfac_y; 
  incr[7] = 0.1767766952966368*(alpha[0]*(3.0*fr[7]-1.732050807568877*fr[3])+alpha[2]*(3.0*fr[2]-1.732050807568877*fr[0]))*dfac_y; 
  incr[8] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[15]-1.0*fr[13])+alpha[0]*(1.732050807568877*fr[12]-1.0*fr[8]))*dfac_y; 
  incr[9] = 0.1767766952966368*(alpha[2]*(3.0*fr[14]-1.732050807568877*fr[10])+alpha[0]*(3.0*fr[9]-1.732050807568877*fr[4]))*dfac_y; 
  incr[10] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[14]-1.0*fr[10])+alpha[2]*(1.732050807568877*fr[9]-1.0*fr[4]))*dfac_y; 
  incr[11] = 0.1767766952966368*(alpha[0]*(3.0*fr[11]-1.732050807568877*fr[6])+alpha[2]*(3.0*fr[5]-1.732050807568877*fr[1]))*dfac_y; 
  incr[12] = 0.1767766952966368*(alpha[2]*(3.0*fr[15]-1.732050807568877*fr[13])+alpha[0]*(3.0*fr[12]-1.732050807568877*fr[8]))*dfac_y; 
  incr[13] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[15]-1.0*fr[13])+alpha[2]*(1.732050807568877*fr[12]-1.0*fr[8]))*dfac_y; 
  incr[14] = 0.1767766952966368*(alpha[0]*(3.0*fr[14]-1.732050807568877*fr[10])+alpha[2]*(3.0*fr[9]-1.732050807568877*fr[4]))*dfac_y; 
  incr[15] = 0.1767766952966368*(alpha[0]*(3.0*fr[15]-1.732050807568877*fr[13])+alpha[2]*(3.0*fr[12]-1.732050807568877*fr[8]))*dfac_y; 
  }
#elif upwindType == QUAD 
  double rCtrlL[8], rCtrlR[8];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrlL[0] = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]+fl[12]+fl[11]+3.0*fl[2])+9.0*(fl[9]+fl[7]+fl[5])))/(36.0*EPSILON-1.732050807568877*(fl[13]+3.0*(fl[4]+fl[3]+fl[1]))+3.0*(fl[10]+fl[8]+fl[6]+3.0*fl[0])); 
  rCtrlL[1] = (3.0*fl[15]+5.196152422706631*(fl[14]-1.0*(fl[12]+fl[11])+3.0*fl[2])+9.0*(fl[5]-1.0*(fl[9]+fl[7])))/(36.0*EPSILON+1.732050807568877*(fl[13]+3.0*(fl[1]-1.0*(fl[4]+fl[3])))+3.0*(fl[10]-1.0*(fl[8]+fl[6]-3.0*fl[0]))); 
  rCtrlL[2] = (3.0*fl[15]-5.196152422706631*(fl[14]-1.0*fl[12]+fl[11]-3.0*fl[2])+9.0*((-1.0*fl[9])+fl[7]-1.0*fl[5]))/(36.0*EPSILON+1.732050807568877*(fl[13]+3.0*((-1.0*fl[4])+fl[3]-1.0*fl[1]))+3.0*((-1.0*fl[10])+fl[8]-1.0*fl[6]+3.0*fl[0])); 
  rCtrlL[3] = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]+fl[12]-1.0*(fl[11]+3.0*fl[2]))+9.0*(fl[9]-1.0*(fl[7]+fl[5]))))/(36.0*EPSILON-1.732050807568877*(fl[13]+3.0*(fl[4]-1.0*(fl[3]+fl[1])))+3.0*((-1.0*(fl[10]+fl[8]))+fl[6]+3.0*fl[0])); 
  rCtrlL[4] = (3.0*fl[15]-5.196152422706631*(fl[14]+fl[12]-1.0*(fl[11]+3.0*fl[2]))+9.0*(fl[9]-1.0*(fl[7]+fl[5])))/(36.0*EPSILON+1.732050807568877*(fl[13]+3.0*(fl[4]-1.0*(fl[3]+fl[1])))+3.0*((-1.0*(fl[10]+fl[8]))+fl[6]+3.0*fl[0])); 
  rCtrlL[5] = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]-1.0*fl[12]+fl[11]-3.0*fl[2])+9.0*((-1.0*fl[9])+fl[7]-1.0*fl[5])))/(36.0*EPSILON-1.732050807568877*(fl[13]+3.0*((-1.0*fl[4])+fl[3]-1.0*fl[1]))+3.0*((-1.0*fl[10])+fl[8]-1.0*fl[6]+3.0*fl[0])); 
  rCtrlL[6] = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]-1.0*(fl[12]+fl[11])+3.0*fl[2])+9.0*(fl[5]-1.0*(fl[9]+fl[7]))))/(36.0*EPSILON-1.732050807568877*(fl[13]+3.0*(fl[1]-1.0*(fl[4]+fl[3])))+3.0*(fl[10]-1.0*(fl[8]+fl[6]-3.0*fl[0]))); 
  rCtrlL[7] = (3.0*fl[15]+5.196152422706631*(fl[14]+fl[12]+fl[11]+3.0*fl[2])+9.0*(fl[9]+fl[7]+fl[5]))/(36.0*EPSILON+1.732050807568877*(fl[13]+3.0*(fl[4]+fl[3]+fl[1]))+3.0*(fl[10]+fl[8]+fl[6]+3.0*fl[0])); 
  rCtrlR[0] = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]+fr[12]+fr[11]+3.0*fr[2])+9.0*(fr[9]+fr[7]+fr[5])))/(36.0*EPSILON-1.732050807568877*(fr[13]+3.0*(fr[4]+fr[3]+fr[1]))+3.0*(fr[10]+fr[8]+fr[6]+3.0*fr[0])); 
  rCtrlR[1] = (3.0*fr[15]+5.196152422706631*(fr[14]-1.0*(fr[12]+fr[11])+3.0*fr[2])+9.0*(fr[5]-1.0*(fr[9]+fr[7])))/(36.0*EPSILON+1.732050807568877*(fr[13]+3.0*(fr[1]-1.0*(fr[4]+fr[3])))+3.0*(fr[10]-1.0*(fr[8]+fr[6]-3.0*fr[0]))); 
  rCtrlR[2] = (3.0*fr[15]-5.196152422706631*(fr[14]-1.0*fr[12]+fr[11]-3.0*fr[2])+9.0*((-1.0*fr[9])+fr[7]-1.0*fr[5]))/(36.0*EPSILON+1.732050807568877*(fr[13]+3.0*((-1.0*fr[4])+fr[3]-1.0*fr[1]))+3.0*((-1.0*fr[10])+fr[8]-1.0*fr[6]+3.0*fr[0])); 
  rCtrlR[3] = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]+fr[12]-1.0*(fr[11]+3.0*fr[2]))+9.0*(fr[9]-1.0*(fr[7]+fr[5]))))/(36.0*EPSILON-1.732050807568877*(fr[13]+3.0*(fr[4]-1.0*(fr[3]+fr[1])))+3.0*((-1.0*(fr[10]+fr[8]))+fr[6]+3.0*fr[0])); 
  rCtrlR[4] = (3.0*fr[15]-5.196152422706631*(fr[14]+fr[12]-1.0*(fr[11]+3.0*fr[2]))+9.0*(fr[9]-1.0*(fr[7]+fr[5])))/(36.0*EPSILON+1.732050807568877*(fr[13]+3.0*(fr[4]-1.0*(fr[3]+fr[1])))+3.0*((-1.0*(fr[10]+fr[8]))+fr[6]+3.0*fr[0])); 
  rCtrlR[5] = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]-1.0*fr[12]+fr[11]-3.0*fr[2])+9.0*((-1.0*fr[9])+fr[7]-1.0*fr[5])))/(36.0*EPSILON-1.732050807568877*(fr[13]+3.0*((-1.0*fr[4])+fr[3]-1.0*fr[1]))+3.0*((-1.0*fr[10])+fr[8]-1.0*fr[6]+3.0*fr[0])); 
  rCtrlR[6] = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]-1.0*(fr[12]+fr[11])+3.0*fr[2])+9.0*(fr[5]-1.0*(fr[9]+fr[7]))))/(36.0*EPSILON-1.732050807568877*(fr[13]+3.0*(fr[1]-1.0*(fr[4]+fr[3])))+3.0*(fr[10]-1.0*(fr[8]+fr[6]-3.0*fr[0]))); 
  rCtrlR[7] = (3.0*fr[15]+5.196152422706631*(fr[14]+fr[12]+fr[11]+3.0*fr[2])+9.0*(fr[9]+fr[7]+fr[5]))/(36.0*EPSILON+1.732050807568877*(fr[13]+3.0*(fr[4]+fr[3]+fr[1]))+3.0*(fr[10]+fr[8]+fr[6]+3.0*fr[0])); 
  double fCtrlL[8], fCtrlR[8];  // fCtrl = anti-limited f evaluated at each control node on y surface 
  // control node [x,vx,vy] = [-1/3,-1/3,-1/3] 
  fCtrlL[0] = -0.02777777777777778*(1.732050807568877*(fl[13]+3.0*(fl[4]+fl[3]+fl[1]))-3.0*(fl[10]+fl[8]+fl[6])-9.0*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.02777777777777778*(1.732050807568877*(fr[13]+3.0*(fr[4]+fr[3]+fr[1]))-3.0*(fr[10]+fr[8]+fr[6])-9.0*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,vx,vy] = [1/3,-1/3,-1/3] 
  fCtrlL[1] = 0.02777777777777778*(1.732050807568877*(fl[13]-3.0*(fl[4]+fl[3])+3.0*fl[1])+3.0*fl[10]-3.0*(fl[8]+fl[6])+9.0*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.02777777777777778*(1.732050807568877*(fr[13]-3.0*(fr[4]+fr[3])+3.0*fr[1])+3.0*fr[10]-3.0*(fr[8]+fr[6])+9.0*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,vx,vy] = [-1/3,1/3,-1/3] 
  fCtrlL[2] = 0.02777777777777778*(1.732050807568877*(fl[13]-3.0*fl[4]+3.0*fl[3]-3.0*fl[1])-3.0*fl[10]+3.0*fl[8]-3.0*fl[6]+9.0*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = 0.02777777777777778*(1.732050807568877*(fr[13]-3.0*fr[4]+3.0*fr[3]-3.0*fr[1])-3.0*fr[10]+3.0*fr[8]-3.0*fr[6]+9.0*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,vx,vy] = [1/3,1/3,-1/3] 
  fCtrlL[3] = -0.02777777777777778*(1.732050807568877*(fl[13]+3.0*fl[4]-3.0*(fl[3]+fl[1]))+3.0*(fl[10]+fl[8])-3.0*fl[6]-9.0*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = -0.02777777777777778*(1.732050807568877*(fr[13]+3.0*fr[4]-3.0*(fr[3]+fr[1]))+3.0*(fr[10]+fr[8])-3.0*fr[6]-9.0*fr[0])*limTheta(rCtrlR[3],-1.0); 
  // control node [x,vx,vy] = [-1/3,-1/3,1/3] 
  fCtrlL[4] = 0.02777777777777778*(1.732050807568877*(fl[13]+3.0*fl[4]-3.0*(fl[3]+fl[1]))-3.0*(fl[10]+fl[8])+3.0*fl[6]+9.0*fl[0])*limTheta(rCtrlL[4],1.0); 
  fCtrlR[4] = 0.02777777777777778*(1.732050807568877*(fr[13]+3.0*fr[4]-3.0*(fr[3]+fr[1]))-3.0*(fr[10]+fr[8])+3.0*fr[6]+9.0*fr[0])*limTheta(rCtrlR[4],-1.0); 
  // control node [x,vx,vy] = [1/3,-1/3,1/3] 
  fCtrlL[5] = -0.02777777777777778*(1.732050807568877*(fl[13]-3.0*fl[4]+3.0*fl[3]-3.0*fl[1])+3.0*fl[10]-3.0*fl[8]+3.0*fl[6]-9.0*fl[0])*limTheta(rCtrlL[5],1.0); 
  fCtrlR[5] = -0.02777777777777778*(1.732050807568877*(fr[13]-3.0*fr[4]+3.0*fr[3]-3.0*fr[1])+3.0*fr[10]-3.0*fr[8]+3.0*fr[6]-9.0*fr[0])*limTheta(rCtrlR[5],-1.0); 
  // control node [x,vx,vy] = [-1/3,1/3,1/3] 
  fCtrlL[6] = -0.02777777777777778*(1.732050807568877*(fl[13]-3.0*(fl[4]+fl[3])+3.0*fl[1])-3.0*fl[10]+3.0*(fl[8]+fl[6])-9.0*fl[0])*limTheta(rCtrlL[6],1.0); 
  fCtrlR[6] = -0.02777777777777778*(1.732050807568877*(fr[13]-3.0*(fr[4]+fr[3])+3.0*fr[1])-3.0*fr[10]+3.0*(fr[8]+fr[6])-9.0*fr[0])*limTheta(rCtrlR[6],-1.0); 
  // control node [x,vx,vy] = [1/3,1/3,1/3] 
  fCtrlL[7] = 0.02777777777777778*(1.732050807568877*(fl[13]+3.0*(fl[4]+fl[3]+fl[1]))+3.0*(fl[10]+fl[8]+fl[6])+9.0*fl[0])*limTheta(rCtrlL[7],1.0); 
  fCtrlR[7] = 0.02777777777777778*(1.732050807568877*(fr[13]+3.0*(fr[4]+fr[3]+fr[1]))+3.0*(fr[10]+fr[8]+fr[6])+9.0*fr[0])*limTheta(rCtrlR[7],-1.0); 
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
  fhatALQuad[0] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[0] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[2] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[2] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[4] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[4] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[6] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[6] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[7] = 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
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
 
  incr[0] = 0.25*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_y; 
  incr[1] = 0.25*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])*dfac_y; 
  incr[2] = -0.4330127018922193*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_y; 
  incr[3] = 0.25*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_y; 
  incr[4] = 0.25*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])*dfac_y; 
  incr[5] = -0.4330127018922193*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])*dfac_y; 
  incr[6] = 0.25*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])*dfac_y; 
  incr[7] = -0.4330127018922193*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_y; 
  incr[8] = 0.25*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])*dfac_y; 
  incr[9] = -0.4330127018922193*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])*dfac_y; 
  incr[10] = 0.25*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])*dfac_y; 
  incr[11] = -0.4330127018922193*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])*dfac_y; 
  incr[12] = -0.4330127018922193*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])*dfac_y; 
  incr[13] = 0.25*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])*dfac_y; 
  incr[14] = -0.4330127018922193*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])*dfac_y; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])*dfac_y; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[16], outrPos[16]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.25 : positivityWeightByDirL[2]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.25 : positivityWeightByDirR[2]/positivityWeightByDirR[0]; 
  outlPos[0] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[1] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[2] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[3] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[4] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[5] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[6] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[7] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[8] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outlPos[9] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outlPos[10] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[11] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[12] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[13] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[14] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outlPos[15] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[0] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[1] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[2] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[3] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[4] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[5] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[6] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[7] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[8] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[9] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[10] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[11] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[12] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[13] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[14] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[15] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  if(outlPos[2] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]+fl[11]-3.0*(fl[4]+fl[3])+3.0*fl[2]-3.0*fl[1])+3.0*fl[10]-3.0*fl[9]+3.0*fl[8]-3.0*fl[7]+3.0*fl[6]-3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]+fr[11]+3.0*(fr[4]+fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[2] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[3] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]-1.0*(fl[12]+fl[11])-3.0*(fl[4]+fl[3])+3.0*(fl[2]+fl[1]))+3.0*fl[10]-3.0*(fl[9]+fl[8]+fl[7]+fl[6])+3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[1] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12]+fr[11])+3.0*(fr[4]+fr[3]+fr[2])-3.0*fr[1])-3.0*(fr[10]+fr[9])+3.0*fr[8]-3.0*fr[7]+3.0*(fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[1] *= limFac; 
  outrPos[3] *= limFac; 
  outrPos[1] *= limFac; 
  if(outlPos[6] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12])+fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2])+3.0*fl[1])+3.0*(fl[10]+fl[9])-3.0*(fl[8]+fl[7])+3.0*(fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[6]); 
  else limFac = 1.0; 
  if(outrPos[4] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]-1.0*fr[12]+fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*(fr[2]+fr[1]))-3.0*fr[10]+3.0*(fr[9]+fr[8])-3.0*(fr[7]+fr[6])+3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[4]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[6] *= limFac; 
  outlPos[4] *= limFac; 
  outrPos[6] *= limFac; 
  outrPos[4] *= limFac; 
  if(outlPos[7] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]-1.0*fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8])-3.0*(fl[7]+fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[5] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]-1.0*fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2]+3.0*fr[1])-3.0*fr[10]+3.0*fr[9]-3.0*(fr[8]+fr[7])+3.0*fr[6]-3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[5]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[5] *= limFac; 
  outrPos[7] *= limFac; 
  outrPos[5] *= limFac; 
  if(outlPos[10] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]-1.0*fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2]+3.0*fl[1])+3.0*fl[10]-3.0*fl[9]+3.0*(fl[8]+fl[7])-3.0*fl[6]+3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[10]); 
  else limFac = 1.0; 
  if(outrPos[8] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]-1.0*fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8])+3.0*(fr[7]+fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[8]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[10] *= limFac; 
  outlPos[8] *= limFac; 
  outrPos[10] *= limFac; 
  outrPos[8] *= limFac; 
  if(outlPos[11] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]-1.0*fl[12]+fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*(fl[2]+fl[1]))+3.0*fl[10]-3.0*(fl[9]+fl[8])+3.0*(fl[7]+fl[6])-3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[11]); 
  else limFac = 1.0; 
  if(outrPos[9] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12])+fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2])+3.0*fr[1])-3.0*(fr[10]+fr[9])+3.0*(fr[8]+fr[7])-3.0*(fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[9]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[11] *= limFac; 
  outlPos[9] *= limFac; 
  outrPos[11] *= limFac; 
  outrPos[9] *= limFac; 
  if(outlPos[14] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12]+fl[11])+3.0*(fl[4]+fl[3]+fl[2])-3.0*fl[1])+3.0*(fl[10]+fl[9])-3.0*fl[8]+3.0*fl[7]-3.0*(fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[14]); 
  else limFac = 1.0; 
  if(outrPos[12] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]-1.0*(fr[12]+fr[11])-3.0*(fr[4]+fr[3])+3.0*(fr[2]+fr[1]))-3.0*fr[10]+3.0*(fr[9]+fr[8]+fr[7]+fr[6])-3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[12]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[14] *= limFac; 
  outlPos[12] *= limFac; 
  outrPos[14] *= limFac; 
  outrPos[12] *= limFac; 
  if(outlPos[15] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]+fl[11]+3.0*(fl[4]+fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[15]); 
  else limFac = 1.0; 
  if(outrPos[13] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]+fr[11]-3.0*(fr[4]+fr[3])+3.0*fr[2]-3.0*fr[1])-3.0*fr[10]+3.0*fr[9]-3.0*fr[8]+3.0*fr[7]-3.0*fr[6]+3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[13]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[15] *= limFac; 
  outlPos[13] *= limFac; 
  outrPos[15] *= limFac; 
  outrPos[13] *= limFac; 
  outr[0] += 0.25*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.4330127018922193*(outrPos[15]-1.0*outrPos[14]+outrPos[13]-1.0*outrPos[12]+outrPos[11]-1.0*outrPos[10]+outrPos[9]-1.0*outrPos[8]+outrPos[7]-1.0*outrPos[6]+outrPos[5]-1.0*outrPos[4]+outrPos[3]-1.0*outrPos[2]+outrPos[1]-1.0*outrPos[0]); 
  outr[2] += 0.4330127018922193*(outrPos[15]+outrPos[14]-1.0*(outrPos[13]+outrPos[12])+outrPos[11]+outrPos[10]-1.0*(outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]-1.0*(outrPos[5]+outrPos[4])+outrPos[3]+outrPos[2]-1.0*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]-1.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]-1.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[4] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]-1.0*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[5] += 0.25*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*(outrPos[8]+outrPos[7])-3.0*(outrPos[6]+outrPos[5])+3.0*(outrPos[4]+outrPos[3])-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 
  outr[6] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*(outrPos[8]+outrPos[7])-3.0*outrPos[6]+3.0*outrPos[5]-3.0*(outrPos[4]+outrPos[3])+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[7] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])-3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[8] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*outrPos[12]+3.0*outrPos[11]-3.0*outrPos[10]+3.0*outrPos[9]-3.0*(outrPos[8]+outrPos[7])+3.0*outrPos[6]-3.0*outrPos[5]+3.0*outrPos[4]-3.0*outrPos[3]+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[9] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12])+3.0*(outrPos[11]+outrPos[10])-3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4])-3.0*(outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[10] += 0.25*(3.0*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12])-3.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])+3.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[11] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*outrPos[12]-3.0*outrPos[11]+3.0*(outrPos[10]+outrPos[9])-3.0*outrPos[8]+3.0*outrPos[7]-3.0*(outrPos[6]+outrPos[5])+3.0*outrPos[4]-3.0*outrPos[3]+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[12] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*outrPos[8]-3.0*outrPos[7]+3.0*(outrPos[6]+outrPos[5])-3.0*(outrPos[4]+outrPos[3])+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[13] += 0.4330127018922193*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*outrPos[8]-3.0*outrPos[7]+3.0*outrPos[6]-3.0*outrPos[5]+3.0*(outrPos[4]+outrPos[3])-3.0*outrPos[2]+3.0*outrPos[1]-3.0*outrPos[0]); 
  outr[14] += 0.4330127018922193*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8])-3.0*(outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])-3.0*(outrPos[1]+outrPos[0])); 
  outr[15] += 0.25*(9.0*outrPos[15]-9.0*(outrPos[14]+outrPos[13])+9.0*outrPos[12]-9.0*outrPos[11]+9.0*(outrPos[10]+outrPos[9])-9.0*(outrPos[8]+outrPos[7])+9.0*(outrPos[6]+outrPos[5])-9.0*outrPos[4]+9.0*outrPos[3]-9.0*(outrPos[2]+outrPos[1])+9.0*outrPos[0]); 

  outl[0] += 0.25*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.4330127018922193*(outlPos[15]-1.0*outlPos[14]+outlPos[13]-1.0*outlPos[12]+outlPos[11]-1.0*outlPos[10]+outlPos[9]-1.0*outlPos[8]+outlPos[7]-1.0*outlPos[6]+outlPos[5]-1.0*outlPos[4]+outlPos[3]-1.0*outlPos[2]+outlPos[1]-1.0*outlPos[0]); 
  outl[2] += 0.4330127018922193*(outlPos[15]+outlPos[14]-1.0*(outlPos[13]+outlPos[12])+outlPos[11]+outlPos[10]-1.0*(outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]-1.0*(outlPos[5]+outlPos[4])+outlPos[3]+outlPos[2]-1.0*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]-1.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]-1.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[4] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]-1.0*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[5] += 0.25*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*(outlPos[8]+outlPos[7])-3.0*(outlPos[6]+outlPos[5])+3.0*(outlPos[4]+outlPos[3])-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  outl[6] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*(outlPos[8]+outlPos[7])-3.0*outlPos[6]+3.0*outlPos[5]-3.0*(outlPos[4]+outlPos[3])+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[7] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])-3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[8] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*outlPos[12]+3.0*outlPos[11]-3.0*outlPos[10]+3.0*outlPos[9]-3.0*(outlPos[8]+outlPos[7])+3.0*outlPos[6]-3.0*outlPos[5]+3.0*outlPos[4]-3.0*outlPos[3]+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[9] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12])+3.0*(outlPos[11]+outlPos[10])-3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4])-3.0*(outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[10] += 0.25*(3.0*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12])-3.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])+3.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[11] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*outlPos[12]-3.0*outlPos[11]+3.0*(outlPos[10]+outlPos[9])-3.0*outlPos[8]+3.0*outlPos[7]-3.0*(outlPos[6]+outlPos[5])+3.0*outlPos[4]-3.0*outlPos[3]+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[12] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*outlPos[8]-3.0*outlPos[7]+3.0*(outlPos[6]+outlPos[5])-3.0*(outlPos[4]+outlPos[3])+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[13] += 0.4330127018922193*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*outlPos[8]-3.0*outlPos[7]+3.0*outlPos[6]-3.0*outlPos[5]+3.0*(outlPos[4]+outlPos[3])-3.0*outlPos[2]+3.0*outlPos[1]-3.0*outlPos[0]); 
  outl[14] += 0.4330127018922193*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8])-3.0*(outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])-3.0*(outlPos[1]+outlPos[0])); 
  outl[15] += 0.25*(9.0*outlPos[15]-9.0*(outlPos[14]+outlPos[13])+9.0*outlPos[12]-9.0*outlPos[11]+9.0*(outlPos[10]+outlPos[9])-9.0*(outlPos[8]+outlPos[7])+9.0*(outlPos[6]+outlPos[5])-9.0*outlPos[4]+9.0*outlPos[3]-9.0*(outlPos[2]+outlPos[1])+9.0*outlPos[0]); 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
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
  alpha[0] = -(1.224744871391589*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*(dfac_v*wv-1.0))/dfac_v; 
  alpha[1] = -(1.224744871391589*BdriftY[0]*Phi[3]*dfac_y*(dfac_v*wv-1.0))/dfac_v; 
  alpha[2] = -(1.224744871391589*BdriftX[0]*Phi[3]*dfac_x*(dfac_v*wv-1.0))/dfac_v; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[7]+alpha[1]*fl[6]+alpha[0]*fl[3])+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[11]+alpha[0]*fl[6])+alpha[2]*fl[5]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[11]+alpha[0]*fl[7])+alpha[1]*fl[5]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.1767766952966368*(3.0*(alpha[2]*fl[7]+alpha[1]*fl[6]+alpha[0]*fl[3])+1.732050807568877*(alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[14]+alpha[1]*fl[13]+alpha[0]*fl[10])+alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[11]+alpha[1]*fl[7]+alpha[2]*fl[6])+alpha[0]*fl[5]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*(alpha[2]*fl[11]+alpha[0]*fl[6])+1.732050807568877*alpha[2]*fl[5]+3.0*alpha[1]*fl[3]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[1]*fl[11]+alpha[0]*fl[7])+1.732050807568877*alpha[1]*fl[5]+3.0*alpha[2]*fl[3]+1.732050807568877*(alpha[0]*fl[2]+fl[0]*alpha[2]))*dfac_v; 
  incr[8] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[15]+alpha[0]*fl[13])+alpha[2]*fl[12]+1.732050807568877*alpha[1]*fl[10]+alpha[0]*fl[8]+alpha[1]*fl[4])*dfac_v; 
  incr[9] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14])+alpha[1]*fl[12]+1.732050807568877*alpha[2]*fl[10]+alpha[0]*fl[9]+alpha[2]*fl[4])*dfac_v; 
  incr[10] = -0.1767766952966368*(3.0*(alpha[2]*fl[14]+alpha[1]*fl[13]+alpha[0]*fl[10])+1.732050807568877*(alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4]))*dfac_v; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[0]*fl[11]+alpha[1]*fl[7]+alpha[2]*fl[6])+1.732050807568877*(alpha[0]*fl[5]+alpha[1]*fl[2]+fl[1]*alpha[2]))*dfac_v; 
  incr[12] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+alpha[0]*fl[12]+alpha[1]*fl[9]+alpha[2]*fl[8])*dfac_v; 
  incr[13] = -0.1767766952966368*(3.0*(alpha[2]*fl[15]+alpha[0]*fl[13])+1.732050807568877*alpha[2]*fl[12]+3.0*alpha[1]*fl[10]+1.732050807568877*(alpha[0]*fl[8]+alpha[1]*fl[4]))*dfac_v; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14])+1.732050807568877*alpha[1]*fl[12]+3.0*alpha[2]*fl[10]+1.732050807568877*(alpha[0]*fl[9]+alpha[2]*fl[4]))*dfac_v; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+1.732050807568877*(alpha[0]*fl[12]+alpha[1]*fl[9]+alpha[2]*fl[8]))*dfac_v; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[7]+alpha[1]*fr[6]+alpha[0]*fr[3])-1.0*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[11]+alpha[0]*fr[6])-1.0*alpha[2]*fr[5]+1.732050807568877*alpha[1]*fr[3]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[11]+alpha[0]*fr[7])-1.0*alpha[1]*fr[5]+1.732050807568877*alpha[2]*fr[3]-1.0*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[2]*fr[7]+alpha[1]*fr[6]+alpha[0]*fr[3])-1.732050807568877*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[14]+alpha[1]*fr[13]+alpha[0]*fr[10])-1.0*(alpha[2]*fr[9]+alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[11]+alpha[1]*fr[7]+alpha[2]*fr[6])-1.0*(alpha[0]*fr[5]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[2]*fr[11]+alpha[0]*fr[6])-1.732050807568877*alpha[2]*fr[5]+3.0*alpha[1]*fr[3]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[1]*fr[11]+alpha[0]*fr[7])-1.732050807568877*alpha[1]*fr[5]+3.0*alpha[2]*fr[3]-1.732050807568877*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[15]+alpha[0]*fr[13])-1.0*alpha[2]*fr[12]+1.732050807568877*alpha[1]*fr[10]-1.0*(alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14])-1.0*alpha[1]*fr[12]+1.732050807568877*alpha[2]*fr[10]-1.0*(alpha[0]*fr[9]+alpha[2]*fr[4]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[2]*fr[14]+alpha[1]*fr[13]+alpha[0]*fr[10])-1.732050807568877*(alpha[2]*fr[9]+alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[0]*fr[11]+alpha[1]*fr[7]+alpha[2]*fr[6])-1.732050807568877*(alpha[0]*fr[5]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.0*(alpha[0]*fr[12]+alpha[1]*fr[9]+alpha[2]*fr[8]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fr[15]+alpha[0]*fr[13])-1.732050807568877*alpha[2]*fr[12]+3.0*alpha[1]*fr[10]-1.732050807568877*(alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14])-1.732050807568877*alpha[1]*fr[12]+3.0*alpha[2]*fr[10]-1.732050807568877*(alpha[0]*fr[9]+alpha[2]*fr[4]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.732050807568877*(alpha[0]*fr[12]+alpha[1]*fr[9]+alpha[2]*fr[8]))*dfac_v; 
  }
#elif upwindType == QUAD 
  double rCtrlL[8], rCtrlR[8];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrlL[0] = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]+fl[13]+fl[11]+3.0*fl[3])+9.0*(fl[10]+fl[7]+fl[6])))/(36.0*EPSILON-1.732050807568877*(fl[12]+3.0*(fl[4]+fl[2]+fl[1]))+3.0*(fl[9]+fl[8]+fl[5]+3.0*fl[0])); 
  rCtrlL[1] = (3.0*fl[15]+5.196152422706631*(fl[14]-1.0*(fl[13]+fl[11])+3.0*fl[3])+9.0*(fl[6]-1.0*(fl[10]+fl[7])))/(36.0*EPSILON+1.732050807568877*(fl[12]+3.0*(fl[1]-1.0*(fl[4]+fl[2])))+3.0*(fl[9]-1.0*(fl[8]+fl[5]-3.0*fl[0]))); 
  rCtrlL[2] = (3.0*fl[15]-5.196152422706631*(fl[14]-1.0*fl[13]+fl[11]-3.0*fl[3])+9.0*((-1.0*fl[10])+fl[7]-1.0*fl[6]))/(36.0*EPSILON+1.732050807568877*(fl[12]+3.0*((-1.0*fl[4])+fl[2]-1.0*fl[1]))+3.0*((-1.0*fl[9])+fl[8]-1.0*fl[5]+3.0*fl[0])); 
  rCtrlL[3] = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]+fl[13]-1.0*(fl[11]+3.0*fl[3]))+9.0*(fl[10]-1.0*(fl[7]+fl[6]))))/(36.0*EPSILON-1.732050807568877*(fl[12]+3.0*(fl[4]-1.0*(fl[2]+fl[1])))+3.0*((-1.0*(fl[9]+fl[8]))+fl[5]+3.0*fl[0])); 
  rCtrlL[4] = (3.0*fl[15]-5.196152422706631*(fl[14]+fl[13]-1.0*(fl[11]+3.0*fl[3]))+9.0*(fl[10]-1.0*(fl[7]+fl[6])))/(36.0*EPSILON+1.732050807568877*(fl[12]+3.0*(fl[4]-1.0*(fl[2]+fl[1])))+3.0*((-1.0*(fl[9]+fl[8]))+fl[5]+3.0*fl[0])); 
  rCtrlL[5] = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]-1.0*fl[13]+fl[11]-3.0*fl[3])+9.0*((-1.0*fl[10])+fl[7]-1.0*fl[6])))/(36.0*EPSILON-1.732050807568877*(fl[12]+3.0*((-1.0*fl[4])+fl[2]-1.0*fl[1]))+3.0*((-1.0*fl[9])+fl[8]-1.0*fl[5]+3.0*fl[0])); 
  rCtrlL[6] = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]-1.0*(fl[13]+fl[11])+3.0*fl[3])+9.0*(fl[6]-1.0*(fl[10]+fl[7]))))/(36.0*EPSILON-1.732050807568877*(fl[12]+3.0*(fl[1]-1.0*(fl[4]+fl[2])))+3.0*(fl[9]-1.0*(fl[8]+fl[5]-3.0*fl[0]))); 
  rCtrlL[7] = (3.0*fl[15]+5.196152422706631*(fl[14]+fl[13]+fl[11]+3.0*fl[3])+9.0*(fl[10]+fl[7]+fl[6]))/(36.0*EPSILON+1.732050807568877*(fl[12]+3.0*(fl[4]+fl[2]+fl[1]))+3.0*(fl[9]+fl[8]+fl[5]+3.0*fl[0])); 
  rCtrlR[0] = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]+fr[13]+fr[11]+3.0*fr[3])+9.0*(fr[10]+fr[7]+fr[6])))/(36.0*EPSILON-1.732050807568877*(fr[12]+3.0*(fr[4]+fr[2]+fr[1]))+3.0*(fr[9]+fr[8]+fr[5]+3.0*fr[0])); 
  rCtrlR[1] = (3.0*fr[15]+5.196152422706631*(fr[14]-1.0*(fr[13]+fr[11])+3.0*fr[3])+9.0*(fr[6]-1.0*(fr[10]+fr[7])))/(36.0*EPSILON+1.732050807568877*(fr[12]+3.0*(fr[1]-1.0*(fr[4]+fr[2])))+3.0*(fr[9]-1.0*(fr[8]+fr[5]-3.0*fr[0]))); 
  rCtrlR[2] = (3.0*fr[15]-5.196152422706631*(fr[14]-1.0*fr[13]+fr[11]-3.0*fr[3])+9.0*((-1.0*fr[10])+fr[7]-1.0*fr[6]))/(36.0*EPSILON+1.732050807568877*(fr[12]+3.0*((-1.0*fr[4])+fr[2]-1.0*fr[1]))+3.0*((-1.0*fr[9])+fr[8]-1.0*fr[5]+3.0*fr[0])); 
  rCtrlR[3] = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]+fr[13]-1.0*(fr[11]+3.0*fr[3]))+9.0*(fr[10]-1.0*(fr[7]+fr[6]))))/(36.0*EPSILON-1.732050807568877*(fr[12]+3.0*(fr[4]-1.0*(fr[2]+fr[1])))+3.0*((-1.0*(fr[9]+fr[8]))+fr[5]+3.0*fr[0])); 
  rCtrlR[4] = (3.0*fr[15]-5.196152422706631*(fr[14]+fr[13]-1.0*(fr[11]+3.0*fr[3]))+9.0*(fr[10]-1.0*(fr[7]+fr[6])))/(36.0*EPSILON+1.732050807568877*(fr[12]+3.0*(fr[4]-1.0*(fr[2]+fr[1])))+3.0*((-1.0*(fr[9]+fr[8]))+fr[5]+3.0*fr[0])); 
  rCtrlR[5] = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]-1.0*fr[13]+fr[11]-3.0*fr[3])+9.0*((-1.0*fr[10])+fr[7]-1.0*fr[6])))/(36.0*EPSILON-1.732050807568877*(fr[12]+3.0*((-1.0*fr[4])+fr[2]-1.0*fr[1]))+3.0*((-1.0*fr[9])+fr[8]-1.0*fr[5]+3.0*fr[0])); 
  rCtrlR[6] = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]-1.0*(fr[13]+fr[11])+3.0*fr[3])+9.0*(fr[6]-1.0*(fr[10]+fr[7]))))/(36.0*EPSILON-1.732050807568877*(fr[12]+3.0*(fr[1]-1.0*(fr[4]+fr[2])))+3.0*(fr[9]-1.0*(fr[8]+fr[5]-3.0*fr[0]))); 
  rCtrlR[7] = (3.0*fr[15]+5.196152422706631*(fr[14]+fr[13]+fr[11]+3.0*fr[3])+9.0*(fr[10]+fr[7]+fr[6]))/(36.0*EPSILON+1.732050807568877*(fr[12]+3.0*(fr[4]+fr[2]+fr[1]))+3.0*(fr[9]+fr[8]+fr[5]+3.0*fr[0])); 
  double fCtrlL[8], fCtrlR[8];  // fCtrl = anti-limited f evaluated at each control node on vx surface 
  // control node [x,y,vy] = [-1/3,-1/3,-1/3] 
  fCtrlL[0] = -0.02777777777777778*(1.732050807568877*(fl[12]+3.0*(fl[4]+fl[2]+fl[1]))-3.0*(fl[9]+fl[8]+fl[5])-9.0*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.02777777777777778*(1.732050807568877*(fr[12]+3.0*(fr[4]+fr[2]+fr[1]))-3.0*(fr[9]+fr[8]+fr[5])-9.0*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,y,vy] = [1/3,-1/3,-1/3] 
  fCtrlL[1] = 0.02777777777777778*(1.732050807568877*(fl[12]-3.0*(fl[4]+fl[2])+3.0*fl[1])+3.0*fl[9]-3.0*(fl[8]+fl[5])+9.0*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.02777777777777778*(1.732050807568877*(fr[12]-3.0*(fr[4]+fr[2])+3.0*fr[1])+3.0*fr[9]-3.0*(fr[8]+fr[5])+9.0*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,y,vy] = [-1/3,1/3,-1/3] 
  fCtrlL[2] = 0.02777777777777778*(1.732050807568877*(fl[12]-3.0*fl[4]+3.0*fl[2]-3.0*fl[1])-3.0*fl[9]+3.0*fl[8]-3.0*fl[5]+9.0*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = 0.02777777777777778*(1.732050807568877*(fr[12]-3.0*fr[4]+3.0*fr[2]-3.0*fr[1])-3.0*fr[9]+3.0*fr[8]-3.0*fr[5]+9.0*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,y,vy] = [1/3,1/3,-1/3] 
  fCtrlL[3] = -0.02777777777777778*(1.732050807568877*(fl[12]+3.0*fl[4]-3.0*(fl[2]+fl[1]))+3.0*(fl[9]+fl[8])-3.0*fl[5]-9.0*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = -0.02777777777777778*(1.732050807568877*(fr[12]+3.0*fr[4]-3.0*(fr[2]+fr[1]))+3.0*(fr[9]+fr[8])-3.0*fr[5]-9.0*fr[0])*limTheta(rCtrlR[3],-1.0); 
  // control node [x,y,vy] = [-1/3,-1/3,1/3] 
  fCtrlL[4] = 0.02777777777777778*(1.732050807568877*(fl[12]+3.0*fl[4]-3.0*(fl[2]+fl[1]))-3.0*(fl[9]+fl[8])+3.0*fl[5]+9.0*fl[0])*limTheta(rCtrlL[4],1.0); 
  fCtrlR[4] = 0.02777777777777778*(1.732050807568877*(fr[12]+3.0*fr[4]-3.0*(fr[2]+fr[1]))-3.0*(fr[9]+fr[8])+3.0*fr[5]+9.0*fr[0])*limTheta(rCtrlR[4],-1.0); 
  // control node [x,y,vy] = [1/3,-1/3,1/3] 
  fCtrlL[5] = -0.02777777777777778*(1.732050807568877*(fl[12]-3.0*fl[4]+3.0*fl[2]-3.0*fl[1])+3.0*fl[9]-3.0*fl[8]+3.0*fl[5]-9.0*fl[0])*limTheta(rCtrlL[5],1.0); 
  fCtrlR[5] = -0.02777777777777778*(1.732050807568877*(fr[12]-3.0*fr[4]+3.0*fr[2]-3.0*fr[1])+3.0*fr[9]-3.0*fr[8]+3.0*fr[5]-9.0*fr[0])*limTheta(rCtrlR[5],-1.0); 
  // control node [x,y,vy] = [-1/3,1/3,1/3] 
  fCtrlL[6] = -0.02777777777777778*(1.732050807568877*(fl[12]-3.0*(fl[4]+fl[2])+3.0*fl[1])-3.0*fl[9]+3.0*(fl[8]+fl[5])-9.0*fl[0])*limTheta(rCtrlL[6],1.0); 
  fCtrlR[6] = -0.02777777777777778*(1.732050807568877*(fr[12]-3.0*(fr[4]+fr[2])+3.0*fr[1])-3.0*fr[9]+3.0*(fr[8]+fr[5])-9.0*fr[0])*limTheta(rCtrlR[6],-1.0); 
  // control node [x,y,vy] = [1/3,1/3,1/3] 
  fCtrlL[7] = 0.02777777777777778*(1.732050807568877*(fl[12]+3.0*(fl[4]+fl[2]+fl[1]))+3.0*(fl[9]+fl[8]+fl[5])+9.0*fl[0])*limTheta(rCtrlL[7],1.0); 
  fCtrlR[7] = 0.02777777777777778*(1.732050807568877*(fr[12]+3.0*(fr[4]+fr[2]+fr[1]))+3.0*(fr[9]+fr[8]+fr[5])+9.0*fr[0])*limTheta(rCtrlR[7],-1.0); 
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
  fhatALQuad[0] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[0] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0])-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[2] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*(alpha[2]+alpha[1]); 
  if(alphaQuad > 0) {
  fhatALQuad[4] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[4] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0])-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[6] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[6] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[7] = 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
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
 
  incr[0] = 0.25*(alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[1] = 0.25*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(alpha[1]*fhatAL[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[3] = -0.4330127018922193*(alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[4] = 0.25*(alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+alpha[0]*fhatAL[3])*dfac_v; 
  incr[5] = 0.25*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 
  incr[6] = -0.4330127018922193*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[7] = -0.4330127018922193*(alpha[1]*fhatAL[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[8] = 0.25*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]+alpha[1]*fhatAL[3])*dfac_v; 
  incr[9] = 0.25*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])*dfac_v; 
  incr[10] = -0.4330127018922193*(alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+alpha[0]*fhatAL[3])*dfac_v; 
  incr[11] = -0.4330127018922193*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 
  incr[12] = 0.25*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5])*dfac_v; 
  incr[13] = -0.4330127018922193*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5]+alpha[1]*fhatAL[3])*dfac_v; 
  incr[14] = -0.4330127018922193*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])*dfac_v; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5])*dfac_v; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[16], outrPos[16]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.25 : positivityWeightByDirL[3]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.25 : positivityWeightByDirR[3]/positivityWeightByDirR[0]; 
  outlPos[0] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[1] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[2] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[3] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[4] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[5] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[6] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[7] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[8] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outlPos[9] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outlPos[10] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[11] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[12] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[13] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[14] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outlPos[15] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[0] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[1] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[2] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[3] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[4] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[5] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[6] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[7] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[8] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[9] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[10] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[11] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[12] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[13] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[14] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[15] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  if(outlPos[4] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]+fl[13]-1.0*fl[12]+fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*(fl[2]+fl[1]))-3.0*fl[10]+3.0*(fl[9]+fl[8])-3.0*(fl[7]+fl[6])+3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[4]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]+fr[11]+3.0*(fr[4]+fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[4] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[4] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[5] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]-1.0*fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2]+3.0*fl[1])-3.0*fl[10]+3.0*fl[9]-3.0*(fl[8]+fl[7])+3.0*fl[6]-3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[5]); 
  else limFac = 1.0; 
  if(outrPos[1] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12]+fr[11])+3.0*(fr[4]+fr[3]+fr[2])-3.0*fr[1])-3.0*(fr[10]+fr[9])+3.0*fr[8]-3.0*fr[7]+3.0*(fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[5] *= limFac; 
  outlPos[1] *= limFac; 
  outrPos[5] *= limFac; 
  outrPos[1] *= limFac; 
  if(outlPos[6] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12])+fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2])+3.0*fl[1])+3.0*(fl[10]+fl[9])-3.0*(fl[8]+fl[7])+3.0*(fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[6]); 
  else limFac = 1.0; 
  if(outrPos[2] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]+fr[11]-3.0*(fr[4]+fr[3])+3.0*fr[2]-3.0*fr[1])+3.0*fr[10]-3.0*fr[9]+3.0*fr[8]-3.0*fr[7]+3.0*fr[6]-3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[2]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[6] *= limFac; 
  outlPos[2] *= limFac; 
  outrPos[6] *= limFac; 
  outrPos[2] *= limFac; 
  if(outlPos[7] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]-1.0*fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8])-3.0*(fl[7]+fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[3] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]+fr[13]-1.0*(fr[12]+fr[11])-3.0*(fr[4]+fr[3])+3.0*(fr[2]+fr[1]))+3.0*fr[10]-3.0*(fr[9]+fr[8]+fr[7]+fr[6])+3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[3]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[3] *= limFac; 
  outrPos[7] *= limFac; 
  outrPos[3] *= limFac; 
  if(outlPos[12] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]+fl[13]-1.0*(fl[12]+fl[11])-3.0*(fl[4]+fl[3])+3.0*(fl[2]+fl[1]))-3.0*fl[10]+3.0*(fl[9]+fl[8]+fl[7]+fl[6])-3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[12]); 
  else limFac = 1.0; 
  if(outrPos[8] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]-1.0*fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8])+3.0*(fr[7]+fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[8]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[12] *= limFac; 
  outlPos[8] *= limFac; 
  outrPos[12] *= limFac; 
  outrPos[8] *= limFac; 
  if(outlPos[13] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]+fl[11]-3.0*(fl[4]+fl[3])+3.0*fl[2]-3.0*fl[1])-3.0*fl[10]+3.0*fl[9]-3.0*fl[8]+3.0*fl[7]-3.0*fl[6]+3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[13]); 
  else limFac = 1.0; 
  if(outrPos[9] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12])+fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2])+3.0*fr[1])-3.0*(fr[10]+fr[9])+3.0*(fr[8]+fr[7])-3.0*(fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[9]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[13] *= limFac; 
  outlPos[9] *= limFac; 
  outrPos[13] *= limFac; 
  outrPos[9] *= limFac; 
  if(outlPos[14] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12]+fl[11])+3.0*(fl[4]+fl[3]+fl[2])-3.0*fl[1])+3.0*(fl[10]+fl[9])-3.0*fl[8]+3.0*fl[7]-3.0*(fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[14]); 
  else limFac = 1.0; 
  if(outrPos[10] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]-1.0*fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2]+3.0*fr[1])+3.0*fr[10]-3.0*fr[9]+3.0*(fr[8]+fr[7])-3.0*fr[6]+3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[10]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[14] *= limFac; 
  outlPos[10] *= limFac; 
  outrPos[14] *= limFac; 
  outrPos[10] *= limFac; 
  if(outlPos[15] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]+fl[11]+3.0*(fl[4]+fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[15]); 
  else limFac = 1.0; 
  if(outrPos[11] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]+fr[13]-1.0*fr[12]+fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*(fr[2]+fr[1]))+3.0*fr[10]-3.0*(fr[9]+fr[8])+3.0*(fr[7]+fr[6])-3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[11]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[15] *= limFac; 
  outlPos[11] *= limFac; 
  outrPos[15] *= limFac; 
  outrPos[11] *= limFac; 
  outr[0] += 0.25*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.4330127018922193*(outrPos[15]-1.0*outrPos[14]+outrPos[13]-1.0*outrPos[12]+outrPos[11]-1.0*outrPos[10]+outrPos[9]-1.0*outrPos[8]+outrPos[7]-1.0*outrPos[6]+outrPos[5]-1.0*outrPos[4]+outrPos[3]-1.0*outrPos[2]+outrPos[1]-1.0*outrPos[0]); 
  outr[2] += 0.4330127018922193*(outrPos[15]+outrPos[14]-1.0*(outrPos[13]+outrPos[12])+outrPos[11]+outrPos[10]-1.0*(outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]-1.0*(outrPos[5]+outrPos[4])+outrPos[3]+outrPos[2]-1.0*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]-1.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]-1.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[4] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]-1.0*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[5] += 0.25*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*(outrPos[8]+outrPos[7])-3.0*(outrPos[6]+outrPos[5])+3.0*(outrPos[4]+outrPos[3])-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 
  outr[6] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*(outrPos[8]+outrPos[7])-3.0*outrPos[6]+3.0*outrPos[5]-3.0*(outrPos[4]+outrPos[3])+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[7] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])-3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[8] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*outrPos[12]+3.0*outrPos[11]-3.0*outrPos[10]+3.0*outrPos[9]-3.0*(outrPos[8]+outrPos[7])+3.0*outrPos[6]-3.0*outrPos[5]+3.0*outrPos[4]-3.0*outrPos[3]+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[9] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12])+3.0*(outrPos[11]+outrPos[10])-3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4])-3.0*(outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[10] += 0.25*(3.0*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12])-3.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])+3.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[11] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*outrPos[12]-3.0*outrPos[11]+3.0*(outrPos[10]+outrPos[9])-3.0*outrPos[8]+3.0*outrPos[7]-3.0*(outrPos[6]+outrPos[5])+3.0*outrPos[4]-3.0*outrPos[3]+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[12] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*outrPos[8]-3.0*outrPos[7]+3.0*(outrPos[6]+outrPos[5])-3.0*(outrPos[4]+outrPos[3])+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[13] += 0.4330127018922193*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*outrPos[8]-3.0*outrPos[7]+3.0*outrPos[6]-3.0*outrPos[5]+3.0*(outrPos[4]+outrPos[3])-3.0*outrPos[2]+3.0*outrPos[1]-3.0*outrPos[0]); 
  outr[14] += 0.4330127018922193*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8])-3.0*(outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])-3.0*(outrPos[1]+outrPos[0])); 
  outr[15] += 0.25*(9.0*outrPos[15]-9.0*(outrPos[14]+outrPos[13])+9.0*outrPos[12]-9.0*outrPos[11]+9.0*(outrPos[10]+outrPos[9])-9.0*(outrPos[8]+outrPos[7])+9.0*(outrPos[6]+outrPos[5])-9.0*outrPos[4]+9.0*outrPos[3]-9.0*(outrPos[2]+outrPos[1])+9.0*outrPos[0]); 

  outl[0] += 0.25*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.4330127018922193*(outlPos[15]-1.0*outlPos[14]+outlPos[13]-1.0*outlPos[12]+outlPos[11]-1.0*outlPos[10]+outlPos[9]-1.0*outlPos[8]+outlPos[7]-1.0*outlPos[6]+outlPos[5]-1.0*outlPos[4]+outlPos[3]-1.0*outlPos[2]+outlPos[1]-1.0*outlPos[0]); 
  outl[2] += 0.4330127018922193*(outlPos[15]+outlPos[14]-1.0*(outlPos[13]+outlPos[12])+outlPos[11]+outlPos[10]-1.0*(outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]-1.0*(outlPos[5]+outlPos[4])+outlPos[3]+outlPos[2]-1.0*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]-1.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]-1.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[4] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]-1.0*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[5] += 0.25*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*(outlPos[8]+outlPos[7])-3.0*(outlPos[6]+outlPos[5])+3.0*(outlPos[4]+outlPos[3])-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  outl[6] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*(outlPos[8]+outlPos[7])-3.0*outlPos[6]+3.0*outlPos[5]-3.0*(outlPos[4]+outlPos[3])+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[7] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])-3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[8] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*outlPos[12]+3.0*outlPos[11]-3.0*outlPos[10]+3.0*outlPos[9]-3.0*(outlPos[8]+outlPos[7])+3.0*outlPos[6]-3.0*outlPos[5]+3.0*outlPos[4]-3.0*outlPos[3]+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[9] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12])+3.0*(outlPos[11]+outlPos[10])-3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4])-3.0*(outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[10] += 0.25*(3.0*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12])-3.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])+3.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[11] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*outlPos[12]-3.0*outlPos[11]+3.0*(outlPos[10]+outlPos[9])-3.0*outlPos[8]+3.0*outlPos[7]-3.0*(outlPos[6]+outlPos[5])+3.0*outlPos[4]-3.0*outlPos[3]+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[12] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*outlPos[8]-3.0*outlPos[7]+3.0*(outlPos[6]+outlPos[5])-3.0*(outlPos[4]+outlPos[3])+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[13] += 0.4330127018922193*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*outlPos[8]-3.0*outlPos[7]+3.0*outlPos[6]-3.0*outlPos[5]+3.0*(outlPos[4]+outlPos[3])-3.0*outlPos[2]+3.0*outlPos[1]-3.0*outlPos[0]); 
  outl[14] += 0.4330127018922193*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8])-3.0*(outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])-3.0*(outlPos[1]+outlPos[0])); 
  outl[15] += 0.25*(9.0*outlPos[15]-9.0*(outlPos[14]+outlPos[13])+9.0*outlPos[12]-9.0*outlPos[11]+9.0*(outlPos[10]+outlPos[9])-9.0*(outlPos[8]+outlPos[7])+9.0*(outlPos[6]+outlPos[5])-9.0*outlPos[4]+9.0*outlPos[3]-9.0*(outlPos[2]+outlPos[1])+9.0*outlPos[0]); 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(0.125*(1.732050807568877*(2.0*BdriftX[1]*m_*wv2+(3.0*BmagInv[1]*Phi[3]+BmagInv[0]*Phi[2])*dfac_y*q_)-2.0*BdriftX[0]*m_*wv2-3.0*(BmagInv[0]*Phi[3]+BmagInv[1]*Phi[2])*dfac_y*q_))/q_; 

  double alpha[8]; 
  alpha[0] = -(0.5*(1.732050807568877*(2.828427124746191*BdriftX[1]*m_*wv2+(4.242640687119286*BmagInv[1]*Phi[3]+1.414213562373095*BmagInv[0]*Phi[2])*dfac_y*q_)-2.828427124746191*BdriftX[0]*m_*wv2-4.242640687119286*(BmagInv[0]*Phi[3]+BmagInv[1]*Phi[2])*dfac_y*q_))/q_; 
  alpha[2] = (0.3333333333333333*(2.449489742783178*BdriftX[0]-4.242640687119286*BdriftX[1])*m_*wv)/(dfac_v*q_); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[6]+fl[3])+alpha[0]*(1.732050807568877*fl[1]+fl[0]))*dfac_x; 
  incr[1] = -0.1767766952966368*(alpha[2]*(3.0*fl[6]+1.732050807568877*fl[3])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[2] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[11]+fl[7])+alpha[0]*(1.732050807568877*fl[5]+fl[2]))*dfac_x; 
  incr[3] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[6]+fl[3])+(1.732050807568877*fl[1]+fl[0])*alpha[2])*dfac_x; 
  incr[4] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[13]+fl[10])+alpha[0]*(1.732050807568877*fl[8]+fl[4]))*dfac_x; 
  incr[5] = -0.1767766952966368*(alpha[2]*(3.0*fl[11]+1.732050807568877*fl[7])+alpha[0]*(3.0*fl[5]+1.732050807568877*fl[2]))*dfac_x; 
  incr[6] = -0.1767766952966368*(alpha[0]*(3.0*fl[6]+1.732050807568877*fl[3])+(3.0*fl[1]+1.732050807568877*fl[0])*alpha[2])*dfac_x; 
  incr[7] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[11]+fl[7])+alpha[2]*(1.732050807568877*fl[5]+fl[2]))*dfac_x; 
  incr[8] = -0.1767766952966368*(alpha[2]*(3.0*fl[13]+1.732050807568877*fl[10])+alpha[0]*(3.0*fl[8]+1.732050807568877*fl[4]))*dfac_x; 
  incr[9] = 0.1767766952966368*(alpha[2]*(1.732050807568877*fl[15]+fl[14])+alpha[0]*(1.732050807568877*fl[12]+fl[9]))*dfac_x; 
  incr[10] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[13]+fl[10])+alpha[2]*(1.732050807568877*fl[8]+fl[4]))*dfac_x; 
  incr[11] = -0.1767766952966368*(alpha[0]*(3.0*fl[11]+1.732050807568877*fl[7])+alpha[2]*(3.0*fl[5]+1.732050807568877*fl[2]))*dfac_x; 
  incr[12] = -0.1767766952966368*(alpha[2]*(3.0*fl[15]+1.732050807568877*fl[14])+alpha[0]*(3.0*fl[12]+1.732050807568877*fl[9]))*dfac_x; 
  incr[13] = -0.1767766952966368*(alpha[0]*(3.0*fl[13]+1.732050807568877*fl[10])+alpha[2]*(3.0*fl[8]+1.732050807568877*fl[4]))*dfac_x; 
  incr[14] = 0.1767766952966368*(alpha[0]*(1.732050807568877*fl[15]+fl[14])+alpha[2]*(1.732050807568877*fl[12]+fl[9]))*dfac_x; 
  incr[15] = -0.1767766952966368*(alpha[0]*(3.0*fl[15]+1.732050807568877*fl[14])+alpha[2]*(3.0*fl[12]+1.732050807568877*fl[9]))*dfac_x; 
  } else { 
  incr[0] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[6]-1.0*fr[3])+alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0]))*dfac_x; 
  incr[1] = 0.1767766952966368*(alpha[2]*(3.0*fr[6]-1.732050807568877*fr[3])+alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0]))*dfac_x; 
  incr[2] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[11]-1.0*fr[7])+alpha[0]*(1.732050807568877*fr[5]-1.0*fr[2]))*dfac_x; 
  incr[3] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[6]-1.0*fr[3])+(1.732050807568877*fr[1]-1.0*fr[0])*alpha[2])*dfac_x; 
  incr[4] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[13]-1.0*fr[10])+alpha[0]*(1.732050807568877*fr[8]-1.0*fr[4]))*dfac_x; 
  incr[5] = 0.1767766952966368*(alpha[2]*(3.0*fr[11]-1.732050807568877*fr[7])+alpha[0]*(3.0*fr[5]-1.732050807568877*fr[2]))*dfac_x; 
  incr[6] = 0.1767766952966368*(alpha[0]*(3.0*fr[6]-1.732050807568877*fr[3])+(3.0*fr[1]-1.732050807568877*fr[0])*alpha[2])*dfac_x; 
  incr[7] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[11]-1.0*fr[7])+alpha[2]*(1.732050807568877*fr[5]-1.0*fr[2]))*dfac_x; 
  incr[8] = 0.1767766952966368*(alpha[2]*(3.0*fr[13]-1.732050807568877*fr[10])+alpha[0]*(3.0*fr[8]-1.732050807568877*fr[4]))*dfac_x; 
  incr[9] = -0.1767766952966368*(alpha[2]*(1.732050807568877*fr[15]-1.0*fr[14])+alpha[0]*(1.732050807568877*fr[12]-1.0*fr[9]))*dfac_x; 
  incr[10] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[13]-1.0*fr[10])+alpha[2]*(1.732050807568877*fr[8]-1.0*fr[4]))*dfac_x; 
  incr[11] = 0.1767766952966368*(alpha[0]*(3.0*fr[11]-1.732050807568877*fr[7])+alpha[2]*(3.0*fr[5]-1.732050807568877*fr[2]))*dfac_x; 
  incr[12] = 0.1767766952966368*(alpha[2]*(3.0*fr[15]-1.732050807568877*fr[14])+alpha[0]*(3.0*fr[12]-1.732050807568877*fr[9]))*dfac_x; 
  incr[13] = 0.1767766952966368*(alpha[0]*(3.0*fr[13]-1.732050807568877*fr[10])+alpha[2]*(3.0*fr[8]-1.732050807568877*fr[4]))*dfac_x; 
  incr[14] = -0.1767766952966368*(alpha[0]*(1.732050807568877*fr[15]-1.0*fr[14])+alpha[2]*(1.732050807568877*fr[12]-1.0*fr[9]))*dfac_x; 
  incr[15] = 0.1767766952966368*(alpha[0]*(3.0*fr[15]-1.732050807568877*fr[14])+alpha[2]*(3.0*fr[12]-1.732050807568877*fr[9]))*dfac_x; 
  }
#elif upwindType == QUAD 
  double rCtrlL[8], rCtrlR[8];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrlL[0] = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[13]+fl[12]+fl[11]+3.0*fl[1])+9.0*(fl[8]+fl[6]+fl[5])))/(36.0*EPSILON-1.732050807568877*(fl[14]+3.0*(fl[4]+fl[3]+fl[2]))+3.0*(fl[10]+fl[9]+fl[7]+3.0*fl[0])); 
  rCtrlL[1] = (3.0*fl[15]+5.196152422706631*(fl[13]-1.0*(fl[12]+fl[11])+3.0*fl[1])+9.0*(fl[5]-1.0*(fl[8]+fl[6])))/(36.0*EPSILON+1.732050807568877*(fl[14]+3.0*(fl[2]-1.0*(fl[4]+fl[3])))+3.0*(fl[10]-1.0*(fl[9]+fl[7]-3.0*fl[0]))); 
  rCtrlL[2] = (3.0*fl[15]-5.196152422706631*(fl[13]-1.0*fl[12]+fl[11]-3.0*fl[1])+9.0*((-1.0*fl[8])+fl[6]-1.0*fl[5]))/(36.0*EPSILON+1.732050807568877*(fl[14]+3.0*((-1.0*fl[4])+fl[3]-1.0*fl[2]))+3.0*((-1.0*fl[10])+fl[9]-1.0*fl[7]+3.0*fl[0])); 
  rCtrlL[3] = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[13]+fl[12]-1.0*(fl[11]+3.0*fl[1]))+9.0*(fl[8]-1.0*(fl[6]+fl[5]))))/(36.0*EPSILON-1.732050807568877*(fl[14]+3.0*(fl[4]-1.0*(fl[3]+fl[2])))+3.0*((-1.0*(fl[10]+fl[9]))+fl[7]+3.0*fl[0])); 
  rCtrlL[4] = (3.0*fl[15]-5.196152422706631*(fl[13]+fl[12]-1.0*(fl[11]+3.0*fl[1]))+9.0*(fl[8]-1.0*(fl[6]+fl[5])))/(36.0*EPSILON+1.732050807568877*(fl[14]+3.0*(fl[4]-1.0*(fl[3]+fl[2])))+3.0*((-1.0*(fl[10]+fl[9]))+fl[7]+3.0*fl[0])); 
  rCtrlL[5] = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[13]-1.0*fl[12]+fl[11]-3.0*fl[1])+9.0*((-1.0*fl[8])+fl[6]-1.0*fl[5])))/(36.0*EPSILON-1.732050807568877*(fl[14]+3.0*((-1.0*fl[4])+fl[3]-1.0*fl[2]))+3.0*((-1.0*fl[10])+fl[9]-1.0*fl[7]+3.0*fl[0])); 
  rCtrlL[6] = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[13]-1.0*(fl[12]+fl[11])+3.0*fl[1])+9.0*(fl[5]-1.0*(fl[8]+fl[6]))))/(36.0*EPSILON-1.732050807568877*(fl[14]+3.0*(fl[2]-1.0*(fl[4]+fl[3])))+3.0*(fl[10]-1.0*(fl[9]+fl[7]-3.0*fl[0]))); 
  rCtrlL[7] = (3.0*fl[15]+5.196152422706631*(fl[13]+fl[12]+fl[11]+3.0*fl[1])+9.0*(fl[8]+fl[6]+fl[5]))/(36.0*EPSILON+1.732050807568877*(fl[14]+3.0*(fl[4]+fl[3]+fl[2]))+3.0*(fl[10]+fl[9]+fl[7]+3.0*fl[0])); 
  rCtrlR[0] = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[13]+fr[12]+fr[11]+3.0*fr[1])+9.0*(fr[8]+fr[6]+fr[5])))/(36.0*EPSILON-1.732050807568877*(fr[14]+3.0*(fr[4]+fr[3]+fr[2]))+3.0*(fr[10]+fr[9]+fr[7]+3.0*fr[0])); 
  rCtrlR[1] = (3.0*fr[15]+5.196152422706631*(fr[13]-1.0*(fr[12]+fr[11])+3.0*fr[1])+9.0*(fr[5]-1.0*(fr[8]+fr[6])))/(36.0*EPSILON+1.732050807568877*(fr[14]+3.0*(fr[2]-1.0*(fr[4]+fr[3])))+3.0*(fr[10]-1.0*(fr[9]+fr[7]-3.0*fr[0]))); 
  rCtrlR[2] = (3.0*fr[15]-5.196152422706631*(fr[13]-1.0*fr[12]+fr[11]-3.0*fr[1])+9.0*((-1.0*fr[8])+fr[6]-1.0*fr[5]))/(36.0*EPSILON+1.732050807568877*(fr[14]+3.0*((-1.0*fr[4])+fr[3]-1.0*fr[2]))+3.0*((-1.0*fr[10])+fr[9]-1.0*fr[7]+3.0*fr[0])); 
  rCtrlR[3] = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[13]+fr[12]-1.0*(fr[11]+3.0*fr[1]))+9.0*(fr[8]-1.0*(fr[6]+fr[5]))))/(36.0*EPSILON-1.732050807568877*(fr[14]+3.0*(fr[4]-1.0*(fr[3]+fr[2])))+3.0*((-1.0*(fr[10]+fr[9]))+fr[7]+3.0*fr[0])); 
  rCtrlR[4] = (3.0*fr[15]-5.196152422706631*(fr[13]+fr[12]-1.0*(fr[11]+3.0*fr[1]))+9.0*(fr[8]-1.0*(fr[6]+fr[5])))/(36.0*EPSILON+1.732050807568877*(fr[14]+3.0*(fr[4]-1.0*(fr[3]+fr[2])))+3.0*((-1.0*(fr[10]+fr[9]))+fr[7]+3.0*fr[0])); 
  rCtrlR[5] = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[13]-1.0*fr[12]+fr[11]-3.0*fr[1])+9.0*((-1.0*fr[8])+fr[6]-1.0*fr[5])))/(36.0*EPSILON-1.732050807568877*(fr[14]+3.0*((-1.0*fr[4])+fr[3]-1.0*fr[2]))+3.0*((-1.0*fr[10])+fr[9]-1.0*fr[7]+3.0*fr[0])); 
  rCtrlR[6] = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[13]-1.0*(fr[12]+fr[11])+3.0*fr[1])+9.0*(fr[5]-1.0*(fr[8]+fr[6]))))/(36.0*EPSILON-1.732050807568877*(fr[14]+3.0*(fr[2]-1.0*(fr[4]+fr[3])))+3.0*(fr[10]-1.0*(fr[9]+fr[7]-3.0*fr[0]))); 
  rCtrlR[7] = (3.0*fr[15]+5.196152422706631*(fr[13]+fr[12]+fr[11]+3.0*fr[1])+9.0*(fr[8]+fr[6]+fr[5]))/(36.0*EPSILON+1.732050807568877*(fr[14]+3.0*(fr[4]+fr[3]+fr[2]))+3.0*(fr[10]+fr[9]+fr[7]+3.0*fr[0])); 
  double fCtrlL[8], fCtrlR[8];  // fCtrl = anti-limited f evaluated at each control node on x surface 
  // control node [y,vx,vy] = [-1/3,-1/3,-1/3] 
  fCtrlL[0] = -0.02777777777777778*(1.732050807568877*(fl[14]+3.0*(fl[4]+fl[3]+fl[2]))-3.0*(fl[10]+fl[9]+fl[7])-9.0*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.02777777777777778*(1.732050807568877*(fr[14]+3.0*(fr[4]+fr[3]+fr[2]))-3.0*(fr[10]+fr[9]+fr[7])-9.0*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [y,vx,vy] = [1/3,-1/3,-1/3] 
  fCtrlL[1] = 0.02777777777777778*(1.732050807568877*(fl[14]-3.0*(fl[4]+fl[3])+3.0*fl[2])+3.0*fl[10]-3.0*(fl[9]+fl[7])+9.0*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.02777777777777778*(1.732050807568877*(fr[14]-3.0*(fr[4]+fr[3])+3.0*fr[2])+3.0*fr[10]-3.0*(fr[9]+fr[7])+9.0*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [y,vx,vy] = [-1/3,1/3,-1/3] 
  fCtrlL[2] = 0.02777777777777778*(1.732050807568877*(fl[14]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2])-3.0*fl[10]+3.0*fl[9]-3.0*fl[7]+9.0*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = 0.02777777777777778*(1.732050807568877*(fr[14]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2])-3.0*fr[10]+3.0*fr[9]-3.0*fr[7]+9.0*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [y,vx,vy] = [1/3,1/3,-1/3] 
  fCtrlL[3] = -0.02777777777777778*(1.732050807568877*(fl[14]+3.0*fl[4]-3.0*(fl[3]+fl[2]))+3.0*(fl[10]+fl[9])-3.0*fl[7]-9.0*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = -0.02777777777777778*(1.732050807568877*(fr[14]+3.0*fr[4]-3.0*(fr[3]+fr[2]))+3.0*(fr[10]+fr[9])-3.0*fr[7]-9.0*fr[0])*limTheta(rCtrlR[3],-1.0); 
  // control node [y,vx,vy] = [-1/3,-1/3,1/3] 
  fCtrlL[4] = 0.02777777777777778*(1.732050807568877*(fl[14]+3.0*fl[4]-3.0*(fl[3]+fl[2]))-3.0*(fl[10]+fl[9])+3.0*fl[7]+9.0*fl[0])*limTheta(rCtrlL[4],1.0); 
  fCtrlR[4] = 0.02777777777777778*(1.732050807568877*(fr[14]+3.0*fr[4]-3.0*(fr[3]+fr[2]))-3.0*(fr[10]+fr[9])+3.0*fr[7]+9.0*fr[0])*limTheta(rCtrlR[4],-1.0); 
  // control node [y,vx,vy] = [1/3,-1/3,1/3] 
  fCtrlL[5] = -0.02777777777777778*(1.732050807568877*(fl[14]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2])+3.0*fl[10]-3.0*fl[9]+3.0*fl[7]-9.0*fl[0])*limTheta(rCtrlL[5],1.0); 
  fCtrlR[5] = -0.02777777777777778*(1.732050807568877*(fr[14]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2])+3.0*fr[10]-3.0*fr[9]+3.0*fr[7]-9.0*fr[0])*limTheta(rCtrlR[5],-1.0); 
  // control node [y,vx,vy] = [-1/3,1/3,1/3] 
  fCtrlL[6] = -0.02777777777777778*(1.732050807568877*(fl[14]-3.0*(fl[4]+fl[3])+3.0*fl[2])-3.0*fl[10]+3.0*(fl[9]+fl[7])-9.0*fl[0])*limTheta(rCtrlL[6],1.0); 
  fCtrlR[6] = -0.02777777777777778*(1.732050807568877*(fr[14]-3.0*(fr[4]+fr[3])+3.0*fr[2])-3.0*fr[10]+3.0*(fr[9]+fr[7])-9.0*fr[0])*limTheta(rCtrlR[6],-1.0); 
  // control node [y,vx,vy] = [1/3,1/3,1/3] 
  fCtrlL[7] = 0.02777777777777778*(1.732050807568877*(fl[14]+3.0*(fl[4]+fl[3]+fl[2]))+3.0*(fl[10]+fl[9]+fl[7])+9.0*fl[0])*limTheta(rCtrlL[7],1.0); 
  fCtrlR[7] = 0.02777777777777778*(1.732050807568877*(fr[14]+3.0*(fr[4]+fr[3]+fr[2]))+3.0*(fr[10]+fr[9]+fr[7])+9.0*fr[0])*limTheta(rCtrlR[7],-1.0); 
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
  fhatALQuad[0] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[0] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[2] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[2] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[4] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[4] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2]; 
  if(alphaQuad > 0) {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[6] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[6] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[2]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[7] = 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
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
 
  incr[0] = 0.25*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_x; 
  incr[1] = -0.4330127018922193*(alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_x; 
  incr[2] = 0.25*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])*dfac_x; 
  incr[3] = 0.25*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_x; 
  incr[4] = 0.25*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])*dfac_x; 
  incr[5] = -0.4330127018922193*(alpha[2]*fhatAL[4]+alpha[0]*fhatAL[1])*dfac_x; 
  incr[6] = -0.4330127018922193*(alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_x; 
  incr[7] = 0.25*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])*dfac_x; 
  incr[8] = -0.4330127018922193*(alpha[2]*fhatAL[6]+alpha[0]*fhatAL[3])*dfac_x; 
  incr[9] = 0.25*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])*dfac_x; 
  incr[10] = 0.25*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])*dfac_x; 
  incr[11] = -0.4330127018922193*(alpha[0]*fhatAL[4]+fhatAL[1]*alpha[2])*dfac_x; 
  incr[12] = -0.4330127018922193*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[5])*dfac_x; 
  incr[13] = -0.4330127018922193*(alpha[0]*fhatAL[6]+alpha[2]*fhatAL[3])*dfac_x; 
  incr[14] = 0.25*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])*dfac_x; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[5])*dfac_x; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[16], outrPos[16]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.25 : positivityWeightByDirL[1]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.25 : positivityWeightByDirR[1]/positivityWeightByDirR[0]; 
  outlPos[0] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[1] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[2] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[3] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[4] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[5] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[6] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[7] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[8] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[9] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[10] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outlPos[11] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outlPos[12] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outlPos[13] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outlPos[14] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[15] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[0] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[1] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[2] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[3] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[4] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[5] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[6] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[7] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[8] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[9] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[10] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[11] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[12] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[13] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[14] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[15] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  if(outlPos[1] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12]+fl[11])+3.0*(fl[4]+fl[3]+fl[2])-3.0*fl[1])-3.0*(fl[10]+fl[9])+3.0*fl[8]-3.0*fl[7]+3.0*(fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[1]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]+fr[11]+3.0*(fr[4]+fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[1] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[1] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[3] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]-1.0*(fl[12]+fl[11])-3.0*(fl[4]+fl[3])+3.0*(fl[2]+fl[1]))+3.0*fl[10]-3.0*(fl[9]+fl[8]+fl[7]+fl[6])+3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[2] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]+fr[11]-3.0*(fr[4]+fr[3])+3.0*fr[2]-3.0*fr[1])+3.0*fr[10]-3.0*fr[9]+3.0*fr[8]-3.0*fr[7]+3.0*fr[6]-3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[2]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[2] *= limFac; 
  outrPos[3] *= limFac; 
  outrPos[2] *= limFac; 
  if(outlPos[5] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]-1.0*fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2]+3.0*fl[1])-3.0*fl[10]+3.0*fl[9]-3.0*(fl[8]+fl[7])+3.0*fl[6]-3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[5]); 
  else limFac = 1.0; 
  if(outrPos[4] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]-1.0*fr[12]+fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*(fr[2]+fr[1]))-3.0*fr[10]+3.0*(fr[9]+fr[8])-3.0*(fr[7]+fr[6])+3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[4]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[5] *= limFac; 
  outlPos[4] *= limFac; 
  outrPos[5] *= limFac; 
  outrPos[4] *= limFac; 
  if(outlPos[7] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]-1.0*fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8])-3.0*(fl[7]+fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[6] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12])+fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2])+3.0*fr[1])+3.0*(fr[10]+fr[9])-3.0*(fr[8]+fr[7])+3.0*(fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[6]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[6] *= limFac; 
  outrPos[7] *= limFac; 
  outrPos[6] *= limFac; 
  if(outlPos[9] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12])+fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2])+3.0*fl[1])-3.0*(fl[10]+fl[9])+3.0*(fl[8]+fl[7])-3.0*(fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[9]); 
  else limFac = 1.0; 
  if(outrPos[8] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]-1.0*fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8])+3.0*(fr[7]+fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[8]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[9] *= limFac; 
  outlPos[8] *= limFac; 
  outrPos[9] *= limFac; 
  outrPos[8] *= limFac; 
  if(outlPos[11] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]-1.0*fl[12]+fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*(fl[2]+fl[1]))+3.0*fl[10]-3.0*(fl[9]+fl[8])+3.0*(fl[7]+fl[6])-3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[11]); 
  else limFac = 1.0; 
  if(outrPos[10] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]-1.0*fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2]+3.0*fr[1])+3.0*fr[10]-3.0*fr[9]+3.0*(fr[8]+fr[7])-3.0*fr[6]+3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[10]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[11] *= limFac; 
  outlPos[10] *= limFac; 
  outrPos[11] *= limFac; 
  outrPos[10] *= limFac; 
  if(outlPos[13] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]+fl[11]-3.0*(fl[4]+fl[3])+3.0*fl[2]-3.0*fl[1])-3.0*fl[10]+3.0*fl[9]-3.0*fl[8]+3.0*fl[7]-3.0*fl[6]+3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[13]); 
  else limFac = 1.0; 
  if(outrPos[12] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]-1.0*(fr[12]+fr[11])-3.0*(fr[4]+fr[3])+3.0*(fr[2]+fr[1]))-3.0*fr[10]+3.0*(fr[9]+fr[8]+fr[7]+fr[6])-3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[12]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[13] *= limFac; 
  outlPos[12] *= limFac; 
  outrPos[13] *= limFac; 
  outrPos[12] *= limFac; 
  if(outlPos[15] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]+fl[11]+3.0*(fl[4]+fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[15]); 
  else limFac = 1.0; 
  if(outrPos[14] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12]+fr[11])+3.0*(fr[4]+fr[3]+fr[2])-3.0*fr[1])+3.0*(fr[10]+fr[9])-3.0*fr[8]+3.0*fr[7]-3.0*(fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[14]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[15] *= limFac; 
  outlPos[14] *= limFac; 
  outrPos[15] *= limFac; 
  outrPos[14] *= limFac; 
  outr[0] += 0.25*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.4330127018922193*(outrPos[15]-1.0*outrPos[14]+outrPos[13]-1.0*outrPos[12]+outrPos[11]-1.0*outrPos[10]+outrPos[9]-1.0*outrPos[8]+outrPos[7]-1.0*outrPos[6]+outrPos[5]-1.0*outrPos[4]+outrPos[3]-1.0*outrPos[2]+outrPos[1]-1.0*outrPos[0]); 
  outr[2] += 0.4330127018922193*(outrPos[15]+outrPos[14]-1.0*(outrPos[13]+outrPos[12])+outrPos[11]+outrPos[10]-1.0*(outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]-1.0*(outrPos[5]+outrPos[4])+outrPos[3]+outrPos[2]-1.0*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]-1.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]-1.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[4] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]-1.0*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[5] += 0.25*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*(outrPos[8]+outrPos[7])-3.0*(outrPos[6]+outrPos[5])+3.0*(outrPos[4]+outrPos[3])-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 
  outr[6] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*(outrPos[8]+outrPos[7])-3.0*outrPos[6]+3.0*outrPos[5]-3.0*(outrPos[4]+outrPos[3])+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[7] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])-3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[8] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*outrPos[12]+3.0*outrPos[11]-3.0*outrPos[10]+3.0*outrPos[9]-3.0*(outrPos[8]+outrPos[7])+3.0*outrPos[6]-3.0*outrPos[5]+3.0*outrPos[4]-3.0*outrPos[3]+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[9] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12])+3.0*(outrPos[11]+outrPos[10])-3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4])-3.0*(outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[10] += 0.25*(3.0*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12])-3.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])+3.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[11] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*outrPos[12]-3.0*outrPos[11]+3.0*(outrPos[10]+outrPos[9])-3.0*outrPos[8]+3.0*outrPos[7]-3.0*(outrPos[6]+outrPos[5])+3.0*outrPos[4]-3.0*outrPos[3]+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[12] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*outrPos[8]-3.0*outrPos[7]+3.0*(outrPos[6]+outrPos[5])-3.0*(outrPos[4]+outrPos[3])+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[13] += 0.4330127018922193*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*outrPos[8]-3.0*outrPos[7]+3.0*outrPos[6]-3.0*outrPos[5]+3.0*(outrPos[4]+outrPos[3])-3.0*outrPos[2]+3.0*outrPos[1]-3.0*outrPos[0]); 
  outr[14] += 0.4330127018922193*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8])-3.0*(outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])-3.0*(outrPos[1]+outrPos[0])); 
  outr[15] += 0.25*(9.0*outrPos[15]-9.0*(outrPos[14]+outrPos[13])+9.0*outrPos[12]-9.0*outrPos[11]+9.0*(outrPos[10]+outrPos[9])-9.0*(outrPos[8]+outrPos[7])+9.0*(outrPos[6]+outrPos[5])-9.0*outrPos[4]+9.0*outrPos[3]-9.0*(outrPos[2]+outrPos[1])+9.0*outrPos[0]); 

  outl[0] += 0.25*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.4330127018922193*(outlPos[15]-1.0*outlPos[14]+outlPos[13]-1.0*outlPos[12]+outlPos[11]-1.0*outlPos[10]+outlPos[9]-1.0*outlPos[8]+outlPos[7]-1.0*outlPos[6]+outlPos[5]-1.0*outlPos[4]+outlPos[3]-1.0*outlPos[2]+outlPos[1]-1.0*outlPos[0]); 
  outl[2] += 0.4330127018922193*(outlPos[15]+outlPos[14]-1.0*(outlPos[13]+outlPos[12])+outlPos[11]+outlPos[10]-1.0*(outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]-1.0*(outlPos[5]+outlPos[4])+outlPos[3]+outlPos[2]-1.0*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]-1.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]-1.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[4] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]-1.0*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[5] += 0.25*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*(outlPos[8]+outlPos[7])-3.0*(outlPos[6]+outlPos[5])+3.0*(outlPos[4]+outlPos[3])-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  outl[6] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*(outlPos[8]+outlPos[7])-3.0*outlPos[6]+3.0*outlPos[5]-3.0*(outlPos[4]+outlPos[3])+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[7] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])-3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[8] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*outlPos[12]+3.0*outlPos[11]-3.0*outlPos[10]+3.0*outlPos[9]-3.0*(outlPos[8]+outlPos[7])+3.0*outlPos[6]-3.0*outlPos[5]+3.0*outlPos[4]-3.0*outlPos[3]+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[9] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12])+3.0*(outlPos[11]+outlPos[10])-3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4])-3.0*(outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[10] += 0.25*(3.0*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12])-3.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])+3.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[11] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*outlPos[12]-3.0*outlPos[11]+3.0*(outlPos[10]+outlPos[9])-3.0*outlPos[8]+3.0*outlPos[7]-3.0*(outlPos[6]+outlPos[5])+3.0*outlPos[4]-3.0*outlPos[3]+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[12] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*outlPos[8]-3.0*outlPos[7]+3.0*(outlPos[6]+outlPos[5])-3.0*(outlPos[4]+outlPos[3])+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[13] += 0.4330127018922193*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*outlPos[8]-3.0*outlPos[7]+3.0*outlPos[6]-3.0*outlPos[5]+3.0*(outlPos[4]+outlPos[3])-3.0*outlPos[2]+3.0*outlPos[1]-3.0*outlPos[0]); 
  outl[14] += 0.4330127018922193*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8])-3.0*(outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])-3.0*(outlPos[1]+outlPos[0])); 
  outl[15] += 0.25*(9.0*outlPos[15]-9.0*(outlPos[14]+outlPos[13])+9.0*outlPos[12]-9.0*outlPos[11]+9.0*(outlPos[10]+outlPos[9])-9.0*(outlPos[8]+outlPos[7])+9.0*(outlPos[6]+outlPos[5])-9.0*outlPos[4]+9.0*outlPos[3]-9.0*(outlPos[2]+outlPos[1])+9.0*outlPos[0]); 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Y_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = (0.125*(2.0*BdriftY[0]*m_*wv2+BmagInv[0]*dfac_x*(1.732050807568877*(Bmag[1]*wm+Phi[1]*q_)-3.0*Phi[3]*q_)))/q_; 

  double alpha[8]; 
  alpha[0] = (0.7071067811865475*(2.0*BdriftY[0]*m_*wv2+BmagInv[0]*dfac_x*(1.732050807568877*(Bmag[1]*wm+Phi[1]*q_)-3.0*Phi[3]*q_)))/q_; 
  alpha[1] = (0.7071067811865475*(2.0*BdriftY[1]*m_*wv2+BmagInv[1]*dfac_x*(1.732050807568877*(Bmag[1]*wm+Phi[1]*q_)-3.0*Phi[3]*q_)))/q_; 
  alpha[2] = (0.8164965809277261*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  alpha[3] = (0.7071067811865475*BmagInv[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alpha[4] = (0.8164965809277261*BdriftY[1]*m_*wv)/(dfac_v*q_); 
  alpha[5] = (0.7071067811865475*Bmag[1]*BmagInv[1]*dfac_x)/(dfac_m*q_); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[5]*fl[12]+alpha[4]*fl[11]+alpha[3]*fl[9])+alpha[5]*fl[8]+1.732050807568877*alpha[2]*fl[7]+alpha[4]*fl[6]+1.732050807568877*alpha[1]*fl[5]+alpha[3]*fl[4]+alpha[2]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_y; 
  incr[1] = 0.1767766952966368*(1.732050807568877*(alpha[3]*fl[12]+alpha[2]*fl[11]+alpha[5]*fl[9])+alpha[3]*fl[8]+1.732050807568877*alpha[4]*fl[7]+alpha[2]*fl[6]+1.732050807568877*alpha[0]*fl[5]+fl[4]*alpha[5]+fl[3]*alpha[4]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_y; 
  incr[2] = -0.1767766952966368*(3.0*(alpha[5]*fl[12]+alpha[4]*fl[11]+alpha[3]*fl[9])+1.732050807568877*alpha[5]*fl[8]+3.0*alpha[2]*fl[7]+1.732050807568877*alpha[4]*fl[6]+3.0*alpha[1]*fl[5]+1.732050807568877*(alpha[3]*fl[4]+alpha[2]*fl[3])+3.0*alpha[0]*fl[2]+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_y; 
  incr[3] = 0.1767766952966368*(1.732050807568877*(alpha[5]*fl[15]+alpha[3]*fl[14])+alpha[5]*fl[13]+1.732050807568877*alpha[1]*fl[11]+alpha[3]*fl[10]+1.732050807568877*alpha[0]*fl[7]+alpha[1]*fl[6]+alpha[4]*(1.732050807568877*fl[5]+fl[1])+alpha[0]*fl[3]+alpha[2]*(1.732050807568877*fl[2]+fl[0]))*dfac_y; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[15]+alpha[2]*fl[14])+alpha[4]*fl[13]+1.732050807568877*alpha[1]*fl[12]+alpha[2]*fl[10]+1.732050807568877*alpha[0]*fl[9]+alpha[1]*fl[8]+alpha[5]*(1.732050807568877*fl[5]+fl[1])+alpha[0]*fl[4]+(1.732050807568877*fl[2]+fl[0])*alpha[3])*dfac_y; 
  incr[5] = -0.1767766952966368*(3.0*(alpha[3]*fl[12]+alpha[2]*fl[11]+alpha[5]*fl[9])+1.732050807568877*alpha[3]*fl[8]+3.0*alpha[4]*fl[7]+1.732050807568877*alpha[2]*fl[6]+3.0*alpha[0]*fl[5]+1.732050807568877*(fl[4]*alpha[5]+fl[3]*alpha[4])+3.0*alpha[1]*fl[2]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_y; 
  incr[6] = 0.1767766952966368*(1.732050807568877*(alpha[3]*fl[15]+alpha[5]*fl[14])+alpha[3]*fl[13]+1.732050807568877*alpha[0]*fl[11]+alpha[5]*fl[10]+1.732050807568877*alpha[1]*fl[7]+alpha[0]*fl[6]+1.732050807568877*alpha[2]*fl[5]+(1.732050807568877*fl[2]+fl[0])*alpha[4]+alpha[1]*fl[3]+fl[1]*alpha[2])*dfac_y; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[5]*fl[15]+alpha[3]*fl[14])+1.732050807568877*alpha[5]*fl[13]+3.0*alpha[1]*fl[11]+1.732050807568877*alpha[3]*fl[10]+3.0*alpha[0]*fl[7]+1.732050807568877*alpha[1]*fl[6]+3.0*alpha[4]*fl[5]+1.732050807568877*(fl[1]*alpha[4]+alpha[0]*fl[3])+alpha[2]*(3.0*fl[2]+1.732050807568877*fl[0]))*dfac_y; 
  incr[8] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[15]+alpha[4]*fl[14])+alpha[2]*fl[13]+1.732050807568877*alpha[0]*fl[12]+alpha[4]*fl[10]+1.732050807568877*alpha[1]*fl[9]+alpha[0]*fl[8]+1.732050807568877*alpha[3]*fl[5]+(1.732050807568877*fl[2]+fl[0])*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3])*dfac_y; 
  incr[9] = -0.1767766952966368*(3.0*(alpha[4]*fl[15]+alpha[2]*fl[14])+1.732050807568877*alpha[4]*fl[13]+3.0*alpha[1]*fl[12]+1.732050807568877*alpha[2]*fl[10]+3.0*alpha[0]*fl[9]+1.732050807568877*alpha[1]*fl[8]+3.0*alpha[5]*fl[5]+1.732050807568877*(fl[1]*alpha[5]+alpha[0]*fl[4])+(3.0*fl[2]+1.732050807568877*fl[0])*alpha[3])*dfac_y; 
  incr[10] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14])+alpha[1]*fl[13]+1.732050807568877*(alpha[4]*fl[12]+alpha[5]*fl[11])+alpha[0]*fl[10]+1.732050807568877*alpha[2]*fl[9]+alpha[4]*fl[8]+1.732050807568877*alpha[3]*fl[7]+alpha[5]*fl[6]+alpha[2]*fl[4]+alpha[3]*fl[3])*dfac_y; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[3]*fl[15]+alpha[5]*fl[14])+1.732050807568877*alpha[3]*fl[13]+3.0*alpha[0]*fl[11]+1.732050807568877*alpha[5]*fl[10]+3.0*alpha[1]*fl[7]+1.732050807568877*alpha[0]*fl[6]+3.0*(alpha[2]*fl[5]+fl[2]*alpha[4])+1.732050807568877*(fl[0]*alpha[4]+alpha[1]*fl[3]+fl[1]*alpha[2]))*dfac_y; 
  incr[12] = -0.1767766952966368*(3.0*(alpha[2]*fl[15]+alpha[4]*fl[14])+1.732050807568877*alpha[2]*fl[13]+3.0*alpha[0]*fl[12]+1.732050807568877*alpha[4]*fl[10]+3.0*alpha[1]*fl[9]+1.732050807568877*alpha[0]*fl[8]+3.0*(alpha[3]*fl[5]+fl[2]*alpha[5])+1.732050807568877*(fl[0]*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3]))*dfac_y; 
  incr[13] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14])+alpha[0]*fl[13]+1.732050807568877*(alpha[2]*fl[12]+alpha[3]*fl[11])+alpha[1]*fl[10]+1.732050807568877*alpha[4]*fl[9]+alpha[2]*fl[8]+1.732050807568877*alpha[5]*fl[7]+alpha[3]*fl[6]+fl[3]*alpha[5]+alpha[4]*fl[4])*dfac_y; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14])+1.732050807568877*alpha[1]*fl[13]+3.0*(alpha[4]*fl[12]+alpha[5]*fl[11])+1.732050807568877*alpha[0]*fl[10]+3.0*alpha[2]*fl[9]+1.732050807568877*alpha[4]*fl[8]+3.0*alpha[3]*fl[7]+1.732050807568877*(alpha[5]*fl[6]+alpha[2]*fl[4]+alpha[3]*fl[3]))*dfac_y; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14])+1.732050807568877*alpha[0]*fl[13]+3.0*(alpha[2]*fl[12]+alpha[3]*fl[11])+1.732050807568877*alpha[1]*fl[10]+3.0*alpha[4]*fl[9]+1.732050807568877*alpha[2]*fl[8]+3.0*alpha[5]*fl[7]+1.732050807568877*(alpha[3]*fl[6]+fl[3]*alpha[5]+alpha[4]*fl[4]))*dfac_y; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fr[12]+alpha[4]*fr[11]+alpha[3]*fr[9])-1.0*alpha[5]*fr[8]+1.732050807568877*alpha[2]*fr[7]-1.0*alpha[4]*fr[6]+1.732050807568877*alpha[1]*fr[5]-1.0*(alpha[3]*fr[4]+alpha[2]*fr[3])+1.732050807568877*alpha[0]*fr[2]-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_y; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fr[12]+alpha[2]*fr[11]+alpha[5]*fr[9])-1.0*alpha[3]*fr[8]+1.732050807568877*alpha[4]*fr[7]-1.0*alpha[2]*fr[6]+1.732050807568877*alpha[0]*fr[5]-1.0*(fr[4]*alpha[5]+fr[3]*alpha[4])+1.732050807568877*alpha[1]*fr[2]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_y; 
  incr[2] = 0.1767766952966368*(3.0*(alpha[5]*fr[12]+alpha[4]*fr[11]+alpha[3]*fr[9])-1.732050807568877*alpha[5]*fr[8]+3.0*alpha[2]*fr[7]-1.732050807568877*alpha[4]*fr[6]+3.0*alpha[1]*fr[5]-1.732050807568877*(alpha[3]*fr[4]+alpha[2]*fr[3])+3.0*alpha[0]*fr[2]-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_y; 
  incr[3] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fr[15]+alpha[3]*fr[14])-1.0*alpha[5]*fr[13]+1.732050807568877*alpha[1]*fr[11]-1.0*alpha[3]*fr[10]+1.732050807568877*alpha[0]*fr[7]-1.0*alpha[1]*fr[6]+1.732050807568877*alpha[4]*fr[5]-1.0*(fr[1]*alpha[4]+alpha[0]*fr[3])+alpha[2]*(1.732050807568877*fr[2]-1.0*fr[0]))*dfac_y; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[15]+alpha[2]*fr[14])-1.0*alpha[4]*fr[13]+1.732050807568877*alpha[1]*fr[12]-1.0*alpha[2]*fr[10]+1.732050807568877*alpha[0]*fr[9]-1.0*alpha[1]*fr[8]+1.732050807568877*alpha[5]*fr[5]-1.0*(fr[1]*alpha[5]+alpha[0]*fr[4])+(1.732050807568877*fr[2]-1.0*fr[0])*alpha[3])*dfac_y; 
  incr[5] = 0.1767766952966368*(3.0*(alpha[3]*fr[12]+alpha[2]*fr[11]+alpha[5]*fr[9])-1.732050807568877*alpha[3]*fr[8]+3.0*alpha[4]*fr[7]-1.732050807568877*alpha[2]*fr[6]+3.0*alpha[0]*fr[5]-1.732050807568877*(fr[4]*alpha[5]+fr[3]*alpha[4])+3.0*alpha[1]*fr[2]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_y; 
  incr[6] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fr[15]+alpha[5]*fr[14])-1.0*alpha[3]*fr[13]+1.732050807568877*alpha[0]*fr[11]-1.0*alpha[5]*fr[10]+1.732050807568877*alpha[1]*fr[7]-1.0*alpha[0]*fr[6]+1.732050807568877*(alpha[2]*fr[5]+fr[2]*alpha[4])-1.0*(fr[0]*alpha[4]+alpha[1]*fr[3]+fr[1]*alpha[2]))*dfac_y; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[5]*fr[15]+alpha[3]*fr[14])-1.732050807568877*alpha[5]*fr[13]+3.0*alpha[1]*fr[11]-1.732050807568877*alpha[3]*fr[10]+3.0*alpha[0]*fr[7]-1.732050807568877*alpha[1]*fr[6]+3.0*alpha[4]*fr[5]-1.732050807568877*(fr[1]*alpha[4]+alpha[0]*fr[3])+alpha[2]*(3.0*fr[2]-1.732050807568877*fr[0]))*dfac_y; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[15]+alpha[4]*fr[14])-1.0*alpha[2]*fr[13]+1.732050807568877*alpha[0]*fr[12]-1.0*alpha[4]*fr[10]+1.732050807568877*alpha[1]*fr[9]-1.0*alpha[0]*fr[8]+1.732050807568877*(alpha[3]*fr[5]+fr[2]*alpha[5])-1.0*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_y; 
  incr[9] = 0.1767766952966368*(3.0*(alpha[4]*fr[15]+alpha[2]*fr[14])-1.732050807568877*alpha[4]*fr[13]+3.0*alpha[1]*fr[12]-1.732050807568877*alpha[2]*fr[10]+3.0*alpha[0]*fr[9]-1.732050807568877*alpha[1]*fr[8]+3.0*alpha[5]*fr[5]-1.732050807568877*(fr[1]*alpha[5]+alpha[0]*fr[4])+(3.0*fr[2]-1.732050807568877*fr[0])*alpha[3])*dfac_y; 
  incr[10] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14])-1.0*alpha[1]*fr[13]+1.732050807568877*(alpha[4]*fr[12]+alpha[5]*fr[11])-1.0*alpha[0]*fr[10]+1.732050807568877*alpha[2]*fr[9]-1.0*alpha[4]*fr[8]+1.732050807568877*alpha[3]*fr[7]-1.0*(alpha[5]*fr[6]+alpha[2]*fr[4]+alpha[3]*fr[3]))*dfac_y; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[3]*fr[15]+alpha[5]*fr[14])-1.732050807568877*alpha[3]*fr[13]+3.0*alpha[0]*fr[11]-1.732050807568877*alpha[5]*fr[10]+3.0*alpha[1]*fr[7]-1.732050807568877*alpha[0]*fr[6]+3.0*(alpha[2]*fr[5]+fr[2]*alpha[4])-1.732050807568877*(fr[0]*alpha[4]+alpha[1]*fr[3]+fr[1]*alpha[2]))*dfac_y; 
  incr[12] = 0.1767766952966368*(3.0*(alpha[2]*fr[15]+alpha[4]*fr[14])-1.732050807568877*alpha[2]*fr[13]+3.0*alpha[0]*fr[12]-1.732050807568877*alpha[4]*fr[10]+3.0*alpha[1]*fr[9]-1.732050807568877*alpha[0]*fr[8]+3.0*(alpha[3]*fr[5]+fr[2]*alpha[5])-1.732050807568877*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_y; 
  incr[13] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14])-1.0*alpha[0]*fr[13]+1.732050807568877*(alpha[2]*fr[12]+alpha[3]*fr[11])-1.0*alpha[1]*fr[10]+1.732050807568877*alpha[4]*fr[9]-1.0*alpha[2]*fr[8]+1.732050807568877*alpha[5]*fr[7]-1.0*(alpha[3]*fr[6]+fr[3]*alpha[5]+alpha[4]*fr[4]))*dfac_y; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14])-1.732050807568877*alpha[1]*fr[13]+3.0*(alpha[4]*fr[12]+alpha[5]*fr[11])-1.732050807568877*alpha[0]*fr[10]+3.0*alpha[2]*fr[9]-1.732050807568877*alpha[4]*fr[8]+3.0*alpha[3]*fr[7]-1.732050807568877*(alpha[5]*fr[6]+alpha[2]*fr[4]+alpha[3]*fr[3]))*dfac_y; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14])-1.732050807568877*alpha[0]*fr[13]+3.0*(alpha[2]*fr[12]+alpha[3]*fr[11])-1.732050807568877*alpha[1]*fr[10]+3.0*alpha[4]*fr[9]-1.732050807568877*alpha[2]*fr[8]+3.0*alpha[5]*fr[7]-1.732050807568877*(alpha[3]*fr[6]+fr[3]*alpha[5]+alpha[4]*fr[4]))*dfac_y; 
  }
#elif upwindType == QUAD 
  double rCtrlL[8], rCtrlR[8];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrlL[0] = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]+fl[12]+fl[11]+3.0*fl[2])+9.0*(fl[9]+fl[7]+fl[5])))/(36.0*EPSILON-1.732050807568877*(fl[13]+3.0*(fl[4]+fl[3]+fl[1]))+3.0*(fl[10]+fl[8]+fl[6]+3.0*fl[0])); 
  rCtrlL[1] = (3.0*fl[15]+5.196152422706631*(fl[14]-1.0*(fl[12]+fl[11])+3.0*fl[2])+9.0*(fl[5]-1.0*(fl[9]+fl[7])))/(36.0*EPSILON+1.732050807568877*(fl[13]+3.0*(fl[1]-1.0*(fl[4]+fl[3])))+3.0*(fl[10]-1.0*(fl[8]+fl[6]-3.0*fl[0]))); 
  rCtrlL[2] = (3.0*fl[15]-5.196152422706631*(fl[14]-1.0*fl[12]+fl[11]-3.0*fl[2])+9.0*((-1.0*fl[9])+fl[7]-1.0*fl[5]))/(36.0*EPSILON+1.732050807568877*(fl[13]+3.0*((-1.0*fl[4])+fl[3]-1.0*fl[1]))+3.0*((-1.0*fl[10])+fl[8]-1.0*fl[6]+3.0*fl[0])); 
  rCtrlL[3] = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]+fl[12]-1.0*(fl[11]+3.0*fl[2]))+9.0*(fl[9]-1.0*(fl[7]+fl[5]))))/(36.0*EPSILON-1.732050807568877*(fl[13]+3.0*(fl[4]-1.0*(fl[3]+fl[1])))+3.0*((-1.0*(fl[10]+fl[8]))+fl[6]+3.0*fl[0])); 
  rCtrlL[4] = (3.0*fl[15]-5.196152422706631*(fl[14]+fl[12]-1.0*(fl[11]+3.0*fl[2]))+9.0*(fl[9]-1.0*(fl[7]+fl[5])))/(36.0*EPSILON+1.732050807568877*(fl[13]+3.0*(fl[4]-1.0*(fl[3]+fl[1])))+3.0*((-1.0*(fl[10]+fl[8]))+fl[6]+3.0*fl[0])); 
  rCtrlL[5] = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]-1.0*fl[12]+fl[11]-3.0*fl[2])+9.0*((-1.0*fl[9])+fl[7]-1.0*fl[5])))/(36.0*EPSILON-1.732050807568877*(fl[13]+3.0*((-1.0*fl[4])+fl[3]-1.0*fl[1]))+3.0*((-1.0*fl[10])+fl[8]-1.0*fl[6]+3.0*fl[0])); 
  rCtrlL[6] = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]-1.0*(fl[12]+fl[11])+3.0*fl[2])+9.0*(fl[5]-1.0*(fl[9]+fl[7]))))/(36.0*EPSILON-1.732050807568877*(fl[13]+3.0*(fl[1]-1.0*(fl[4]+fl[3])))+3.0*(fl[10]-1.0*(fl[8]+fl[6]-3.0*fl[0]))); 
  rCtrlL[7] = (3.0*fl[15]+5.196152422706631*(fl[14]+fl[12]+fl[11]+3.0*fl[2])+9.0*(fl[9]+fl[7]+fl[5]))/(36.0*EPSILON+1.732050807568877*(fl[13]+3.0*(fl[4]+fl[3]+fl[1]))+3.0*(fl[10]+fl[8]+fl[6]+3.0*fl[0])); 
  rCtrlR[0] = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]+fr[12]+fr[11]+3.0*fr[2])+9.0*(fr[9]+fr[7]+fr[5])))/(36.0*EPSILON-1.732050807568877*(fr[13]+3.0*(fr[4]+fr[3]+fr[1]))+3.0*(fr[10]+fr[8]+fr[6]+3.0*fr[0])); 
  rCtrlR[1] = (3.0*fr[15]+5.196152422706631*(fr[14]-1.0*(fr[12]+fr[11])+3.0*fr[2])+9.0*(fr[5]-1.0*(fr[9]+fr[7])))/(36.0*EPSILON+1.732050807568877*(fr[13]+3.0*(fr[1]-1.0*(fr[4]+fr[3])))+3.0*(fr[10]-1.0*(fr[8]+fr[6]-3.0*fr[0]))); 
  rCtrlR[2] = (3.0*fr[15]-5.196152422706631*(fr[14]-1.0*fr[12]+fr[11]-3.0*fr[2])+9.0*((-1.0*fr[9])+fr[7]-1.0*fr[5]))/(36.0*EPSILON+1.732050807568877*(fr[13]+3.0*((-1.0*fr[4])+fr[3]-1.0*fr[1]))+3.0*((-1.0*fr[10])+fr[8]-1.0*fr[6]+3.0*fr[0])); 
  rCtrlR[3] = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]+fr[12]-1.0*(fr[11]+3.0*fr[2]))+9.0*(fr[9]-1.0*(fr[7]+fr[5]))))/(36.0*EPSILON-1.732050807568877*(fr[13]+3.0*(fr[4]-1.0*(fr[3]+fr[1])))+3.0*((-1.0*(fr[10]+fr[8]))+fr[6]+3.0*fr[0])); 
  rCtrlR[4] = (3.0*fr[15]-5.196152422706631*(fr[14]+fr[12]-1.0*(fr[11]+3.0*fr[2]))+9.0*(fr[9]-1.0*(fr[7]+fr[5])))/(36.0*EPSILON+1.732050807568877*(fr[13]+3.0*(fr[4]-1.0*(fr[3]+fr[1])))+3.0*((-1.0*(fr[10]+fr[8]))+fr[6]+3.0*fr[0])); 
  rCtrlR[5] = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]-1.0*fr[12]+fr[11]-3.0*fr[2])+9.0*((-1.0*fr[9])+fr[7]-1.0*fr[5])))/(36.0*EPSILON-1.732050807568877*(fr[13]+3.0*((-1.0*fr[4])+fr[3]-1.0*fr[1]))+3.0*((-1.0*fr[10])+fr[8]-1.0*fr[6]+3.0*fr[0])); 
  rCtrlR[6] = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]-1.0*(fr[12]+fr[11])+3.0*fr[2])+9.0*(fr[5]-1.0*(fr[9]+fr[7]))))/(36.0*EPSILON-1.732050807568877*(fr[13]+3.0*(fr[1]-1.0*(fr[4]+fr[3])))+3.0*(fr[10]-1.0*(fr[8]+fr[6]-3.0*fr[0]))); 
  rCtrlR[7] = (3.0*fr[15]+5.196152422706631*(fr[14]+fr[12]+fr[11]+3.0*fr[2])+9.0*(fr[9]+fr[7]+fr[5]))/(36.0*EPSILON+1.732050807568877*(fr[13]+3.0*(fr[4]+fr[3]+fr[1]))+3.0*(fr[10]+fr[8]+fr[6]+3.0*fr[0])); 
  double fCtrlL[8], fCtrlR[8];  // fCtrl = anti-limited f evaluated at each control node on y surface 
  // control node [x,vx,vy] = [-1/3,-1/3,-1/3] 
  fCtrlL[0] = -0.02777777777777778*(1.732050807568877*(fl[13]+3.0*(fl[4]+fl[3]+fl[1]))-3.0*(fl[10]+fl[8]+fl[6])-9.0*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.02777777777777778*(1.732050807568877*(fr[13]+3.0*(fr[4]+fr[3]+fr[1]))-3.0*(fr[10]+fr[8]+fr[6])-9.0*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,vx,vy] = [1/3,-1/3,-1/3] 
  fCtrlL[1] = 0.02777777777777778*(1.732050807568877*(fl[13]-3.0*(fl[4]+fl[3])+3.0*fl[1])+3.0*fl[10]-3.0*(fl[8]+fl[6])+9.0*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.02777777777777778*(1.732050807568877*(fr[13]-3.0*(fr[4]+fr[3])+3.0*fr[1])+3.0*fr[10]-3.0*(fr[8]+fr[6])+9.0*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,vx,vy] = [-1/3,1/3,-1/3] 
  fCtrlL[2] = 0.02777777777777778*(1.732050807568877*(fl[13]-3.0*fl[4]+3.0*fl[3]-3.0*fl[1])-3.0*fl[10]+3.0*fl[8]-3.0*fl[6]+9.0*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = 0.02777777777777778*(1.732050807568877*(fr[13]-3.0*fr[4]+3.0*fr[3]-3.0*fr[1])-3.0*fr[10]+3.0*fr[8]-3.0*fr[6]+9.0*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,vx,vy] = [1/3,1/3,-1/3] 
  fCtrlL[3] = -0.02777777777777778*(1.732050807568877*(fl[13]+3.0*fl[4]-3.0*(fl[3]+fl[1]))+3.0*(fl[10]+fl[8])-3.0*fl[6]-9.0*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = -0.02777777777777778*(1.732050807568877*(fr[13]+3.0*fr[4]-3.0*(fr[3]+fr[1]))+3.0*(fr[10]+fr[8])-3.0*fr[6]-9.0*fr[0])*limTheta(rCtrlR[3],-1.0); 
  // control node [x,vx,vy] = [-1/3,-1/3,1/3] 
  fCtrlL[4] = 0.02777777777777778*(1.732050807568877*(fl[13]+3.0*fl[4]-3.0*(fl[3]+fl[1]))-3.0*(fl[10]+fl[8])+3.0*fl[6]+9.0*fl[0])*limTheta(rCtrlL[4],1.0); 
  fCtrlR[4] = 0.02777777777777778*(1.732050807568877*(fr[13]+3.0*fr[4]-3.0*(fr[3]+fr[1]))-3.0*(fr[10]+fr[8])+3.0*fr[6]+9.0*fr[0])*limTheta(rCtrlR[4],-1.0); 
  // control node [x,vx,vy] = [1/3,-1/3,1/3] 
  fCtrlL[5] = -0.02777777777777778*(1.732050807568877*(fl[13]-3.0*fl[4]+3.0*fl[3]-3.0*fl[1])+3.0*fl[10]-3.0*fl[8]+3.0*fl[6]-9.0*fl[0])*limTheta(rCtrlL[5],1.0); 
  fCtrlR[5] = -0.02777777777777778*(1.732050807568877*(fr[13]-3.0*fr[4]+3.0*fr[3]-3.0*fr[1])+3.0*fr[10]-3.0*fr[8]+3.0*fr[6]-9.0*fr[0])*limTheta(rCtrlR[5],-1.0); 
  // control node [x,vx,vy] = [-1/3,1/3,1/3] 
  fCtrlL[6] = -0.02777777777777778*(1.732050807568877*(fl[13]-3.0*(fl[4]+fl[3])+3.0*fl[1])-3.0*fl[10]+3.0*(fl[8]+fl[6])-9.0*fl[0])*limTheta(rCtrlL[6],1.0); 
  fCtrlR[6] = -0.02777777777777778*(1.732050807568877*(fr[13]-3.0*(fr[4]+fr[3])+3.0*fr[1])-3.0*fr[10]+3.0*(fr[8]+fr[6])-9.0*fr[0])*limTheta(rCtrlR[6],-1.0); 
  // control node [x,vx,vy] = [1/3,1/3,1/3] 
  fCtrlL[7] = 0.02777777777777778*(1.732050807568877*(fl[13]+3.0*(fl[4]+fl[3]+fl[1]))+3.0*(fl[10]+fl[8]+fl[6])+9.0*fl[0])*limTheta(rCtrlL[7],1.0); 
  fCtrlR[7] = 0.02777777777777778*(1.732050807568877*(fr[13]+3.0*(fr[4]+fr[3]+fr[1]))+3.0*(fr[10]+fr[8]+fr[6])+9.0*fr[0])*limTheta(rCtrlR[7],-1.0); 
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
  fhatALQuad[0] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[0] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0])-0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2]); 
  if(alphaQuad > 0) {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[5]-0.3535533905932737*(alpha[4]+alpha[3])+0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[2] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = (-0.3535533905932737*alpha[5])+0.3535533905932737*alpha[4]-0.3535533905932737*alpha[3]+0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = (-0.3535533905932737*alpha[5])+0.3535533905932737*(alpha[4]+alpha[3])-0.3535533905932737*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[4] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[4] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*alpha[5]-0.3535533905932737*alpha[4]+0.3535533905932737*alpha[3]-0.3535533905932737*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = (-0.3535533905932737*(alpha[5]+alpha[4]))+0.3535533905932737*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[6] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[6] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[7] = 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
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
 
  incr[0] = 0.25*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_y; 
  incr[1] = 0.25*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_y; 
  incr[2] = -0.4330127018922193*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_y; 
  incr[3] = 0.25*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_y; 
  incr[4] = 0.25*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_y; 
  incr[5] = -0.4330127018922193*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_y; 
  incr[6] = 0.25*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_y; 
  incr[7] = -0.4330127018922193*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_y; 
  incr[8] = 0.25*(alpha[2]*fhatAL[7]+alpha[4]*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_y; 
  incr[9] = -0.4330127018922193*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_y; 
  incr[10] = 0.25*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_y; 
  incr[11] = -0.4330127018922193*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_y; 
  incr[12] = -0.4330127018922193*(alpha[2]*fhatAL[7]+alpha[4]*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_y; 
  incr[13] = 0.25*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])*dfac_y; 
  incr[14] = -0.4330127018922193*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_y; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])*dfac_y; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[16], outrPos[16]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.25 : positivityWeightByDirL[2]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.25 : positivityWeightByDirR[2]/positivityWeightByDirR[0]; 
  outlPos[0] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[1] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[2] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[3] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[4] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[5] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[6] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[7] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[8] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outlPos[9] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outlPos[10] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[11] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[12] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[13] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[14] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outlPos[15] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[0] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[1] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[2] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[3] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[4] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[5] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[6] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[7] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[8] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[9] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[10] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[11] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[12] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[13] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[14] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[15] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  if(outlPos[2] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]+fl[11]-3.0*(fl[4]+fl[3])+3.0*fl[2]-3.0*fl[1])+3.0*fl[10]-3.0*fl[9]+3.0*fl[8]-3.0*fl[7]+3.0*fl[6]-3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]+fr[11]+3.0*(fr[4]+fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[2] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[3] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]-1.0*(fl[12]+fl[11])-3.0*(fl[4]+fl[3])+3.0*(fl[2]+fl[1]))+3.0*fl[10]-3.0*(fl[9]+fl[8]+fl[7]+fl[6])+3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[1] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12]+fr[11])+3.0*(fr[4]+fr[3]+fr[2])-3.0*fr[1])-3.0*(fr[10]+fr[9])+3.0*fr[8]-3.0*fr[7]+3.0*(fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[1] *= limFac; 
  outrPos[3] *= limFac; 
  outrPos[1] *= limFac; 
  if(outlPos[6] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12])+fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2])+3.0*fl[1])+3.0*(fl[10]+fl[9])-3.0*(fl[8]+fl[7])+3.0*(fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[6]); 
  else limFac = 1.0; 
  if(outrPos[4] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]-1.0*fr[12]+fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*(fr[2]+fr[1]))-3.0*fr[10]+3.0*(fr[9]+fr[8])-3.0*(fr[7]+fr[6])+3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[4]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[6] *= limFac; 
  outlPos[4] *= limFac; 
  outrPos[6] *= limFac; 
  outrPos[4] *= limFac; 
  if(outlPos[7] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]-1.0*fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8])-3.0*(fl[7]+fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[5] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]-1.0*fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2]+3.0*fr[1])-3.0*fr[10]+3.0*fr[9]-3.0*(fr[8]+fr[7])+3.0*fr[6]-3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[5]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[5] *= limFac; 
  outrPos[7] *= limFac; 
  outrPos[5] *= limFac; 
  if(outlPos[10] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]-1.0*fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2]+3.0*fl[1])+3.0*fl[10]-3.0*fl[9]+3.0*(fl[8]+fl[7])-3.0*fl[6]+3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[10]); 
  else limFac = 1.0; 
  if(outrPos[8] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]-1.0*fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8])+3.0*(fr[7]+fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[8]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[10] *= limFac; 
  outlPos[8] *= limFac; 
  outrPos[10] *= limFac; 
  outrPos[8] *= limFac; 
  if(outlPos[11] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]-1.0*fl[12]+fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*(fl[2]+fl[1]))+3.0*fl[10]-3.0*(fl[9]+fl[8])+3.0*(fl[7]+fl[6])-3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[11]); 
  else limFac = 1.0; 
  if(outrPos[9] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12])+fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2])+3.0*fr[1])-3.0*(fr[10]+fr[9])+3.0*(fr[8]+fr[7])-3.0*(fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[9]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[11] *= limFac; 
  outlPos[9] *= limFac; 
  outrPos[11] *= limFac; 
  outrPos[9] *= limFac; 
  if(outlPos[14] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12]+fl[11])+3.0*(fl[4]+fl[3]+fl[2])-3.0*fl[1])+3.0*(fl[10]+fl[9])-3.0*fl[8]+3.0*fl[7]-3.0*(fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[14]); 
  else limFac = 1.0; 
  if(outrPos[12] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]-1.0*(fr[12]+fr[11])-3.0*(fr[4]+fr[3])+3.0*(fr[2]+fr[1]))-3.0*fr[10]+3.0*(fr[9]+fr[8]+fr[7]+fr[6])-3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[12]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[14] *= limFac; 
  outlPos[12] *= limFac; 
  outrPos[14] *= limFac; 
  outrPos[12] *= limFac; 
  if(outlPos[15] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]+fl[11]+3.0*(fl[4]+fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[15]); 
  else limFac = 1.0; 
  if(outrPos[13] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]+fr[11]-3.0*(fr[4]+fr[3])+3.0*fr[2]-3.0*fr[1])-3.0*fr[10]+3.0*fr[9]-3.0*fr[8]+3.0*fr[7]-3.0*fr[6]+3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[13]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[15] *= limFac; 
  outlPos[13] *= limFac; 
  outrPos[15] *= limFac; 
  outrPos[13] *= limFac; 
  outr[0] += 0.25*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.4330127018922193*(outrPos[15]-1.0*outrPos[14]+outrPos[13]-1.0*outrPos[12]+outrPos[11]-1.0*outrPos[10]+outrPos[9]-1.0*outrPos[8]+outrPos[7]-1.0*outrPos[6]+outrPos[5]-1.0*outrPos[4]+outrPos[3]-1.0*outrPos[2]+outrPos[1]-1.0*outrPos[0]); 
  outr[2] += 0.4330127018922193*(outrPos[15]+outrPos[14]-1.0*(outrPos[13]+outrPos[12])+outrPos[11]+outrPos[10]-1.0*(outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]-1.0*(outrPos[5]+outrPos[4])+outrPos[3]+outrPos[2]-1.0*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]-1.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]-1.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[4] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]-1.0*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[5] += 0.25*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*(outrPos[8]+outrPos[7])-3.0*(outrPos[6]+outrPos[5])+3.0*(outrPos[4]+outrPos[3])-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 
  outr[6] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*(outrPos[8]+outrPos[7])-3.0*outrPos[6]+3.0*outrPos[5]-3.0*(outrPos[4]+outrPos[3])+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[7] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])-3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[8] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*outrPos[12]+3.0*outrPos[11]-3.0*outrPos[10]+3.0*outrPos[9]-3.0*(outrPos[8]+outrPos[7])+3.0*outrPos[6]-3.0*outrPos[5]+3.0*outrPos[4]-3.0*outrPos[3]+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[9] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12])+3.0*(outrPos[11]+outrPos[10])-3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4])-3.0*(outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[10] += 0.25*(3.0*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12])-3.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])+3.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[11] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*outrPos[12]-3.0*outrPos[11]+3.0*(outrPos[10]+outrPos[9])-3.0*outrPos[8]+3.0*outrPos[7]-3.0*(outrPos[6]+outrPos[5])+3.0*outrPos[4]-3.0*outrPos[3]+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[12] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*outrPos[8]-3.0*outrPos[7]+3.0*(outrPos[6]+outrPos[5])-3.0*(outrPos[4]+outrPos[3])+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[13] += 0.4330127018922193*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*outrPos[8]-3.0*outrPos[7]+3.0*outrPos[6]-3.0*outrPos[5]+3.0*(outrPos[4]+outrPos[3])-3.0*outrPos[2]+3.0*outrPos[1]-3.0*outrPos[0]); 
  outr[14] += 0.4330127018922193*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8])-3.0*(outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])-3.0*(outrPos[1]+outrPos[0])); 
  outr[15] += 0.25*(9.0*outrPos[15]-9.0*(outrPos[14]+outrPos[13])+9.0*outrPos[12]-9.0*outrPos[11]+9.0*(outrPos[10]+outrPos[9])-9.0*(outrPos[8]+outrPos[7])+9.0*(outrPos[6]+outrPos[5])-9.0*outrPos[4]+9.0*outrPos[3]-9.0*(outrPos[2]+outrPos[1])+9.0*outrPos[0]); 

  outl[0] += 0.25*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.4330127018922193*(outlPos[15]-1.0*outlPos[14]+outlPos[13]-1.0*outlPos[12]+outlPos[11]-1.0*outlPos[10]+outlPos[9]-1.0*outlPos[8]+outlPos[7]-1.0*outlPos[6]+outlPos[5]-1.0*outlPos[4]+outlPos[3]-1.0*outlPos[2]+outlPos[1]-1.0*outlPos[0]); 
  outl[2] += 0.4330127018922193*(outlPos[15]+outlPos[14]-1.0*(outlPos[13]+outlPos[12])+outlPos[11]+outlPos[10]-1.0*(outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]-1.0*(outlPos[5]+outlPos[4])+outlPos[3]+outlPos[2]-1.0*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]-1.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]-1.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[4] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]-1.0*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[5] += 0.25*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*(outlPos[8]+outlPos[7])-3.0*(outlPos[6]+outlPos[5])+3.0*(outlPos[4]+outlPos[3])-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  outl[6] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*(outlPos[8]+outlPos[7])-3.0*outlPos[6]+3.0*outlPos[5]-3.0*(outlPos[4]+outlPos[3])+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[7] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])-3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[8] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*outlPos[12]+3.0*outlPos[11]-3.0*outlPos[10]+3.0*outlPos[9]-3.0*(outlPos[8]+outlPos[7])+3.0*outlPos[6]-3.0*outlPos[5]+3.0*outlPos[4]-3.0*outlPos[3]+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[9] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12])+3.0*(outlPos[11]+outlPos[10])-3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4])-3.0*(outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[10] += 0.25*(3.0*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12])-3.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])+3.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[11] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*outlPos[12]-3.0*outlPos[11]+3.0*(outlPos[10]+outlPos[9])-3.0*outlPos[8]+3.0*outlPos[7]-3.0*(outlPos[6]+outlPos[5])+3.0*outlPos[4]-3.0*outlPos[3]+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[12] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*outlPos[8]-3.0*outlPos[7]+3.0*(outlPos[6]+outlPos[5])-3.0*(outlPos[4]+outlPos[3])+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[13] += 0.4330127018922193*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*outlPos[8]-3.0*outlPos[7]+3.0*outlPos[6]-3.0*outlPos[5]+3.0*(outlPos[4]+outlPos[3])-3.0*outlPos[2]+3.0*outlPos[1]-3.0*outlPos[0]); 
  outl[14] += 0.4330127018922193*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8])-3.0*(outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])-3.0*(outlPos[1]+outlPos[0])); 
  outl[15] += 0.25*(9.0*outlPos[15]-9.0*(outlPos[14]+outlPos[13])+9.0*outlPos[12]-9.0*outlPos[11]+9.0*(outlPos[10]+outlPos[9])-9.0*(outlPos[8]+outlPos[7])+9.0*(outlPos[6]+outlPos[5])-9.0*outlPos[4]+9.0*outlPos[3]-9.0*(outlPos[2]+outlPos[1])+9.0*outlPos[0]); 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, 
                        const double *fl, const double *fr, double *outl, double *outr) 
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
  alpha[0] = -(1.224744871391589*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*(dfac_v*wv-1.0))/(dfac_v*q_); 
  alpha[1] = -(1.224744871391589*(BdriftX[1]*Bmag[1]*dfac_x*wm+((BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*q_)*(dfac_v*wv-1.0))/(dfac_v*q_); 
  alpha[2] = -(1.224744871391589*BdriftX[0]*Phi[3]*dfac_x*(dfac_v*wv-1.0))/dfac_v; 
  alpha[3] = -(0.7071067811865475*BdriftX[0]*Bmag[1]*dfac_x*(dfac_v*wv-1.0))/(dfac_m*dfac_v*q_); 
  alpha[4] = -(1.224744871391589*BdriftX[1]*Phi[3]*dfac_x*(dfac_v*wv-1.0))/dfac_v; 
  alpha[5] = -(0.7071067811865475*BdriftX[1]*Bmag[1]*dfac_x*(dfac_v*wv-1.0))/(dfac_m*dfac_v*q_); 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[5]*fl[13]+alpha[4]*fl[11]+alpha[3]*fl[10])+alpha[5]*fl[8]+1.732050807568877*(alpha[2]*fl[7]+alpha[1]*fl[6])+alpha[4]*fl[5]+alpha[3]*fl[4]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*(alpha[3]*fl[13]+alpha[2]*fl[11]+alpha[5]*fl[10])+alpha[3]*fl[8]+1.732050807568877*(alpha[4]*fl[7]+alpha[0]*fl[6])+alpha[2]*fl[5]+fl[4]*alpha[5]+fl[2]*alpha[4]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(1.732050807568877*(alpha[5]*fl[15]+alpha[3]*fl[14])+alpha[5]*fl[12]+1.732050807568877*alpha[1]*fl[11]+alpha[3]*fl[9]+1.732050807568877*(alpha[0]*fl[7]+alpha[4]*fl[6])+alpha[1]*fl[5]+fl[1]*alpha[4]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.1767766952966368*(3.0*(alpha[5]*fl[13]+alpha[4]*fl[11]+alpha[3]*fl[10])+1.732050807568877*alpha[5]*fl[8]+3.0*(alpha[2]*fl[7]+alpha[1]*fl[6])+1.732050807568877*(alpha[4]*fl[5]+alpha[3]*fl[4])+3.0*alpha[0]*fl[3]+1.732050807568877*(alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+alpha[4]*fl[12]+1.732050807568877*alpha[0]*fl[10]+alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[5]*(1.732050807568877*fl[6]+fl[1])+alpha[0]*fl[4]+alpha[3]*(1.732050807568877*fl[3]+fl[0]))*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*(alpha[3]*fl[15]+alpha[5]*fl[14])+alpha[3]*fl[12]+1.732050807568877*alpha[0]*fl[11]+alpha[5]*fl[9]+1.732050807568877*(alpha[1]*fl[7]+alpha[2]*fl[6])+alpha[0]*fl[5]+(1.732050807568877*fl[3]+fl[0])*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*(alpha[3]*fl[13]+alpha[2]*fl[11]+alpha[5]*fl[10])+1.732050807568877*alpha[3]*fl[8]+3.0*(alpha[4]*fl[7]+alpha[0]*fl[6])+1.732050807568877*(alpha[2]*fl[5]+fl[4]*alpha[5]+fl[2]*alpha[4])+3.0*alpha[1]*fl[3]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[5]*fl[15]+alpha[3]*fl[14])+1.732050807568877*alpha[5]*fl[12]+3.0*alpha[1]*fl[11]+1.732050807568877*alpha[3]*fl[9]+3.0*(alpha[0]*fl[7]+alpha[4]*fl[6])+1.732050807568877*(alpha[1]*fl[5]+fl[1]*alpha[4])+3.0*alpha[2]*fl[3]+1.732050807568877*(alpha[0]*fl[2]+fl[0]*alpha[2]))*dfac_v; 
  incr[8] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+alpha[2]*fl[12]+1.732050807568877*alpha[1]*fl[10]+alpha[4]*fl[9]+alpha[0]*fl[8]+1.732050807568877*alpha[3]*fl[6]+(1.732050807568877*fl[3]+fl[0])*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3])*dfac_v; 
  incr[9] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+alpha[1]*fl[12]+1.732050807568877*(alpha[5]*fl[11]+alpha[2]*fl[10])+alpha[0]*fl[9]+alpha[4]*fl[8]+1.732050807568877*alpha[3]*fl[7]+alpha[5]*fl[5]+alpha[2]*fl[4]+fl[2]*alpha[3])*dfac_v; 
  incr[10] = -0.1767766952966368*(3.0*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+1.732050807568877*alpha[4]*fl[12]+3.0*alpha[0]*fl[10]+1.732050807568877*(alpha[2]*fl[9]+alpha[1]*fl[8])+3.0*alpha[5]*fl[6]+1.732050807568877*(fl[1]*alpha[5]+alpha[0]*fl[4])+alpha[3]*(3.0*fl[3]+1.732050807568877*fl[0]))*dfac_v; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[3]*fl[15]+alpha[5]*fl[14])+1.732050807568877*alpha[3]*fl[12]+3.0*alpha[0]*fl[11]+1.732050807568877*alpha[5]*fl[9]+3.0*(alpha[1]*fl[7]+alpha[2]*fl[6])+1.732050807568877*alpha[0]*fl[5]+3.0*fl[3]*alpha[4]+1.732050807568877*(fl[0]*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2]))*dfac_v; 
  incr[12] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+alpha[0]*fl[12]+1.732050807568877*(alpha[3]*fl[11]+alpha[4]*fl[10])+alpha[1]*fl[9]+alpha[2]*fl[8]+1.732050807568877*alpha[5]*fl[7]+alpha[3]*fl[5]+fl[2]*alpha[5]+alpha[4]*fl[4])*dfac_v; 
  incr[13] = -0.1767766952966368*(3.0*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+1.732050807568877*alpha[2]*fl[12]+3.0*alpha[1]*fl[10]+1.732050807568877*(alpha[4]*fl[9]+alpha[0]*fl[8])+3.0*(alpha[3]*fl[6]+fl[3]*alpha[5])+1.732050807568877*(fl[0]*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3]))*dfac_v; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+1.732050807568877*alpha[1]*fl[12]+3.0*(alpha[5]*fl[11]+alpha[2]*fl[10])+1.732050807568877*(alpha[0]*fl[9]+alpha[4]*fl[8])+3.0*alpha[3]*fl[7]+1.732050807568877*(alpha[5]*fl[5]+alpha[2]*fl[4]+fl[2]*alpha[3]))*dfac_v; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+1.732050807568877*alpha[0]*fl[12]+3.0*(alpha[3]*fl[11]+alpha[4]*fl[10])+1.732050807568877*(alpha[1]*fl[9]+alpha[2]*fl[8])+3.0*alpha[5]*fl[7]+1.732050807568877*(alpha[3]*fl[5]+fl[2]*alpha[5]+alpha[4]*fl[4]))*dfac_v; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fr[13]+alpha[4]*fr[11]+alpha[3]*fr[10])-1.0*alpha[5]*fr[8]+1.732050807568877*(alpha[2]*fr[7]+alpha[1]*fr[6])-1.0*(alpha[4]*fr[5]+alpha[3]*fr[4])+1.732050807568877*alpha[0]*fr[3]-1.0*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fr[13]+alpha[2]*fr[11]+alpha[5]*fr[10])-1.0*alpha[3]*fr[8]+1.732050807568877*(alpha[4]*fr[7]+alpha[0]*fr[6])-1.0*(alpha[2]*fr[5]+fr[4]*alpha[5]+fr[2]*alpha[4])+1.732050807568877*alpha[1]*fr[3]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fr[15]+alpha[3]*fr[14])-1.0*alpha[5]*fr[12]+1.732050807568877*alpha[1]*fr[11]-1.0*alpha[3]*fr[9]+1.732050807568877*(alpha[0]*fr[7]+alpha[4]*fr[6])-1.0*(alpha[1]*fr[5]+fr[1]*alpha[4])+1.732050807568877*alpha[2]*fr[3]-1.0*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[5]*fr[13]+alpha[4]*fr[11]+alpha[3]*fr[10])-1.732050807568877*alpha[5]*fr[8]+3.0*(alpha[2]*fr[7]+alpha[1]*fr[6])-1.732050807568877*(alpha[4]*fr[5]+alpha[3]*fr[4])+3.0*alpha[0]*fr[3]-1.732050807568877*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.0*alpha[4]*fr[12]+1.732050807568877*alpha[0]*fr[10]-1.0*(alpha[2]*fr[9]+alpha[1]*fr[8])+1.732050807568877*alpha[5]*fr[6]-1.0*(fr[1]*alpha[5]+alpha[0]*fr[4])+alpha[3]*(1.732050807568877*fr[3]-1.0*fr[0]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fr[15]+alpha[5]*fr[14])-1.0*alpha[3]*fr[12]+1.732050807568877*alpha[0]*fr[11]-1.0*alpha[5]*fr[9]+1.732050807568877*(alpha[1]*fr[7]+alpha[2]*fr[6])-1.0*alpha[0]*fr[5]+1.732050807568877*fr[3]*alpha[4]-1.0*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[3]*fr[13]+alpha[2]*fr[11]+alpha[5]*fr[10])-1.732050807568877*alpha[3]*fr[8]+3.0*(alpha[4]*fr[7]+alpha[0]*fr[6])-1.732050807568877*(alpha[2]*fr[5]+fr[4]*alpha[5]+fr[2]*alpha[4])+3.0*alpha[1]*fr[3]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[5]*fr[15]+alpha[3]*fr[14])-1.732050807568877*alpha[5]*fr[12]+3.0*alpha[1]*fr[11]-1.732050807568877*alpha[3]*fr[9]+3.0*(alpha[0]*fr[7]+alpha[4]*fr[6])-1.732050807568877*(alpha[1]*fr[5]+fr[1]*alpha[4])+3.0*alpha[2]*fr[3]-1.732050807568877*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.0*alpha[2]*fr[12]+1.732050807568877*alpha[1]*fr[10]-1.0*(alpha[4]*fr[9]+alpha[0]*fr[8])+1.732050807568877*(alpha[3]*fr[6]+fr[3]*alpha[5])-1.0*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.0*alpha[1]*fr[12]+1.732050807568877*(alpha[5]*fr[11]+alpha[2]*fr[10])-1.0*(alpha[0]*fr[9]+alpha[4]*fr[8])+1.732050807568877*alpha[3]*fr[7]-1.0*(alpha[5]*fr[5]+alpha[2]*fr[4]+fr[2]*alpha[3]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.732050807568877*alpha[4]*fr[12]+3.0*alpha[0]*fr[10]-1.732050807568877*(alpha[2]*fr[9]+alpha[1]*fr[8])+3.0*alpha[5]*fr[6]-1.732050807568877*(fr[1]*alpha[5]+alpha[0]*fr[4])+alpha[3]*(3.0*fr[3]-1.732050807568877*fr[0]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[3]*fr[15]+alpha[5]*fr[14])-1.732050807568877*alpha[3]*fr[12]+3.0*alpha[0]*fr[11]-1.732050807568877*alpha[5]*fr[9]+3.0*(alpha[1]*fr[7]+alpha[2]*fr[6])-1.732050807568877*alpha[0]*fr[5]+3.0*fr[3]*alpha[4]-1.732050807568877*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.0*alpha[0]*fr[12]+1.732050807568877*(alpha[3]*fr[11]+alpha[4]*fr[10])-1.0*(alpha[1]*fr[9]+alpha[2]*fr[8])+1.732050807568877*alpha[5]*fr[7]-1.0*(alpha[3]*fr[5]+fr[2]*alpha[5]+alpha[4]*fr[4]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.732050807568877*alpha[2]*fr[12]+3.0*alpha[1]*fr[10]-1.732050807568877*(alpha[4]*fr[9]+alpha[0]*fr[8])+3.0*(alpha[3]*fr[6]+fr[3]*alpha[5])-1.732050807568877*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.732050807568877*alpha[1]*fr[12]+3.0*(alpha[5]*fr[11]+alpha[2]*fr[10])-1.732050807568877*(alpha[0]*fr[9]+alpha[4]*fr[8])+3.0*alpha[3]*fr[7]-1.732050807568877*(alpha[5]*fr[5]+alpha[2]*fr[4]+fr[2]*alpha[3]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.732050807568877*alpha[0]*fr[12]+3.0*(alpha[3]*fr[11]+alpha[4]*fr[10])-1.732050807568877*(alpha[1]*fr[9]+alpha[2]*fr[8])+3.0*alpha[5]*fr[7]-1.732050807568877*(alpha[3]*fr[5]+fr[2]*alpha[5]+alpha[4]*fr[4]))*dfac_v; 
  }
#elif upwindType == QUAD 
  double rCtrlL[8], rCtrlR[8];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrlL[0] = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]+fl[13]+fl[11]+3.0*fl[3])+9.0*(fl[10]+fl[7]+fl[6])))/(36.0*EPSILON-1.732050807568877*(fl[12]+3.0*(fl[4]+fl[2]+fl[1]))+3.0*(fl[9]+fl[8]+fl[5]+3.0*fl[0])); 
  rCtrlL[1] = (3.0*fl[15]+5.196152422706631*(fl[14]-1.0*(fl[13]+fl[11])+3.0*fl[3])+9.0*(fl[6]-1.0*(fl[10]+fl[7])))/(36.0*EPSILON+1.732050807568877*(fl[12]+3.0*(fl[1]-1.0*(fl[4]+fl[2])))+3.0*(fl[9]-1.0*(fl[8]+fl[5]-3.0*fl[0]))); 
  rCtrlL[2] = (3.0*fl[15]-5.196152422706631*(fl[14]-1.0*fl[13]+fl[11]-3.0*fl[3])+9.0*((-1.0*fl[10])+fl[7]-1.0*fl[6]))/(36.0*EPSILON+1.732050807568877*(fl[12]+3.0*((-1.0*fl[4])+fl[2]-1.0*fl[1]))+3.0*((-1.0*fl[9])+fl[8]-1.0*fl[5]+3.0*fl[0])); 
  rCtrlL[3] = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]+fl[13]-1.0*(fl[11]+3.0*fl[3]))+9.0*(fl[10]-1.0*(fl[7]+fl[6]))))/(36.0*EPSILON-1.732050807568877*(fl[12]+3.0*(fl[4]-1.0*(fl[2]+fl[1])))+3.0*((-1.0*(fl[9]+fl[8]))+fl[5]+3.0*fl[0])); 
  rCtrlL[4] = (3.0*fl[15]-5.196152422706631*(fl[14]+fl[13]-1.0*(fl[11]+3.0*fl[3]))+9.0*(fl[10]-1.0*(fl[7]+fl[6])))/(36.0*EPSILON+1.732050807568877*(fl[12]+3.0*(fl[4]-1.0*(fl[2]+fl[1])))+3.0*((-1.0*(fl[9]+fl[8]))+fl[5]+3.0*fl[0])); 
  rCtrlL[5] = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]-1.0*fl[13]+fl[11]-3.0*fl[3])+9.0*((-1.0*fl[10])+fl[7]-1.0*fl[6])))/(36.0*EPSILON-1.732050807568877*(fl[12]+3.0*((-1.0*fl[4])+fl[2]-1.0*fl[1]))+3.0*((-1.0*fl[9])+fl[8]-1.0*fl[5]+3.0*fl[0])); 
  rCtrlL[6] = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]-1.0*(fl[13]+fl[11])+3.0*fl[3])+9.0*(fl[6]-1.0*(fl[10]+fl[7]))))/(36.0*EPSILON-1.732050807568877*(fl[12]+3.0*(fl[1]-1.0*(fl[4]+fl[2])))+3.0*(fl[9]-1.0*(fl[8]+fl[5]-3.0*fl[0]))); 
  rCtrlL[7] = (3.0*fl[15]+5.196152422706631*(fl[14]+fl[13]+fl[11]+3.0*fl[3])+9.0*(fl[10]+fl[7]+fl[6]))/(36.0*EPSILON+1.732050807568877*(fl[12]+3.0*(fl[4]+fl[2]+fl[1]))+3.0*(fl[9]+fl[8]+fl[5]+3.0*fl[0])); 
  rCtrlR[0] = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]+fr[13]+fr[11]+3.0*fr[3])+9.0*(fr[10]+fr[7]+fr[6])))/(36.0*EPSILON-1.732050807568877*(fr[12]+3.0*(fr[4]+fr[2]+fr[1]))+3.0*(fr[9]+fr[8]+fr[5]+3.0*fr[0])); 
  rCtrlR[1] = (3.0*fr[15]+5.196152422706631*(fr[14]-1.0*(fr[13]+fr[11])+3.0*fr[3])+9.0*(fr[6]-1.0*(fr[10]+fr[7])))/(36.0*EPSILON+1.732050807568877*(fr[12]+3.0*(fr[1]-1.0*(fr[4]+fr[2])))+3.0*(fr[9]-1.0*(fr[8]+fr[5]-3.0*fr[0]))); 
  rCtrlR[2] = (3.0*fr[15]-5.196152422706631*(fr[14]-1.0*fr[13]+fr[11]-3.0*fr[3])+9.0*((-1.0*fr[10])+fr[7]-1.0*fr[6]))/(36.0*EPSILON+1.732050807568877*(fr[12]+3.0*((-1.0*fr[4])+fr[2]-1.0*fr[1]))+3.0*((-1.0*fr[9])+fr[8]-1.0*fr[5]+3.0*fr[0])); 
  rCtrlR[3] = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]+fr[13]-1.0*(fr[11]+3.0*fr[3]))+9.0*(fr[10]-1.0*(fr[7]+fr[6]))))/(36.0*EPSILON-1.732050807568877*(fr[12]+3.0*(fr[4]-1.0*(fr[2]+fr[1])))+3.0*((-1.0*(fr[9]+fr[8]))+fr[5]+3.0*fr[0])); 
  rCtrlR[4] = (3.0*fr[15]-5.196152422706631*(fr[14]+fr[13]-1.0*(fr[11]+3.0*fr[3]))+9.0*(fr[10]-1.0*(fr[7]+fr[6])))/(36.0*EPSILON+1.732050807568877*(fr[12]+3.0*(fr[4]-1.0*(fr[2]+fr[1])))+3.0*((-1.0*(fr[9]+fr[8]))+fr[5]+3.0*fr[0])); 
  rCtrlR[5] = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]-1.0*fr[13]+fr[11]-3.0*fr[3])+9.0*((-1.0*fr[10])+fr[7]-1.0*fr[6])))/(36.0*EPSILON-1.732050807568877*(fr[12]+3.0*((-1.0*fr[4])+fr[2]-1.0*fr[1]))+3.0*((-1.0*fr[9])+fr[8]-1.0*fr[5]+3.0*fr[0])); 
  rCtrlR[6] = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]-1.0*(fr[13]+fr[11])+3.0*fr[3])+9.0*(fr[6]-1.0*(fr[10]+fr[7]))))/(36.0*EPSILON-1.732050807568877*(fr[12]+3.0*(fr[1]-1.0*(fr[4]+fr[2])))+3.0*(fr[9]-1.0*(fr[8]+fr[5]-3.0*fr[0]))); 
  rCtrlR[7] = (3.0*fr[15]+5.196152422706631*(fr[14]+fr[13]+fr[11]+3.0*fr[3])+9.0*(fr[10]+fr[7]+fr[6]))/(36.0*EPSILON+1.732050807568877*(fr[12]+3.0*(fr[4]+fr[2]+fr[1]))+3.0*(fr[9]+fr[8]+fr[5]+3.0*fr[0])); 
  double fCtrlL[8], fCtrlR[8];  // fCtrl = anti-limited f evaluated at each control node on vx surface 
  // control node [x,y,vy] = [-1/3,-1/3,-1/3] 
  fCtrlL[0] = -0.02777777777777778*(1.732050807568877*(fl[12]+3.0*(fl[4]+fl[2]+fl[1]))-3.0*(fl[9]+fl[8]+fl[5])-9.0*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.02777777777777778*(1.732050807568877*(fr[12]+3.0*(fr[4]+fr[2]+fr[1]))-3.0*(fr[9]+fr[8]+fr[5])-9.0*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,y,vy] = [1/3,-1/3,-1/3] 
  fCtrlL[1] = 0.02777777777777778*(1.732050807568877*(fl[12]-3.0*(fl[4]+fl[2])+3.0*fl[1])+3.0*fl[9]-3.0*(fl[8]+fl[5])+9.0*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.02777777777777778*(1.732050807568877*(fr[12]-3.0*(fr[4]+fr[2])+3.0*fr[1])+3.0*fr[9]-3.0*(fr[8]+fr[5])+9.0*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,y,vy] = [-1/3,1/3,-1/3] 
  fCtrlL[2] = 0.02777777777777778*(1.732050807568877*(fl[12]-3.0*fl[4]+3.0*fl[2]-3.0*fl[1])-3.0*fl[9]+3.0*fl[8]-3.0*fl[5]+9.0*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = 0.02777777777777778*(1.732050807568877*(fr[12]-3.0*fr[4]+3.0*fr[2]-3.0*fr[1])-3.0*fr[9]+3.0*fr[8]-3.0*fr[5]+9.0*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,y,vy] = [1/3,1/3,-1/3] 
  fCtrlL[3] = -0.02777777777777778*(1.732050807568877*(fl[12]+3.0*fl[4]-3.0*(fl[2]+fl[1]))+3.0*(fl[9]+fl[8])-3.0*fl[5]-9.0*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = -0.02777777777777778*(1.732050807568877*(fr[12]+3.0*fr[4]-3.0*(fr[2]+fr[1]))+3.0*(fr[9]+fr[8])-3.0*fr[5]-9.0*fr[0])*limTheta(rCtrlR[3],-1.0); 
  // control node [x,y,vy] = [-1/3,-1/3,1/3] 
  fCtrlL[4] = 0.02777777777777778*(1.732050807568877*(fl[12]+3.0*fl[4]-3.0*(fl[2]+fl[1]))-3.0*(fl[9]+fl[8])+3.0*fl[5]+9.0*fl[0])*limTheta(rCtrlL[4],1.0); 
  fCtrlR[4] = 0.02777777777777778*(1.732050807568877*(fr[12]+3.0*fr[4]-3.0*(fr[2]+fr[1]))-3.0*(fr[9]+fr[8])+3.0*fr[5]+9.0*fr[0])*limTheta(rCtrlR[4],-1.0); 
  // control node [x,y,vy] = [1/3,-1/3,1/3] 
  fCtrlL[5] = -0.02777777777777778*(1.732050807568877*(fl[12]-3.0*fl[4]+3.0*fl[2]-3.0*fl[1])+3.0*fl[9]-3.0*fl[8]+3.0*fl[5]-9.0*fl[0])*limTheta(rCtrlL[5],1.0); 
  fCtrlR[5] = -0.02777777777777778*(1.732050807568877*(fr[12]-3.0*fr[4]+3.0*fr[2]-3.0*fr[1])+3.0*fr[9]-3.0*fr[8]+3.0*fr[5]-9.0*fr[0])*limTheta(rCtrlR[5],-1.0); 
  // control node [x,y,vy] = [-1/3,1/3,1/3] 
  fCtrlL[6] = -0.02777777777777778*(1.732050807568877*(fl[12]-3.0*(fl[4]+fl[2])+3.0*fl[1])-3.0*fl[9]+3.0*(fl[8]+fl[5])-9.0*fl[0])*limTheta(rCtrlL[6],1.0); 
  fCtrlR[6] = -0.02777777777777778*(1.732050807568877*(fr[12]-3.0*(fr[4]+fr[2])+3.0*fr[1])-3.0*fr[9]+3.0*(fr[8]+fr[5])-9.0*fr[0])*limTheta(rCtrlR[6],-1.0); 
  // control node [x,y,vy] = [1/3,1/3,1/3] 
  fCtrlL[7] = 0.02777777777777778*(1.732050807568877*(fl[12]+3.0*(fl[4]+fl[2]+fl[1]))+3.0*(fl[9]+fl[8]+fl[5])+9.0*fl[0])*limTheta(rCtrlL[7],1.0); 
  fCtrlR[7] = 0.02777777777777778*(1.732050807568877*(fr[12]+3.0*(fr[4]+fr[2]+fr[1]))+3.0*(fr[9]+fr[8]+fr[5])+9.0*fr[0])*limTheta(rCtrlR[7],-1.0); 
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
  fhatALQuad[0] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5]+fL_AL[4])+1.414213562373095*(fL_AL[3]+fL_AL[2]+fL_AL[1])-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[0] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5]+fR_AL[4])+1.414213562373095*(fR_AL[3]+fR_AL[2]+fR_AL[1])-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0])-0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2]); 
  if(alphaQuad > 0) {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*(fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2])+1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = 0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*(fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2])+1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.3535533905932737*alpha[5]-0.3535533905932737*(alpha[4]+alpha[3])+0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*fL_AL[5]-1.414213562373095*(fL_AL[4]+fL_AL[3])+1.414213562373095*fL_AL[2]-1.414213562373095*fL_AL[1]+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[2] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*fR_AL[5]-1.414213562373095*(fR_AL[4]+fR_AL[3])+1.414213562373095*fR_AL[2]-1.414213562373095*fR_AL[1]+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = (-0.3535533905932737*alpha[5])+0.3535533905932737*alpha[4]-0.3535533905932737*alpha[3]+0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6]+fL_AL[5])-1.414213562373095*fL_AL[4]+1.414213562373095*fL_AL[3]-1.414213562373095*(fL_AL[2]+fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[3] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6]+fR_AL[5])-1.414213562373095*fR_AL[4]+1.414213562373095*fR_AL[3]-1.414213562373095*(fR_AL[2]+fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = (-0.3535533905932737*alpha[5])+0.3535533905932737*(alpha[4]+alpha[3])-0.3535533905932737*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[4] = 0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*(fL_AL[6]+fL_AL[5])+1.414213562373095*(fL_AL[4]+fL_AL[3])-1.414213562373095*(fL_AL[2]+fL_AL[1])+1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[4] = 0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*(fR_AL[6]+fR_AL[5])+1.414213562373095*(fR_AL[4]+fR_AL[3])-1.414213562373095*(fR_AL[2]+fR_AL[1])+1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*alpha[5]-0.3535533905932737*alpha[4]+0.3535533905932737*alpha[3]-0.3535533905932737*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fL_AL[7]+fL_AL[6])-1.414213562373095*fL_AL[5]+1.414213562373095*fL_AL[4]-1.414213562373095*fL_AL[3]+1.414213562373095*fL_AL[2]-1.414213562373095*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[5] = -0.25*(1.414213562373095*(fR_AL[7]+fR_AL[6])-1.414213562373095*fR_AL[5]+1.414213562373095*fR_AL[4]-1.414213562373095*fR_AL[3]+1.414213562373095*fR_AL[2]-1.414213562373095*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = (-0.3535533905932737*(alpha[5]+alpha[4]))+0.3535533905932737*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[6] = -0.25*(1.414213562373095*fL_AL[7]-1.414213562373095*fL_AL[6]+1.414213562373095*(fL_AL[5]+fL_AL[4])-1.414213562373095*(fL_AL[3]+fL_AL[2])+1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[6] = -0.25*(1.414213562373095*fR_AL[7]-1.414213562373095*fR_AL[6]+1.414213562373095*(fR_AL[5]+fR_AL[4])-1.414213562373095*(fR_AL[3]+fR_AL[2])+1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = 0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[7] = 0.3535533905932738*(fL_AL[7]+fL_AL[6]+fL_AL[5]+fL_AL[4]+fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[7] = 0.3535533905932738*(fR_AL[7]+fR_AL[6]+fR_AL[5]+fR_AL[4]+fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
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
 
  incr[0] = 0.25*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[1] = 0.25*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[3] = -0.4330127018922193*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[4] = 0.25*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_v; 
  incr[5] = 0.25*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 
  incr[6] = -0.4330127018922193*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[7] = -0.4330127018922193*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[8] = 0.25*(alpha[2]*fhatAL[7]+alpha[4]*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_v; 
  incr[9] = 0.25*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_v; 
  incr[10] = -0.4330127018922193*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_v; 
  incr[11] = -0.4330127018922193*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 
  incr[12] = 0.25*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])*dfac_v; 
  incr[13] = -0.4330127018922193*(alpha[2]*fhatAL[7]+alpha[4]*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_v; 
  incr[14] = -0.4330127018922193*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_v; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])*dfac_v; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[16], outrPos[16]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.25 : positivityWeightByDirL[3]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.25 : positivityWeightByDirR[3]/positivityWeightByDirR[0]; 
  outlPos[0] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[1] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[2] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[3] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[4] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[5] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outlPos[6] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outlPos[7] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outlPos[8] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outlPos[9] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outlPos[10] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[11] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[12] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[13] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outlPos[14] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outlPos[15] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[0] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[1] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*incr[8]-3.0*incr[7]+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[2] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*incr[8]-3.0*incr[7]+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[3] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8]+incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[4] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8])-3.0*(incr[7]+incr[6])+3.0*incr[5]+9.0*incr[0]); 
  outrPos[5] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*(incr[8]+incr[7])+3.0*incr[6]-3.0*incr[5]+9.0*incr[0]); 
  outrPos[6] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*(incr[8]+incr[7])+3.0*(incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[7] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8])-3.0*(incr[7]+incr[6]+incr[5])-9.0*incr[0]); 
  outrPos[8] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]+incr[12]-1.0*incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2]+incr[1]))-3.0*(incr[10]+incr[9]+incr[8])+3.0*(incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[9] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12])+incr[11]+3.0*incr[4]-3.0*(incr[3]+incr[2])+3.0*incr[1])-3.0*(incr[10]+incr[9])+3.0*(incr[8]+incr[7])-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[10] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]-1.0*incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*incr[2]+3.0*incr[1])+3.0*incr[10]-3.0*incr[9]+3.0*(incr[8]+incr[7])-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[11] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]-1.0*incr[12]+incr[11]-3.0*incr[4]+3.0*incr[3]-3.0*(incr[2]+incr[1]))+3.0*incr[10]-3.0*(incr[9]+incr[8])+3.0*(incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[12] = -0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]+incr[13]-1.0*(incr[12]+incr[11])-3.0*(incr[4]+incr[3])+3.0*(incr[2]+incr[1]))-3.0*incr[10]+3.0*(incr[9]+incr[8]+incr[7]+incr[6])-3.0*incr[5]-9.0*incr[0]); 
  outrPos[13] = -0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]-1.0*incr[13]+incr[12]+incr[11]-3.0*(incr[4]+incr[3])+3.0*incr[2]-3.0*incr[1])-3.0*incr[10]+3.0*incr[9]-3.0*incr[8]+3.0*incr[7]-3.0*incr[6]+3.0*incr[5]-9.0*incr[0]); 
  outrPos[14] = 0.02777777777777778*((-1.0*incr[15])+1.732050807568877*(incr[14]-1.0*(incr[13]+incr[12]+incr[11])+3.0*(incr[4]+incr[3]+incr[2])-3.0*incr[1])+3.0*(incr[10]+incr[9])-3.0*incr[8]+3.0*incr[7]-3.0*(incr[6]+incr[5])+9.0*incr[0]); 
  outrPos[15] = 0.02777777777777778*(incr[15]+1.732050807568877*(incr[14]+incr[13]+incr[12]+incr[11]+3.0*(incr[4]+incr[3]+incr[2]+incr[1]))+3.0*(incr[10]+incr[9]+incr[8]+incr[7]+incr[6]+incr[5])+9.0*incr[0]); 
  if(outlPos[4] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]+fl[13]-1.0*fl[12]+fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*(fl[2]+fl[1]))-3.0*fl[10]+3.0*(fl[9]+fl[8])-3.0*(fl[7]+fl[6])+3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[4]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]+fr[11]+3.0*(fr[4]+fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[4] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[4] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[5] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]-1.0*fl[11]-3.0*fl[4]+3.0*fl[3]-3.0*fl[2]+3.0*fl[1])-3.0*fl[10]+3.0*fl[9]-3.0*(fl[8]+fl[7])+3.0*fl[6]-3.0*fl[5]+9.0*fl[0]))/dtApprox/outlPos[5]); 
  else limFac = 1.0; 
  if(outrPos[1] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12]+fr[11])+3.0*(fr[4]+fr[3]+fr[2])-3.0*fr[1])-3.0*(fr[10]+fr[9])+3.0*fr[8]-3.0*fr[7]+3.0*(fr[6]+fr[5])-9.0*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[5] *= limFac; 
  outlPos[1] *= limFac; 
  outrPos[5] *= limFac; 
  outrPos[1] *= limFac; 
  if(outlPos[6] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12])+fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2])+3.0*fl[1])+3.0*(fl[10]+fl[9])-3.0*(fl[8]+fl[7])+3.0*(fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[6]); 
  else limFac = 1.0; 
  if(outrPos[2] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]+fr[11]-3.0*(fr[4]+fr[3])+3.0*fr[2]-3.0*fr[1])+3.0*fr[10]-3.0*fr[9]+3.0*fr[8]-3.0*fr[7]+3.0*fr[6]-3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[2]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[6] *= limFac; 
  outlPos[2] *= limFac; 
  outrPos[6] *= limFac; 
  outrPos[2] *= limFac; 
  if(outlPos[7] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]-1.0*fl[11]+3.0*fl[4]-3.0*(fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8])-3.0*(fl[7]+fl[6]+fl[5])-9.0*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[3] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]+fr[13]-1.0*(fr[12]+fr[11])-3.0*(fr[4]+fr[3])+3.0*(fr[2]+fr[1]))+3.0*fr[10]-3.0*(fr[9]+fr[8]+fr[7]+fr[6])+3.0*fr[5]+9.0*fr[0]))/dtApprox/outrPos[3]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[3] *= limFac; 
  outrPos[7] *= limFac; 
  outrPos[3] *= limFac; 
  if(outlPos[12] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]+fl[13]-1.0*(fl[12]+fl[11])-3.0*(fl[4]+fl[3])+3.0*(fl[2]+fl[1]))-3.0*fl[10]+3.0*(fl[9]+fl[8]+fl[7]+fl[6])-3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[12]); 
  else limFac = 1.0; 
  if(outrPos[8] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]+fr[13]+fr[12]-1.0*fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2]+fr[1]))-3.0*(fr[10]+fr[9]+fr[8])+3.0*(fr[7]+fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[8]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[12] *= limFac; 
  outlPos[8] *= limFac; 
  outrPos[12] *= limFac; 
  outrPos[8] *= limFac; 
  if(outlPos[13] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]-1.0*fl[13]+fl[12]+fl[11]-3.0*(fl[4]+fl[3])+3.0*fl[2]-3.0*fl[1])-3.0*fl[10]+3.0*fl[9]-3.0*fl[8]+3.0*fl[7]-3.0*fl[6]+3.0*fl[5]-9.0*fl[0]))/dtApprox/outlPos[13]); 
  else limFac = 1.0; 
  if(outrPos[9] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]-1.0*(fr[13]+fr[12])+fr[11]+3.0*fr[4]-3.0*(fr[3]+fr[2])+3.0*fr[1])-3.0*(fr[10]+fr[9])+3.0*(fr[8]+fr[7])-3.0*(fr[6]+fr[5])+9.0*fr[0]))/dtApprox/outrPos[9]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[13] *= limFac; 
  outlPos[9] *= limFac; 
  outrPos[13] *= limFac; 
  outrPos[9] *= limFac; 
  if(outlPos[14] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*((-1.0*fl[15])+1.732050807568877*(fl[14]-1.0*(fl[13]+fl[12]+fl[11])+3.0*(fl[4]+fl[3]+fl[2])-3.0*fl[1])+3.0*(fl[10]+fl[9])-3.0*fl[8]+3.0*fl[7]-3.0*(fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[14]); 
  else limFac = 1.0; 
  if(outrPos[10] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*((-1.0*fr[15])+1.732050807568877*(fr[14]-1.0*fr[13]+fr[12]-1.0*fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*fr[2]+3.0*fr[1])+3.0*fr[10]-3.0*fr[9]+3.0*(fr[8]+fr[7])-3.0*fr[6]+3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[10]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[14] *= limFac; 
  outlPos[10] *= limFac; 
  outrPos[14] *= limFac; 
  outrPos[10] *= limFac; 
  if(outlPos[15] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]+fl[11]+3.0*(fl[4]+fl[3]+fl[2]+fl[1]))+3.0*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+9.0*fl[0]))/dtApprox/outlPos[15]); 
  else limFac = 1.0; 
  if(outrPos[11] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]+fr[13]-1.0*fr[12]+fr[11]-3.0*fr[4]+3.0*fr[3]-3.0*(fr[2]+fr[1]))+3.0*fr[10]-3.0*(fr[9]+fr[8])+3.0*(fr[7]+fr[6])-3.0*fr[5]-9.0*fr[0]))/dtApprox/outrPos[11]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[15] *= limFac; 
  outlPos[11] *= limFac; 
  outrPos[15] *= limFac; 
  outrPos[11] *= limFac; 
  outr[0] += 0.25*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.4330127018922193*(outrPos[15]-1.0*outrPos[14]+outrPos[13]-1.0*outrPos[12]+outrPos[11]-1.0*outrPos[10]+outrPos[9]-1.0*outrPos[8]+outrPos[7]-1.0*outrPos[6]+outrPos[5]-1.0*outrPos[4]+outrPos[3]-1.0*outrPos[2]+outrPos[1]-1.0*outrPos[0]); 
  outr[2] += 0.4330127018922193*(outrPos[15]+outrPos[14]-1.0*(outrPos[13]+outrPos[12])+outrPos[11]+outrPos[10]-1.0*(outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]-1.0*(outrPos[5]+outrPos[4])+outrPos[3]+outrPos[2]-1.0*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]-1.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8])+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]-1.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[4] += 0.4330127018922193*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]-1.0*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[5] += 0.25*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*(outrPos[8]+outrPos[7])-3.0*(outrPos[6]+outrPos[5])+3.0*(outrPos[4]+outrPos[3])-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 
  outr[6] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*(outrPos[8]+outrPos[7])-3.0*outrPos[6]+3.0*outrPos[5]-3.0*(outrPos[4]+outrPos[3])+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[7] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])-3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[8] += 0.25*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*outrPos[12]+3.0*outrPos[11]-3.0*outrPos[10]+3.0*outrPos[9]-3.0*(outrPos[8]+outrPos[7])+3.0*outrPos[6]-3.0*outrPos[5]+3.0*outrPos[4]-3.0*outrPos[3]+3.0*outrPos[2]-3.0*outrPos[1]+3.0*outrPos[0]); 
  outr[9] += 0.25*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12])+3.0*(outrPos[11]+outrPos[10])-3.0*(outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4])-3.0*(outrPos[3]+outrPos[2])+3.0*(outrPos[1]+outrPos[0])); 
  outr[10] += 0.25*(3.0*(outrPos[15]+outrPos[14]+outrPos[13]+outrPos[12])-3.0*(outrPos[11]+outrPos[10]+outrPos[9]+outrPos[8]+outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])+3.0*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0])); 
  outr[11] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*outrPos[12]-3.0*outrPos[11]+3.0*(outrPos[10]+outrPos[9])-3.0*outrPos[8]+3.0*outrPos[7]-3.0*(outrPos[6]+outrPos[5])+3.0*outrPos[4]-3.0*outrPos[3]+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[12] += 0.4330127018922193*(3.0*outrPos[15]-3.0*(outrPos[14]+outrPos[13])+3.0*(outrPos[12]+outrPos[11])-3.0*(outrPos[10]+outrPos[9])+3.0*outrPos[8]-3.0*outrPos[7]+3.0*(outrPos[6]+outrPos[5])-3.0*(outrPos[4]+outrPos[3])+3.0*(outrPos[2]+outrPos[1])-3.0*outrPos[0]); 
  outr[13] += 0.4330127018922193*(3.0*outrPos[15]-3.0*outrPos[14]+3.0*outrPos[13]-3.0*(outrPos[12]+outrPos[11])+3.0*outrPos[10]-3.0*outrPos[9]+3.0*outrPos[8]-3.0*outrPos[7]+3.0*outrPos[6]-3.0*outrPos[5]+3.0*(outrPos[4]+outrPos[3])-3.0*outrPos[2]+3.0*outrPos[1]-3.0*outrPos[0]); 
  outr[14] += 0.4330127018922193*(3.0*(outrPos[15]+outrPos[14])-3.0*(outrPos[13]+outrPos[12]+outrPos[11]+outrPos[10])+3.0*(outrPos[9]+outrPos[8])-3.0*(outrPos[7]+outrPos[6])+3.0*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])-3.0*(outrPos[1]+outrPos[0])); 
  outr[15] += 0.25*(9.0*outrPos[15]-9.0*(outrPos[14]+outrPos[13])+9.0*outrPos[12]-9.0*outrPos[11]+9.0*(outrPos[10]+outrPos[9])-9.0*(outrPos[8]+outrPos[7])+9.0*(outrPos[6]+outrPos[5])-9.0*outrPos[4]+9.0*outrPos[3]-9.0*(outrPos[2]+outrPos[1])+9.0*outrPos[0]); 

  outl[0] += 0.25*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.4330127018922193*(outlPos[15]-1.0*outlPos[14]+outlPos[13]-1.0*outlPos[12]+outlPos[11]-1.0*outlPos[10]+outlPos[9]-1.0*outlPos[8]+outlPos[7]-1.0*outlPos[6]+outlPos[5]-1.0*outlPos[4]+outlPos[3]-1.0*outlPos[2]+outlPos[1]-1.0*outlPos[0]); 
  outl[2] += 0.4330127018922193*(outlPos[15]+outlPos[14]-1.0*(outlPos[13]+outlPos[12])+outlPos[11]+outlPos[10]-1.0*(outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]-1.0*(outlPos[5]+outlPos[4])+outlPos[3]+outlPos[2]-1.0*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]-1.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8])+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]-1.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[4] += 0.4330127018922193*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]-1.0*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[5] += 0.25*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*(outlPos[8]+outlPos[7])-3.0*(outlPos[6]+outlPos[5])+3.0*(outlPos[4]+outlPos[3])-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  outl[6] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*(outlPos[8]+outlPos[7])-3.0*outlPos[6]+3.0*outlPos[5]-3.0*(outlPos[4]+outlPos[3])+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[7] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])-3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[8] += 0.25*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*outlPos[12]+3.0*outlPos[11]-3.0*outlPos[10]+3.0*outlPos[9]-3.0*(outlPos[8]+outlPos[7])+3.0*outlPos[6]-3.0*outlPos[5]+3.0*outlPos[4]-3.0*outlPos[3]+3.0*outlPos[2]-3.0*outlPos[1]+3.0*outlPos[0]); 
  outl[9] += 0.25*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12])+3.0*(outlPos[11]+outlPos[10])-3.0*(outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4])-3.0*(outlPos[3]+outlPos[2])+3.0*(outlPos[1]+outlPos[0])); 
  outl[10] += 0.25*(3.0*(outlPos[15]+outlPos[14]+outlPos[13]+outlPos[12])-3.0*(outlPos[11]+outlPos[10]+outlPos[9]+outlPos[8]+outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])+3.0*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0])); 
  outl[11] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*outlPos[12]-3.0*outlPos[11]+3.0*(outlPos[10]+outlPos[9])-3.0*outlPos[8]+3.0*outlPos[7]-3.0*(outlPos[6]+outlPos[5])+3.0*outlPos[4]-3.0*outlPos[3]+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[12] += 0.4330127018922193*(3.0*outlPos[15]-3.0*(outlPos[14]+outlPos[13])+3.0*(outlPos[12]+outlPos[11])-3.0*(outlPos[10]+outlPos[9])+3.0*outlPos[8]-3.0*outlPos[7]+3.0*(outlPos[6]+outlPos[5])-3.0*(outlPos[4]+outlPos[3])+3.0*(outlPos[2]+outlPos[1])-3.0*outlPos[0]); 
  outl[13] += 0.4330127018922193*(3.0*outlPos[15]-3.0*outlPos[14]+3.0*outlPos[13]-3.0*(outlPos[12]+outlPos[11])+3.0*outlPos[10]-3.0*outlPos[9]+3.0*outlPos[8]-3.0*outlPos[7]+3.0*outlPos[6]-3.0*outlPos[5]+3.0*(outlPos[4]+outlPos[3])-3.0*outlPos[2]+3.0*outlPos[1]-3.0*outlPos[0]); 
  outl[14] += 0.4330127018922193*(3.0*(outlPos[15]+outlPos[14])-3.0*(outlPos[13]+outlPos[12]+outlPos[11]+outlPos[10])+3.0*(outlPos[9]+outlPos[8])-3.0*(outlPos[7]+outlPos[6])+3.0*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])-3.0*(outlPos[1]+outlPos[0])); 
  outl[15] += 0.25*(9.0*outlPos[15]-9.0*(outlPos[14]+outlPos[13])+9.0*outlPos[12]-9.0*outlPos[11]+9.0*(outlPos[10]+outlPos[9])-9.0*(outlPos[8]+outlPos[7])+9.0*(outlPos[6]+outlPos[5])-9.0*outlPos[4]+9.0*outlPos[3]-9.0*(outlPos[2]+outlPos[1])+9.0*outlPos[0]); 
  return std::abs(alpha0); 
} 
