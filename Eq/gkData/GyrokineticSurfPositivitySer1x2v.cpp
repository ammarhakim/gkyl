#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = 1.414213562373095*wv; 

  double alpha[8]; 
  alpha[0] = 2.828427124746191*wv; 
  if (alpha0>0) { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[1] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rVal[5] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rVal[6] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[7] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cfl); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cfl); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[1]+fhatALVal[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[1]+1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[4]+fhatALVal[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[5]+fhatALVal[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[4]+1.732050807568877*fhatALVal[2])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[5]+1.732050807568877*fhatALVal[3])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[7]+fhatALVal[6])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[7]+1.732050807568877*fhatALVal[6])*dfac_x; 

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
  } else { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[1] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rVal[5] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rVal[6] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[7] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cfl); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cfl); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[1]-1.0*fhatALVal[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[1]-1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[4]-1.0*fhatALVal[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[5]-1.0*fhatALVal[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[4]-1.732050807568877*fhatALVal[2])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[5]-1.732050807568877*fhatALVal[3])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[7]-1.0*fhatALVal[6])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[7]-1.732050807568877*fhatALVal[6])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.732050807568877*Phi[1]*dfac_x*q_)/m_; 

  double alpha[8]; 
  alpha[0] = -(3.464101615137754*Phi[1]*dfac_x*q_)/m_; 
  if (alpha0>0) { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[2] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rVal[5] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[6] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rVal[7] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cfl); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cfl); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[2]+fhatALVal[0])*dfac_v; 
  incr[1] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[4]+fhatALVal[1])*dfac_v; 
  incr[2] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[2]+1.732050807568877*fhatALVal[0])*dfac_v; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[6]+fhatALVal[3])*dfac_v; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[4]+1.732050807568877*fhatALVal[1])*dfac_v; 
  incr[5] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[7]+fhatALVal[5])*dfac_v; 
  incr[6] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[6]+1.732050807568877*fhatALVal[3])*dfac_v; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[7]+1.732050807568877*fhatALVal[5])*dfac_v; 

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
  } else { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[2] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rVal[5] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[6] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rVal[7] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cfl); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cfl); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[2]-1.0*fhatALVal[0])*dfac_v; 
  incr[1] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[4]-1.0*fhatALVal[1])*dfac_v; 
  incr[2] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[2]-1.732050807568877*fhatALVal[0])*dfac_v; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[6]-1.0*fhatALVal[3])*dfac_v; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[4]-1.732050807568877*fhatALVal[1])*dfac_v; 
  incr[5] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[7]-1.0*fhatALVal[5])*dfac_v; 
  incr[6] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[6]-1.732050807568877*fhatALVal[3])*dfac_v; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[7]-1.732050807568877*fhatALVal[5])*dfac_v; 

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
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dAparcfl, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = 1.414213562373095*wv; 

  double alpha[8]; 
  alpha[0] = 2.828427124746191*wv; 
  if (alpha0>0) { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[1] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rVal[5] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rVal[6] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[7] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cfl); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cfl); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[1]+fhatALVal[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[1]+1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[4]+fhatALVal[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[5]+fhatALVal[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[4]+1.732050807568877*fhatALVal[2])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[5]+1.732050807568877*fhatALVal[3])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[7]+fhatALVal[6])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[7]+1.732050807568877*fhatALVal[6])*dfac_x; 

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
  } else { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[1] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rVal[5] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rVal[6] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[7] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cfl); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cfl); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[1]-1.0*fhatALVal[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[1]-1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[4]-1.0*fhatALVal[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[5]-1.0*fhatALVal[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[4]-1.732050807568877*fhatALVal[2])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[5]-1.732050807568877*fhatALVal[3])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[7]-1.0*fhatALVal[6])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[7]-1.732050807568877*fhatALVal[6])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dAparcfl, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.0*(1.732050807568877*Phi[1]*dfac_x+dAparcfl[0])*q_)/m_; 

  double alpha[8]; 
  alpha[0] = (-(3.464101615137754*Phi[1]*dfac_x*q_)/m_)-(2.0*dAparcfl[0]*q_)/m_; 
  alpha[1] = -(2.0*dAparcfl[1]*q_)/m_; 
  if (alpha0>0) { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[2] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rVal[5] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[6] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rVal[7] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cfl); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cfl); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.1767766952966368*(1.732050807568877*alpha[1]*fhatALVal[4]+1.732050807568877*alpha[0]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*alpha[0]*fhatALVal[4]+1.732050807568877*alpha[1]*fhatALVal[2]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.1767766952966368*(3.0*alpha[1]*fhatALVal[4]+3.0*alpha[0]*fhatALVal[2]+1.732050807568877*alpha[1]*fhatALVal[1]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.1767766952966368*(1.732050807568877*alpha[1]*fhatALVal[7]+1.732050807568877*alpha[0]*fhatALVal[6]+alpha[1]*fhatALVal[5]+alpha[0]*fhatALVal[3])*dfac_v; 
  incr[4] = -0.1767766952966368*(3.0*alpha[0]*fhatALVal[4]+3.0*alpha[1]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[1]+1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*alpha[0]*fhatALVal[7]+1.732050807568877*alpha[1]*fhatALVal[6]+alpha[0]*fhatALVal[5]+alpha[1]*fhatALVal[3])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*alpha[1]*fhatALVal[7]+3.0*alpha[0]*fhatALVal[6]+1.732050807568877*alpha[1]*fhatALVal[5]+1.732050807568877*alpha[0]*fhatALVal[3])*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*alpha[0]*fhatALVal[7]+3.0*alpha[1]*fhatALVal[6]+1.732050807568877*alpha[0]*fhatALVal[5]+1.732050807568877*alpha[1]*fhatALVal[3])*dfac_v; 

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
  } else { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[2] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rVal[5] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[6] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rVal[7] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cfl); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cfl); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = -0.1767766952966368*(1.732050807568877*alpha[1]*fhatALVal[4]+1.732050807568877*alpha[0]*fhatALVal[2]-1.0*alpha[1]*fhatALVal[1]-1.0*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*alpha[0]*fhatALVal[4]+1.732050807568877*alpha[1]*fhatALVal[2]-1.0*alpha[0]*fhatALVal[1]-1.0*fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(3.0*alpha[1]*fhatALVal[4]+3.0*alpha[0]*fhatALVal[2]-1.732050807568877*alpha[1]*fhatALVal[1]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = -0.1767766952966368*(1.732050807568877*alpha[1]*fhatALVal[7]+1.732050807568877*alpha[0]*fhatALVal[6]-1.0*alpha[1]*fhatALVal[5]-1.0*alpha[0]*fhatALVal[3])*dfac_v; 
  incr[4] = 0.1767766952966368*(3.0*alpha[0]*fhatALVal[4]+3.0*alpha[1]*fhatALVal[2]-1.732050807568877*alpha[0]*fhatALVal[1]-1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*alpha[0]*fhatALVal[7]+1.732050807568877*alpha[1]*fhatALVal[6]-1.0*alpha[0]*fhatALVal[5]-1.0*alpha[1]*fhatALVal[3])*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*alpha[1]*fhatALVal[7]+3.0*alpha[0]*fhatALVal[6]-1.732050807568877*alpha[1]*fhatALVal[5]-1.732050807568877*alpha[0]*fhatALVal[3])*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*alpha[0]*fhatALVal[7]+3.0*alpha[1]*fhatALVal[6]-1.732050807568877*alpha[0]*fhatALVal[5]-1.732050807568877*alpha[1]*fhatALVal[3])*dfac_v; 

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
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = 1.414213562373095*wv; 

  double alpha[8]; 
  alpha[0] = 2.828427124746191*wv; 
  if (alpha0>0) { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[1] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rVal[5] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rVal[6] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[7] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cfl); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cfl); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[1]+fhatALVal[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[1]+1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[4]+fhatALVal[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[5]+fhatALVal[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[4]+1.732050807568877*fhatALVal[2])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[5]+1.732050807568877*fhatALVal[3])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[7]+fhatALVal[6])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[7]+1.732050807568877*fhatALVal[6])*dfac_x; 

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
  } else { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[1] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rVal[5] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rVal[6] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[7] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cfl); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cfl); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[1]-1.0*fhatALVal[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[1]-1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[4]-1.0*fhatALVal[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[5]-1.0*fhatALVal[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[4]-1.732050807568877*fhatALVal[2])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[5]-1.732050807568877*fhatALVal[3])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[7]-1.0*fhatALVal[6])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[7]-1.732050807568877*fhatALVal[6])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.0*(1.732050807568877*Bmag[1]*dfac_x*wm+1.732050807568877*Phi[1]*dfac_x*q_))/m_; 

  double alpha[8]; 
  alpha[0] = (-(3.464101615137754*Bmag[1]*dfac_x*wm)/m_)-(3.464101615137754*Phi[1]*dfac_x*q_)/m_; 
  alpha[3] = -(2.0*Bmag[1]*dfac_x)/(dfac_m*m_); 
  if (alpha0>0) { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[2] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rVal[5] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[6] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rVal[7] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cfl); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cfl); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.1767766952966368*(1.732050807568877*alpha[3]*fhatALVal[6]+alpha[3]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*alpha[3]*fhatALVal[7]+alpha[3]*fhatALVal[5]+1.732050807568877*alpha[0]*fhatALVal[4]+alpha[0]*fhatALVal[1])*dfac_v; 
  incr[2] = -0.1767766952966368*(3.0*alpha[3]*fhatALVal[6]+1.732050807568877*alpha[3]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.1767766952966368*(1.732050807568877*alpha[0]*fhatALVal[6]+alpha[0]*fhatALVal[3]+1.732050807568877*fhatALVal[2]*alpha[3]+fhatALVal[0]*alpha[3])*dfac_v; 
  incr[4] = -0.1767766952966368*(3.0*alpha[3]*fhatALVal[7]+1.732050807568877*alpha[3]*fhatALVal[5]+3.0*alpha[0]*fhatALVal[4]+1.732050807568877*alpha[0]*fhatALVal[1])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*alpha[0]*fhatALVal[7]+alpha[0]*fhatALVal[5]+1.732050807568877*alpha[3]*fhatALVal[4]+fhatALVal[1]*alpha[3])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*alpha[0]*fhatALVal[6]+1.732050807568877*alpha[0]*fhatALVal[3]+3.0*fhatALVal[2]*alpha[3]+1.732050807568877*fhatALVal[0]*alpha[3])*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*alpha[0]*fhatALVal[7]+1.732050807568877*alpha[0]*fhatALVal[5]+3.0*alpha[3]*fhatALVal[4]+1.732050807568877*fhatALVal[1]*alpha[3])*dfac_v; 

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
  } else { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[2] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rVal[5] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[6] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rVal[7] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cfl); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cfl); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = -0.1767766952966368*(1.732050807568877*alpha[3]*fhatALVal[6]-1.0*alpha[3]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]-1.0*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*alpha[3]*fhatALVal[7]-1.0*alpha[3]*fhatALVal[5]+1.732050807568877*alpha[0]*fhatALVal[4]-1.0*alpha[0]*fhatALVal[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(3.0*alpha[3]*fhatALVal[6]-1.732050807568877*alpha[3]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = -0.1767766952966368*(1.732050807568877*alpha[0]*fhatALVal[6]-1.0*alpha[0]*fhatALVal[3]+1.732050807568877*fhatALVal[2]*alpha[3]-1.0*fhatALVal[0]*alpha[3])*dfac_v; 
  incr[4] = 0.1767766952966368*(3.0*alpha[3]*fhatALVal[7]-1.732050807568877*alpha[3]*fhatALVal[5]+3.0*alpha[0]*fhatALVal[4]-1.732050807568877*alpha[0]*fhatALVal[1])*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*alpha[0]*fhatALVal[7]-1.0*alpha[0]*fhatALVal[5]+1.732050807568877*alpha[3]*fhatALVal[4]-1.0*fhatALVal[1]*alpha[3])*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*alpha[0]*fhatALVal[6]-1.732050807568877*alpha[0]*fhatALVal[3]+3.0*fhatALVal[2]*alpha[3]-1.732050807568877*fhatALVal[0]*alpha[3])*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*alpha[0]*fhatALVal[7]-1.732050807568877*alpha[0]*fhatALVal[5]+3.0*alpha[3]*fhatALVal[4]-1.732050807568877*fhatALVal[1]*alpha[3])*dfac_v; 

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
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dAparcfl, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = 1.414213562373095*wv; 

  double alpha[8]; 
  alpha[0] = 2.828427124746191*wv; 
  if (alpha0>0) { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[1] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rVal[5] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rVal[6] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[7] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cfl); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cfl); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[1]+fhatALVal[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[1]+1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[4]+fhatALVal[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[5]+fhatALVal[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[4]+1.732050807568877*fhatALVal[2])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[5]+1.732050807568877*fhatALVal[3])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[7]+fhatALVal[6])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fhatALVal[7]+1.732050807568877*fhatALVal[6])*dfac_x; 

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
  } else { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[1] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rVal[5] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rVal[6] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[7] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cfl); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cfl); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[1]-1.0*fhatALVal[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[1]-1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[4]-1.0*fhatALVal[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[5]-1.0*fhatALVal[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[4]-1.732050807568877*fhatALVal[2])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[5]-1.732050807568877*fhatALVal[3])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(1.732050807568877*fhatALVal[7]-1.0*fhatALVal[6])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fhatALVal[7]-1.732050807568877*fhatALVal[6])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dAparcfl, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.0*(1.732050807568877*Bmag[1]*dfac_x*wm+(1.732050807568877*Phi[1]*dfac_x+dAparcfl[0])*q_))/m_; 

  double alpha[8]; 
  alpha[0] = (-(3.464101615137754*Bmag[1]*dfac_x*wm)/m_)-(3.464101615137754*Phi[1]*dfac_x*q_)/m_-(2.0*dAparcfl[0]*q_)/m_; 
  alpha[1] = -(2.0*dAparcfl[1]*q_)/m_; 
  alpha[3] = -(2.0*Bmag[1]*dfac_x)/(dfac_m*m_); 
  if (alpha0>0) { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[2] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rVal[5] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[6] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rVal[7] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cfl); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cfl); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.1767766952966368*(1.732050807568877*alpha[3]*fhatALVal[6]+1.732050807568877*alpha[1]*fhatALVal[4]+alpha[3]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*alpha[3]*fhatALVal[7]+alpha[3]*fhatALVal[5]+1.732050807568877*alpha[0]*fhatALVal[4]+1.732050807568877*alpha[1]*fhatALVal[2]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.1767766952966368*(3.0*alpha[3]*fhatALVal[6]+3.0*alpha[1]*fhatALVal[4]+1.732050807568877*alpha[3]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]+1.732050807568877*alpha[1]*fhatALVal[1]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.1767766952966368*(1.732050807568877*alpha[1]*fhatALVal[7]+1.732050807568877*alpha[0]*fhatALVal[6]+alpha[1]*fhatALVal[5]+alpha[0]*fhatALVal[3]+1.732050807568877*fhatALVal[2]*alpha[3]+fhatALVal[0]*alpha[3])*dfac_v; 
  incr[4] = -0.1767766952966368*(3.0*alpha[3]*fhatALVal[7]+1.732050807568877*alpha[3]*fhatALVal[5]+3.0*alpha[0]*fhatALVal[4]+3.0*alpha[1]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[1]+1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*alpha[0]*fhatALVal[7]+1.732050807568877*alpha[1]*fhatALVal[6]+alpha[0]*fhatALVal[5]+1.732050807568877*alpha[3]*fhatALVal[4]+alpha[1]*fhatALVal[3]+fhatALVal[1]*alpha[3])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*alpha[1]*fhatALVal[7]+3.0*alpha[0]*fhatALVal[6]+1.732050807568877*alpha[1]*fhatALVal[5]+1.732050807568877*alpha[0]*fhatALVal[3]+3.0*fhatALVal[2]*alpha[3]+1.732050807568877*fhatALVal[0]*alpha[3])*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*alpha[0]*fhatALVal[7]+3.0*alpha[1]*fhatALVal[6]+1.732050807568877*alpha[0]*fhatALVal[5]+3.0*alpha[3]*fhatALVal[4]+1.732050807568877*alpha[1]*fhatALVal[3]+1.732050807568877*fhatALVal[1]*alpha[3])*dfac_v; 

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
  } else { 
  double rVal[8];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[2] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[3] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[4] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rVal[5] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[6] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rVal[7] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cfl); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cfl); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cfl); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cfl); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cfl); 
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = -0.1767766952966368*(1.732050807568877*alpha[3]*fhatALVal[6]+1.732050807568877*alpha[1]*fhatALVal[4]-1.0*alpha[3]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]-1.0*alpha[1]*fhatALVal[1]-1.0*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*alpha[3]*fhatALVal[7]-1.0*alpha[3]*fhatALVal[5]+1.732050807568877*alpha[0]*fhatALVal[4]+1.732050807568877*alpha[1]*fhatALVal[2]-1.0*alpha[0]*fhatALVal[1]-1.0*fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(3.0*alpha[3]*fhatALVal[6]+3.0*alpha[1]*fhatALVal[4]-1.732050807568877*alpha[3]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]-1.732050807568877*alpha[1]*fhatALVal[1]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = -0.1767766952966368*(1.732050807568877*alpha[1]*fhatALVal[7]+1.732050807568877*alpha[0]*fhatALVal[6]-1.0*alpha[1]*fhatALVal[5]-1.0*alpha[0]*fhatALVal[3]+1.732050807568877*fhatALVal[2]*alpha[3]-1.0*fhatALVal[0]*alpha[3])*dfac_v; 
  incr[4] = 0.1767766952966368*(3.0*alpha[3]*fhatALVal[7]-1.732050807568877*alpha[3]*fhatALVal[5]+3.0*alpha[0]*fhatALVal[4]+3.0*alpha[1]*fhatALVal[2]-1.732050807568877*alpha[0]*fhatALVal[1]-1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*alpha[0]*fhatALVal[7]+1.732050807568877*alpha[1]*fhatALVal[6]-1.0*alpha[0]*fhatALVal[5]+1.732050807568877*alpha[3]*fhatALVal[4]-1.0*alpha[1]*fhatALVal[3]-1.0*fhatALVal[1]*alpha[3])*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*alpha[1]*fhatALVal[7]+3.0*alpha[0]*fhatALVal[6]-1.732050807568877*alpha[1]*fhatALVal[5]-1.732050807568877*alpha[0]*fhatALVal[3]+3.0*fhatALVal[2]*alpha[3]-1.732050807568877*fhatALVal[0]*alpha[3])*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*alpha[0]*fhatALVal[7]+3.0*alpha[1]*fhatALVal[6]-1.732050807568877*alpha[0]*fhatALVal[5]+3.0*alpha[3]*fhatALVal[4]-1.732050807568877*alpha[1]*fhatALVal[3]-1.732050807568877*fhatALVal[1]*alpha[3])*dfac_v; 

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
  } 
  return std::abs(alpha0); 
} 
