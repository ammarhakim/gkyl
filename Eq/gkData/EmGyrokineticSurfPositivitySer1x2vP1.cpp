#include <GyrokineticModDecl.h> 
double EmGyrokineticSurfPositivity1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.1767766952966368*(1.732050807568877*BmagInv[0]*geoY[0]*Apar[1]*dfac_x+2.0*Gradpar[0])*wv; 

  double alpha[4]; 
  alpha[0] = 0.5*(2.449489742783178*BmagInv[0]*geoY[0]*Apar[1]*dfac_x+2.828427124746191*Gradpar[0])*wv; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rVal[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[3],1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[1] = -0.6123724356957944*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*fhatALVal[1]*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*fhatALVal[2]*dfac_x; 
  incr[4] = -0.6123724356957944*alpha[0]*fhatALVal[1]*dfac_x; 
  incr[5] = -0.6123724356957944*alpha[0]*fhatALVal[2]*dfac_x; 
  incr[6] = 0.3535533905932737*alpha[0]*fhatALVal[3]*dfac_x; 
  incr[7] = -0.6123724356957944*alpha[0]*fhatALVal[3]*dfac_x; 

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
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rVal[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[3],-1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[1] = -0.6123724356957944*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*fhatALVal[1]*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*fhatALVal[2]*dfac_x; 
  incr[4] = -0.6123724356957944*alpha[0]*fhatALVal[1]*dfac_x; 
  incr[5] = -0.6123724356957944*alpha[0]*fhatALVal[2]*dfac_x; 
  incr[6] = 0.3535533905932737*alpha[0]*fhatALVal[3]*dfac_x; 
  incr[7] = -0.6123724356957944*alpha[0]*fhatALVal[3]*dfac_x; 

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
double EmGyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.125*(Phi[1]*dfac_x*(3.0*BmagInv[0]*geoY[0]*Apar[1]*dfac_x+3.464101615137754*Gradpar[0])+2.828427124746191*dApardt[0])*q_)/m_; 

  double alpha[4]; 
  alpha[0] = -(0.5*(Phi[1]*(3.0*BmagInv[0]*geoY[0]*Apar[1]*dfac_x2+3.464101615137754*Gradpar[0]*dfac_x)+2.828427124746191*dApardt[0])*q_)/m_; 
  alpha[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rVal[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[3],1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.6123724356957944*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.3535533905932737*(alpha[1]*fhatALVal[3]+alpha[0]*fhatALVal[2])*dfac_v; 
  incr[4] = -0.6123724356957944*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[5] = 0.3535533905932737*(alpha[0]*fhatALVal[3]+alpha[1]*fhatALVal[2])*dfac_v; 
  incr[6] = -0.6123724356957944*(alpha[1]*fhatALVal[3]+alpha[0]*fhatALVal[2])*dfac_v; 
  incr[7] = -0.6123724356957944*(alpha[0]*fhatALVal[3]+alpha[1]*fhatALVal[2])*dfac_v; 

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
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rVal[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[3],-1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.6123724356957944*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.3535533905932737*(alpha[1]*fhatALVal[3]+alpha[0]*fhatALVal[2])*dfac_v; 
  incr[4] = -0.6123724356957944*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[5] = 0.3535533905932737*(alpha[0]*fhatALVal[3]+alpha[1]*fhatALVal[2])*dfac_v; 
  incr[6] = -0.6123724356957944*(alpha[1]*fhatALVal[3]+alpha[0]*fhatALVal[2])*dfac_v; 
  incr[7] = -0.6123724356957944*(alpha[0]*fhatALVal[3]+alpha[1]*fhatALVal[2])*dfac_v; 

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
double EmGyrokineticSurfPositivity1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.125*wv*(1.414213562373095*Bmag[1]*((5.196152422706631*BmagInv[1]-3.0*BmagInv[0])*geoY[1]+geoY[0]*(1.732050807568877*BmagInv[0]-3.0*BmagInv[1]))*dfac_x*m_*wv+(((((5.196152422706631*Apar[0]-9.0*Apar[1])*Bmag[1]-7.348469228349534*Apar[1])*BmagInv[1]+BmagInv[0]*((5.196152422706631*Apar[1]-3.0*Apar[0])*Bmag[1]+4.242640687119286*Apar[1]))*geoY[1]+geoY[0]*(((5.196152422706631*Apar[1]-3.0*Apar[0])*Bmag[1]+4.242640687119286*Apar[1])*BmagInv[1]+BmagInv[0]*((1.732050807568877*Apar[0]-3.0*Apar[1])*Bmag[1]-2.449489742783178*Apar[1])))*dfac_x+2.828427124746191*(1.732050807568877*Gradpar[1]-1.0*Gradpar[0]))*q_))/q_; 

  double alpha[4]; 
  alpha[0] = -(0.5*(Bmag[1]*((7.348469228349534*BmagInv[1]-4.242640687119286*BmagInv[0])*geoY[1]+geoY[0]*(2.449489742783178*BmagInv[0]-4.242640687119286*BmagInv[1]))*dfac_x*m_*wv2+(((((5.196152422706631*Apar[0]-9.0*Apar[1])*Bmag[1]-7.348469228349534*Apar[1])*BmagInv[1]+BmagInv[0]*((5.196152422706631*Apar[1]-3.0*Apar[0])*Bmag[1]+4.242640687119286*Apar[1]))*geoY[1]+geoY[0]*(((5.196152422706631*Apar[1]-3.0*Apar[0])*Bmag[1]+4.242640687119286*Apar[1])*BmagInv[1]+BmagInv[0]*((1.732050807568877*Apar[0]-3.0*Apar[1])*Bmag[1]-2.449489742783178*Apar[1])))*dfac_x+4.898979485566357*Gradpar[1]-2.828427124746191*Gradpar[0])*q_*wv))/q_; 
  alpha[1] = -(0.5*Bmag[1]*((4.242640687119286*BmagInv[1]-2.449489742783178*BmagInv[0])*geoY[1]+geoY[0]*(1.414213562373095*BmagInv[0]-2.449489742783178*BmagInv[1]))*dfac_x*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[5]-1.0*(5.196152422706631*fl[4]+9.0*fl[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[2]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[5])-9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[6])+3.0*(fl[3]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rVal[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]+fl[4])+9.0*fl[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[3],1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = -0.6123724356957944*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = 0.3535533905932737*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 
  incr[3] = 0.3535533905932737*(alpha[1]*fhatALVal[3]+alpha[0]*fhatALVal[2])*dfac_x; 
  incr[4] = -0.6123724356957944*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 
  incr[5] = -0.6123724356957944*(alpha[1]*fhatALVal[3]+alpha[0]*fhatALVal[2])*dfac_x; 
  incr[6] = 0.3535533905932737*(alpha[0]*fhatALVal[3]+alpha[1]*fhatALVal[2])*dfac_x; 
  incr[7] = -0.6123724356957944*(alpha[0]*fhatALVal[3]+alpha[1]*fhatALVal[2])*dfac_x; 

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
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[5]-1.0*(5.196152422706631*fr[4]+9.0*fr[1])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[2]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[5])-9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[6])+3.0*(fr[3]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rVal[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]+fr[4])+9.0*fr[1]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[3],-1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = -0.6123724356957944*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = 0.3535533905932737*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 
  incr[3] = 0.3535533905932737*(alpha[1]*fhatALVal[3]+alpha[0]*fhatALVal[2])*dfac_x; 
  incr[4] = -0.6123724356957944*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 
  incr[5] = -0.6123724356957944*(alpha[1]*fhatALVal[3]+alpha[0]*fhatALVal[2])*dfac_x; 
  incr[6] = 0.3535533905932737*(alpha[0]*fhatALVal[3]+alpha[1]*fhatALVal[2])*dfac_x; 
  incr[7] = -0.6123724356957944*(alpha[0]*fhatALVal[3]+alpha[1]*fhatALVal[2])*dfac_x; 

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
double EmGyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.0625*(Bmag[1]*dfac_x*(6.0*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_v*dfac_x*m_*(Bmag[1]*wm+Phi[1]*q_)*wv+(dfac_v*(3.0*(((1.414213562373095*Apar[0]*Bmag[1]-2.0*Apar[1])*BmagInv[1]+1.414213562373095*BmagInv[0]*Apar[1]*Bmag[1])*geoY[1]+geoY[0]*(1.414213562373095*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(1.414213562373095*Apar[0]*Bmag[1]-2.0*Apar[1])))*dfac_x-6.928203230275509*Gradpar[0])*q_-6.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x*m_)*wm)+q_*(dfac_v*(3.0*Phi[1]*(((1.414213562373095*Apar[0]*Bmag[1]-2.0*Apar[1])*BmagInv[1]+1.414213562373095*BmagInv[0]*Apar[1]*Bmag[1])*geoY[1]+geoY[0]*(1.414213562373095*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(1.414213562373095*Apar[0]*Bmag[1]-2.0*Apar[1])))*dfac_x2-6.928203230275509*Gradpar[0]*Phi[1]*dfac_x-5.656854249492382*dApardt[0])*q_-6.0*Bmag[1]*Phi[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_)))/(dfac_v*m_*q_); 

  double alpha[4]; 
  alpha[0] = (0.25*(Bmag[1]*(6.0*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_v*dfac_x2*m_*(Bmag[1]*wm+Phi[1]*q_)*wv+(dfac_v*((((4.242640687119286*Apar[0]*Bmag[1]-6.0*Apar[1])*BmagInv[1]+4.242640687119286*BmagInv[0]*Apar[1]*Bmag[1])*geoY[1]+geoY[0]*(4.242640687119286*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(4.242640687119286*Apar[0]*Bmag[1]-6.0*Apar[1])))*dfac_x2-6.928203230275509*Gradpar[0]*dfac_x)*q_-6.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_)*wm-6.0*Phi[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_*q_)+dfac_v*(Phi[1]*((((4.242640687119286*Apar[0]*Bmag[1]-6.0*Apar[1])*BmagInv[1]+4.242640687119286*BmagInv[0]*Apar[1]*Bmag[1])*geoY[1]+geoY[0]*(4.242640687119286*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(4.242640687119286*Apar[0]*Bmag[1]-6.0*Apar[1])))*dfac_x2-6.928203230275509*Gradpar[0]*dfac_x)-5.656854249492382*dApardt[0])*q2))/(dfac_v*m_*q_); 
  alpha[1] = (0.05*(Bmag[1]*(30.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_v*dfac_x2*m_*(Bmag[1]*wm+Phi[1]*q_)*wv+(dfac_v*(((38.18376618407357*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(21.21320343559643*Apar[0]*Bmag[1]-30.0*Apar[1]))*geoY[1]+geoY[0]*((21.21320343559643*Apar[0]*Bmag[1]-30.0*Apar[1])*BmagInv[1]+21.21320343559643*BmagInv[0]*Apar[1]*Bmag[1]))*dfac_x2-34.64101615137754*Gradpar[1]*dfac_x)*q_-30.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_)*wm-30.0*Phi[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_*q_)+dfac_v*(Phi[1]*(((38.18376618407357*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(21.21320343559643*Apar[0]*Bmag[1]-30.0*Apar[1]))*geoY[1]+geoY[0]*((21.21320343559643*Apar[0]*Bmag[1]-30.0*Apar[1])*BmagInv[1]+21.21320343559643*BmagInv[0]*Apar[1]*Bmag[1]))*dfac_x2-34.64101615137754*Gradpar[1]*dfac_x)-28.28427124746191*dApardt[1])*q2))/(dfac_v*m_*q_); 
  alpha[2] = (0.25*Bmag[1]*(dfac_v*(3.464101615137754*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_*wv+((((2.449489742783178*Apar[0]*Bmag[1]-3.464101615137754*Apar[1])*BmagInv[1]+2.449489742783178*BmagInv[0]*Apar[1]*Bmag[1])*geoY[1]+geoY[0]*(2.449489742783178*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(2.449489742783178*Apar[0]*Bmag[1]-3.464101615137754*Apar[1])))*dfac_x2-4.0*Gradpar[0]*dfac_x)*q_)-3.464101615137754*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_))/(dfac_m*dfac_v*m_*q_); 
  alpha[3] = (0.05*Bmag[1]*(dfac_v*(17.32050807568877*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_*wv+(((22.0454076850486*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(12.24744871391589*Apar[0]*Bmag[1]-17.32050807568877*Apar[1]))*geoY[1]+geoY[0]*((12.24744871391589*Apar[0]*Bmag[1]-17.32050807568877*Apar[1])*BmagInv[1]+12.24744871391589*BmagInv[0]*Apar[1]*Bmag[1]))*dfac_x2-20.0*Gradpar[1]*dfac_x)*q_)-17.32050807568877*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_))/(dfac_m*dfac_v*m_*q_); 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rVal[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[3],1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*(alpha[3]*fhatALVal[3]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(alpha[2]*fhatALVal[3]+fhatALVal[2]*alpha[3]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.6123724356957944*(alpha[3]*fhatALVal[3]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.3535533905932737*(alpha[1]*fhatALVal[3]+fhatALVal[1]*alpha[3]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_v; 
  incr[4] = -0.6123724356957944*(alpha[2]*fhatALVal[3]+fhatALVal[2]*alpha[3]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[5] = 0.3535533905932737*(alpha[0]*fhatALVal[3]+fhatALVal[0]*alpha[3]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_v; 
  incr[6] = -0.6123724356957944*(alpha[1]*fhatALVal[3]+fhatALVal[1]*alpha[3]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_v; 
  incr[7] = -0.6123724356957944*(alpha[0]*fhatALVal[3]+fhatALVal[0]*alpha[3]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_v; 

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
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rVal[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rVal[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rVal[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[3],-1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*(alpha[3]*fhatALVal[3]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(alpha[2]*fhatALVal[3]+fhatALVal[2]*alpha[3]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.6123724356957944*(alpha[3]*fhatALVal[3]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.3535533905932737*(alpha[1]*fhatALVal[3]+fhatALVal[1]*alpha[3]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_v; 
  incr[4] = -0.6123724356957944*(alpha[2]*fhatALVal[3]+fhatALVal[2]*alpha[3]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[5] = 0.3535533905932737*(alpha[0]*fhatALVal[3]+fhatALVal[0]*alpha[3]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_v; 
  incr[6] = -0.6123724356957944*(alpha[1]*fhatALVal[3]+fhatALVal[1]*alpha[3]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_v; 
  incr[7] = -0.6123724356957944*(alpha[0]*fhatALVal[3]+fhatALVal[0]*alpha[3]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_v; 

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
