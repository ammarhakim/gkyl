#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cflVal); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cflVal); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cflVal); 
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
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cflVal); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cflVal); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cflVal); 
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
double GyrokineticSurf1x2vSer_X_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  Ghat[0] = 1.118033988749895*fr[7]*wv+1.118033988749895*fl[7]*wv-0.8660254037844386*fr[1]*wv+0.8660254037844386*fl[1]*wv+0.5*fr[0]*wv+0.5*fl[0]*wv+(0.6454972243679028*fr[11])/dfac_v+(0.6454972243679028*fl[11])/dfac_v-(0.5*fr[4])/dfac_v+(0.5*fl[4])/dfac_v+(0.2886751345948129*fr[2])/dfac_v+(0.2886751345948129*fl[2])/dfac_v-1.118033988749895*fr[7]*amax+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = 1.118033988749895*fr[11]*wv+1.118033988749895*fl[11]*wv-0.8660254037844386*fr[4]*wv+0.8660254037844386*fl[4]*wv+0.5*fr[2]*wv+0.5*fl[2]*wv-(0.447213595499958*fr[12])/dfac_v+(0.447213595499958*fl[12])/dfac_v+(0.2581988897471612*fr[8])/dfac_v+(0.2581988897471612*fl[8])/dfac_v+(0.6454972243679029*fr[7])/dfac_v+(0.6454972243679029*fl[7])/dfac_v-(0.5*fr[1])/dfac_v+(0.5*fl[1])/dfac_v+(0.2886751345948129*fr[0])/dfac_v+(0.2886751345948129*fl[0])/dfac_v-1.118033988749895*fr[11]*amax+1.118033988749895*fl[11]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[3] = 1.118033988749895*fr[13]*wv+1.118033988749895*fl[13]*wv-0.8660254037844386*fr[5]*wv+0.8660254037844386*fl[5]*wv+0.5*fr[3]*wv+0.5*fl[3]*wv+(0.6454972243679029*fr[17])/dfac_v+(0.6454972243679029*fl[17])/dfac_v-(0.5*fr[10])/dfac_v+(0.5*fl[10])/dfac_v+(0.2886751345948129*fr[6])/dfac_v+(0.2886751345948129*fl[6])/dfac_v-1.118033988749895*fr[13]*amax+1.118033988749895*fl[13]*amax+0.8660254037844386*fr[5]*amax+0.8660254037844386*fl[5]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[6] = 1.118033988749895*fr[17]*wv+1.118033988749895*fl[17]*wv-0.8660254037844386*fr[10]*wv+0.8660254037844386*fl[10]*wv+0.5*fr[6]*wv+0.5*fl[6]*wv-(0.4472135954999579*fr[18])/dfac_v+(0.4472135954999579*fl[18])/dfac_v+(0.2581988897471611*fr[14])/dfac_v+(0.2581988897471611*fl[14])/dfac_v+(0.6454972243679028*fr[13])/dfac_v+(0.6454972243679028*fl[13])/dfac_v-(0.5*fr[5])/dfac_v+(0.5*fl[5])/dfac_v+(0.2886751345948129*fr[3])/dfac_v+(0.2886751345948129*fl[3])/dfac_v-1.118033988749895*fr[17]*amax+1.118033988749895*fl[17]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[6]*amax+0.5*fl[6]*amax; 
  Ghat[8] = (-0.8660254037844387*fr[12]*wv)+0.8660254037844387*fl[12]*wv+0.5*fr[8]*wv+0.5*fl[8]*wv+(0.5773502691896257*fr[11])/dfac_v+(0.5773502691896257*fl[11])/dfac_v-(0.4472135954999579*fr[4])/dfac_v+(0.4472135954999579*fl[4])/dfac_v+(0.2581988897471612*fr[2])/dfac_v+(0.2581988897471612*fl[2])/dfac_v+0.8660254037844387*fr[12]*amax+0.8660254037844387*fl[12]*amax-0.5*fr[8]*amax+0.5*fl[8]*amax; 
  Ghat[9] = (-0.8660254037844387*fr[15]*wv)+0.8660254037844387*fl[15]*wv+0.5*fr[9]*wv+0.5*fl[9]*wv-(0.5*fr[19])/dfac_v+(0.5*fl[19])/dfac_v+(0.2886751345948129*fr[16])/dfac_v+(0.2886751345948129*fl[16])/dfac_v+0.8660254037844387*fr[15]*amax+0.8660254037844387*fl[15]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[14] = (-0.8660254037844387*fr[18]*wv)+0.8660254037844387*fl[18]*wv+0.5*fr[14]*wv+0.5*fl[14]*wv+(0.5773502691896257*fr[17])/dfac_v+(0.5773502691896257*fl[17])/dfac_v-(0.447213595499958*fr[10])/dfac_v+(0.447213595499958*fl[10])/dfac_v+(0.2581988897471611*fr[6])/dfac_v+(0.2581988897471611*fl[6])/dfac_v+0.8660254037844387*fr[18]*amax+0.8660254037844387*fl[18]*amax-0.5*fr[14]*amax+0.5*fl[14]*amax; 
  Ghat[16] = (-0.8660254037844387*fr[19]*wv)+0.8660254037844387*fl[19]*wv+0.5*fr[16]*wv+0.5*fl[16]*wv-(0.5*fr[15])/dfac_v+(0.5*fl[15])/dfac_v+(0.2886751345948129*fr[9])/dfac_v+(0.2886751345948129*fl[9])/dfac_v+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[16]*amax+0.5*fl[16]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = 0.5*Ghat[3]*dfac_x; 
  incr[4] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[5] = -0.8660254037844386*Ghat[3]*dfac_x; 
  incr[6] = 0.5*Ghat[6]*dfac_x; 
  incr[7] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[8] = 0.5*Ghat[8]*dfac_x; 
  incr[9] = 0.5*Ghat[9]*dfac_x; 
  incr[10] = -0.8660254037844386*Ghat[6]*dfac_x; 
  incr[11] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[12] = -0.8660254037844387*Ghat[8]*dfac_x; 
  incr[13] = 1.118033988749895*Ghat[3]*dfac_x; 
  incr[14] = 0.5*Ghat[14]*dfac_x; 
  incr[15] = -0.8660254037844387*Ghat[9]*dfac_x; 
  incr[16] = 0.5*Ghat[16]*dfac_x; 
  incr[17] = 1.118033988749895*Ghat[6]*dfac_x; 
  incr[18] = -0.8660254037844387*Ghat[14]*dfac_x; 
  incr[19] = -0.8660254037844387*Ghat[16]*dfac_x; 

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
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cflVal); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cflVal); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cflVal); 
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
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cflVal); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cflVal); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cflVal); 
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
double GyrokineticSurf1x2vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.732050807568877*Phi[1]*dfac_x*q_)/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  Ghat[0] = (-(3.061862178478972*Phi[2]*fr[12]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[8]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[8]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-1.118033988749895*fr[8]*amax+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Phi[1]*fr[12]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[12]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[11]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[11]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[8]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[8]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[7]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[4]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-1.118033988749895*fr[12]*amax+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[3] = (-(3.061862178478972*Phi[2]*fr[18]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[14]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[14]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[10]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[3]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[3]*dfac_x*q_)/m_-1.118033988749895*fr[14]*amax+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[5] = (-(1.369306393762915*Phi[1]*fr[18]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[18]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[17]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[17]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[14]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[14]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[13]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[13]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[10]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[10]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[3]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[3]*dfac_x*q_)/m_-1.118033988749895*fr[18]*amax+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 
  Ghat[7] = (-(2.738612787525831*Phi[2]*fr[12]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[12]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[11]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[11]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[7]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[4]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_+0.8660254037844387*fr[11]*amax+0.8660254037844387*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax; 
  Ghat[9] = (2.371708245126284*Phi[2]*fr[19]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[19]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[16]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[16]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[9]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[9]*dfac_x*q_)/m_+0.8660254037844387*fr[16]*amax+0.8660254037844387*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[13] = (-(2.738612787525831*Phi[2]*fr[18]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[18]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[17]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[17]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[13]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[13]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[10]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[5]*dfac_x*q_)/m_+0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax; 
  Ghat[15] = (1.060660171779821*Phi[1]*fr[19]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[19]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[16]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[16]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[9]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[9]*dfac_x*q_)/m_+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = 0.5*Ghat[3]*dfac_v; 
  incr[4] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[5] = 0.5*Ghat[5]*dfac_v; 
  incr[6] = -0.8660254037844386*Ghat[3]*dfac_v; 
  incr[7] = 0.5*Ghat[7]*dfac_v; 
  incr[8] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[9] = 0.5*Ghat[9]*dfac_v; 
  incr[10] = -0.8660254037844386*Ghat[5]*dfac_v; 
  incr[11] = -0.8660254037844387*Ghat[7]*dfac_v; 
  incr[12] = 1.118033988749895*Ghat[1]*dfac_v; 
  incr[13] = 0.5*Ghat[13]*dfac_v; 
  incr[14] = 1.118033988749895*Ghat[3]*dfac_v; 
  incr[15] = 0.5*Ghat[15]*dfac_v; 
  incr[16] = -0.8660254037844387*Ghat[9]*dfac_v; 
  incr[17] = -0.8660254037844387*Ghat[13]*dfac_v; 
  incr[18] = 1.118033988749895*Ghat[5]*dfac_v; 
  incr[19] = -0.8660254037844387*Ghat[15]*dfac_v; 

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
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cflVal); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cflVal); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cflVal); 
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
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cflVal); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cflVal); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cflVal); 
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
double EmGyrokineticSurf1x2vSer_X_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  Ghat[0] = 1.118033988749895*fr[7]*wv+1.118033988749895*fl[7]*wv-0.8660254037844386*fr[1]*wv+0.8660254037844386*fl[1]*wv+0.5*fr[0]*wv+0.5*fl[0]*wv+(0.6454972243679028*fr[11])/dfac_v+(0.6454972243679028*fl[11])/dfac_v-(0.5*fr[4])/dfac_v+(0.5*fl[4])/dfac_v+(0.2886751345948129*fr[2])/dfac_v+(0.2886751345948129*fl[2])/dfac_v-1.118033988749895*fr[7]*amax+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = 1.118033988749895*fr[11]*wv+1.118033988749895*fl[11]*wv-0.8660254037844386*fr[4]*wv+0.8660254037844386*fl[4]*wv+0.5*fr[2]*wv+0.5*fl[2]*wv-(0.447213595499958*fr[12])/dfac_v+(0.447213595499958*fl[12])/dfac_v+(0.2581988897471612*fr[8])/dfac_v+(0.2581988897471612*fl[8])/dfac_v+(0.6454972243679029*fr[7])/dfac_v+(0.6454972243679029*fl[7])/dfac_v-(0.5*fr[1])/dfac_v+(0.5*fl[1])/dfac_v+(0.2886751345948129*fr[0])/dfac_v+(0.2886751345948129*fl[0])/dfac_v-1.118033988749895*fr[11]*amax+1.118033988749895*fl[11]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[3] = 1.118033988749895*fr[13]*wv+1.118033988749895*fl[13]*wv-0.8660254037844386*fr[5]*wv+0.8660254037844386*fl[5]*wv+0.5*fr[3]*wv+0.5*fl[3]*wv+(0.6454972243679029*fr[17])/dfac_v+(0.6454972243679029*fl[17])/dfac_v-(0.5*fr[10])/dfac_v+(0.5*fl[10])/dfac_v+(0.2886751345948129*fr[6])/dfac_v+(0.2886751345948129*fl[6])/dfac_v-1.118033988749895*fr[13]*amax+1.118033988749895*fl[13]*amax+0.8660254037844386*fr[5]*amax+0.8660254037844386*fl[5]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[6] = 1.118033988749895*fr[17]*wv+1.118033988749895*fl[17]*wv-0.8660254037844386*fr[10]*wv+0.8660254037844386*fl[10]*wv+0.5*fr[6]*wv+0.5*fl[6]*wv-(0.4472135954999579*fr[18])/dfac_v+(0.4472135954999579*fl[18])/dfac_v+(0.2581988897471611*fr[14])/dfac_v+(0.2581988897471611*fl[14])/dfac_v+(0.6454972243679028*fr[13])/dfac_v+(0.6454972243679028*fl[13])/dfac_v-(0.5*fr[5])/dfac_v+(0.5*fl[5])/dfac_v+(0.2886751345948129*fr[3])/dfac_v+(0.2886751345948129*fl[3])/dfac_v-1.118033988749895*fr[17]*amax+1.118033988749895*fl[17]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[6]*amax+0.5*fl[6]*amax; 
  Ghat[8] = (-0.8660254037844387*fr[12]*wv)+0.8660254037844387*fl[12]*wv+0.5*fr[8]*wv+0.5*fl[8]*wv+(0.5773502691896257*fr[11])/dfac_v+(0.5773502691896257*fl[11])/dfac_v-(0.4472135954999579*fr[4])/dfac_v+(0.4472135954999579*fl[4])/dfac_v+(0.2581988897471612*fr[2])/dfac_v+(0.2581988897471612*fl[2])/dfac_v+0.8660254037844387*fr[12]*amax+0.8660254037844387*fl[12]*amax-0.5*fr[8]*amax+0.5*fl[8]*amax; 
  Ghat[9] = (-0.8660254037844387*fr[15]*wv)+0.8660254037844387*fl[15]*wv+0.5*fr[9]*wv+0.5*fl[9]*wv-(0.5*fr[19])/dfac_v+(0.5*fl[19])/dfac_v+(0.2886751345948129*fr[16])/dfac_v+(0.2886751345948129*fl[16])/dfac_v+0.8660254037844387*fr[15]*amax+0.8660254037844387*fl[15]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[14] = (-0.8660254037844387*fr[18]*wv)+0.8660254037844387*fl[18]*wv+0.5*fr[14]*wv+0.5*fl[14]*wv+(0.5773502691896257*fr[17])/dfac_v+(0.5773502691896257*fl[17])/dfac_v-(0.447213595499958*fr[10])/dfac_v+(0.447213595499958*fl[10])/dfac_v+(0.2581988897471611*fr[6])/dfac_v+(0.2581988897471611*fl[6])/dfac_v+0.8660254037844387*fr[18]*amax+0.8660254037844387*fl[18]*amax-0.5*fr[14]*amax+0.5*fl[14]*amax; 
  Ghat[16] = (-0.8660254037844387*fr[19]*wv)+0.8660254037844387*fl[19]*wv+0.5*fr[16]*wv+0.5*fl[16]*wv-(0.5*fr[15])/dfac_v+(0.5*fl[15])/dfac_v+(0.2886751345948129*fr[9])/dfac_v+(0.2886751345948129*fl[9])/dfac_v+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[16]*amax+0.5*fl[16]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = 0.5*Ghat[3]*dfac_x; 
  incr[4] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[5] = -0.8660254037844386*Ghat[3]*dfac_x; 
  incr[6] = 0.5*Ghat[6]*dfac_x; 
  incr[7] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[8] = 0.5*Ghat[8]*dfac_x; 
  incr[9] = 0.5*Ghat[9]*dfac_x; 
  incr[10] = -0.8660254037844386*Ghat[6]*dfac_x; 
  incr[11] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[12] = -0.8660254037844387*Ghat[8]*dfac_x; 
  incr[13] = 1.118033988749895*Ghat[3]*dfac_x; 
  incr[14] = 0.5*Ghat[14]*dfac_x; 
  incr[15] = -0.8660254037844387*Ghat[9]*dfac_x; 
  incr[16] = 0.5*Ghat[16]*dfac_x; 
  incr[17] = 1.118033988749895*Ghat[6]*dfac_x; 
  incr[18] = -0.8660254037844387*Ghat[14]*dfac_x; 
  incr[19] = -0.8660254037844387*Ghat[16]*dfac_x; 

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
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.0*(1.732050807568877*Phi[1]*dfac_x+dApardt[0])*q_)/m_; 

  double alpha[8]; 
  alpha[0] = (-(3.464101615137754*Phi[1]*dfac_x*q_)/m_)-(2.0*dApardt[0]*q_)/m_; 
  alpha[1] = -(2.0*dApardt[1]*q_)/m_; 
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
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cflVal); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cflVal); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cflVal); 
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
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cflVal); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cflVal); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cflVal); 
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
double EmGyrokineticSurf1x2vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.0*(1.732050807568877*Phi[1]*dfac_x+dApardt[0])*q_)/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  Ghat[0] = (-(3.061862178478972*Phi[2]*fr[12]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[8]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[8]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-(0.7905694150420948*dApardt[1]*fr[12]*q_)/m_-(0.7905694150420948*dApardt[1]*fl[12]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[11]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[0]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[0]*fl[8]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[7]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[7]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[4]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[2]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[0]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[0]*q_)/m_-1.118033988749895*fr[8]*amax+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Phi[1]*fr[12]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[12]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[11]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[11]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[8]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[8]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[7]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[4]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[2]*fr[12]*q_)/m_-(0.7905694150420948*dApardt[0]*fr[12]*q_)/m_-(0.7071067811865475*dApardt[2]*fl[12]*q_)/m_-(0.7905694150420948*dApardt[0]*fl[12]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[11]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[1]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[1]*fl[8]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[7]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[7]*q_)/m_+(0.5477225575051661*dApardt[2]*fr[4]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[4]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[4]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[2]*q_)/m_-(0.3162277660168379*fr[1]*dApardt[2]*q_)/m_-(0.3162277660168379*fl[1]*dApardt[2]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[1]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[1]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[1]*q_)/m_-1.118033988749895*fr[12]*amax+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[3] = (-(3.061862178478972*Phi[2]*fr[18]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[14]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[14]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[10]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[3]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[3]*dfac_x*q_)/m_-(0.7905694150420947*dApardt[1]*fr[18]*q_)/m_-(0.7905694150420947*dApardt[1]*fl[18]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[17]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[17]*q_)/m_-(0.7905694150420948*dApardt[0]*fr[14]*q_)/m_-(0.7905694150420948*dApardt[0]*fl[14]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[13]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[13]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[10]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[6]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[5]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[3]*q_)/m_-1.118033988749895*fr[14]*amax+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[5] = (-(1.369306393762915*Phi[1]*fr[18]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[18]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[17]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[17]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[14]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[14]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[13]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[13]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[10]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[10]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[3]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[3]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[2]*fr[18]*q_)/m_-(0.7905694150420947*dApardt[0]*fr[18]*q_)/m_-(0.7071067811865475*dApardt[2]*fl[18]*q_)/m_-(0.7905694150420947*dApardt[0]*fl[18]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[17]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[17]*q_)/m_-(0.7905694150420948*dApardt[1]*fr[14]*q_)/m_-(0.7905694150420948*dApardt[1]*fl[14]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[13]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[13]*q_)/m_+(0.5477225575051661*dApardt[2]*fr[10]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[10]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[10]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[6]*q_)/m_-(0.3162277660168379*dApardt[2]*fr[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[5]*q_)/m_-(0.3162277660168379*dApardt[2]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[3]*q_)/m_-1.118033988749895*fr[18]*amax+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 
  Ghat[7] = (-(2.738612787525831*Phi[2]*fr[12]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[12]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[11]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[11]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[7]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[4]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[1]*fr[12]*q_)/m_-(0.7071067811865475*dApardt[1]*fl[12]*q_)/m_+(0.3912303982179757*dApardt[2]*fr[11]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[11]*q_)/m_-(0.3912303982179757*dApardt[2]*fl[11]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[2]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[2]*fl[8]*q_)/m_-(0.2258769757263128*dApardt[2]*fr[7]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[7]*q_)/m_-(0.2258769757263128*dApardt[2]*fl[7]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[7]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[4]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[2]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[2]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[2]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[1]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[1]*q_)/m_+0.8660254037844387*fr[11]*amax+0.8660254037844387*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax; 
  Ghat[9] = (2.371708245126284*Phi[2]*fr[19]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[19]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[16]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[16]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[9]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[9]*dfac_x*q_)/m_+(0.6123724356957944*dApardt[1]*fr[19]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[19]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[16]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[16]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[15]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[9]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[9]*q_)/m_+0.8660254037844387*fr[16]*amax+0.8660254037844387*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[13] = (-(2.738612787525831*Phi[2]*fr[18]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[18]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[17]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[17]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[13]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[13]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[10]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[1]*fr[18]*q_)/m_-(0.7071067811865475*dApardt[1]*fl[18]*q_)/m_+(0.3912303982179757*dApardt[2]*fr[17]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[17]*q_)/m_-(0.3912303982179757*dApardt[2]*fl[17]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[17]*q_)/m_-(0.7905694150420947*dApardt[2]*fr[14]*q_)/m_-(0.7905694150420947*dApardt[2]*fl[14]*q_)/m_-(0.2258769757263128*dApardt[2]*fr[13]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[13]*q_)/m_-(0.2258769757263128*dApardt[2]*fl[13]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[13]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[10]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[6]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[5]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[3]*q_)/m_+0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax; 
  Ghat[15] = (1.060660171779821*Phi[1]*fr[19]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[19]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[16]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[16]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[9]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[9]*dfac_x*q_)/m_+(0.5477225575051661*dApardt[2]*fr[19]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[19]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[19]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[19]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[16]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[16]*q_)/m_-(0.3162277660168379*dApardt[2]*fr[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[15]*q_)/m_-(0.3162277660168379*dApardt[2]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[9]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[9]*q_)/m_+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = 0.5*Ghat[3]*dfac_v; 
  incr[4] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[5] = 0.5*Ghat[5]*dfac_v; 
  incr[6] = -0.8660254037844386*Ghat[3]*dfac_v; 
  incr[7] = 0.5*Ghat[7]*dfac_v; 
  incr[8] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[9] = 0.5*Ghat[9]*dfac_v; 
  incr[10] = -0.8660254037844386*Ghat[5]*dfac_v; 
  incr[11] = -0.8660254037844387*Ghat[7]*dfac_v; 
  incr[12] = 1.118033988749895*Ghat[1]*dfac_v; 
  incr[13] = 0.5*Ghat[13]*dfac_v; 
  incr[14] = 1.118033988749895*Ghat[3]*dfac_v; 
  incr[15] = 0.5*Ghat[15]*dfac_v; 
  incr[16] = -0.8660254037844387*Ghat[9]*dfac_v; 
  incr[17] = -0.8660254037844387*Ghat[13]*dfac_v; 
  incr[18] = 1.118033988749895*Ghat[5]*dfac_v; 
  incr[19] = -0.8660254037844387*Ghat[15]*dfac_v; 

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
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cflVal); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cflVal); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cflVal); 
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
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cflVal); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cflVal); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cflVal); 
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
double GyrokineticSurf1x2vSer_X_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  Ghat[0] = 1.118033988749895*fr[7]*wv+1.118033988749895*fl[7]*wv-0.8660254037844386*fr[1]*wv+0.8660254037844386*fl[1]*wv+0.5*fr[0]*wv+0.5*fl[0]*wv+(0.6454972243679028*fr[11])/dfac_v+(0.6454972243679028*fl[11])/dfac_v-(0.5*fr[4])/dfac_v+(0.5*fl[4])/dfac_v+(0.2886751345948129*fr[2])/dfac_v+(0.2886751345948129*fl[2])/dfac_v-1.118033988749895*fr[7]*amax+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = 1.118033988749895*fr[11]*wv+1.118033988749895*fl[11]*wv-0.8660254037844386*fr[4]*wv+0.8660254037844386*fl[4]*wv+0.5*fr[2]*wv+0.5*fl[2]*wv-(0.447213595499958*fr[12])/dfac_v+(0.447213595499958*fl[12])/dfac_v+(0.2581988897471612*fr[8])/dfac_v+(0.2581988897471612*fl[8])/dfac_v+(0.6454972243679029*fr[7])/dfac_v+(0.6454972243679029*fl[7])/dfac_v-(0.5*fr[1])/dfac_v+(0.5*fl[1])/dfac_v+(0.2886751345948129*fr[0])/dfac_v+(0.2886751345948129*fl[0])/dfac_v-1.118033988749895*fr[11]*amax+1.118033988749895*fl[11]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[3] = 1.118033988749895*fr[13]*wv+1.118033988749895*fl[13]*wv-0.8660254037844386*fr[5]*wv+0.8660254037844386*fl[5]*wv+0.5*fr[3]*wv+0.5*fl[3]*wv+(0.6454972243679029*fr[17])/dfac_v+(0.6454972243679029*fl[17])/dfac_v-(0.5*fr[10])/dfac_v+(0.5*fl[10])/dfac_v+(0.2886751345948129*fr[6])/dfac_v+(0.2886751345948129*fl[6])/dfac_v-1.118033988749895*fr[13]*amax+1.118033988749895*fl[13]*amax+0.8660254037844386*fr[5]*amax+0.8660254037844386*fl[5]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[6] = 1.118033988749895*fr[17]*wv+1.118033988749895*fl[17]*wv-0.8660254037844386*fr[10]*wv+0.8660254037844386*fl[10]*wv+0.5*fr[6]*wv+0.5*fl[6]*wv-(0.4472135954999579*fr[18])/dfac_v+(0.4472135954999579*fl[18])/dfac_v+(0.2581988897471611*fr[14])/dfac_v+(0.2581988897471611*fl[14])/dfac_v+(0.6454972243679028*fr[13])/dfac_v+(0.6454972243679028*fl[13])/dfac_v-(0.5*fr[5])/dfac_v+(0.5*fl[5])/dfac_v+(0.2886751345948129*fr[3])/dfac_v+(0.2886751345948129*fl[3])/dfac_v-1.118033988749895*fr[17]*amax+1.118033988749895*fl[17]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[6]*amax+0.5*fl[6]*amax; 
  Ghat[8] = (-0.8660254037844387*fr[12]*wv)+0.8660254037844387*fl[12]*wv+0.5*fr[8]*wv+0.5*fl[8]*wv+(0.5773502691896257*fr[11])/dfac_v+(0.5773502691896257*fl[11])/dfac_v-(0.4472135954999579*fr[4])/dfac_v+(0.4472135954999579*fl[4])/dfac_v+(0.2581988897471612*fr[2])/dfac_v+(0.2581988897471612*fl[2])/dfac_v+0.8660254037844387*fr[12]*amax+0.8660254037844387*fl[12]*amax-0.5*fr[8]*amax+0.5*fl[8]*amax; 
  Ghat[9] = (-0.8660254037844387*fr[15]*wv)+0.8660254037844387*fl[15]*wv+0.5*fr[9]*wv+0.5*fl[9]*wv-(0.5*fr[19])/dfac_v+(0.5*fl[19])/dfac_v+(0.2886751345948129*fr[16])/dfac_v+(0.2886751345948129*fl[16])/dfac_v+0.8660254037844387*fr[15]*amax+0.8660254037844387*fl[15]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[14] = (-0.8660254037844387*fr[18]*wv)+0.8660254037844387*fl[18]*wv+0.5*fr[14]*wv+0.5*fl[14]*wv+(0.5773502691896257*fr[17])/dfac_v+(0.5773502691896257*fl[17])/dfac_v-(0.447213595499958*fr[10])/dfac_v+(0.447213595499958*fl[10])/dfac_v+(0.2581988897471611*fr[6])/dfac_v+(0.2581988897471611*fl[6])/dfac_v+0.8660254037844387*fr[18]*amax+0.8660254037844387*fl[18]*amax-0.5*fr[14]*amax+0.5*fl[14]*amax; 
  Ghat[16] = (-0.8660254037844387*fr[19]*wv)+0.8660254037844387*fl[19]*wv+0.5*fr[16]*wv+0.5*fl[16]*wv-(0.5*fr[15])/dfac_v+(0.5*fl[15])/dfac_v+(0.2886751345948129*fr[9])/dfac_v+(0.2886751345948129*fl[9])/dfac_v+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[16]*amax+0.5*fl[16]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = 0.5*Ghat[3]*dfac_x; 
  incr[4] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[5] = -0.8660254037844386*Ghat[3]*dfac_x; 
  incr[6] = 0.5*Ghat[6]*dfac_x; 
  incr[7] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[8] = 0.5*Ghat[8]*dfac_x; 
  incr[9] = 0.5*Ghat[9]*dfac_x; 
  incr[10] = -0.8660254037844386*Ghat[6]*dfac_x; 
  incr[11] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[12] = -0.8660254037844387*Ghat[8]*dfac_x; 
  incr[13] = 1.118033988749895*Ghat[3]*dfac_x; 
  incr[14] = 0.5*Ghat[14]*dfac_x; 
  incr[15] = -0.8660254037844387*Ghat[9]*dfac_x; 
  incr[16] = 0.5*Ghat[16]*dfac_x; 
  incr[17] = 1.118033988749895*Ghat[6]*dfac_x; 
  incr[18] = -0.8660254037844387*Ghat[14]*dfac_x; 
  incr[19] = -0.8660254037844387*Ghat[16]*dfac_x; 

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
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cflVal); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cflVal); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cflVal); 
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
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cflVal); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cflVal); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cflVal); 
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
double GyrokineticSurf1x2vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.0*(1.732050807568877*Bmag[1]*dfac_x*wm+1.732050807568877*Phi[1]*dfac_x*q_))/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  Ghat[0] = (-(3.061862178478972*Bmag[2]*fr[12]*dfac_x*wm)/m_)-(3.061862178478972*Bmag[2]*fl[12]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fr[8]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fl[8]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[4]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[4]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[2]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[2]*dfac_x*wm)/m_-(1.369306393762915*fr[1]*Bmag[2]*dfac_x*wm)/m_-(1.369306393762915*fl[1]*Bmag[2]*dfac_x*wm)/m_-(0.6123724356957944*fr[0]*Bmag[1]*dfac_x*wm)/m_-(0.6123724356957944*fl[0]*Bmag[1]*dfac_x*wm)/m_-(3.061862178478972*Phi[2]*fr[12]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[8]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[8]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-(1.767766952966368*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fr[14]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fl[14]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[6]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[6]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[3]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[3]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[8]*amax+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Bmag[1]*fr[12]*dfac_x*wm)/m_)-(1.369306393762915*Bmag[1]*fl[12]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[11]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[11]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fr[8]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fl[8]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[7]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[7]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[4]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[4]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[2]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[2]*dfac_x*wm)/m_-(1.369306393762915*fr[0]*Bmag[2]*dfac_x*wm)/m_-(1.369306393762915*fl[0]*Bmag[2]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[1]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[1]*dfac_x*wm)/m_-(1.369306393762915*Phi[1]*fr[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[12]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[11]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[11]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[8]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[8]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[7]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[4]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-(0.7905694150420947*Bmag[1]*fr[18]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fl[18]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[17]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[17]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fr[14]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fl[14]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[13]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[13]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[10]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[10]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[6]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[6]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[5]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[5]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fr[3]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fl[3]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[12]*amax+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[3] = (-(3.061862178478972*Bmag[2]*fr[18]*dfac_x*wm)/m_)-(3.061862178478972*Bmag[2]*fl[18]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fr[14]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fl[14]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[10]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[10]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[6]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[6]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[5]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[5]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[3]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[3]*dfac_x*wm)/m_-(3.061862178478972*Phi[2]*fr[18]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[14]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[14]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[10]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[3]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[3]*dfac_x*q_)/m_+(1.224744871391589*Bmag[2]*fr[19]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[19]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[16]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[16]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[15]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[15]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fr[12]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fl[12]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[9]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[9]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fr[8]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fl[8]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[4]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[4]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[2]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fr[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fl[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*fr[0]*Bmag[1]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*fl[0]*Bmag[1]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[14]*amax+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[5] = (-(1.369306393762915*Bmag[1]*fr[18]*dfac_x*wm)/m_)-(1.369306393762915*Bmag[1]*fl[18]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[17]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[17]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fr[14]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fl[14]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[13]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[13]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[10]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[10]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[6]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[6]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[5]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[5]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[3]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[3]*dfac_x*wm)/m_-(1.369306393762915*Phi[1]*fr[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[18]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[17]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[17]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[14]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[14]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[13]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[13]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[10]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[10]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[3]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[3]*dfac_x*q_)/m_+(0.5477225575051661*Bmag[1]*fr[19]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[19]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[16]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[16]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[15]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[15]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fr[12]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fl[12]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[11]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[11]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[9]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[9]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fr[8]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fl[8]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[7]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[7]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[4]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[4]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[2]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fr[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fl[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[1]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[1]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[18]*amax+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 
  Ghat[7] = (-(2.738612787525831*Bmag[2]*fr[12]*dfac_x*wm)/m_)-(2.738612787525831*Bmag[2]*fl[12]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[11]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[11]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[7]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[7]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[4]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[4]*dfac_x*wm)/m_-(1.224744871391589*fr[1]*Bmag[2]*dfac_x*wm)/m_-(1.224744871391589*fl[1]*Bmag[2]*dfac_x*wm)/m_-(2.738612787525831*Phi[2]*fr[12]*dfac_x*q_)/m_-(2.738612787525831*Phi[2]*fl[12]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[11]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[11]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[7]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[4]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_-(1.58113883008419*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[17]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[17]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[13]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[13]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[11]*amax+0.8660254037844387*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax; 
  Ghat[9] = (2.371708245126284*Bmag[2]*fr[19]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[19]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[16]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[16]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[15]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[15]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[9]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[9]*dfac_x*wm)/m_+(2.371708245126284*Phi[2]*fr[19]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[19]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[16]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[16]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[9]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[9]*dfac_x*q_)/m_-(1.58113883008419*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fr[14]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fl[14]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[6]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[6]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[3]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[3]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[16]*amax+0.8660254037844387*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[13] = (-(2.738612787525831*Bmag[2]*fr[18]*dfac_x*wm)/m_)-(2.738612787525831*Bmag[2]*fl[18]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[17]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[17]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[13]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[13]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[10]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[10]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[5]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[5]*dfac_x*wm)/m_-(2.738612787525831*Phi[2]*fr[18]*dfac_x*q_)/m_-(2.738612787525831*Phi[2]*fl[18]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[17]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[17]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[13]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[13]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[10]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[5]*dfac_x*q_)/m_+(1.095445115010332*Bmag[2]*fr[19]*dfac_x)/(dfac_m*m_)-(1.095445115010332*Bmag[2]*fl[19]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fr[15]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fl[15]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fr[12]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[12]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[11]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[11]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[7]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[7]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[4]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[4]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*fr[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*fl[1]*Bmag[2]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax; 
  Ghat[15] = (1.060660171779821*Bmag[1]*fr[19]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[19]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[16]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[16]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[15]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[15]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[9]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[9]*dfac_x*wm)/m_+(1.060660171779821*Phi[1]*fr[19]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[19]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[16]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[16]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[9]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[9]*dfac_x*q_)/m_-(0.7071067811865475*Bmag[1]*fr[18]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fl[18]*dfac_x)/(dfac_m*m_)+(1.095445115010332*Bmag[2]*fr[17]*dfac_x)/(dfac_m*m_)-(1.095445115010332*Bmag[2]*fl[17]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fr[14]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[14]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fr[13]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fl[13]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[10]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[10]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[6]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[6]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[5]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[3]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[3]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = 0.5*Ghat[3]*dfac_v; 
  incr[4] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[5] = 0.5*Ghat[5]*dfac_v; 
  incr[6] = -0.8660254037844386*Ghat[3]*dfac_v; 
  incr[7] = 0.5*Ghat[7]*dfac_v; 
  incr[8] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[9] = 0.5*Ghat[9]*dfac_v; 
  incr[10] = -0.8660254037844386*Ghat[5]*dfac_v; 
  incr[11] = -0.8660254037844387*Ghat[7]*dfac_v; 
  incr[12] = 1.118033988749895*Ghat[1]*dfac_v; 
  incr[13] = 0.5*Ghat[13]*dfac_v; 
  incr[14] = 1.118033988749895*Ghat[3]*dfac_v; 
  incr[15] = 0.5*Ghat[15]*dfac_v; 
  incr[16] = -0.8660254037844387*Ghat[9]*dfac_v; 
  incr[17] = -0.8660254037844387*Ghat[13]*dfac_v; 
  incr[18] = 1.118033988749895*Ghat[5]*dfac_v; 
  incr[19] = -0.8660254037844387*Ghat[15]*dfac_v; 

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
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fl[6]-3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*fl[3]-1.0*(3.0*fl[2]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cflVal); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[2]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cflVal); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[6]+3.0*(fl[3]+fl[2])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cflVal); 
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
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.06804138174397717*(1.732050807568877*fr[6]-3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*fr[3]-1.0*(3.0*fr[2]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cflVal); 
  fqVal[5] = -0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[2]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cflVal); 
  fqVal[6] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[6]+3.0*(fr[3]+fr[2])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cflVal); 
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
double EmGyrokineticSurf1x2vSer_X_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  Ghat[0] = 1.118033988749895*fr[7]*wv+1.118033988749895*fl[7]*wv-0.8660254037844386*fr[1]*wv+0.8660254037844386*fl[1]*wv+0.5*fr[0]*wv+0.5*fl[0]*wv+(0.6454972243679028*fr[11])/dfac_v+(0.6454972243679028*fl[11])/dfac_v-(0.5*fr[4])/dfac_v+(0.5*fl[4])/dfac_v+(0.2886751345948129*fr[2])/dfac_v+(0.2886751345948129*fl[2])/dfac_v-1.118033988749895*fr[7]*amax+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = 1.118033988749895*fr[11]*wv+1.118033988749895*fl[11]*wv-0.8660254037844386*fr[4]*wv+0.8660254037844386*fl[4]*wv+0.5*fr[2]*wv+0.5*fl[2]*wv-(0.447213595499958*fr[12])/dfac_v+(0.447213595499958*fl[12])/dfac_v+(0.2581988897471612*fr[8])/dfac_v+(0.2581988897471612*fl[8])/dfac_v+(0.6454972243679029*fr[7])/dfac_v+(0.6454972243679029*fl[7])/dfac_v-(0.5*fr[1])/dfac_v+(0.5*fl[1])/dfac_v+(0.2886751345948129*fr[0])/dfac_v+(0.2886751345948129*fl[0])/dfac_v-1.118033988749895*fr[11]*amax+1.118033988749895*fl[11]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[3] = 1.118033988749895*fr[13]*wv+1.118033988749895*fl[13]*wv-0.8660254037844386*fr[5]*wv+0.8660254037844386*fl[5]*wv+0.5*fr[3]*wv+0.5*fl[3]*wv+(0.6454972243679029*fr[17])/dfac_v+(0.6454972243679029*fl[17])/dfac_v-(0.5*fr[10])/dfac_v+(0.5*fl[10])/dfac_v+(0.2886751345948129*fr[6])/dfac_v+(0.2886751345948129*fl[6])/dfac_v-1.118033988749895*fr[13]*amax+1.118033988749895*fl[13]*amax+0.8660254037844386*fr[5]*amax+0.8660254037844386*fl[5]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[6] = 1.118033988749895*fr[17]*wv+1.118033988749895*fl[17]*wv-0.8660254037844386*fr[10]*wv+0.8660254037844386*fl[10]*wv+0.5*fr[6]*wv+0.5*fl[6]*wv-(0.4472135954999579*fr[18])/dfac_v+(0.4472135954999579*fl[18])/dfac_v+(0.2581988897471611*fr[14])/dfac_v+(0.2581988897471611*fl[14])/dfac_v+(0.6454972243679028*fr[13])/dfac_v+(0.6454972243679028*fl[13])/dfac_v-(0.5*fr[5])/dfac_v+(0.5*fl[5])/dfac_v+(0.2886751345948129*fr[3])/dfac_v+(0.2886751345948129*fl[3])/dfac_v-1.118033988749895*fr[17]*amax+1.118033988749895*fl[17]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[6]*amax+0.5*fl[6]*amax; 
  Ghat[8] = (-0.8660254037844387*fr[12]*wv)+0.8660254037844387*fl[12]*wv+0.5*fr[8]*wv+0.5*fl[8]*wv+(0.5773502691896257*fr[11])/dfac_v+(0.5773502691896257*fl[11])/dfac_v-(0.4472135954999579*fr[4])/dfac_v+(0.4472135954999579*fl[4])/dfac_v+(0.2581988897471612*fr[2])/dfac_v+(0.2581988897471612*fl[2])/dfac_v+0.8660254037844387*fr[12]*amax+0.8660254037844387*fl[12]*amax-0.5*fr[8]*amax+0.5*fl[8]*amax; 
  Ghat[9] = (-0.8660254037844387*fr[15]*wv)+0.8660254037844387*fl[15]*wv+0.5*fr[9]*wv+0.5*fl[9]*wv-(0.5*fr[19])/dfac_v+(0.5*fl[19])/dfac_v+(0.2886751345948129*fr[16])/dfac_v+(0.2886751345948129*fl[16])/dfac_v+0.8660254037844387*fr[15]*amax+0.8660254037844387*fl[15]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[14] = (-0.8660254037844387*fr[18]*wv)+0.8660254037844387*fl[18]*wv+0.5*fr[14]*wv+0.5*fl[14]*wv+(0.5773502691896257*fr[17])/dfac_v+(0.5773502691896257*fl[17])/dfac_v-(0.447213595499958*fr[10])/dfac_v+(0.447213595499958*fl[10])/dfac_v+(0.2581988897471611*fr[6])/dfac_v+(0.2581988897471611*fl[6])/dfac_v+0.8660254037844387*fr[18]*amax+0.8660254037844387*fl[18]*amax-0.5*fr[14]*amax+0.5*fl[14]*amax; 
  Ghat[16] = (-0.8660254037844387*fr[19]*wv)+0.8660254037844387*fl[19]*wv+0.5*fr[16]*wv+0.5*fl[16]*wv-(0.5*fr[15])/dfac_v+(0.5*fl[15])/dfac_v+(0.2886751345948129*fr[9])/dfac_v+(0.2886751345948129*fl[9])/dfac_v+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[16]*amax+0.5*fl[16]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = 0.5*Ghat[3]*dfac_x; 
  incr[4] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[5] = -0.8660254037844386*Ghat[3]*dfac_x; 
  incr[6] = 0.5*Ghat[6]*dfac_x; 
  incr[7] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[8] = 0.5*Ghat[8]*dfac_x; 
  incr[9] = 0.5*Ghat[9]*dfac_x; 
  incr[10] = -0.8660254037844386*Ghat[6]*dfac_x; 
  incr[11] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[12] = -0.8660254037844387*Ghat[8]*dfac_x; 
  incr[13] = 1.118033988749895*Ghat[3]*dfac_x; 
  incr[14] = 0.5*Ghat[14]*dfac_x; 
  incr[15] = -0.8660254037844387*Ghat[9]*dfac_x; 
  incr[16] = 0.5*Ghat[16]*dfac_x; 
  incr[17] = 1.118033988749895*Ghat[6]*dfac_x; 
  incr[18] = -0.8660254037844387*Ghat[14]*dfac_x; 
  incr[19] = -0.8660254037844387*Ghat[16]*dfac_x; 

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
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.0*(1.732050807568877*Bmag[1]*dfac_x*wm+(1.732050807568877*Phi[1]*dfac_x+dApardt[0])*q_))/m_; 

  double alpha[8]; 
  alpha[0] = (-(3.464101615137754*Bmag[1]*dfac_x*wm)/m_)-(3.464101615137754*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[0]*q_)/m_; 
  alpha[1] = -(2.0*dApardt[1]*q_)/m_; 
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
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*fl[3]-1.0*(3.0*fl[1]+5.196152422706631*fl[0]))*limTheta(rVal[3],1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[4],1.0,cflVal); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[5],1.0,cflVal); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[1]-1.0*fl[3])-5.196152422706631*fl[0])*limTheta(rVal[6],1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])*limTheta(rVal[7],1.0,cflVal); 
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
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.06804138174397717*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*fr[3]-1.0*(3.0*fr[1]+5.196152422706631*fr[0]))*limTheta(rVal[3],-1.0,cflVal); 
  fqVal[4] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[4],-1.0,cflVal); 
  fqVal[5] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[5],-1.0,cflVal); 
  fqVal[6] = -0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[1]-1.0*fr[3])-5.196152422706631*fr[0])*limTheta(rVal[6],-1.0,cflVal); 
  fqVal[7] = 0.06804138174397717*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])*limTheta(rVal[7],-1.0,cflVal); 
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
double EmGyrokineticSurf1x2vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.0*(1.732050807568877*Bmag[1]*dfac_x*wm+(1.732050807568877*Phi[1]*dfac_x+dApardt[0])*q_))/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  Ghat[0] = (-(3.061862178478972*Bmag[2]*fr[12]*dfac_x*wm)/m_)-(3.061862178478972*Bmag[2]*fl[12]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fr[8]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fl[8]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[4]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[4]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[2]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[2]*dfac_x*wm)/m_-(1.369306393762915*fr[1]*Bmag[2]*dfac_x*wm)/m_-(1.369306393762915*fl[1]*Bmag[2]*dfac_x*wm)/m_-(0.6123724356957944*fr[0]*Bmag[1]*dfac_x*wm)/m_-(0.6123724356957944*fl[0]*Bmag[1]*dfac_x*wm)/m_-(3.061862178478972*Phi[2]*fr[12]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[8]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[8]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-(0.7905694150420948*dApardt[1]*fr[12]*q_)/m_-(0.7905694150420948*dApardt[1]*fl[12]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[11]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[0]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[0]*fl[8]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[7]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[7]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[4]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[2]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[0]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[0]*q_)/m_-(1.767766952966368*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fr[14]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fl[14]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[6]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[6]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[3]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[3]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[8]*amax+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Bmag[1]*fr[12]*dfac_x*wm)/m_)-(1.369306393762915*Bmag[1]*fl[12]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[11]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[11]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fr[8]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fl[8]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[7]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[7]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[4]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[4]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[2]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[2]*dfac_x*wm)/m_-(1.369306393762915*fr[0]*Bmag[2]*dfac_x*wm)/m_-(1.369306393762915*fl[0]*Bmag[2]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[1]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[1]*dfac_x*wm)/m_-(1.369306393762915*Phi[1]*fr[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[12]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[11]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[11]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[8]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[8]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[7]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[4]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[2]*fr[12]*q_)/m_-(0.7905694150420948*dApardt[0]*fr[12]*q_)/m_-(0.7071067811865475*dApardt[2]*fl[12]*q_)/m_-(0.7905694150420948*dApardt[0]*fl[12]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[11]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[1]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[1]*fl[8]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[7]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[7]*q_)/m_+(0.5477225575051661*dApardt[2]*fr[4]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[4]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[4]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[2]*q_)/m_-(0.3162277660168379*fr[1]*dApardt[2]*q_)/m_-(0.3162277660168379*fl[1]*dApardt[2]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[1]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[1]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[1]*q_)/m_-(0.7905694150420947*Bmag[1]*fr[18]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fl[18]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[17]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[17]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fr[14]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fl[14]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[13]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[13]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[10]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[10]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[6]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[6]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[5]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[5]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fr[3]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fl[3]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[12]*amax+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[3] = (-(3.061862178478972*Bmag[2]*fr[18]*dfac_x*wm)/m_)-(3.061862178478972*Bmag[2]*fl[18]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fr[14]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fl[14]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[10]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[10]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[6]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[6]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[5]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[5]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[3]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[3]*dfac_x*wm)/m_-(3.061862178478972*Phi[2]*fr[18]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[14]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[14]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[10]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[3]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[3]*dfac_x*q_)/m_-(0.7905694150420947*dApardt[1]*fr[18]*q_)/m_-(0.7905694150420947*dApardt[1]*fl[18]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[17]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[17]*q_)/m_-(0.7905694150420948*dApardt[0]*fr[14]*q_)/m_-(0.7905694150420948*dApardt[0]*fl[14]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[13]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[13]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[10]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[6]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[5]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[3]*q_)/m_+(1.224744871391589*Bmag[2]*fr[19]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[19]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[16]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[16]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[15]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[15]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fr[12]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fl[12]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[9]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[9]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fr[8]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fl[8]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[4]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[4]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[2]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fr[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fl[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*fr[0]*Bmag[1]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*fl[0]*Bmag[1]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[14]*amax+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[5] = (-(1.369306393762915*Bmag[1]*fr[18]*dfac_x*wm)/m_)-(1.369306393762915*Bmag[1]*fl[18]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[17]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[17]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fr[14]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fl[14]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[13]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[13]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[10]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[10]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[6]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[6]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[5]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[5]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[3]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[3]*dfac_x*wm)/m_-(1.369306393762915*Phi[1]*fr[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[18]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[17]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[17]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[14]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[14]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[13]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[13]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[10]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[10]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[3]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[3]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[2]*fr[18]*q_)/m_-(0.7905694150420947*dApardt[0]*fr[18]*q_)/m_-(0.7071067811865475*dApardt[2]*fl[18]*q_)/m_-(0.7905694150420947*dApardt[0]*fl[18]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[17]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[17]*q_)/m_-(0.7905694150420948*dApardt[1]*fr[14]*q_)/m_-(0.7905694150420948*dApardt[1]*fl[14]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[13]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[13]*q_)/m_+(0.5477225575051661*dApardt[2]*fr[10]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[10]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[10]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[6]*q_)/m_-(0.3162277660168379*dApardt[2]*fr[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[5]*q_)/m_-(0.3162277660168379*dApardt[2]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[3]*q_)/m_+(0.5477225575051661*Bmag[1]*fr[19]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[19]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[16]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[16]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[15]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[15]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fr[12]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fl[12]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[11]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[11]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[9]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[9]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fr[8]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fl[8]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[7]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[7]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[4]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[4]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[2]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fr[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fl[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[1]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[1]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[18]*amax+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 
  Ghat[7] = (-(2.738612787525831*Bmag[2]*fr[12]*dfac_x*wm)/m_)-(2.738612787525831*Bmag[2]*fl[12]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[11]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[11]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[7]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[7]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[4]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[4]*dfac_x*wm)/m_-(1.224744871391589*fr[1]*Bmag[2]*dfac_x*wm)/m_-(1.224744871391589*fl[1]*Bmag[2]*dfac_x*wm)/m_-(2.738612787525831*Phi[2]*fr[12]*dfac_x*q_)/m_-(2.738612787525831*Phi[2]*fl[12]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[11]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[11]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[7]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[4]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[1]*fr[12]*q_)/m_-(0.7071067811865475*dApardt[1]*fl[12]*q_)/m_+(0.3912303982179757*dApardt[2]*fr[11]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[11]*q_)/m_-(0.3912303982179757*dApardt[2]*fl[11]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[2]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[2]*fl[8]*q_)/m_-(0.2258769757263128*dApardt[2]*fr[7]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[7]*q_)/m_-(0.2258769757263128*dApardt[2]*fl[7]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[7]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[4]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[2]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[2]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[2]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[1]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[1]*q_)/m_-(1.58113883008419*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[17]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[17]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[13]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[13]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[11]*amax+0.8660254037844387*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax; 
  Ghat[9] = (2.371708245126284*Bmag[2]*fr[19]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[19]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[16]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[16]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[15]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[15]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[9]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[9]*dfac_x*wm)/m_+(2.371708245126284*Phi[2]*fr[19]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[19]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[16]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[16]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[9]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[9]*dfac_x*q_)/m_+(0.6123724356957944*dApardt[1]*fr[19]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[19]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[16]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[16]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[15]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[9]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[9]*q_)/m_-(1.58113883008419*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fr[14]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fl[14]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[6]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[6]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[3]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[3]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[16]*amax+0.8660254037844387*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[13] = (-(2.738612787525831*Bmag[2]*fr[18]*dfac_x*wm)/m_)-(2.738612787525831*Bmag[2]*fl[18]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[17]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[17]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[13]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[13]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[10]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[10]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[5]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[5]*dfac_x*wm)/m_-(2.738612787525831*Phi[2]*fr[18]*dfac_x*q_)/m_-(2.738612787525831*Phi[2]*fl[18]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[17]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[17]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[13]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[13]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[10]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[1]*fr[18]*q_)/m_-(0.7071067811865475*dApardt[1]*fl[18]*q_)/m_+(0.3912303982179757*dApardt[2]*fr[17]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[17]*q_)/m_-(0.3912303982179757*dApardt[2]*fl[17]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[17]*q_)/m_-(0.7905694150420947*dApardt[2]*fr[14]*q_)/m_-(0.7905694150420947*dApardt[2]*fl[14]*q_)/m_-(0.2258769757263128*dApardt[2]*fr[13]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[13]*q_)/m_-(0.2258769757263128*dApardt[2]*fl[13]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[13]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[10]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[6]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[5]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[3]*q_)/m_+(1.095445115010332*Bmag[2]*fr[19]*dfac_x)/(dfac_m*m_)-(1.095445115010332*Bmag[2]*fl[19]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fr[15]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fl[15]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fr[12]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[12]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[11]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[11]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[7]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[7]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[4]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[4]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*fr[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*fl[1]*Bmag[2]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax; 
  Ghat[15] = (1.060660171779821*Bmag[1]*fr[19]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[19]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[16]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[16]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[15]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[15]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[9]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[9]*dfac_x*wm)/m_+(1.060660171779821*Phi[1]*fr[19]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[19]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[16]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[16]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[9]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[9]*dfac_x*q_)/m_+(0.5477225575051661*dApardt[2]*fr[19]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[19]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[19]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[19]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[16]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[16]*q_)/m_-(0.3162277660168379*dApardt[2]*fr[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[15]*q_)/m_-(0.3162277660168379*dApardt[2]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[9]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[9]*q_)/m_-(0.7071067811865475*Bmag[1]*fr[18]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fl[18]*dfac_x)/(dfac_m*m_)+(1.095445115010332*Bmag[2]*fr[17]*dfac_x)/(dfac_m*m_)-(1.095445115010332*Bmag[2]*fl[17]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fr[14]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[14]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fr[13]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fl[13]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[10]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[10]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[6]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[6]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[5]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[3]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[3]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = 0.5*Ghat[3]*dfac_v; 
  incr[4] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[5] = 0.5*Ghat[5]*dfac_v; 
  incr[6] = -0.8660254037844386*Ghat[3]*dfac_v; 
  incr[7] = 0.5*Ghat[7]*dfac_v; 
  incr[8] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[9] = 0.5*Ghat[9]*dfac_v; 
  incr[10] = -0.8660254037844386*Ghat[5]*dfac_v; 
  incr[11] = -0.8660254037844387*Ghat[7]*dfac_v; 
  incr[12] = 1.118033988749895*Ghat[1]*dfac_v; 
  incr[13] = 0.5*Ghat[13]*dfac_v; 
  incr[14] = 1.118033988749895*Ghat[3]*dfac_v; 
  incr[15] = 0.5*Ghat[15]*dfac_v; 
  incr[16] = -0.8660254037844387*Ghat[9]*dfac_v; 
  incr[17] = -0.8660254037844387*Ghat[13]*dfac_v; 
  incr[18] = 1.118033988749895*Ghat[5]*dfac_v; 
  incr[19] = -0.8660254037844387*Ghat[15]*dfac_v; 

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
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
