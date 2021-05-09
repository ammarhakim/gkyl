#include <VmLBOModDecl.h> 
double VmLBOconstNuSurfNonUniform1x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[4]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[4]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 

  double fjump[4]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 

  double Ghat[4]; 
  for(unsigned short int i=0; i<4; ++i){ 
    Ghat[i]=0.0; 
  }; 

  const double dxvl1R2 = pow(dxvl[1],2);
  const double dxvl1R3 = pow(dxvl[1],3);
  const double dxvl1R4 = pow(dxvl[1],4);
  const double dxvr1R2 = pow(dxvr[1],2);
  const double dxvr1R3 = pow(dxvr[1],3);
  const double dxvr1R4 = pow(dxvr[1],4);

  Ghat[0] = (-(1.0*((17.32050807568877*nuVtSqSum[0]*dxvl1R2*dxvr1R2+3.464101615137754*nuVtSqSum[0]*dxvl1R3*dxvr[1]-3.464101615137754*nuVtSqSum[0]*dxvl1R4)*fr[2]+((-3.464101615137754*nuVtSqSum[0]*dxvr1R4)+3.464101615137754*nuVtSqSum[0]*dxvl[1]*dxvr1R3+17.32050807568877*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[2]+(18.0*dxvl1R2*dxvr1R2*fl[1]-18.0*dxvl1R2*dxvr1R2*fr[1])*nuVtSqSum[1]+(18.0*fl[0]-18.0*fr[0])*nuVtSqSum[0]*dxvl1R2*dxvr1R2))/(1.414213562373095*dxvl[1]*dxvr1R4+4.242640687119286*dxvl1R2*dxvr1R3+4.242640687119286*dxvl1R3*dxvr1R2+1.414213562373095*dxvl1R4*dxvr[1]))-0.5*(1.732050807568877*fjump[2]+fjump[0])+0.3535533905932737*(1.732050807568877*alphaDrag[0]*favg[2]+alphaDrag[1]*favg[1]+alphaDrag[0]*favg[0]); 
  Ghat[1] = (-(1.0*((12.24744871391589*dxvl1R2*dxvr1R2+2.449489742783178*dxvl1R3*dxvr[1]-2.449489742783178*dxvl1R4)*nuVtSqSum[1]*fr[2]+((-2.449489742783178*dxvr1R4)+2.449489742783178*dxvl[1]*dxvr1R3+12.24744871391589*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[2]+(12.72792206135786*fl[0]-12.72792206135786*fr[0])*dxvl1R2*dxvr1R2*nuVtSqSum[1]-12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[1]+12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[1]))/(dxvl[1]*dxvr1R4+3.0*dxvl1R2*dxvr1R3+3.0*dxvl1R3*dxvr1R2+dxvl1R4*dxvr[1]))+0.3535533905932737*(1.732050807568877*alphaDrag[1]*favg[2]+alphaDrag[0]*favg[1]+favg[0]*alphaDrag[1])-0.5*fjump[1]; 
  Ghat[3] = (18.0*nuVtSqSum[0]*dxvl[1]*dxvr[1]*fr[3]-18.0*nuVtSqSum[0]*dxvl[1]*dxvr[1]*fl[3])/(1.414213562373095*dxvr1R3+4.242640687119286*dxvl[1]*dxvr1R2+4.242640687119286*dxvl1R2*dxvr[1]+1.414213562373095*dxvl1R3)-0.5*fjump[3]+0.3535533905932737*alphaDrag[0]*favg[3]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 

  double incr2[4]; 

  incr2[2] = -(1.0*((5.0*nuVtSqSum[0]*dxvl1R2*dxvr[1]+3.0*nuVtSqSum[0]*dxvl1R3)*fr[2]+((-3.0*nuVtSqSum[0]*dxvr1R3)-5.0*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[2]+(((-5.196152422706631*dxvl1R2*dxvr[1])-1.732050807568877*dxvl1R3)*fr[1]+((-1.732050807568877*dxvr1R3)-5.196152422706631*dxvl[1]*dxvr1R2)*fl[1])*nuVtSqSum[1]-1.732050807568877*fl[0]*nuVtSqSum[0]*dxvr1R3-5.196152422706631*fl[0]*nuVtSqSum[0]*dxvl[1]*dxvr1R2-5.196152422706631*fr[0]*nuVtSqSum[0]*dxvl1R2*dxvr[1]-1.732050807568877*fr[0]*nuVtSqSum[0]*dxvl1R3))/(2.828427124746191*dxvr1R3+8.485281374238571*dxvl[1]*dxvr1R2+8.485281374238571*dxvl1R2*dxvr[1]+2.828427124746191*dxvl1R3); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 

  return std::abs(wl[1]-(0.7071067811865475*sumNuUx[0])/nuSum); 
} 
double VmLBOconstNuSurfNonUniform1x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[4]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUy = &nuUSum[2]; 

  double favg[4]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 

  double fjump[4]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[2]*nuSum+0.7071067811865475*dxvl[2]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 

  double Ghat[4]; 
  for(unsigned short int i=0; i<4; ++i){ 
    Ghat[i]=0.0; 
  }; 

  const double dxvl2R2 = pow(dxvl[2],2);
  const double dxvl2R3 = pow(dxvl[2],3);
  const double dxvl2R4 = pow(dxvl[2],4);
  const double dxvr2R2 = pow(dxvr[2],2);
  const double dxvr2R3 = pow(dxvr[2],3);
  const double dxvr2R4 = pow(dxvr[2],4);

  Ghat[0] = (-(1.0*((17.32050807568877*nuVtSqSum[0]*dxvl2R2*dxvr2R2+3.464101615137754*nuVtSqSum[0]*dxvl2R3*dxvr[2]-3.464101615137754*nuVtSqSum[0]*dxvl2R4)*fr[3]+((-3.464101615137754*nuVtSqSum[0]*dxvr2R4)+3.464101615137754*nuVtSqSum[0]*dxvl[2]*dxvr2R3+17.32050807568877*nuVtSqSum[0]*dxvl2R2*dxvr2R2)*fl[3]+((18.0*fl[1]-18.0*fr[1])*nuVtSqSum[1]+(18.0*fl[0]-18.0*fr[0])*nuVtSqSum[0])*dxvl2R2*dxvr2R2))/(1.414213562373095*dxvl[2]*dxvr2R4+4.242640687119286*dxvl2R2*dxvr2R3+4.242640687119286*dxvl2R3*dxvr2R2+1.414213562373095*dxvl2R4*dxvr[2]))-0.5*(1.732050807568877*fjump[3]+fjump[0])+0.3535533905932737*(1.732050807568877*alphaDrag[0]*favg[3]+alphaDrag[1]*favg[1]+alphaDrag[0]*favg[0]); 
  Ghat[1] = (-(1.0*((12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[1]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R4)*fr[3]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R4)+2.449489742783178*nuVtSqSum[1]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2)*fl[3]+((12.72792206135786*fl[0]-12.72792206135786*fr[0])*nuVtSqSum[1]-12.72792206135786*nuVtSqSum[0]*fr[1]+12.72792206135786*nuVtSqSum[0]*fl[1])*dxvl2R2*dxvr2R2))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))+0.3535533905932737*(1.732050807568877*alphaDrag[1]*favg[3]+alphaDrag[0]*favg[1]+favg[0]*alphaDrag[1])-0.5*fjump[1]; 
  Ghat[2] = (18.0*nuVtSqSum[0]*dxvl[2]*dxvr[2]*fr[2]-18.0*nuVtSqSum[0]*dxvl[2]*dxvr[2]*fl[2])/(1.414213562373095*dxvr2R3+4.242640687119286*dxvl[2]*dxvr2R2+4.242640687119286*dxvl2R2*dxvr[2]+1.414213562373095*dxvl2R3)-0.5*fjump[2]+0.3535533905932737*alphaDrag[0]*favg[2]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 

  double incr2[4]; 

  incr2[3] = -(1.0*((5.0*nuVtSqSum[0]*dxvl2R2*dxvr[2]+3.0*nuVtSqSum[0]*dxvl2R3)*fr[3]+((-3.0*nuVtSqSum[0]*dxvr2R3)-5.0*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[3]+((-1.732050807568877*fl[1]*nuVtSqSum[1])-1.732050807568877*fl[0]*nuVtSqSum[0])*dxvr2R3+((-5.196152422706631*fl[1]*nuVtSqSum[1])-5.196152422706631*fl[0]*nuVtSqSum[0])*dxvl[2]*dxvr2R2+((-5.196152422706631*fr[1]*nuVtSqSum[1])-5.196152422706631*fr[0]*nuVtSqSum[0])*dxvl2R2*dxvr[2]+((-1.732050807568877*fr[1]*nuVtSqSum[1])-1.732050807568877*fr[0]*nuVtSqSum[0])*dxvl2R3))/(2.828427124746191*dxvr2R3+8.485281374238571*dxvl[2]*dxvr2R2+8.485281374238571*dxvl2R2*dxvr[2]+2.828427124746191*dxvl2R3); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 

  return std::abs(wl[2]-(0.7071067811865475*sumNuUy[0])/nuSum); 
} 
