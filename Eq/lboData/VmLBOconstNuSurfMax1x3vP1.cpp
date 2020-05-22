#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[5]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 

  double fjump[5]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 

  double Ghat[5]; 
  for(unsigned short int i=0; i<5; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.530931089239486*nuVtSqSum[0]*fr[2])-1.530931089239486*nuVtSqSum[0]*fl[2]+1.590990257669731*fr[1]*nuVtSqSum[1]-1.590990257669731*fl[1]*nuVtSqSum[1]+1.590990257669731*fr[0]*nuVtSqSum[0]-1.590990257669731*fl[0]*nuVtSqSum[0])*rdv-0.8660254037844386*fjump[2]+alphaDrag[0]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = ((-1.530931089239486*nuVtSqSum[1]*fr[2])-1.530931089239486*nuVtSqSum[1]*fl[2]+1.590990257669731*fr[0]*nuVtSqSum[1]-1.590990257669731*fl[0]*nuVtSqSum[1]+1.590990257669731*nuVtSqSum[0]*fr[1]-1.590990257669731*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[1]+0.3535533905932737*alphaDrag[0]*favg[1]; 
  Ghat[3] = (1.590990257669731*nuVtSqSum[0]*fr[3]-1.590990257669731*nuVtSqSum[0]*fl[3])*rdv-0.5*fjump[3]+0.3535533905932737*alphaDrag[0]*favg[3]; 
  Ghat[4] = (1.590990257669731*nuVtSqSum[0]*fr[4]-1.590990257669731*nuVtSqSum[0]*fl[4])*rdv-0.5*fjump[4]+0.3535533905932737*alphaDrag[0]*favg[4]; 

  double incr1[5]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = -0.5*Ghat[4]; 

  double incr2[5]; 
  incr2[2] = nuVtSqSum[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*(fr[1]+fl[1])*nuVtSqSum[1]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 

  return std::abs(wl[1]-(0.7071067811865475*sumNuUx[0])/nuSum); 
} 
double VmLBOconstNuSurf1x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUy = &nuUSum[2]; 

  double favg[5]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 

  double fjump[5]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[2]*nuSum+0.7071067811865475*dxvl[2]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 

  double Ghat[5]; 
  for(unsigned short int i=0; i<5; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.530931089239486*nuVtSqSum[0]*fr[3])-1.530931089239486*nuVtSqSum[0]*fl[3]+1.590990257669731*fr[1]*nuVtSqSum[1]-1.590990257669731*fl[1]*nuVtSqSum[1]+1.590990257669731*fr[0]*nuVtSqSum[0]-1.590990257669731*fl[0]*nuVtSqSum[0])*rdv-0.8660254037844386*fjump[3]+alphaDrag[0]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = ((-1.530931089239486*nuVtSqSum[1]*fr[3])-1.530931089239486*nuVtSqSum[1]*fl[3]+1.590990257669731*fr[0]*nuVtSqSum[1]-1.590990257669731*fl[0]*nuVtSqSum[1]+1.590990257669731*nuVtSqSum[0]*fr[1]-1.590990257669731*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.5*fjump[1]+0.3535533905932737*alphaDrag[0]*favg[1]; 
  Ghat[2] = (1.590990257669731*nuVtSqSum[0]*fr[2]-1.590990257669731*nuVtSqSum[0]*fl[2])*rdv-0.5*fjump[2]+0.3535533905932737*alphaDrag[0]*favg[2]; 
  Ghat[4] = (1.590990257669731*nuVtSqSum[0]*fr[4]-1.590990257669731*nuVtSqSum[0]*fl[4])*rdv-0.5*fjump[4]+0.3535533905932737*alphaDrag[0]*favg[4]; 

  double incr1[5]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 

  double incr2[5]; 
  incr2[3] = nuVtSqSum[0]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*(fr[1]+fl[1])*nuVtSqSum[1]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 

  return std::abs(wl[2]-(0.7071067811865475*sumNuUy[0])/nuSum); 
} 
double VmLBOconstNuSurf1x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[3]; 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  const double *sumNuUz = &nuUSum[4]; 

  double favg[5]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 

  double fjump[5]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(-1*fr[4])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[3]*nuSum+0.7071067811865475*dxvl[3]*nuSum-1.0*sumNuUz[0]; 
  alphaDrag[1] = -1.0*sumNuUz[1]; 

  double Ghat[5]; 
  for(unsigned short int i=0; i<5; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.530931089239486*nuVtSqSum[0]*fr[4])-1.530931089239486*nuVtSqSum[0]*fl[4]+1.590990257669731*fr[1]*nuVtSqSum[1]-1.590990257669731*fl[1]*nuVtSqSum[1]+1.590990257669731*fr[0]*nuVtSqSum[0]-1.590990257669731*fl[0]*nuVtSqSum[0])*rdv-0.8660254037844386*fjump[4]+alphaDrag[0]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = ((-1.530931089239486*nuVtSqSum[1]*fr[4])-1.530931089239486*nuVtSqSum[1]*fl[4]+1.590990257669731*fr[0]*nuVtSqSum[1]-1.590990257669731*fl[0]*nuVtSqSum[1]+1.590990257669731*nuVtSqSum[0]*fr[1]-1.590990257669731*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[0])-0.5*fjump[1]+0.3535533905932737*alphaDrag[0]*favg[1]; 
  Ghat[2] = (1.590990257669731*nuVtSqSum[0]*fr[2]-1.590990257669731*nuVtSqSum[0]*fl[2])*rdv-0.5*fjump[2]+0.3535533905932737*alphaDrag[0]*favg[2]; 
  Ghat[3] = (1.590990257669731*nuVtSqSum[0]*fr[3]-1.590990257669731*nuVtSqSum[0]*fl[3])*rdv-0.5*fjump[3]+0.3535533905932737*alphaDrag[0]*favg[3]; 

  double incr1[5]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[0]; 

  double incr2[5]; 
  incr2[4] = nuVtSqSum[0]*((-0.3535533905932737*fr[4])+0.3535533905932737*fl[4]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*(fr[1]+fl[1])*nuVtSqSum[1]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 

  return std::abs(wl[3]-(0.7071067811865475*sumNuUz[0])/nuSum); 
} 
