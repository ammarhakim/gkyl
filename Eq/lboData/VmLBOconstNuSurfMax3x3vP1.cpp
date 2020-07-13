#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf3x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:          Cell-center coordinates. 
  // dxv[6]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[12]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[3]; 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[7]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 

  double fjump[7]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 

  double alphaDrag[4]; 
  alphaDrag[0] = 2.828427124746191*wl[3]*nuSum+1.414213562373095*dxvl[3]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 
  alphaDrag[3] = -1.0*sumNuUx[3]; 

  double Ghat[7]; 
  for(unsigned short int i=0; i<7; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-0.7654655446197428*nuVtSqSum[0]*fr[4])-0.7654655446197428*nuVtSqSum[0]*fl[4]+0.7954951288348656*fr[3]*nuVtSqSum[3]-0.7954951288348656*fl[3]*nuVtSqSum[3]+0.7954951288348656*fr[2]*nuVtSqSum[2]-0.7954951288348656*fl[2]*nuVtSqSum[2]+0.7954951288348656*fr[1]*nuVtSqSum[1]-0.7954951288348656*fl[1]*nuVtSqSum[1]+0.7954951288348656*fr[0]*nuVtSqSum[0]-0.7954951288348656*fl[0]*nuVtSqSum[0])*rdv-0.8660254037844386*fjump[4]+alphaDrag[0]*(0.3061862178478971*favg[4]+0.1767766952966368*favg[0])+0.1767766952966368*alphaDrag[3]*favg[3]+0.1767766952966368*alphaDrag[2]*favg[2]+0.1767766952966368*alphaDrag[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = ((-0.7654655446197428*nuVtSqSum[1]*fr[4])-0.7654655446197428*nuVtSqSum[1]*fl[4]+0.7954951288348656*fr[0]*nuVtSqSum[1]-0.7954951288348656*fl[0]*nuVtSqSum[1]+0.7954951288348656*nuVtSqSum[0]*fr[1]-0.7954951288348656*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.3061862178478971*favg[4]+0.1767766952966368*favg[0])-0.5*fjump[1]+0.1767766952966368*alphaDrag[0]*favg[1]; 
  Ghat[2] = ((-0.7654655446197428*nuVtSqSum[2]*fr[4])-0.7654655446197428*nuVtSqSum[2]*fl[4]+0.7954951288348656*fr[0]*nuVtSqSum[2]-0.7954951288348656*fl[0]*nuVtSqSum[2]+0.7954951288348656*nuVtSqSum[0]*fr[2]-0.7954951288348656*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[2]*(0.3061862178478971*favg[4]+0.1767766952966368*favg[0])-0.5*fjump[2]+0.1767766952966368*alphaDrag[0]*favg[2]; 
  Ghat[3] = ((-0.7654655446197428*nuVtSqSum[3]*fr[4])-0.7654655446197428*nuVtSqSum[3]*fl[4]+0.7954951288348656*fr[0]*nuVtSqSum[3]-0.7954951288348656*fl[0]*nuVtSqSum[3]+0.7954951288348656*nuVtSqSum[0]*fr[3]-0.7954951288348656*nuVtSqSum[0]*fl[3])*rdv+alphaDrag[3]*(0.3061862178478971*favg[4]+0.1767766952966368*favg[0])-0.5*fjump[3]+0.1767766952966368*alphaDrag[0]*favg[3]; 
  Ghat[5] = (0.7954951288348656*nuVtSqSum[0]*fr[5]-0.7954951288348656*nuVtSqSum[0]*fl[5])*rdv-0.5*fjump[5]+0.1767766952966368*alphaDrag[0]*favg[5]; 
  Ghat[6] = (0.7954951288348656*nuVtSqSum[0]*fr[6]-0.7954951288348656*nuVtSqSum[0]*fl[6])*rdv-0.5*fjump[6]+0.1767766952966368*alphaDrag[0]*favg[6]; 

  double incr1[7]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[0]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = -0.5*Ghat[6]; 

  double incr2[7]; 
  incr2[4] = nuVtSqSum[0]*((-0.1767766952966368*fr[4])+0.1767766952966368*fl[4]+0.1530931089239486*(fr[0]+fl[0]))+0.1530931089239486*((fr[3]+fl[3])*nuVtSqSum[3]+(fr[2]+fl[2])*nuVtSqSum[2]+(fr[1]+fl[1])*nuVtSqSum[1]); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 

  return std::abs(wl[3]-(0.3535533905932737*sumNuUx[0])/nuSum); 
} 
double VmLBOconstNuSurf3x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:          Cell-center coordinates. 
  // dxv[6]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[12]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[4]; 
  double rdv2L = 2.0/dxvl[4]; 
  double rdv2R = 2.0/dxvr[4]; 
  double rdvSq4L = 4.0/(dxvl[4]*dxvl[4]); 
  double rdvSq4R = 4.0/(dxvr[4]*dxvr[4]); 

  const double *sumNuUy = &nuUSum[4]; 

  double favg[7]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 

  double fjump[7]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 

  double alphaDrag[4]; 
  alphaDrag[0] = 2.828427124746191*wl[4]*nuSum+1.414213562373095*dxvl[4]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 
  alphaDrag[3] = -1.0*sumNuUy[3]; 

  double Ghat[7]; 
  for(unsigned short int i=0; i<7; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-0.7654655446197428*nuVtSqSum[0]*fr[5])-0.7654655446197428*nuVtSqSum[0]*fl[5]+0.7954951288348656*fr[3]*nuVtSqSum[3]-0.7954951288348656*fl[3]*nuVtSqSum[3]+0.7954951288348656*fr[2]*nuVtSqSum[2]-0.7954951288348656*fl[2]*nuVtSqSum[2]+0.7954951288348656*fr[1]*nuVtSqSum[1]-0.7954951288348656*fl[1]*nuVtSqSum[1]+0.7954951288348656*fr[0]*nuVtSqSum[0]-0.7954951288348656*fl[0]*nuVtSqSum[0])*rdv-0.8660254037844386*fjump[5]+alphaDrag[0]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[0])+0.1767766952966368*alphaDrag[3]*favg[3]+0.1767766952966368*alphaDrag[2]*favg[2]+0.1767766952966368*alphaDrag[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = ((-0.7654655446197428*nuVtSqSum[1]*fr[5])-0.7654655446197428*nuVtSqSum[1]*fl[5]+0.7954951288348656*fr[0]*nuVtSqSum[1]-0.7954951288348656*fl[0]*nuVtSqSum[1]+0.7954951288348656*nuVtSqSum[0]*fr[1]-0.7954951288348656*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[0])-0.5*fjump[1]+0.1767766952966368*alphaDrag[0]*favg[1]; 
  Ghat[2] = ((-0.7654655446197428*nuVtSqSum[2]*fr[5])-0.7654655446197428*nuVtSqSum[2]*fl[5]+0.7954951288348656*fr[0]*nuVtSqSum[2]-0.7954951288348656*fl[0]*nuVtSqSum[2]+0.7954951288348656*nuVtSqSum[0]*fr[2]-0.7954951288348656*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[2]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[0])-0.5*fjump[2]+0.1767766952966368*alphaDrag[0]*favg[2]; 
  Ghat[3] = ((-0.7654655446197428*nuVtSqSum[3]*fr[5])-0.7654655446197428*nuVtSqSum[3]*fl[5]+0.7954951288348656*fr[0]*nuVtSqSum[3]-0.7954951288348656*fl[0]*nuVtSqSum[3]+0.7954951288348656*nuVtSqSum[0]*fr[3]-0.7954951288348656*nuVtSqSum[0]*fl[3])*rdv+alphaDrag[3]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[0])-0.5*fjump[3]+0.1767766952966368*alphaDrag[0]*favg[3]; 
  Ghat[4] = (0.7954951288348656*nuVtSqSum[0]*fr[4]-0.7954951288348656*nuVtSqSum[0]*fl[4])*rdv-0.5*fjump[4]+0.1767766952966368*alphaDrag[0]*favg[4]; 
  Ghat[6] = (0.7954951288348656*nuVtSqSum[0]*fr[6]-0.7954951288348656*nuVtSqSum[0]*fl[6])*rdv-0.5*fjump[6]+0.1767766952966368*alphaDrag[0]*favg[6]; 

  double incr1[7]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[0]; 
  incr1[6] = -0.5*Ghat[6]; 

  double incr2[7]; 
  incr2[5] = nuVtSqSum[0]*((-0.1767766952966368*fr[5])+0.1767766952966368*fl[5]+0.1530931089239486*(fr[0]+fl[0]))+0.1530931089239486*((fr[3]+fl[3])*nuVtSqSum[3]+(fr[2]+fl[2])*nuVtSqSum[2]+(fr[1]+fl[1])*nuVtSqSum[1]); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 

  return std::abs(wl[4]-(0.3535533905932737*sumNuUy[0])/nuSum); 
} 
double VmLBOconstNuSurf3x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:          Cell-center coordinates. 
  // dxv[6]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[12]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[5]; 
  double rdv2L = 2.0/dxvl[5]; 
  double rdv2R = 2.0/dxvr[5]; 
  double rdvSq4L = 4.0/(dxvl[5]*dxvl[5]); 
  double rdvSq4R = 4.0/(dxvr[5]*dxvr[5]); 

  const double *sumNuUz = &nuUSum[8]; 

  double favg[7]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 

  double fjump[7]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(-1*fr[6])); 

  double alphaDrag[4]; 
  alphaDrag[0] = 2.828427124746191*wl[5]*nuSum+1.414213562373095*dxvl[5]*nuSum-1.0*sumNuUz[0]; 
  alphaDrag[1] = -1.0*sumNuUz[1]; 
  alphaDrag[2] = -1.0*sumNuUz[2]; 
  alphaDrag[3] = -1.0*sumNuUz[3]; 

  double Ghat[7]; 
  for(unsigned short int i=0; i<7; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-0.7654655446197428*nuVtSqSum[0]*fr[6])-0.7654655446197428*nuVtSqSum[0]*fl[6]+0.7954951288348656*fr[3]*nuVtSqSum[3]-0.7954951288348656*fl[3]*nuVtSqSum[3]+0.7954951288348656*fr[2]*nuVtSqSum[2]-0.7954951288348656*fl[2]*nuVtSqSum[2]+0.7954951288348656*fr[1]*nuVtSqSum[1]-0.7954951288348656*fl[1]*nuVtSqSum[1]+0.7954951288348656*fr[0]*nuVtSqSum[0]-0.7954951288348656*fl[0]*nuVtSqSum[0])*rdv-0.8660254037844386*fjump[6]+alphaDrag[0]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[0])+0.1767766952966368*alphaDrag[3]*favg[3]+0.1767766952966368*alphaDrag[2]*favg[2]+0.1767766952966368*alphaDrag[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = ((-0.7654655446197428*nuVtSqSum[1]*fr[6])-0.7654655446197428*nuVtSqSum[1]*fl[6]+0.7954951288348656*fr[0]*nuVtSqSum[1]-0.7954951288348656*fl[0]*nuVtSqSum[1]+0.7954951288348656*nuVtSqSum[0]*fr[1]-0.7954951288348656*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[0])-0.5*fjump[1]+0.1767766952966368*alphaDrag[0]*favg[1]; 
  Ghat[2] = ((-0.7654655446197428*nuVtSqSum[2]*fr[6])-0.7654655446197428*nuVtSqSum[2]*fl[6]+0.7954951288348656*fr[0]*nuVtSqSum[2]-0.7954951288348656*fl[0]*nuVtSqSum[2]+0.7954951288348656*nuVtSqSum[0]*fr[2]-0.7954951288348656*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[2]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[0])-0.5*fjump[2]+0.1767766952966368*alphaDrag[0]*favg[2]; 
  Ghat[3] = ((-0.7654655446197428*nuVtSqSum[3]*fr[6])-0.7654655446197428*nuVtSqSum[3]*fl[6]+0.7954951288348656*fr[0]*nuVtSqSum[3]-0.7954951288348656*fl[0]*nuVtSqSum[3]+0.7954951288348656*nuVtSqSum[0]*fr[3]-0.7954951288348656*nuVtSqSum[0]*fl[3])*rdv+alphaDrag[3]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[0])-0.5*fjump[3]+0.1767766952966368*alphaDrag[0]*favg[3]; 
  Ghat[4] = (0.7954951288348656*nuVtSqSum[0]*fr[4]-0.7954951288348656*nuVtSqSum[0]*fl[4])*rdv-0.5*fjump[4]+0.1767766952966368*alphaDrag[0]*favg[4]; 
  Ghat[5] = (0.7954951288348656*nuVtSqSum[0]*fr[5]-0.7954951288348656*nuVtSqSum[0]*fl[5])*rdv-0.5*fjump[5]+0.1767766952966368*alphaDrag[0]*favg[5]; 

  double incr1[7]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[0]; 

  double incr2[7]; 
  incr2[6] = nuVtSqSum[0]*((-0.1767766952966368*fr[6])+0.1767766952966368*fl[6]+0.1530931089239486*(fr[0]+fl[0]))+0.1530931089239486*((fr[3]+fl[3])*nuVtSqSum[3]+(fr[2]+fl[2])*nuVtSqSum[2]+(fr[1]+fl[1])*nuVtSqSum[1]); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 

  return std::abs(wl[5]-(0.3535533905932737*sumNuUz[0])/nuSum); 
} 
