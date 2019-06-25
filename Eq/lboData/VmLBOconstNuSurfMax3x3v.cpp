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

  Ghat[0] = ((-0.7654655446197428*nuVtSqSum[0]*fr[4])-0.7654655446197428*nuVtSqSum[0]*fl[4]+0.7954951288348656*fr[3]*nuVtSqSum[3]-0.7954951288348656*fl[3]*nuVtSqSum[3]+0.7954951288348656*fr[2]*nuVtSqSum[2]-0.7954951288348656*fl[2]*nuVtSqSum[2]+0.7954951288348656*fr[1]*nuVtSqSum[1]-0.7954951288348656*fl[1]*nuVtSqSum[1]+0.7954951288348656*fr[0]*nuVtSqSum[0]-0.7954951288348656*fl[0]*nuVtSqSum[0])*rdv-0.8660254037844386*fjump[4]+alphaDrag[0]*(0.3061862178478971*favg[4]+0.1767766952966368*favg[0])+0.1767766952966368*favg[3]*alphaDrag[3]+0.1767766952966368*favg[2]*alphaDrag[2]+0.1767766952966368*favg[1]*alphaDrag[1]-0.5*fjump[0]; 
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
double VmLBOconstNuSurf3x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:          Cell-center coordinates. 
  // dxv[6]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[30]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[10]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[3]; 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[28]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = -1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = -1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = -1*fr[20]+fl[20]; 
  favg[21] = 1*fr[21]+fl[21]; 
  favg[22] = 1*fr[22]+fl[22]; 
  favg[23] = 1*fr[23]+fl[23]; 
  favg[24] = 1*fr[24]+fl[24]; 
  favg[25] = 1*fr[25]+fl[25]; 
  favg[26] = 1*fr[26]+fl[26]; 
  favg[27] = 1*fr[27]+fl[27]; 

  double fjump[28]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(-1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = nuSum*vMuMidMax*(fl[16]-(-1*fr[16])); 
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(1*fr[19])); 
  fjump[20] = nuSum*vMuMidMax*(fl[20]-(-1*fr[20])); 
  fjump[21] = nuSum*vMuMidMax*(fl[21]-(1*fr[21])); 
  fjump[22] = nuSum*vMuMidMax*(fl[22]-(1*fr[22])); 
  fjump[23] = nuSum*vMuMidMax*(fl[23]-(1*fr[23])); 
  fjump[24] = nuSum*vMuMidMax*(fl[24]-(1*fr[24])); 
  fjump[25] = nuSum*vMuMidMax*(fl[25]-(1*fr[25])); 
  fjump[26] = nuSum*vMuMidMax*(fl[26]-(1*fr[26])); 
  fjump[27] = nuSum*vMuMidMax*(fl[27]-(1*fr[27])); 

  double alphaDrag[10]; 
  alphaDrag[0] = 2.828427124746191*wl[3]*nuSum+1.414213562373095*dxvl[3]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 
  alphaDrag[3] = -1.0*sumNuUx[3]; 
  alphaDrag[4] = -1.0*sumNuUx[4]; 
  alphaDrag[5] = -1.0*sumNuUx[5]; 
  alphaDrag[6] = -1.0*sumNuUx[6]; 
  alphaDrag[7] = -1.0*sumNuUx[7]; 
  alphaDrag[8] = -1.0*sumNuUx[8]; 
  alphaDrag[9] = -1.0*sumNuUx[9]; 

  double Ghat[28]; 
  for(unsigned short int i=0; i<28; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (0.9486832980505137*nuVtSqSum[0]*fr[25]-0.9486832980505137*nuVtSqSum[0]*fl[25]+1.325825214724776*nuVtSqSum[9]*fr[24]-1.325825214724776*nuVtSqSum[9]*fl[24]+1.325825214724776*nuVtSqSum[8]*fr[23]-1.325825214724776*nuVtSqSum[8]*fl[23]+1.325825214724776*nuVtSqSum[7]*fr[22]-1.325825214724776*nuVtSqSum[7]*fl[22]-1.684024198163434*nuVtSqSum[3]*fr[12]-1.684024198163434*nuVtSqSum[3]*fl[12]-1.684024198163434*nuVtSqSum[2]*fr[11]-1.684024198163434*nuVtSqSum[2]*fl[11]-1.684024198163434*nuVtSqSum[1]*fr[10]-1.684024198163434*nuVtSqSum[1]*fl[10]+1.325825214724776*nuVtSqSum[6]*fr[9]-1.325825214724776*nuVtSqSum[6]*fl[9]+1.325825214724776*nuVtSqSum[5]*fr[8]-1.325825214724776*nuVtSqSum[5]*fl[8]+1.325825214724776*nuVtSqSum[4]*fr[7]-1.325825214724776*nuVtSqSum[4]*fl[7]-1.684024198163434*nuVtSqSum[0]*fr[4]-1.684024198163434*nuVtSqSum[0]*fl[4]+1.325825214724776*fr[3]*nuVtSqSum[3]-1.325825214724776*fl[3]*nuVtSqSum[3]+1.325825214724776*fr[2]*nuVtSqSum[2]-1.325825214724776*fl[2]*nuVtSqSum[2]+1.325825214724776*fr[1]*nuVtSqSum[1]-1.325825214724776*fl[1]*nuVtSqSum[1]+1.325825214724776*fr[0]*nuVtSqSum[0]-1.325825214724776*fl[0]*nuVtSqSum[0])*rdv-1.118033988749895*fjump[25]+alphaDrag[0]*(0.3952847075210473*favg[25]+0.3061862178478971*favg[4]+0.1767766952966368*favg[0])+0.1767766952966368*alphaDrag[9]*favg[24]+0.1767766952966368*alphaDrag[8]*favg[23]+0.1767766952966368*alphaDrag[7]*favg[22]+alphaDrag[3]*(0.3061862178478971*favg[12]+0.1767766952966368*favg[3])+alphaDrag[2]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[2])+alphaDrag[1]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[6]*favg[9]+0.1767766952966368*alphaDrag[5]*favg[8]+0.1767766952966368*alphaDrag[4]*favg[7]-0.8660254037844386*fjump[4]-0.5*fjump[0]; 
  Ghat[1] = (0.9486832980505137*nuVtSqSum[1]*fr[25]-0.9486832980505137*nuVtSqSum[1]*fl[25]+1.185854122563142*nuVtSqSum[1]*fr[22]-1.185854122563142*nuVtSqSum[1]*fl[22]-1.684024198163434*nuVtSqSum[5]*fr[12]-1.684024198163434*nuVtSqSum[5]*fl[12]-1.684024198163434*nuVtSqSum[4]*fr[11]-1.684024198163434*nuVtSqSum[4]*fl[11]-1.506237033139206*nuVtSqSum[7]*fr[10]-1.684024198163434*nuVtSqSum[0]*fr[10]-1.506237033139206*nuVtSqSum[7]*fl[10]-1.684024198163434*nuVtSqSum[0]*fl[10]+1.325825214724776*nuVtSqSum[3]*fr[8]-1.325825214724776*nuVtSqSum[3]*fl[8]+1.185854122563142*fr[1]*nuVtSqSum[7]-1.185854122563142*fl[1]*nuVtSqSum[7]+1.325825214724776*nuVtSqSum[2]*fr[7]-1.325825214724776*nuVtSqSum[2]*fl[7]+1.325825214724776*fr[3]*nuVtSqSum[5]-1.325825214724776*fl[3]*nuVtSqSum[5]+1.325825214724776*fr[2]*nuVtSqSum[4]-1.325825214724776*fl[2]*nuVtSqSum[4]-1.684024198163434*nuVtSqSum[1]*fr[4]-1.684024198163434*nuVtSqSum[1]*fl[4]+1.325825214724776*fr[0]*nuVtSqSum[1]-1.325825214724776*fl[0]*nuVtSqSum[1]+1.325825214724776*nuVtSqSum[0]*fr[1]-1.325825214724776*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.3952847075210473*favg[25]+0.1581138830084189*favg[22]+0.3061862178478971*favg[4]+0.1767766952966368*favg[0])+alphaDrag[5]*(0.3061862178478971*favg[12]+0.1767766952966368*favg[3])+alphaDrag[4]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[2])-0.8660254037844386*fjump[10]+alphaDrag[0]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[1])+alphaDrag[7]*(0.273861278752583*favg[10]+0.1581138830084189*favg[1])+0.1767766952966368*alphaDrag[3]*favg[8]+0.1767766952966368*alphaDrag[2]*favg[7]-0.5*fjump[1]; 
  Ghat[2] = (0.9486832980505137*nuVtSqSum[2]*fr[25]-0.9486832980505137*nuVtSqSum[2]*fl[25]+1.185854122563142*nuVtSqSum[2]*fr[23]-1.185854122563142*nuVtSqSum[2]*fl[23]-1.684024198163434*nuVtSqSum[6]*fr[12]-1.684024198163434*nuVtSqSum[6]*fl[12]-1.506237033139206*nuVtSqSum[8]*fr[11]-1.684024198163434*nuVtSqSum[0]*fr[11]-1.506237033139206*nuVtSqSum[8]*fl[11]-1.684024198163434*nuVtSqSum[0]*fl[11]-1.684024198163434*nuVtSqSum[4]*fr[10]-1.684024198163434*nuVtSqSum[4]*fl[10]+1.325825214724776*nuVtSqSum[3]*fr[9]-1.325825214724776*nuVtSqSum[3]*fl[9]+1.185854122563142*fr[2]*nuVtSqSum[8]-1.185854122563142*fl[2]*nuVtSqSum[8]+1.325825214724776*nuVtSqSum[1]*fr[7]-1.325825214724776*nuVtSqSum[1]*fl[7]+1.325825214724776*fr[3]*nuVtSqSum[6]-1.325825214724776*fl[3]*nuVtSqSum[6]+1.325825214724776*fr[1]*nuVtSqSum[4]-1.325825214724776*fl[1]*nuVtSqSum[4]-1.684024198163434*nuVtSqSum[2]*fr[4]-1.684024198163434*nuVtSqSum[2]*fl[4]+1.325825214724776*fr[0]*nuVtSqSum[2]-1.325825214724776*fl[0]*nuVtSqSum[2]+1.325825214724776*nuVtSqSum[0]*fr[2]-1.325825214724776*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[2]*(0.3952847075210473*favg[25]+0.1581138830084189*favg[23]+0.3061862178478971*favg[4]+0.1767766952966368*favg[0])+alphaDrag[6]*(0.3061862178478971*favg[12]+0.1767766952966368*favg[3])-0.8660254037844386*fjump[11]+alphaDrag[0]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[2])+alphaDrag[8]*(0.273861278752583*favg[11]+0.1581138830084189*favg[2])+alphaDrag[4]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[3]*favg[9]+0.1767766952966368*alphaDrag[1]*favg[7]-0.5*fjump[2]; 
  Ghat[3] = (0.9486832980505137*nuVtSqSum[3]*fr[25]-0.9486832980505137*nuVtSqSum[3]*fl[25]+1.185854122563142*nuVtSqSum[3]*fr[24]-1.185854122563142*nuVtSqSum[3]*fl[24]-1.506237033139206*nuVtSqSum[9]*fr[12]-1.684024198163434*nuVtSqSum[0]*fr[12]-1.506237033139206*nuVtSqSum[9]*fl[12]-1.684024198163434*nuVtSqSum[0]*fl[12]-1.684024198163434*nuVtSqSum[6]*fr[11]-1.684024198163434*nuVtSqSum[6]*fl[11]-1.684024198163434*nuVtSqSum[5]*fr[10]-1.684024198163434*nuVtSqSum[5]*fl[10]+1.185854122563142*fr[3]*nuVtSqSum[9]-1.185854122563142*fl[3]*nuVtSqSum[9]+1.325825214724776*nuVtSqSum[2]*fr[9]-1.325825214724776*nuVtSqSum[2]*fl[9]+1.325825214724776*nuVtSqSum[1]*fr[8]-1.325825214724776*nuVtSqSum[1]*fl[8]+1.325825214724776*fr[2]*nuVtSqSum[6]-1.325825214724776*fl[2]*nuVtSqSum[6]+1.325825214724776*fr[1]*nuVtSqSum[5]-1.325825214724776*fl[1]*nuVtSqSum[5]-1.684024198163434*nuVtSqSum[3]*fr[4]-1.684024198163434*nuVtSqSum[3]*fl[4]+1.325825214724776*fr[0]*nuVtSqSum[3]-1.325825214724776*fl[0]*nuVtSqSum[3]+1.325825214724776*nuVtSqSum[0]*fr[3]-1.325825214724776*nuVtSqSum[0]*fl[3])*rdv+alphaDrag[3]*(0.3952847075210473*favg[25]+0.1581138830084189*favg[24]+0.3061862178478971*favg[4]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[12]+alphaDrag[0]*(0.3061862178478971*favg[12]+0.1767766952966368*favg[3])+alphaDrag[9]*(0.273861278752583*favg[12]+0.1581138830084189*favg[3])+alphaDrag[6]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[2])+alphaDrag[5]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[2]*favg[9]+0.1767766952966368*alphaDrag[1]*favg[8]-0.5*fjump[3]; 
  Ghat[5] = ((-1.684024198163434*nuVtSqSum[0]*fr[16])-1.684024198163434*nuVtSqSum[0]*fl[16]+1.325825214724776*nuVtSqSum[3]*fr[15]-1.325825214724776*nuVtSqSum[3]*fl[15]+1.325825214724776*nuVtSqSum[2]*fr[14]-1.325825214724776*nuVtSqSum[2]*fl[14]+1.325825214724776*nuVtSqSum[1]*fr[13]-1.325825214724776*nuVtSqSum[1]*fl[13]+1.325825214724776*nuVtSqSum[0]*fr[5]-1.325825214724776*nuVtSqSum[0]*fl[5])*rdv-0.8660254037844386*fjump[16]+alphaDrag[0]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[5])+0.1767766952966368*alphaDrag[3]*favg[15]+0.1767766952966368*alphaDrag[2]*favg[14]+0.1767766952966368*alphaDrag[1]*favg[13]-0.5*fjump[5]; 
  Ghat[6] = ((-1.684024198163434*nuVtSqSum[0]*fr[20])-1.684024198163434*nuVtSqSum[0]*fl[20]+1.325825214724776*nuVtSqSum[3]*fr[19]-1.325825214724776*nuVtSqSum[3]*fl[19]+1.325825214724776*nuVtSqSum[2]*fr[18]-1.325825214724776*nuVtSqSum[2]*fl[18]+1.325825214724776*nuVtSqSum[1]*fr[17]-1.325825214724776*nuVtSqSum[1]*fl[17]+1.325825214724776*nuVtSqSum[0]*fr[6]-1.325825214724776*nuVtSqSum[0]*fl[6])*rdv-0.8660254037844386*fjump[20]+alphaDrag[0]*(0.3061862178478971*favg[20]+0.1767766952966368*favg[6])+0.1767766952966368*alphaDrag[3]*favg[19]+0.1767766952966368*alphaDrag[2]*favg[18]+0.1767766952966368*alphaDrag[1]*favg[17]-0.5*fjump[6]; 
  Ghat[7] = (0.9486832980505137*nuVtSqSum[4]*fr[25]-0.9486832980505137*nuVtSqSum[4]*fl[25]+1.185854122563142*nuVtSqSum[4]*fr[23]-1.185854122563142*nuVtSqSum[4]*fl[23]+1.185854122563142*nuVtSqSum[4]*fr[22]-1.185854122563142*nuVtSqSum[4]*fl[22]-1.684024198163434*nuVtSqSum[1]*fr[11]-1.684024198163434*nuVtSqSum[1]*fl[11]-1.684024198163434*nuVtSqSum[2]*fr[10]-1.684024198163434*nuVtSqSum[2]*fl[10]+1.325825214724776*nuVtSqSum[5]*fr[9]-1.325825214724776*nuVtSqSum[5]*fl[9]+1.185854122563142*fr[7]*nuVtSqSum[8]-1.185854122563142*fl[7]*nuVtSqSum[8]+1.325825214724776*nuVtSqSum[6]*fr[8]-1.325825214724776*nuVtSqSum[6]*fl[8]+1.185854122563142*fr[7]*nuVtSqSum[7]-1.185854122563142*fl[7]*nuVtSqSum[7]+1.325825214724776*nuVtSqSum[0]*fr[7]-1.325825214724776*nuVtSqSum[0]*fl[7]-1.684024198163434*fr[4]*nuVtSqSum[4]-1.684024198163434*fl[4]*nuVtSqSum[4]+1.325825214724776*fr[0]*nuVtSqSum[4]-1.325825214724776*fl[0]*nuVtSqSum[4]+1.325825214724776*fr[1]*nuVtSqSum[2]-1.325825214724776*fl[1]*nuVtSqSum[2]+1.325825214724776*nuVtSqSum[1]*fr[2]-1.325825214724776*nuVtSqSum[1]*fl[2])*rdv+alphaDrag[4]*(0.3952847075210473*favg[25]+0.1581138830084189*favg[23]+0.1581138830084189*favg[22]+0.3061862178478971*favg[4]+0.1767766952966368*favg[0])+alphaDrag[1]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[2])+alphaDrag[2]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[5]*favg[9]+0.1581138830084189*favg[7]*alphaDrag[8]+0.1767766952966368*alphaDrag[6]*favg[8]+0.1581138830084189*favg[7]*alphaDrag[7]-0.5*fjump[7]+0.1767766952966368*alphaDrag[0]*favg[7]; 
  Ghat[8] = (0.9486832980505137*nuVtSqSum[5]*fr[25]-0.9486832980505137*nuVtSqSum[5]*fl[25]+1.185854122563142*nuVtSqSum[5]*fr[24]-1.185854122563142*nuVtSqSum[5]*fl[24]+1.185854122563142*nuVtSqSum[5]*fr[22]-1.185854122563142*nuVtSqSum[5]*fl[22]-1.684024198163434*nuVtSqSum[1]*fr[12]-1.684024198163434*nuVtSqSum[1]*fl[12]-1.684024198163434*nuVtSqSum[3]*fr[10]-1.684024198163434*nuVtSqSum[3]*fl[10]+1.185854122563142*fr[8]*nuVtSqSum[9]-1.185854122563142*fl[8]*nuVtSqSum[9]+1.325825214724776*nuVtSqSum[4]*fr[9]-1.325825214724776*nuVtSqSum[4]*fl[9]+1.185854122563142*nuVtSqSum[7]*fr[8]+1.325825214724776*nuVtSqSum[0]*fr[8]-1.185854122563142*nuVtSqSum[7]*fl[8]-1.325825214724776*nuVtSqSum[0]*fl[8]+1.325825214724776*nuVtSqSum[6]*fr[7]-1.325825214724776*nuVtSqSum[6]*fl[7]-1.684024198163434*fr[4]*nuVtSqSum[5]-1.684024198163434*fl[4]*nuVtSqSum[5]+1.325825214724776*fr[0]*nuVtSqSum[5]-1.325825214724776*fl[0]*nuVtSqSum[5]+1.325825214724776*fr[1]*nuVtSqSum[3]-1.325825214724776*fl[1]*nuVtSqSum[3]+1.325825214724776*nuVtSqSum[1]*fr[3]-1.325825214724776*nuVtSqSum[1]*fl[3])*rdv+alphaDrag[5]*(0.3952847075210473*favg[25]+0.1581138830084189*favg[24]+0.1581138830084189*favg[22]+0.3061862178478971*favg[4]+0.1767766952966368*favg[0])+alphaDrag[1]*(0.3061862178478971*favg[12]+0.1767766952966368*favg[3])+alphaDrag[3]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[1])+0.1581138830084189*favg[8]*alphaDrag[9]+0.1767766952966368*alphaDrag[4]*favg[9]-0.5*fjump[8]+0.1581138830084189*alphaDrag[7]*favg[8]+0.1767766952966368*alphaDrag[0]*favg[8]+0.1767766952966368*alphaDrag[6]*favg[7]; 
  Ghat[9] = (0.9486832980505137*nuVtSqSum[6]*fr[25]-0.9486832980505137*nuVtSqSum[6]*fl[25]+1.185854122563142*nuVtSqSum[6]*fr[24]-1.185854122563142*nuVtSqSum[6]*fl[24]+1.185854122563142*nuVtSqSum[6]*fr[23]-1.185854122563142*nuVtSqSum[6]*fl[23]-1.684024198163434*nuVtSqSum[2]*fr[12]-1.684024198163434*nuVtSqSum[2]*fl[12]-1.684024198163434*nuVtSqSum[3]*fr[11]-1.684024198163434*nuVtSqSum[3]*fl[11]+1.185854122563142*fr[9]*nuVtSqSum[9]-1.185854122563142*fl[9]*nuVtSqSum[9]+1.185854122563142*nuVtSqSum[8]*fr[9]+1.325825214724776*nuVtSqSum[0]*fr[9]-1.185854122563142*nuVtSqSum[8]*fl[9]-1.325825214724776*nuVtSqSum[0]*fl[9]+1.325825214724776*nuVtSqSum[4]*fr[8]-1.325825214724776*nuVtSqSum[4]*fl[8]+1.325825214724776*nuVtSqSum[5]*fr[7]-1.325825214724776*nuVtSqSum[5]*fl[7]-1.684024198163434*fr[4]*nuVtSqSum[6]-1.684024198163434*fl[4]*nuVtSqSum[6]+1.325825214724776*fr[0]*nuVtSqSum[6]-1.325825214724776*fl[0]*nuVtSqSum[6]+1.325825214724776*fr[2]*nuVtSqSum[3]-1.325825214724776*fl[2]*nuVtSqSum[3]+1.325825214724776*nuVtSqSum[2]*fr[3]-1.325825214724776*nuVtSqSum[2]*fl[3])*rdv+alphaDrag[6]*(0.3952847075210473*favg[25]+0.1581138830084189*favg[24]+0.1581138830084189*favg[23]+0.3061862178478971*favg[4]+0.1767766952966368*favg[0])+alphaDrag[2]*(0.3061862178478971*favg[12]+0.1767766952966368*favg[3])+alphaDrag[3]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[2])+0.1581138830084189*favg[9]*alphaDrag[9]-0.5*fjump[9]+0.1581138830084189*alphaDrag[8]*favg[9]+0.1767766952966368*alphaDrag[0]*favg[9]+0.1767766952966368*alphaDrag[4]*favg[8]+0.1767766952966368*alphaDrag[5]*favg[7]; 
  Ghat[13] = ((-1.684024198163434*nuVtSqSum[1]*fr[16])-1.684024198163434*nuVtSqSum[1]*fl[16]+1.325825214724776*nuVtSqSum[5]*fr[15]-1.325825214724776*nuVtSqSum[5]*fl[15]+1.325825214724776*nuVtSqSum[4]*fr[14]-1.325825214724776*nuVtSqSum[4]*fl[14]+1.185854122563142*nuVtSqSum[7]*fr[13]+1.325825214724776*nuVtSqSum[0]*fr[13]-1.185854122563142*nuVtSqSum[7]*fl[13]-1.325825214724776*nuVtSqSum[0]*fl[13]+1.325825214724776*nuVtSqSum[1]*fr[5]-1.325825214724776*nuVtSqSum[1]*fl[5])*rdv+alphaDrag[1]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[5])+0.1767766952966368*alphaDrag[5]*favg[15]+0.1767766952966368*alphaDrag[4]*favg[14]-0.5*fjump[13]+0.1581138830084189*alphaDrag[7]*favg[13]+0.1767766952966368*alphaDrag[0]*favg[13]; 
  Ghat[14] = ((-1.684024198163434*nuVtSqSum[2]*fr[16])-1.684024198163434*nuVtSqSum[2]*fl[16]+1.325825214724776*nuVtSqSum[6]*fr[15]-1.325825214724776*nuVtSqSum[6]*fl[15]+1.185854122563142*nuVtSqSum[8]*fr[14]+1.325825214724776*nuVtSqSum[0]*fr[14]-1.185854122563142*nuVtSqSum[8]*fl[14]-1.325825214724776*nuVtSqSum[0]*fl[14]+1.325825214724776*nuVtSqSum[4]*fr[13]-1.325825214724776*nuVtSqSum[4]*fl[13]+1.325825214724776*nuVtSqSum[2]*fr[5]-1.325825214724776*nuVtSqSum[2]*fl[5])*rdv+alphaDrag[2]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[5])+0.1767766952966368*alphaDrag[6]*favg[15]-0.5*fjump[14]+0.1581138830084189*alphaDrag[8]*favg[14]+0.1767766952966368*alphaDrag[0]*favg[14]+0.1767766952966368*alphaDrag[4]*favg[13]; 
  Ghat[15] = ((-1.684024198163434*nuVtSqSum[3]*fr[16])-1.684024198163434*nuVtSqSum[3]*fl[16]+1.185854122563142*nuVtSqSum[9]*fr[15]+1.325825214724776*nuVtSqSum[0]*fr[15]-1.185854122563142*nuVtSqSum[9]*fl[15]-1.325825214724776*nuVtSqSum[0]*fl[15]+1.325825214724776*nuVtSqSum[6]*fr[14]-1.325825214724776*nuVtSqSum[6]*fl[14]+1.325825214724776*nuVtSqSum[5]*fr[13]-1.325825214724776*nuVtSqSum[5]*fl[13]+1.325825214724776*nuVtSqSum[3]*fr[5]-1.325825214724776*nuVtSqSum[3]*fl[5])*rdv+alphaDrag[3]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[5])-0.5*fjump[15]+0.1581138830084189*alphaDrag[9]*favg[15]+0.1767766952966368*alphaDrag[0]*favg[15]+0.1767766952966368*alphaDrag[6]*favg[14]+0.1767766952966368*alphaDrag[5]*favg[13]; 
  Ghat[17] = ((-1.684024198163434*nuVtSqSum[1]*fr[20])-1.684024198163434*nuVtSqSum[1]*fl[20]+1.325825214724776*nuVtSqSum[5]*fr[19]-1.325825214724776*nuVtSqSum[5]*fl[19]+1.325825214724776*nuVtSqSum[4]*fr[18]-1.325825214724776*nuVtSqSum[4]*fl[18]+1.185854122563142*nuVtSqSum[7]*fr[17]+1.325825214724776*nuVtSqSum[0]*fr[17]-1.185854122563142*nuVtSqSum[7]*fl[17]-1.325825214724776*nuVtSqSum[0]*fl[17]+1.325825214724776*nuVtSqSum[1]*fr[6]-1.325825214724776*nuVtSqSum[1]*fl[6])*rdv+alphaDrag[1]*(0.3061862178478971*favg[20]+0.1767766952966368*favg[6])+0.1767766952966368*alphaDrag[5]*favg[19]+0.1767766952966368*alphaDrag[4]*favg[18]-0.5*fjump[17]+0.1581138830084189*alphaDrag[7]*favg[17]+0.1767766952966368*alphaDrag[0]*favg[17]; 
  Ghat[18] = ((-1.684024198163434*nuVtSqSum[2]*fr[20])-1.684024198163434*nuVtSqSum[2]*fl[20]+1.325825214724776*nuVtSqSum[6]*fr[19]-1.325825214724776*nuVtSqSum[6]*fl[19]+1.185854122563142*nuVtSqSum[8]*fr[18]+1.325825214724776*nuVtSqSum[0]*fr[18]-1.185854122563142*nuVtSqSum[8]*fl[18]-1.325825214724776*nuVtSqSum[0]*fl[18]+1.325825214724776*nuVtSqSum[4]*fr[17]-1.325825214724776*nuVtSqSum[4]*fl[17]+1.325825214724776*nuVtSqSum[2]*fr[6]-1.325825214724776*nuVtSqSum[2]*fl[6])*rdv+alphaDrag[2]*(0.3061862178478971*favg[20]+0.1767766952966368*favg[6])+0.1767766952966368*alphaDrag[6]*favg[19]-0.5*fjump[18]+0.1581138830084189*alphaDrag[8]*favg[18]+0.1767766952966368*alphaDrag[0]*favg[18]+0.1767766952966368*alphaDrag[4]*favg[17]; 
  Ghat[19] = ((-1.684024198163434*nuVtSqSum[3]*fr[20])-1.684024198163434*nuVtSqSum[3]*fl[20]+1.185854122563142*nuVtSqSum[9]*fr[19]+1.325825214724776*nuVtSqSum[0]*fr[19]-1.185854122563142*nuVtSqSum[9]*fl[19]-1.325825214724776*nuVtSqSum[0]*fl[19]+1.325825214724776*nuVtSqSum[6]*fr[18]-1.325825214724776*nuVtSqSum[6]*fl[18]+1.325825214724776*nuVtSqSum[5]*fr[17]-1.325825214724776*nuVtSqSum[5]*fl[17]+1.325825214724776*nuVtSqSum[3]*fr[6]-1.325825214724776*nuVtSqSum[3]*fl[6])*rdv+alphaDrag[3]*(0.3061862178478971*favg[20]+0.1767766952966368*favg[6])-0.5*fjump[19]+0.1581138830084189*alphaDrag[9]*favg[19]+0.1767766952966368*alphaDrag[0]*favg[19]+0.1767766952966368*alphaDrag[6]*favg[18]+0.1767766952966368*alphaDrag[5]*favg[17]; 
  Ghat[21] = (1.325825214724776*nuVtSqSum[0]*fr[21]-1.325825214724776*nuVtSqSum[0]*fl[21])*rdv-0.5*fjump[21]+0.1767766952966368*alphaDrag[0]*favg[21]; 
  Ghat[22] = (0.9486832980505137*nuVtSqSum[7]*fr[25]-0.9486832980505137*nuVtSqSum[7]*fl[25]+0.8470386589736728*nuVtSqSum[7]*fr[22]+1.325825214724776*nuVtSqSum[0]*fr[22]-0.8470386589736728*nuVtSqSum[7]*fl[22]-1.325825214724776*nuVtSqSum[0]*fl[22]-1.506237033139206*nuVtSqSum[1]*fr[10]-1.506237033139206*nuVtSqSum[1]*fl[10]+1.185854122563142*nuVtSqSum[5]*fr[8]-1.185854122563142*nuVtSqSum[5]*fl[8]-1.684024198163434*fr[4]*nuVtSqSum[7]-1.684024198163434*fl[4]*nuVtSqSum[7]+1.325825214724776*fr[0]*nuVtSqSum[7]-1.325825214724776*fl[0]*nuVtSqSum[7]+1.185854122563142*nuVtSqSum[4]*fr[7]-1.185854122563142*nuVtSqSum[4]*fl[7]+1.185854122563142*fr[1]*nuVtSqSum[1]-1.185854122563142*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[7]*(0.3952847075210473*favg[25]+0.1129384878631564*favg[22]+0.3061862178478971*favg[4]+0.1767766952966368*favg[0])-0.5*fjump[22]+0.1767766952966368*alphaDrag[0]*favg[22]+alphaDrag[1]*(0.273861278752583*favg[10]+0.1581138830084189*favg[1])+0.1581138830084189*alphaDrag[5]*favg[8]+0.1581138830084189*alphaDrag[4]*favg[7]; 
  Ghat[23] = (0.9486832980505137*nuVtSqSum[8]*fr[25]-0.9486832980505137*nuVtSqSum[8]*fl[25]+0.8470386589736728*nuVtSqSum[8]*fr[23]+1.325825214724776*nuVtSqSum[0]*fr[23]-0.8470386589736728*nuVtSqSum[8]*fl[23]-1.325825214724776*nuVtSqSum[0]*fl[23]-1.506237033139206*nuVtSqSum[2]*fr[11]-1.506237033139206*nuVtSqSum[2]*fl[11]+1.185854122563142*nuVtSqSum[6]*fr[9]-1.185854122563142*nuVtSqSum[6]*fl[9]-1.684024198163434*fr[4]*nuVtSqSum[8]-1.684024198163434*fl[4]*nuVtSqSum[8]+1.325825214724776*fr[0]*nuVtSqSum[8]-1.325825214724776*fl[0]*nuVtSqSum[8]+1.185854122563142*nuVtSqSum[4]*fr[7]-1.185854122563142*nuVtSqSum[4]*fl[7]+1.185854122563142*fr[2]*nuVtSqSum[2]-1.185854122563142*fl[2]*nuVtSqSum[2])*rdv+alphaDrag[8]*(0.3952847075210473*favg[25]+0.1129384878631564*favg[23]+0.3061862178478971*favg[4]+0.1767766952966368*favg[0])-0.5*fjump[23]+0.1767766952966368*alphaDrag[0]*favg[23]+alphaDrag[2]*(0.273861278752583*favg[11]+0.1581138830084189*favg[2])+0.1581138830084189*alphaDrag[6]*favg[9]+0.1581138830084189*alphaDrag[4]*favg[7]; 
  Ghat[24] = (0.9486832980505137*nuVtSqSum[9]*fr[25]-0.9486832980505137*nuVtSqSum[9]*fl[25]+0.8470386589736728*nuVtSqSum[9]*fr[24]+1.325825214724776*nuVtSqSum[0]*fr[24]-0.8470386589736728*nuVtSqSum[9]*fl[24]-1.325825214724776*nuVtSqSum[0]*fl[24]-1.506237033139206*nuVtSqSum[3]*fr[12]-1.506237033139206*nuVtSqSum[3]*fl[12]-1.684024198163434*fr[4]*nuVtSqSum[9]-1.684024198163434*fl[4]*nuVtSqSum[9]+1.325825214724776*fr[0]*nuVtSqSum[9]-1.325825214724776*fl[0]*nuVtSqSum[9]+1.185854122563142*nuVtSqSum[6]*fr[9]-1.185854122563142*nuVtSqSum[6]*fl[9]+1.185854122563142*nuVtSqSum[5]*fr[8]-1.185854122563142*nuVtSqSum[5]*fl[8]+1.185854122563142*fr[3]*nuVtSqSum[3]-1.185854122563142*fl[3]*nuVtSqSum[3])*rdv+alphaDrag[9]*(0.3952847075210473*favg[25]+0.1129384878631564*favg[24]+0.3061862178478971*favg[4]+0.1767766952966368*favg[0])-0.5*fjump[24]+0.1767766952966368*alphaDrag[0]*favg[24]+alphaDrag[3]*(0.273861278752583*favg[12]+0.1581138830084189*favg[3])+0.1581138830084189*alphaDrag[6]*favg[9]+0.1581138830084189*alphaDrag[5]*favg[8]; 
  Ghat[26] = (1.325825214724776*nuVtSqSum[0]*fr[26]-1.325825214724776*nuVtSqSum[0]*fl[26])*rdv-0.5*fjump[26]+0.1767766952966368*alphaDrag[0]*favg[26]; 
  Ghat[27] = (1.325825214724776*nuVtSqSum[0]*fr[27]-1.325825214724776*nuVtSqSum[0]*fl[27])*rdv-0.5*fjump[27]+0.1767766952966368*alphaDrag[0]*favg[27]; 

  double incr1[28]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[0]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = -0.5*Ghat[6]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = -0.5*Ghat[9]; 
  incr1[10] = 0.8660254037844386*Ghat[1]; 
  incr1[11] = 0.8660254037844386*Ghat[2]; 
  incr1[12] = 0.8660254037844386*Ghat[3]; 
  incr1[13] = -0.5*Ghat[13]; 
  incr1[14] = -0.5*Ghat[14]; 
  incr1[15] = -0.5*Ghat[15]; 
  incr1[16] = 0.8660254037844386*Ghat[5]; 
  incr1[17] = -0.5*Ghat[17]; 
  incr1[18] = -0.5*Ghat[18]; 
  incr1[19] = -0.5*Ghat[19]; 
  incr1[20] = 0.8660254037844386*Ghat[6]; 
  incr1[21] = -0.5*Ghat[21]; 
  incr1[22] = -0.5*Ghat[22]; 
  incr1[23] = -0.5*Ghat[23]; 
  incr1[24] = -0.5*Ghat[24]; 
  incr1[25] = -1.118033988749895*Ghat[0]; 
  incr1[26] = -0.5*Ghat[26]; 
  incr1[27] = -0.5*Ghat[27]; 

  double incr2[28]; 
  incr2[4] = nuVtSqSum[0]*(0.1497678868178187*(fr[25]+fl[25])-0.215446597392776*fr[4]+0.215446597392776*fl[4]+0.1530931089239486*(fr[0]+fl[0]))+0.1530931089239486*(nuVtSqSum[9]*(fr[24]+fl[24])+nuVtSqSum[8]*(fr[23]+fl[23])+nuVtSqSum[7]*(fr[22]+fl[22]))+nuVtSqSum[3]*((-0.215446597392776*fr[12])+0.215446597392776*fl[12]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[2]*((-0.215446597392776*fr[11])+0.215446597392776*fl[11]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.215446597392776*fr[10])+0.215446597392776*fl[10]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[6]*(fr[9]+fl[9])+nuVtSqSum[5]*(fr[8]+fl[8])+nuVtSqSum[4]*(fr[7]+fl[7])); 
  incr2[10] = nuVtSqSum[1]*(0.1497678868178187*(fr[25]+fl[25])+0.1369306393762915*(fr[22]+fl[22])-0.215446597392776*fr[4]+0.215446597392776*fl[4]+0.1530931089239486*(fr[0]+fl[0]))+nuVtSqSum[5]*((-0.215446597392776*fr[12])+0.215446597392776*fl[12]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[4]*((-0.215446597392776*fr[11])+0.215446597392776*fl[11]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[7]*((-0.1927012949165104*fr[10])+0.1927012949165104*fl[10]+0.1369306393762915*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.215446597392776*fr[10])+0.215446597392776*fl[10]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[3]*(fr[8]+fl[8])+nuVtSqSum[2]*(fr[7]+fl[7])); 
  incr2[11] = nuVtSqSum[2]*(0.1497678868178187*(fr[25]+fl[25])+0.1369306393762915*(fr[23]+fl[23])-0.215446597392776*fr[4]+0.215446597392776*fl[4]+0.1530931089239486*(fr[0]+fl[0]))+nuVtSqSum[6]*((-0.215446597392776*fr[12])+0.215446597392776*fl[12]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[8]*((-0.1927012949165104*fr[11])+0.1927012949165104*fl[11]+0.1369306393762915*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.215446597392776*fr[11])+0.215446597392776*fl[11]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[4]*((-0.215446597392776*fr[10])+0.215446597392776*fl[10]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[3]*(fr[9]+fl[9])+nuVtSqSum[1]*(fr[7]+fl[7])); 
  incr2[12] = nuVtSqSum[3]*(0.1497678868178187*(fr[25]+fl[25])+0.1369306393762915*(fr[24]+fl[24])-0.215446597392776*fr[4]+0.215446597392776*fl[4]+0.1530931089239486*(fr[0]+fl[0]))+nuVtSqSum[9]*((-0.1927012949165104*fr[12])+0.1927012949165104*fl[12]+0.1369306393762915*(fr[3]+fl[3]))+nuVtSqSum[0]*((-0.215446597392776*fr[12])+0.215446597392776*fl[12]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[6]*((-0.215446597392776*fr[11])+0.215446597392776*fl[11]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[5]*((-0.215446597392776*fr[10])+0.215446597392776*fl[10]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[2]*(fr[9]+fl[9])+nuVtSqSum[1]*(fr[8]+fl[8])); 
  incr2[16] = nuVtSqSum[0]*((-0.215446597392776*fr[16])+0.215446597392776*fl[16]+0.1530931089239486*(fr[5]+fl[5]))+0.1530931089239486*(nuVtSqSum[3]*(fr[15]+fl[15])+nuVtSqSum[2]*(fr[14]+fl[14])+nuVtSqSum[1]*(fr[13]+fl[13])); 
  incr2[20] = nuVtSqSum[0]*((-0.215446597392776*fr[20])+0.215446597392776*fl[20]+0.1530931089239486*(fr[6]+fl[6]))+0.1530931089239486*(nuVtSqSum[3]*(fr[19]+fl[19])+nuVtSqSum[2]*(fr[18]+fl[18])+nuVtSqSum[1]*(fr[17]+fl[17])); 
  incr2[25] = nuVtSqSum[0]*((-0.5800485314420891*(fr[25]+fl[25]))+0.8344210836992756*fr[4]-0.8344210836992756*fl[4]-0.592927061281571*(fr[0]+fl[0]))-0.592927061281571*(nuVtSqSum[9]*(fr[24]+fl[24])+nuVtSqSum[8]*(fr[23]+fl[23])+nuVtSqSum[7]*(fr[22]+fl[22]))+nuVtSqSum[3]*(0.8344210836992756*fr[12]-0.8344210836992756*fl[12]-0.592927061281571*(fr[3]+fl[3]))+nuVtSqSum[2]*(0.8344210836992756*fr[11]-0.8344210836992756*fl[11]-0.592927061281571*(fr[2]+fl[2]))+nuVtSqSum[1]*(0.8344210836992756*fr[10]-0.8344210836992756*fl[10]-0.592927061281571*(fr[1]+fl[1]))-0.592927061281571*(nuVtSqSum[6]*(fr[9]+fl[9])+nuVtSqSum[5]*(fr[8]+fl[8])+nuVtSqSum[4]*(fr[7]+fl[7])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr1[9]*rdv2R; 
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr2[11]*rdvSq4R+incr1[11]*rdv2R; 
  outr[12] += incr2[12]*rdvSq4R+incr1[12]*rdv2R; 
  outr[13] += incr1[13]*rdv2R; 
  outr[14] += incr1[14]*rdv2R; 
  outr[15] += incr1[15]*rdv2R; 
  outr[16] += incr2[16]*rdvSq4R+incr1[16]*rdv2R; 
  outr[17] += incr1[17]*rdv2R; 
  outr[18] += incr1[18]*rdv2R; 
  outr[19] += incr1[19]*rdv2R; 
  outr[20] += incr2[20]*rdvSq4R+incr1[20]*rdv2R; 
  outr[21] += incr1[21]*rdv2R; 
  outr[22] += incr1[22]*rdv2R; 
  outr[23] += incr1[23]*rdv2R; 
  outr[24] += incr1[24]*rdv2R; 
  outr[25] += incr2[25]*rdvSq4R+incr1[25]*rdv2R; 
  outr[26] += incr1[26]*rdv2R; 
  outr[27] += incr1[27]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += -1.0*incr1[9]*rdv2L; 
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += incr1[11]*rdv2L-1.0*incr2[11]*rdvSq4L; 
  outl[12] += incr1[12]*rdv2L-1.0*incr2[12]*rdvSq4L; 
  outl[13] += -1.0*incr1[13]*rdv2L; 
  outl[14] += -1.0*incr1[14]*rdv2L; 
  outl[15] += -1.0*incr1[15]*rdv2L; 
  outl[16] += incr1[16]*rdv2L-1.0*incr2[16]*rdvSq4L; 
  outl[17] += -1.0*incr1[17]*rdv2L; 
  outl[18] += -1.0*incr1[18]*rdv2L; 
  outl[19] += -1.0*incr1[19]*rdv2L; 
  outl[20] += incr1[20]*rdv2L-1.0*incr2[20]*rdvSq4L; 
  outl[21] += -1.0*incr1[21]*rdv2L; 
  outl[22] += -1.0*incr1[22]*rdv2L; 
  outl[23] += -1.0*incr1[23]*rdv2L; 
  outl[24] += -1.0*incr1[24]*rdv2L; 
  outl[25] += incr2[25]*rdvSq4L-1.0*incr1[25]*rdv2L; 
  outl[26] += -1.0*incr1[26]*rdv2L; 
  outl[27] += -1.0*incr1[27]*rdv2L; 

  return std::abs((0.3952847075210473*sumNuUx[9])/nuSum+(0.3952847075210473*sumNuUx[8])/nuSum+(0.3952847075210473*sumNuUx[7])/nuSum-(0.3535533905932737*sumNuUx[0])/nuSum+wl[3]); 
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

  Ghat[0] = ((-0.7654655446197428*nuVtSqSum[0]*fr[5])-0.7654655446197428*nuVtSqSum[0]*fl[5]+0.7954951288348656*fr[3]*nuVtSqSum[3]-0.7954951288348656*fl[3]*nuVtSqSum[3]+0.7954951288348656*fr[2]*nuVtSqSum[2]-0.7954951288348656*fl[2]*nuVtSqSum[2]+0.7954951288348656*fr[1]*nuVtSqSum[1]-0.7954951288348656*fl[1]*nuVtSqSum[1]+0.7954951288348656*fr[0]*nuVtSqSum[0]-0.7954951288348656*fl[0]*nuVtSqSum[0])*rdv-0.8660254037844386*fjump[5]+alphaDrag[0]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[0])+0.1767766952966368*favg[3]*alphaDrag[3]+0.1767766952966368*favg[2]*alphaDrag[2]+0.1767766952966368*favg[1]*alphaDrag[1]-0.5*fjump[0]; 
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
double VmLBOconstNuSurf3x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:          Cell-center coordinates. 
  // dxv[6]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[30]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[10]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[4]; 
  double rdv2L = 2.0/dxvl[4]; 
  double rdv2R = 2.0/dxvr[4]; 
  double rdvSq4L = 4.0/(dxvl[4]*dxvl[4]); 
  double rdvSq4R = 4.0/(dxvr[4]*dxvr[4]); 

  const double *sumNuUy = &nuUSum[10]; 

  double favg[28]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 
  favg[16] = -1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  favg[21] = -1*fr[21]+fl[21]; 
  favg[22] = 1*fr[22]+fl[22]; 
  favg[23] = 1*fr[23]+fl[23]; 
  favg[24] = 1*fr[24]+fl[24]; 
  favg[25] = 1*fr[25]+fl[25]; 
  favg[26] = 1*fr[26]+fl[26]; 
  favg[27] = 1*fr[27]+fl[27]; 

  double fjump[28]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(-1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(-1*fr[15])); 
  fjump[16] = nuSum*vMuMidMax*(fl[16]-(-1*fr[16])); 
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(1*fr[19])); 
  fjump[20] = nuSum*vMuMidMax*(fl[20]-(1*fr[20])); 
  fjump[21] = nuSum*vMuMidMax*(fl[21]-(-1*fr[21])); 
  fjump[22] = nuSum*vMuMidMax*(fl[22]-(1*fr[22])); 
  fjump[23] = nuSum*vMuMidMax*(fl[23]-(1*fr[23])); 
  fjump[24] = nuSum*vMuMidMax*(fl[24]-(1*fr[24])); 
  fjump[25] = nuSum*vMuMidMax*(fl[25]-(1*fr[25])); 
  fjump[26] = nuSum*vMuMidMax*(fl[26]-(1*fr[26])); 
  fjump[27] = nuSum*vMuMidMax*(fl[27]-(1*fr[27])); 

  double alphaDrag[10]; 
  alphaDrag[0] = 2.828427124746191*wl[4]*nuSum+1.414213562373095*dxvl[4]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 
  alphaDrag[3] = -1.0*sumNuUy[3]; 
  alphaDrag[4] = -1.0*sumNuUy[4]; 
  alphaDrag[5] = -1.0*sumNuUy[5]; 
  alphaDrag[6] = -1.0*sumNuUy[6]; 
  alphaDrag[7] = -1.0*sumNuUy[7]; 
  alphaDrag[8] = -1.0*sumNuUy[8]; 
  alphaDrag[9] = -1.0*sumNuUy[9]; 

  double Ghat[28]; 
  for(unsigned short int i=0; i<28; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (0.9486832980505137*nuVtSqSum[0]*fr[26]-0.9486832980505137*nuVtSqSum[0]*fl[26]+1.325825214724776*nuVtSqSum[9]*fr[24]-1.325825214724776*nuVtSqSum[9]*fl[24]+1.325825214724776*nuVtSqSum[8]*fr[23]-1.325825214724776*nuVtSqSum[8]*fl[23]+1.325825214724776*nuVtSqSum[7]*fr[22]-1.325825214724776*nuVtSqSum[7]*fl[22]-1.684024198163434*nuVtSqSum[3]*fr[15]-1.684024198163434*nuVtSqSum[3]*fl[15]-1.684024198163434*nuVtSqSum[2]*fr[14]-1.684024198163434*nuVtSqSum[2]*fl[14]-1.684024198163434*nuVtSqSum[1]*fr[13]-1.684024198163434*nuVtSqSum[1]*fl[13]+1.325825214724776*nuVtSqSum[6]*fr[9]-1.325825214724776*nuVtSqSum[6]*fl[9]+1.325825214724776*nuVtSqSum[5]*fr[8]-1.325825214724776*nuVtSqSum[5]*fl[8]+1.325825214724776*nuVtSqSum[4]*fr[7]-1.325825214724776*nuVtSqSum[4]*fl[7]-1.684024198163434*nuVtSqSum[0]*fr[5]-1.684024198163434*nuVtSqSum[0]*fl[5]+1.325825214724776*fr[3]*nuVtSqSum[3]-1.325825214724776*fl[3]*nuVtSqSum[3]+1.325825214724776*fr[2]*nuVtSqSum[2]-1.325825214724776*fl[2]*nuVtSqSum[2]+1.325825214724776*fr[1]*nuVtSqSum[1]-1.325825214724776*fl[1]*nuVtSqSum[1]+1.325825214724776*fr[0]*nuVtSqSum[0]-1.325825214724776*fl[0]*nuVtSqSum[0])*rdv-1.118033988749895*fjump[26]+alphaDrag[0]*(0.3952847075210473*favg[26]+0.3061862178478971*favg[5]+0.1767766952966368*favg[0])+0.1767766952966368*alphaDrag[9]*favg[24]+0.1767766952966368*alphaDrag[8]*favg[23]+0.1767766952966368*alphaDrag[7]*favg[22]+alphaDrag[3]*(0.3061862178478971*favg[15]+0.1767766952966368*favg[3])+alphaDrag[2]*(0.3061862178478971*favg[14]+0.1767766952966368*favg[2])+alphaDrag[1]*(0.3061862178478971*favg[13]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[6]*favg[9]+0.1767766952966368*alphaDrag[5]*favg[8]+0.1767766952966368*alphaDrag[4]*favg[7]-0.8660254037844386*fjump[5]-0.5*fjump[0]; 
  Ghat[1] = (0.9486832980505137*nuVtSqSum[1]*fr[26]-0.9486832980505137*nuVtSqSum[1]*fl[26]+1.185854122563142*nuVtSqSum[1]*fr[22]-1.185854122563142*nuVtSqSum[1]*fl[22]-1.684024198163434*nuVtSqSum[5]*fr[15]-1.684024198163434*nuVtSqSum[5]*fl[15]-1.684024198163434*nuVtSqSum[4]*fr[14]-1.684024198163434*nuVtSqSum[4]*fl[14]-1.506237033139206*nuVtSqSum[7]*fr[13]-1.684024198163434*nuVtSqSum[0]*fr[13]-1.506237033139206*nuVtSqSum[7]*fl[13]-1.684024198163434*nuVtSqSum[0]*fl[13]+1.325825214724776*nuVtSqSum[3]*fr[8]-1.325825214724776*nuVtSqSum[3]*fl[8]+1.185854122563142*fr[1]*nuVtSqSum[7]-1.185854122563142*fl[1]*nuVtSqSum[7]+1.325825214724776*nuVtSqSum[2]*fr[7]-1.325825214724776*nuVtSqSum[2]*fl[7]+1.325825214724776*fr[3]*nuVtSqSum[5]-1.325825214724776*fl[3]*nuVtSqSum[5]-1.684024198163434*nuVtSqSum[1]*fr[5]-1.684024198163434*nuVtSqSum[1]*fl[5]+1.325825214724776*fr[2]*nuVtSqSum[4]-1.325825214724776*fl[2]*nuVtSqSum[4]+1.325825214724776*fr[0]*nuVtSqSum[1]-1.325825214724776*fl[0]*nuVtSqSum[1]+1.325825214724776*nuVtSqSum[0]*fr[1]-1.325825214724776*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.3952847075210473*favg[26]+0.1581138830084189*favg[22]+0.3061862178478971*favg[5]+0.1767766952966368*favg[0])+alphaDrag[5]*(0.3061862178478971*favg[15]+0.1767766952966368*favg[3])+alphaDrag[4]*(0.3061862178478971*favg[14]+0.1767766952966368*favg[2])-0.8660254037844386*fjump[13]+alphaDrag[0]*(0.3061862178478971*favg[13]+0.1767766952966368*favg[1])+alphaDrag[7]*(0.273861278752583*favg[13]+0.1581138830084189*favg[1])+0.1767766952966368*alphaDrag[3]*favg[8]+0.1767766952966368*alphaDrag[2]*favg[7]-0.5*fjump[1]; 
  Ghat[2] = (0.9486832980505137*nuVtSqSum[2]*fr[26]-0.9486832980505137*nuVtSqSum[2]*fl[26]+1.185854122563142*nuVtSqSum[2]*fr[23]-1.185854122563142*nuVtSqSum[2]*fl[23]-1.684024198163434*nuVtSqSum[6]*fr[15]-1.684024198163434*nuVtSqSum[6]*fl[15]-1.506237033139206*nuVtSqSum[8]*fr[14]-1.684024198163434*nuVtSqSum[0]*fr[14]-1.506237033139206*nuVtSqSum[8]*fl[14]-1.684024198163434*nuVtSqSum[0]*fl[14]-1.684024198163434*nuVtSqSum[4]*fr[13]-1.684024198163434*nuVtSqSum[4]*fl[13]+1.325825214724776*nuVtSqSum[3]*fr[9]-1.325825214724776*nuVtSqSum[3]*fl[9]+1.185854122563142*fr[2]*nuVtSqSum[8]-1.185854122563142*fl[2]*nuVtSqSum[8]+1.325825214724776*nuVtSqSum[1]*fr[7]-1.325825214724776*nuVtSqSum[1]*fl[7]+1.325825214724776*fr[3]*nuVtSqSum[6]-1.325825214724776*fl[3]*nuVtSqSum[6]-1.684024198163434*nuVtSqSum[2]*fr[5]-1.684024198163434*nuVtSqSum[2]*fl[5]+1.325825214724776*fr[1]*nuVtSqSum[4]-1.325825214724776*fl[1]*nuVtSqSum[4]+1.325825214724776*fr[0]*nuVtSqSum[2]-1.325825214724776*fl[0]*nuVtSqSum[2]+1.325825214724776*nuVtSqSum[0]*fr[2]-1.325825214724776*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[2]*(0.3952847075210473*favg[26]+0.1581138830084189*favg[23]+0.3061862178478971*favg[5]+0.1767766952966368*favg[0])+alphaDrag[6]*(0.3061862178478971*favg[15]+0.1767766952966368*favg[3])-0.8660254037844386*fjump[14]+alphaDrag[0]*(0.3061862178478971*favg[14]+0.1767766952966368*favg[2])+alphaDrag[8]*(0.273861278752583*favg[14]+0.1581138830084189*favg[2])+alphaDrag[4]*(0.3061862178478971*favg[13]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[3]*favg[9]+0.1767766952966368*alphaDrag[1]*favg[7]-0.5*fjump[2]; 
  Ghat[3] = (0.9486832980505137*nuVtSqSum[3]*fr[26]-0.9486832980505137*nuVtSqSum[3]*fl[26]+1.185854122563142*nuVtSqSum[3]*fr[24]-1.185854122563142*nuVtSqSum[3]*fl[24]-1.506237033139206*nuVtSqSum[9]*fr[15]-1.684024198163434*nuVtSqSum[0]*fr[15]-1.506237033139206*nuVtSqSum[9]*fl[15]-1.684024198163434*nuVtSqSum[0]*fl[15]-1.684024198163434*nuVtSqSum[6]*fr[14]-1.684024198163434*nuVtSqSum[6]*fl[14]-1.684024198163434*nuVtSqSum[5]*fr[13]-1.684024198163434*nuVtSqSum[5]*fl[13]+1.185854122563142*fr[3]*nuVtSqSum[9]-1.185854122563142*fl[3]*nuVtSqSum[9]+1.325825214724776*nuVtSqSum[2]*fr[9]-1.325825214724776*nuVtSqSum[2]*fl[9]+1.325825214724776*nuVtSqSum[1]*fr[8]-1.325825214724776*nuVtSqSum[1]*fl[8]+1.325825214724776*fr[2]*nuVtSqSum[6]-1.325825214724776*fl[2]*nuVtSqSum[6]+1.325825214724776*fr[1]*nuVtSqSum[5]-1.325825214724776*fl[1]*nuVtSqSum[5]-1.684024198163434*nuVtSqSum[3]*fr[5]-1.684024198163434*nuVtSqSum[3]*fl[5]+1.325825214724776*fr[0]*nuVtSqSum[3]-1.325825214724776*fl[0]*nuVtSqSum[3]+1.325825214724776*nuVtSqSum[0]*fr[3]-1.325825214724776*nuVtSqSum[0]*fl[3])*rdv+alphaDrag[3]*(0.3952847075210473*favg[26]+0.1581138830084189*favg[24]+0.3061862178478971*favg[5]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[15]+alphaDrag[0]*(0.3061862178478971*favg[15]+0.1767766952966368*favg[3])+alphaDrag[9]*(0.273861278752583*favg[15]+0.1581138830084189*favg[3])+alphaDrag[6]*(0.3061862178478971*favg[14]+0.1767766952966368*favg[2])+alphaDrag[5]*(0.3061862178478971*favg[13]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[2]*favg[9]+0.1767766952966368*alphaDrag[1]*favg[8]-0.5*fjump[3]; 
  Ghat[4] = ((-1.684024198163434*nuVtSqSum[0]*fr[16])-1.684024198163434*nuVtSqSum[0]*fl[16]+1.325825214724776*nuVtSqSum[3]*fr[12]-1.325825214724776*nuVtSqSum[3]*fl[12]+1.325825214724776*nuVtSqSum[2]*fr[11]-1.325825214724776*nuVtSqSum[2]*fl[11]+1.325825214724776*nuVtSqSum[1]*fr[10]-1.325825214724776*nuVtSqSum[1]*fl[10]+1.325825214724776*nuVtSqSum[0]*fr[4]-1.325825214724776*nuVtSqSum[0]*fl[4])*rdv-0.8660254037844386*fjump[16]+alphaDrag[0]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[4])+0.1767766952966368*alphaDrag[3]*favg[12]+0.1767766952966368*alphaDrag[2]*favg[11]+0.1767766952966368*alphaDrag[1]*favg[10]-0.5*fjump[4]; 
  Ghat[6] = ((-1.684024198163434*nuVtSqSum[0]*fr[21])-1.684024198163434*nuVtSqSum[0]*fl[21]+1.325825214724776*nuVtSqSum[3]*fr[19]-1.325825214724776*nuVtSqSum[3]*fl[19]+1.325825214724776*nuVtSqSum[2]*fr[18]-1.325825214724776*nuVtSqSum[2]*fl[18]+1.325825214724776*nuVtSqSum[1]*fr[17]-1.325825214724776*nuVtSqSum[1]*fl[17]+1.325825214724776*nuVtSqSum[0]*fr[6]-1.325825214724776*nuVtSqSum[0]*fl[6])*rdv-0.8660254037844386*fjump[21]+alphaDrag[0]*(0.3061862178478971*favg[21]+0.1767766952966368*favg[6])+0.1767766952966368*alphaDrag[3]*favg[19]+0.1767766952966368*alphaDrag[2]*favg[18]+0.1767766952966368*alphaDrag[1]*favg[17]-0.5*fjump[6]; 
  Ghat[7] = (0.9486832980505137*nuVtSqSum[4]*fr[26]-0.9486832980505137*nuVtSqSum[4]*fl[26]+1.185854122563142*nuVtSqSum[4]*fr[23]-1.185854122563142*nuVtSqSum[4]*fl[23]+1.185854122563142*nuVtSqSum[4]*fr[22]-1.185854122563142*nuVtSqSum[4]*fl[22]-1.684024198163434*nuVtSqSum[1]*fr[14]-1.684024198163434*nuVtSqSum[1]*fl[14]-1.684024198163434*nuVtSqSum[2]*fr[13]-1.684024198163434*nuVtSqSum[2]*fl[13]+1.325825214724776*nuVtSqSum[5]*fr[9]-1.325825214724776*nuVtSqSum[5]*fl[9]+1.185854122563142*fr[7]*nuVtSqSum[8]-1.185854122563142*fl[7]*nuVtSqSum[8]+1.325825214724776*nuVtSqSum[6]*fr[8]-1.325825214724776*nuVtSqSum[6]*fl[8]+1.185854122563142*fr[7]*nuVtSqSum[7]-1.185854122563142*fl[7]*nuVtSqSum[7]+1.325825214724776*nuVtSqSum[0]*fr[7]-1.325825214724776*nuVtSqSum[0]*fl[7]-1.684024198163434*nuVtSqSum[4]*fr[5]-1.684024198163434*nuVtSqSum[4]*fl[5]+1.325825214724776*fr[0]*nuVtSqSum[4]-1.325825214724776*fl[0]*nuVtSqSum[4]+1.325825214724776*fr[1]*nuVtSqSum[2]-1.325825214724776*fl[1]*nuVtSqSum[2]+1.325825214724776*nuVtSqSum[1]*fr[2]-1.325825214724776*nuVtSqSum[1]*fl[2])*rdv+alphaDrag[4]*(0.3952847075210473*favg[26]+0.1581138830084189*favg[23]+0.1581138830084189*favg[22]+0.3061862178478971*favg[5]+0.1767766952966368*favg[0])+alphaDrag[1]*(0.3061862178478971*favg[14]+0.1767766952966368*favg[2])+alphaDrag[2]*(0.3061862178478971*favg[13]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[5]*favg[9]+0.1581138830084189*favg[7]*alphaDrag[8]+0.1767766952966368*alphaDrag[6]*favg[8]+0.1581138830084189*favg[7]*alphaDrag[7]-0.5*fjump[7]+0.1767766952966368*alphaDrag[0]*favg[7]; 
  Ghat[8] = (0.9486832980505137*nuVtSqSum[5]*fr[26]-0.9486832980505137*nuVtSqSum[5]*fl[26]+1.185854122563142*nuVtSqSum[5]*fr[24]-1.185854122563142*nuVtSqSum[5]*fl[24]+1.185854122563142*nuVtSqSum[5]*fr[22]-1.185854122563142*nuVtSqSum[5]*fl[22]-1.684024198163434*nuVtSqSum[1]*fr[15]-1.684024198163434*nuVtSqSum[1]*fl[15]-1.684024198163434*nuVtSqSum[3]*fr[13]-1.684024198163434*nuVtSqSum[3]*fl[13]+1.185854122563142*fr[8]*nuVtSqSum[9]-1.185854122563142*fl[8]*nuVtSqSum[9]+1.325825214724776*nuVtSqSum[4]*fr[9]-1.325825214724776*nuVtSqSum[4]*fl[9]+1.185854122563142*nuVtSqSum[7]*fr[8]+1.325825214724776*nuVtSqSum[0]*fr[8]-1.185854122563142*nuVtSqSum[7]*fl[8]-1.325825214724776*nuVtSqSum[0]*fl[8]+1.325825214724776*nuVtSqSum[6]*fr[7]-1.325825214724776*nuVtSqSum[6]*fl[7]-1.684024198163434*fr[5]*nuVtSqSum[5]-1.684024198163434*fl[5]*nuVtSqSum[5]+1.325825214724776*fr[0]*nuVtSqSum[5]-1.325825214724776*fl[0]*nuVtSqSum[5]+1.325825214724776*fr[1]*nuVtSqSum[3]-1.325825214724776*fl[1]*nuVtSqSum[3]+1.325825214724776*nuVtSqSum[1]*fr[3]-1.325825214724776*nuVtSqSum[1]*fl[3])*rdv+alphaDrag[5]*(0.3952847075210473*favg[26]+0.1581138830084189*favg[24]+0.1581138830084189*favg[22]+0.3061862178478971*favg[5]+0.1767766952966368*favg[0])+alphaDrag[1]*(0.3061862178478971*favg[15]+0.1767766952966368*favg[3])+alphaDrag[3]*(0.3061862178478971*favg[13]+0.1767766952966368*favg[1])+0.1581138830084189*favg[8]*alphaDrag[9]+0.1767766952966368*alphaDrag[4]*favg[9]-0.5*fjump[8]+0.1581138830084189*alphaDrag[7]*favg[8]+0.1767766952966368*alphaDrag[0]*favg[8]+0.1767766952966368*alphaDrag[6]*favg[7]; 
  Ghat[9] = (0.9486832980505137*nuVtSqSum[6]*fr[26]-0.9486832980505137*nuVtSqSum[6]*fl[26]+1.185854122563142*nuVtSqSum[6]*fr[24]-1.185854122563142*nuVtSqSum[6]*fl[24]+1.185854122563142*nuVtSqSum[6]*fr[23]-1.185854122563142*nuVtSqSum[6]*fl[23]-1.684024198163434*nuVtSqSum[2]*fr[15]-1.684024198163434*nuVtSqSum[2]*fl[15]-1.684024198163434*nuVtSqSum[3]*fr[14]-1.684024198163434*nuVtSqSum[3]*fl[14]+1.185854122563142*fr[9]*nuVtSqSum[9]-1.185854122563142*fl[9]*nuVtSqSum[9]+1.185854122563142*nuVtSqSum[8]*fr[9]+1.325825214724776*nuVtSqSum[0]*fr[9]-1.185854122563142*nuVtSqSum[8]*fl[9]-1.325825214724776*nuVtSqSum[0]*fl[9]+1.325825214724776*nuVtSqSum[4]*fr[8]-1.325825214724776*nuVtSqSum[4]*fl[8]+1.325825214724776*nuVtSqSum[5]*fr[7]-1.325825214724776*nuVtSqSum[5]*fl[7]-1.684024198163434*fr[5]*nuVtSqSum[6]-1.684024198163434*fl[5]*nuVtSqSum[6]+1.325825214724776*fr[0]*nuVtSqSum[6]-1.325825214724776*fl[0]*nuVtSqSum[6]+1.325825214724776*fr[2]*nuVtSqSum[3]-1.325825214724776*fl[2]*nuVtSqSum[3]+1.325825214724776*nuVtSqSum[2]*fr[3]-1.325825214724776*nuVtSqSum[2]*fl[3])*rdv+alphaDrag[6]*(0.3952847075210473*favg[26]+0.1581138830084189*favg[24]+0.1581138830084189*favg[23]+0.3061862178478971*favg[5]+0.1767766952966368*favg[0])+alphaDrag[2]*(0.3061862178478971*favg[15]+0.1767766952966368*favg[3])+alphaDrag[3]*(0.3061862178478971*favg[14]+0.1767766952966368*favg[2])+0.1581138830084189*favg[9]*alphaDrag[9]-0.5*fjump[9]+0.1581138830084189*alphaDrag[8]*favg[9]+0.1767766952966368*alphaDrag[0]*favg[9]+0.1767766952966368*alphaDrag[4]*favg[8]+0.1767766952966368*alphaDrag[5]*favg[7]; 
  Ghat[10] = ((-1.684024198163434*nuVtSqSum[1]*fr[16])-1.684024198163434*nuVtSqSum[1]*fl[16]+1.325825214724776*nuVtSqSum[5]*fr[12]-1.325825214724776*nuVtSqSum[5]*fl[12]+1.325825214724776*nuVtSqSum[4]*fr[11]-1.325825214724776*nuVtSqSum[4]*fl[11]+1.185854122563142*nuVtSqSum[7]*fr[10]+1.325825214724776*nuVtSqSum[0]*fr[10]-1.185854122563142*nuVtSqSum[7]*fl[10]-1.325825214724776*nuVtSqSum[0]*fl[10]+1.325825214724776*nuVtSqSum[1]*fr[4]-1.325825214724776*nuVtSqSum[1]*fl[4])*rdv+alphaDrag[1]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[4])+0.1767766952966368*alphaDrag[5]*favg[12]+0.1767766952966368*alphaDrag[4]*favg[11]-0.5*fjump[10]+0.1581138830084189*alphaDrag[7]*favg[10]+0.1767766952966368*alphaDrag[0]*favg[10]; 
  Ghat[11] = ((-1.684024198163434*nuVtSqSum[2]*fr[16])-1.684024198163434*nuVtSqSum[2]*fl[16]+1.325825214724776*nuVtSqSum[6]*fr[12]-1.325825214724776*nuVtSqSum[6]*fl[12]+1.185854122563142*nuVtSqSum[8]*fr[11]+1.325825214724776*nuVtSqSum[0]*fr[11]-1.185854122563142*nuVtSqSum[8]*fl[11]-1.325825214724776*nuVtSqSum[0]*fl[11]+1.325825214724776*nuVtSqSum[4]*fr[10]-1.325825214724776*nuVtSqSum[4]*fl[10]+1.325825214724776*nuVtSqSum[2]*fr[4]-1.325825214724776*nuVtSqSum[2]*fl[4])*rdv+alphaDrag[2]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[4])+0.1767766952966368*alphaDrag[6]*favg[12]-0.5*fjump[11]+0.1581138830084189*alphaDrag[8]*favg[11]+0.1767766952966368*alphaDrag[0]*favg[11]+0.1767766952966368*alphaDrag[4]*favg[10]; 
  Ghat[12] = ((-1.684024198163434*nuVtSqSum[3]*fr[16])-1.684024198163434*nuVtSqSum[3]*fl[16]+1.185854122563142*nuVtSqSum[9]*fr[12]+1.325825214724776*nuVtSqSum[0]*fr[12]-1.185854122563142*nuVtSqSum[9]*fl[12]-1.325825214724776*nuVtSqSum[0]*fl[12]+1.325825214724776*nuVtSqSum[6]*fr[11]-1.325825214724776*nuVtSqSum[6]*fl[11]+1.325825214724776*nuVtSqSum[5]*fr[10]-1.325825214724776*nuVtSqSum[5]*fl[10]+1.325825214724776*nuVtSqSum[3]*fr[4]-1.325825214724776*nuVtSqSum[3]*fl[4])*rdv+alphaDrag[3]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[4])-0.5*fjump[12]+0.1581138830084189*alphaDrag[9]*favg[12]+0.1767766952966368*alphaDrag[0]*favg[12]+0.1767766952966368*alphaDrag[6]*favg[11]+0.1767766952966368*alphaDrag[5]*favg[10]; 
  Ghat[17] = ((-1.684024198163434*nuVtSqSum[1]*fr[21])-1.684024198163434*nuVtSqSum[1]*fl[21]+1.325825214724776*nuVtSqSum[5]*fr[19]-1.325825214724776*nuVtSqSum[5]*fl[19]+1.325825214724776*nuVtSqSum[4]*fr[18]-1.325825214724776*nuVtSqSum[4]*fl[18]+1.185854122563142*nuVtSqSum[7]*fr[17]+1.325825214724776*nuVtSqSum[0]*fr[17]-1.185854122563142*nuVtSqSum[7]*fl[17]-1.325825214724776*nuVtSqSum[0]*fl[17]+1.325825214724776*nuVtSqSum[1]*fr[6]-1.325825214724776*nuVtSqSum[1]*fl[6])*rdv+alphaDrag[1]*(0.3061862178478971*favg[21]+0.1767766952966368*favg[6])+0.1767766952966368*alphaDrag[5]*favg[19]+0.1767766952966368*alphaDrag[4]*favg[18]-0.5*fjump[17]+0.1581138830084189*alphaDrag[7]*favg[17]+0.1767766952966368*alphaDrag[0]*favg[17]; 
  Ghat[18] = ((-1.684024198163434*nuVtSqSum[2]*fr[21])-1.684024198163434*nuVtSqSum[2]*fl[21]+1.325825214724776*nuVtSqSum[6]*fr[19]-1.325825214724776*nuVtSqSum[6]*fl[19]+1.185854122563142*nuVtSqSum[8]*fr[18]+1.325825214724776*nuVtSqSum[0]*fr[18]-1.185854122563142*nuVtSqSum[8]*fl[18]-1.325825214724776*nuVtSqSum[0]*fl[18]+1.325825214724776*nuVtSqSum[4]*fr[17]-1.325825214724776*nuVtSqSum[4]*fl[17]+1.325825214724776*nuVtSqSum[2]*fr[6]-1.325825214724776*nuVtSqSum[2]*fl[6])*rdv+alphaDrag[2]*(0.3061862178478971*favg[21]+0.1767766952966368*favg[6])+0.1767766952966368*alphaDrag[6]*favg[19]-0.5*fjump[18]+0.1581138830084189*alphaDrag[8]*favg[18]+0.1767766952966368*alphaDrag[0]*favg[18]+0.1767766952966368*alphaDrag[4]*favg[17]; 
  Ghat[19] = ((-1.684024198163434*nuVtSqSum[3]*fr[21])-1.684024198163434*nuVtSqSum[3]*fl[21]+1.185854122563142*nuVtSqSum[9]*fr[19]+1.325825214724776*nuVtSqSum[0]*fr[19]-1.185854122563142*nuVtSqSum[9]*fl[19]-1.325825214724776*nuVtSqSum[0]*fl[19]+1.325825214724776*nuVtSqSum[6]*fr[18]-1.325825214724776*nuVtSqSum[6]*fl[18]+1.325825214724776*nuVtSqSum[5]*fr[17]-1.325825214724776*nuVtSqSum[5]*fl[17]+1.325825214724776*nuVtSqSum[3]*fr[6]-1.325825214724776*nuVtSqSum[3]*fl[6])*rdv+alphaDrag[3]*(0.3061862178478971*favg[21]+0.1767766952966368*favg[6])-0.5*fjump[19]+0.1581138830084189*alphaDrag[9]*favg[19]+0.1767766952966368*alphaDrag[0]*favg[19]+0.1767766952966368*alphaDrag[6]*favg[18]+0.1767766952966368*alphaDrag[5]*favg[17]; 
  Ghat[20] = (1.325825214724776*nuVtSqSum[0]*fr[20]-1.325825214724776*nuVtSqSum[0]*fl[20])*rdv-0.5*fjump[20]+0.1767766952966368*alphaDrag[0]*favg[20]; 
  Ghat[22] = (0.9486832980505137*nuVtSqSum[7]*fr[26]-0.9486832980505137*nuVtSqSum[7]*fl[26]+0.8470386589736728*nuVtSqSum[7]*fr[22]+1.325825214724776*nuVtSqSum[0]*fr[22]-0.8470386589736728*nuVtSqSum[7]*fl[22]-1.325825214724776*nuVtSqSum[0]*fl[22]-1.506237033139206*nuVtSqSum[1]*fr[13]-1.506237033139206*nuVtSqSum[1]*fl[13]+1.185854122563142*nuVtSqSum[5]*fr[8]-1.185854122563142*nuVtSqSum[5]*fl[8]-1.684024198163434*fr[5]*nuVtSqSum[7]-1.684024198163434*fl[5]*nuVtSqSum[7]+1.325825214724776*fr[0]*nuVtSqSum[7]-1.325825214724776*fl[0]*nuVtSqSum[7]+1.185854122563142*nuVtSqSum[4]*fr[7]-1.185854122563142*nuVtSqSum[4]*fl[7]+1.185854122563142*fr[1]*nuVtSqSum[1]-1.185854122563142*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[7]*(0.3952847075210473*favg[26]+0.1129384878631564*favg[22]+0.3061862178478971*favg[5]+0.1767766952966368*favg[0])-0.5*fjump[22]+0.1767766952966368*alphaDrag[0]*favg[22]+alphaDrag[1]*(0.273861278752583*favg[13]+0.1581138830084189*favg[1])+0.1581138830084189*alphaDrag[5]*favg[8]+0.1581138830084189*alphaDrag[4]*favg[7]; 
  Ghat[23] = (0.9486832980505137*nuVtSqSum[8]*fr[26]-0.9486832980505137*nuVtSqSum[8]*fl[26]+0.8470386589736728*nuVtSqSum[8]*fr[23]+1.325825214724776*nuVtSqSum[0]*fr[23]-0.8470386589736728*nuVtSqSum[8]*fl[23]-1.325825214724776*nuVtSqSum[0]*fl[23]-1.506237033139206*nuVtSqSum[2]*fr[14]-1.506237033139206*nuVtSqSum[2]*fl[14]+1.185854122563142*nuVtSqSum[6]*fr[9]-1.185854122563142*nuVtSqSum[6]*fl[9]-1.684024198163434*fr[5]*nuVtSqSum[8]-1.684024198163434*fl[5]*nuVtSqSum[8]+1.325825214724776*fr[0]*nuVtSqSum[8]-1.325825214724776*fl[0]*nuVtSqSum[8]+1.185854122563142*nuVtSqSum[4]*fr[7]-1.185854122563142*nuVtSqSum[4]*fl[7]+1.185854122563142*fr[2]*nuVtSqSum[2]-1.185854122563142*fl[2]*nuVtSqSum[2])*rdv+alphaDrag[8]*(0.3952847075210473*favg[26]+0.1129384878631564*favg[23]+0.3061862178478971*favg[5]+0.1767766952966368*favg[0])-0.5*fjump[23]+0.1767766952966368*alphaDrag[0]*favg[23]+alphaDrag[2]*(0.273861278752583*favg[14]+0.1581138830084189*favg[2])+0.1581138830084189*alphaDrag[6]*favg[9]+0.1581138830084189*alphaDrag[4]*favg[7]; 
  Ghat[24] = (0.9486832980505137*nuVtSqSum[9]*fr[26]-0.9486832980505137*nuVtSqSum[9]*fl[26]+0.8470386589736728*nuVtSqSum[9]*fr[24]+1.325825214724776*nuVtSqSum[0]*fr[24]-0.8470386589736728*nuVtSqSum[9]*fl[24]-1.325825214724776*nuVtSqSum[0]*fl[24]-1.506237033139206*nuVtSqSum[3]*fr[15]-1.506237033139206*nuVtSqSum[3]*fl[15]-1.684024198163434*fr[5]*nuVtSqSum[9]-1.684024198163434*fl[5]*nuVtSqSum[9]+1.325825214724776*fr[0]*nuVtSqSum[9]-1.325825214724776*fl[0]*nuVtSqSum[9]+1.185854122563142*nuVtSqSum[6]*fr[9]-1.185854122563142*nuVtSqSum[6]*fl[9]+1.185854122563142*nuVtSqSum[5]*fr[8]-1.185854122563142*nuVtSqSum[5]*fl[8]+1.185854122563142*fr[3]*nuVtSqSum[3]-1.185854122563142*fl[3]*nuVtSqSum[3])*rdv+alphaDrag[9]*(0.3952847075210473*favg[26]+0.1129384878631564*favg[24]+0.3061862178478971*favg[5]+0.1767766952966368*favg[0])-0.5*fjump[24]+0.1767766952966368*alphaDrag[0]*favg[24]+alphaDrag[3]*(0.273861278752583*favg[15]+0.1581138830084189*favg[3])+0.1581138830084189*alphaDrag[6]*favg[9]+0.1581138830084189*alphaDrag[5]*favg[8]; 
  Ghat[25] = (1.325825214724776*nuVtSqSum[0]*fr[25]-1.325825214724776*nuVtSqSum[0]*fl[25])*rdv-0.5*fjump[25]+0.1767766952966368*alphaDrag[0]*favg[25]; 
  Ghat[27] = (1.325825214724776*nuVtSqSum[0]*fr[27]-1.325825214724776*nuVtSqSum[0]*fl[27])*rdv-0.5*fjump[27]+0.1767766952966368*alphaDrag[0]*favg[27]; 

  double incr1[28]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[0]; 
  incr1[6] = -0.5*Ghat[6]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = -0.5*Ghat[9]; 
  incr1[10] = -0.5*Ghat[10]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = 0.8660254037844386*Ghat[1]; 
  incr1[14] = 0.8660254037844386*Ghat[2]; 
  incr1[15] = 0.8660254037844386*Ghat[3]; 
  incr1[16] = 0.8660254037844386*Ghat[4]; 
  incr1[17] = -0.5*Ghat[17]; 
  incr1[18] = -0.5*Ghat[18]; 
  incr1[19] = -0.5*Ghat[19]; 
  incr1[20] = -0.5*Ghat[20]; 
  incr1[21] = 0.8660254037844386*Ghat[6]; 
  incr1[22] = -0.5*Ghat[22]; 
  incr1[23] = -0.5*Ghat[23]; 
  incr1[24] = -0.5*Ghat[24]; 
  incr1[25] = -0.5*Ghat[25]; 
  incr1[26] = -1.118033988749895*Ghat[0]; 
  incr1[27] = -0.5*Ghat[27]; 

  double incr2[28]; 
  incr2[5] = nuVtSqSum[0]*(0.1497678868178187*(fr[26]+fl[26])-0.215446597392776*fr[5]+0.215446597392776*fl[5]+0.1530931089239486*(fr[0]+fl[0]))+0.1530931089239486*(nuVtSqSum[9]*(fr[24]+fl[24])+nuVtSqSum[8]*(fr[23]+fl[23])+nuVtSqSum[7]*(fr[22]+fl[22]))+nuVtSqSum[3]*((-0.215446597392776*fr[15])+0.215446597392776*fl[15]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[2]*((-0.215446597392776*fr[14])+0.215446597392776*fl[14]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.215446597392776*fr[13])+0.215446597392776*fl[13]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[6]*(fr[9]+fl[9])+nuVtSqSum[5]*(fr[8]+fl[8])+nuVtSqSum[4]*(fr[7]+fl[7])); 
  incr2[13] = nuVtSqSum[1]*(0.1497678868178187*(fr[26]+fl[26])+0.1369306393762915*(fr[22]+fl[22])-0.215446597392776*fr[5]+0.215446597392776*fl[5]+0.1530931089239486*(fr[0]+fl[0]))+nuVtSqSum[5]*((-0.215446597392776*fr[15])+0.215446597392776*fl[15]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[4]*((-0.215446597392776*fr[14])+0.215446597392776*fl[14]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[7]*((-0.1927012949165104*fr[13])+0.1927012949165104*fl[13]+0.1369306393762915*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.215446597392776*fr[13])+0.215446597392776*fl[13]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[3]*(fr[8]+fl[8])+nuVtSqSum[2]*(fr[7]+fl[7])); 
  incr2[14] = nuVtSqSum[2]*(0.1497678868178187*(fr[26]+fl[26])+0.1369306393762915*(fr[23]+fl[23])-0.215446597392776*fr[5]+0.215446597392776*fl[5]+0.1530931089239486*(fr[0]+fl[0]))+nuVtSqSum[6]*((-0.215446597392776*fr[15])+0.215446597392776*fl[15]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[8]*((-0.1927012949165104*fr[14])+0.1927012949165104*fl[14]+0.1369306393762915*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.215446597392776*fr[14])+0.215446597392776*fl[14]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[4]*((-0.215446597392776*fr[13])+0.215446597392776*fl[13]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[3]*(fr[9]+fl[9])+nuVtSqSum[1]*(fr[7]+fl[7])); 
  incr2[15] = nuVtSqSum[3]*(0.1497678868178187*(fr[26]+fl[26])+0.1369306393762915*(fr[24]+fl[24])-0.215446597392776*fr[5]+0.215446597392776*fl[5]+0.1530931089239486*(fr[0]+fl[0]))+nuVtSqSum[9]*((-0.1927012949165104*fr[15])+0.1927012949165104*fl[15]+0.1369306393762915*(fr[3]+fl[3]))+nuVtSqSum[0]*((-0.215446597392776*fr[15])+0.215446597392776*fl[15]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[6]*((-0.215446597392776*fr[14])+0.215446597392776*fl[14]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[5]*((-0.215446597392776*fr[13])+0.215446597392776*fl[13]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[2]*(fr[9]+fl[9])+nuVtSqSum[1]*(fr[8]+fl[8])); 
  incr2[16] = nuVtSqSum[0]*((-0.215446597392776*fr[16])+0.215446597392776*fl[16]+0.1530931089239486*(fr[4]+fl[4]))+0.1530931089239486*(nuVtSqSum[3]*(fr[12]+fl[12])+nuVtSqSum[2]*(fr[11]+fl[11])+nuVtSqSum[1]*(fr[10]+fl[10])); 
  incr2[21] = nuVtSqSum[0]*((-0.215446597392776*fr[21])+0.215446597392776*fl[21]+0.1530931089239486*(fr[6]+fl[6]))+0.1530931089239486*(nuVtSqSum[3]*(fr[19]+fl[19])+nuVtSqSum[2]*(fr[18]+fl[18])+nuVtSqSum[1]*(fr[17]+fl[17])); 
  incr2[26] = nuVtSqSum[0]*((-0.5800485314420891*(fr[26]+fl[26]))+0.8344210836992756*fr[5]-0.8344210836992756*fl[5]-0.592927061281571*(fr[0]+fl[0]))-0.592927061281571*(nuVtSqSum[9]*(fr[24]+fl[24])+nuVtSqSum[8]*(fr[23]+fl[23])+nuVtSqSum[7]*(fr[22]+fl[22]))+nuVtSqSum[3]*(0.8344210836992756*fr[15]-0.8344210836992756*fl[15]-0.592927061281571*(fr[3]+fl[3]))+nuVtSqSum[2]*(0.8344210836992756*fr[14]-0.8344210836992756*fl[14]-0.592927061281571*(fr[2]+fl[2]))+nuVtSqSum[1]*(0.8344210836992756*fr[13]-0.8344210836992756*fl[13]-0.592927061281571*(fr[1]+fl[1]))-0.592927061281571*(nuVtSqSum[6]*(fr[9]+fl[9])+nuVtSqSum[5]*(fr[8]+fl[8])+nuVtSqSum[4]*(fr[7]+fl[7])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr1[9]*rdv2R; 
  outr[10] += incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 
  outr[16] += incr2[16]*rdvSq4R+incr1[16]*rdv2R; 
  outr[17] += incr1[17]*rdv2R; 
  outr[18] += incr1[18]*rdv2R; 
  outr[19] += incr1[19]*rdv2R; 
  outr[20] += incr1[20]*rdv2R; 
  outr[21] += incr2[21]*rdvSq4R+incr1[21]*rdv2R; 
  outr[22] += incr1[22]*rdv2R; 
  outr[23] += incr1[23]*rdv2R; 
  outr[24] += incr1[24]*rdv2R; 
  outr[25] += incr1[25]*rdv2R; 
  outr[26] += incr2[26]*rdvSq4R+incr1[26]*rdv2R; 
  outr[27] += incr1[27]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += -1.0*incr1[9]*rdv2L; 
  outl[10] += -1.0*incr1[10]*rdv2L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += incr1[13]*rdv2L-1.0*incr2[13]*rdvSq4L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += incr1[15]*rdv2L-1.0*incr2[15]*rdvSq4L; 
  outl[16] += incr1[16]*rdv2L-1.0*incr2[16]*rdvSq4L; 
  outl[17] += -1.0*incr1[17]*rdv2L; 
  outl[18] += -1.0*incr1[18]*rdv2L; 
  outl[19] += -1.0*incr1[19]*rdv2L; 
  outl[20] += -1.0*incr1[20]*rdv2L; 
  outl[21] += incr1[21]*rdv2L-1.0*incr2[21]*rdvSq4L; 
  outl[22] += -1.0*incr1[22]*rdv2L; 
  outl[23] += -1.0*incr1[23]*rdv2L; 
  outl[24] += -1.0*incr1[24]*rdv2L; 
  outl[25] += -1.0*incr1[25]*rdv2L; 
  outl[26] += incr2[26]*rdvSq4L-1.0*incr1[26]*rdv2L; 
  outl[27] += -1.0*incr1[27]*rdv2L; 

  return std::abs((0.3952847075210473*sumNuUy[9])/nuSum+(0.3952847075210473*sumNuUy[8])/nuSum+(0.3952847075210473*sumNuUy[7])/nuSum-(0.3535533905932737*sumNuUy[0])/nuSum+wl[4]); 
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

  Ghat[0] = ((-0.7654655446197428*nuVtSqSum[0]*fr[6])-0.7654655446197428*nuVtSqSum[0]*fl[6]+0.7954951288348656*fr[3]*nuVtSqSum[3]-0.7954951288348656*fl[3]*nuVtSqSum[3]+0.7954951288348656*fr[2]*nuVtSqSum[2]-0.7954951288348656*fl[2]*nuVtSqSum[2]+0.7954951288348656*fr[1]*nuVtSqSum[1]-0.7954951288348656*fl[1]*nuVtSqSum[1]+0.7954951288348656*fr[0]*nuVtSqSum[0]-0.7954951288348656*fl[0]*nuVtSqSum[0])*rdv-0.8660254037844386*fjump[6]+alphaDrag[0]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[0])+0.1767766952966368*favg[3]*alphaDrag[3]+0.1767766952966368*favg[2]*alphaDrag[2]+0.1767766952966368*favg[1]*alphaDrag[1]-0.5*fjump[0]; 
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
double VmLBOconstNuSurf3x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:          Cell-center coordinates. 
  // dxv[6]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[30]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[10]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[5]; 
  double rdv2L = 2.0/dxvl[5]; 
  double rdv2R = 2.0/dxvr[5]; 
  double rdvSq4L = 4.0/(dxvl[5]*dxvl[5]); 
  double rdvSq4R = 4.0/(dxvr[5]*dxvr[5]); 

  const double *sumNuUz = &nuUSum[20]; 

  double favg[28]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = -1*fr[17]+fl[17]; 
  favg[18] = -1*fr[18]+fl[18]; 
  favg[19] = -1*fr[19]+fl[19]; 
  favg[20] = -1*fr[20]+fl[20]; 
  favg[21] = -1*fr[21]+fl[21]; 
  favg[22] = 1*fr[22]+fl[22]; 
  favg[23] = 1*fr[23]+fl[23]; 
  favg[24] = 1*fr[24]+fl[24]; 
  favg[25] = 1*fr[25]+fl[25]; 
  favg[26] = 1*fr[26]+fl[26]; 
  favg[27] = 1*fr[27]+fl[27]; 

  double fjump[28]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = nuSum*vMuMidMax*(fl[16]-(1*fr[16])); 
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(-1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(-1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(-1*fr[19])); 
  fjump[20] = nuSum*vMuMidMax*(fl[20]-(-1*fr[20])); 
  fjump[21] = nuSum*vMuMidMax*(fl[21]-(-1*fr[21])); 
  fjump[22] = nuSum*vMuMidMax*(fl[22]-(1*fr[22])); 
  fjump[23] = nuSum*vMuMidMax*(fl[23]-(1*fr[23])); 
  fjump[24] = nuSum*vMuMidMax*(fl[24]-(1*fr[24])); 
  fjump[25] = nuSum*vMuMidMax*(fl[25]-(1*fr[25])); 
  fjump[26] = nuSum*vMuMidMax*(fl[26]-(1*fr[26])); 
  fjump[27] = nuSum*vMuMidMax*(fl[27]-(1*fr[27])); 

  double alphaDrag[10]; 
  alphaDrag[0] = 2.828427124746191*wl[5]*nuSum+1.414213562373095*dxvl[5]*nuSum-1.0*sumNuUz[0]; 
  alphaDrag[1] = -1.0*sumNuUz[1]; 
  alphaDrag[2] = -1.0*sumNuUz[2]; 
  alphaDrag[3] = -1.0*sumNuUz[3]; 
  alphaDrag[4] = -1.0*sumNuUz[4]; 
  alphaDrag[5] = -1.0*sumNuUz[5]; 
  alphaDrag[6] = -1.0*sumNuUz[6]; 
  alphaDrag[7] = -1.0*sumNuUz[7]; 
  alphaDrag[8] = -1.0*sumNuUz[8]; 
  alphaDrag[9] = -1.0*sumNuUz[9]; 

  double Ghat[28]; 
  for(unsigned short int i=0; i<28; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (0.9486832980505137*nuVtSqSum[0]*fr[27]-0.9486832980505137*nuVtSqSum[0]*fl[27]+1.325825214724776*nuVtSqSum[9]*fr[24]-1.325825214724776*nuVtSqSum[9]*fl[24]+1.325825214724776*nuVtSqSum[8]*fr[23]-1.325825214724776*nuVtSqSum[8]*fl[23]+1.325825214724776*nuVtSqSum[7]*fr[22]-1.325825214724776*nuVtSqSum[7]*fl[22]-1.684024198163434*nuVtSqSum[3]*fr[19]-1.684024198163434*nuVtSqSum[3]*fl[19]-1.684024198163434*nuVtSqSum[2]*fr[18]-1.684024198163434*nuVtSqSum[2]*fl[18]-1.684024198163434*nuVtSqSum[1]*fr[17]-1.684024198163434*nuVtSqSum[1]*fl[17]+1.325825214724776*nuVtSqSum[6]*fr[9]-1.325825214724776*nuVtSqSum[6]*fl[9]+1.325825214724776*nuVtSqSum[5]*fr[8]-1.325825214724776*nuVtSqSum[5]*fl[8]+1.325825214724776*nuVtSqSum[4]*fr[7]-1.325825214724776*nuVtSqSum[4]*fl[7]-1.684024198163434*nuVtSqSum[0]*fr[6]-1.684024198163434*nuVtSqSum[0]*fl[6]+1.325825214724776*fr[3]*nuVtSqSum[3]-1.325825214724776*fl[3]*nuVtSqSum[3]+1.325825214724776*fr[2]*nuVtSqSum[2]-1.325825214724776*fl[2]*nuVtSqSum[2]+1.325825214724776*fr[1]*nuVtSqSum[1]-1.325825214724776*fl[1]*nuVtSqSum[1]+1.325825214724776*fr[0]*nuVtSqSum[0]-1.325825214724776*fl[0]*nuVtSqSum[0])*rdv-1.118033988749895*fjump[27]+alphaDrag[0]*(0.3952847075210473*favg[27]+0.3061862178478971*favg[6]+0.1767766952966368*favg[0])+0.1767766952966368*alphaDrag[9]*favg[24]+0.1767766952966368*alphaDrag[8]*favg[23]+0.1767766952966368*alphaDrag[7]*favg[22]+alphaDrag[3]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[3])+alphaDrag[2]*(0.3061862178478971*favg[18]+0.1767766952966368*favg[2])+alphaDrag[1]*(0.3061862178478971*favg[17]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[6]*favg[9]+0.1767766952966368*alphaDrag[5]*favg[8]+0.1767766952966368*alphaDrag[4]*favg[7]-0.8660254037844386*fjump[6]-0.5*fjump[0]; 
  Ghat[1] = (0.9486832980505137*nuVtSqSum[1]*fr[27]-0.9486832980505137*nuVtSqSum[1]*fl[27]+1.185854122563142*nuVtSqSum[1]*fr[22]-1.185854122563142*nuVtSqSum[1]*fl[22]-1.684024198163434*nuVtSqSum[5]*fr[19]-1.684024198163434*nuVtSqSum[5]*fl[19]-1.684024198163434*nuVtSqSum[4]*fr[18]-1.684024198163434*nuVtSqSum[4]*fl[18]-1.506237033139206*nuVtSqSum[7]*fr[17]-1.684024198163434*nuVtSqSum[0]*fr[17]-1.506237033139206*nuVtSqSum[7]*fl[17]-1.684024198163434*nuVtSqSum[0]*fl[17]+1.325825214724776*nuVtSqSum[3]*fr[8]-1.325825214724776*nuVtSqSum[3]*fl[8]+1.185854122563142*fr[1]*nuVtSqSum[7]-1.185854122563142*fl[1]*nuVtSqSum[7]+1.325825214724776*nuVtSqSum[2]*fr[7]-1.325825214724776*nuVtSqSum[2]*fl[7]-1.684024198163434*nuVtSqSum[1]*fr[6]-1.684024198163434*nuVtSqSum[1]*fl[6]+1.325825214724776*fr[3]*nuVtSqSum[5]-1.325825214724776*fl[3]*nuVtSqSum[5]+1.325825214724776*fr[2]*nuVtSqSum[4]-1.325825214724776*fl[2]*nuVtSqSum[4]+1.325825214724776*fr[0]*nuVtSqSum[1]-1.325825214724776*fl[0]*nuVtSqSum[1]+1.325825214724776*nuVtSqSum[0]*fr[1]-1.325825214724776*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.3952847075210473*favg[27]+0.1581138830084189*favg[22]+0.3061862178478971*favg[6]+0.1767766952966368*favg[0])+alphaDrag[5]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[3])+alphaDrag[4]*(0.3061862178478971*favg[18]+0.1767766952966368*favg[2])-0.8660254037844386*fjump[17]+alphaDrag[0]*(0.3061862178478971*favg[17]+0.1767766952966368*favg[1])+alphaDrag[7]*(0.273861278752583*favg[17]+0.1581138830084189*favg[1])+0.1767766952966368*alphaDrag[3]*favg[8]+0.1767766952966368*alphaDrag[2]*favg[7]-0.5*fjump[1]; 
  Ghat[2] = (0.9486832980505137*nuVtSqSum[2]*fr[27]-0.9486832980505137*nuVtSqSum[2]*fl[27]+1.185854122563142*nuVtSqSum[2]*fr[23]-1.185854122563142*nuVtSqSum[2]*fl[23]-1.684024198163434*nuVtSqSum[6]*fr[19]-1.684024198163434*nuVtSqSum[6]*fl[19]-1.506237033139206*nuVtSqSum[8]*fr[18]-1.684024198163434*nuVtSqSum[0]*fr[18]-1.506237033139206*nuVtSqSum[8]*fl[18]-1.684024198163434*nuVtSqSum[0]*fl[18]-1.684024198163434*nuVtSqSum[4]*fr[17]-1.684024198163434*nuVtSqSum[4]*fl[17]+1.325825214724776*nuVtSqSum[3]*fr[9]-1.325825214724776*nuVtSqSum[3]*fl[9]+1.185854122563142*fr[2]*nuVtSqSum[8]-1.185854122563142*fl[2]*nuVtSqSum[8]+1.325825214724776*nuVtSqSum[1]*fr[7]-1.325825214724776*nuVtSqSum[1]*fl[7]+1.325825214724776*fr[3]*nuVtSqSum[6]-1.325825214724776*fl[3]*nuVtSqSum[6]-1.684024198163434*nuVtSqSum[2]*fr[6]-1.684024198163434*nuVtSqSum[2]*fl[6]+1.325825214724776*fr[1]*nuVtSqSum[4]-1.325825214724776*fl[1]*nuVtSqSum[4]+1.325825214724776*fr[0]*nuVtSqSum[2]-1.325825214724776*fl[0]*nuVtSqSum[2]+1.325825214724776*nuVtSqSum[0]*fr[2]-1.325825214724776*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[2]*(0.3952847075210473*favg[27]+0.1581138830084189*favg[23]+0.3061862178478971*favg[6]+0.1767766952966368*favg[0])+alphaDrag[6]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[3])-0.8660254037844386*fjump[18]+alphaDrag[0]*(0.3061862178478971*favg[18]+0.1767766952966368*favg[2])+alphaDrag[8]*(0.273861278752583*favg[18]+0.1581138830084189*favg[2])+alphaDrag[4]*(0.3061862178478971*favg[17]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[3]*favg[9]+0.1767766952966368*alphaDrag[1]*favg[7]-0.5*fjump[2]; 
  Ghat[3] = (0.9486832980505137*nuVtSqSum[3]*fr[27]-0.9486832980505137*nuVtSqSum[3]*fl[27]+1.185854122563142*nuVtSqSum[3]*fr[24]-1.185854122563142*nuVtSqSum[3]*fl[24]-1.506237033139206*nuVtSqSum[9]*fr[19]-1.684024198163434*nuVtSqSum[0]*fr[19]-1.506237033139206*nuVtSqSum[9]*fl[19]-1.684024198163434*nuVtSqSum[0]*fl[19]-1.684024198163434*nuVtSqSum[6]*fr[18]-1.684024198163434*nuVtSqSum[6]*fl[18]-1.684024198163434*nuVtSqSum[5]*fr[17]-1.684024198163434*nuVtSqSum[5]*fl[17]+1.185854122563142*fr[3]*nuVtSqSum[9]-1.185854122563142*fl[3]*nuVtSqSum[9]+1.325825214724776*nuVtSqSum[2]*fr[9]-1.325825214724776*nuVtSqSum[2]*fl[9]+1.325825214724776*nuVtSqSum[1]*fr[8]-1.325825214724776*nuVtSqSum[1]*fl[8]+1.325825214724776*fr[2]*nuVtSqSum[6]-1.325825214724776*fl[2]*nuVtSqSum[6]-1.684024198163434*nuVtSqSum[3]*fr[6]-1.684024198163434*nuVtSqSum[3]*fl[6]+1.325825214724776*fr[1]*nuVtSqSum[5]-1.325825214724776*fl[1]*nuVtSqSum[5]+1.325825214724776*fr[0]*nuVtSqSum[3]-1.325825214724776*fl[0]*nuVtSqSum[3]+1.325825214724776*nuVtSqSum[0]*fr[3]-1.325825214724776*nuVtSqSum[0]*fl[3])*rdv+alphaDrag[3]*(0.3952847075210473*favg[27]+0.1581138830084189*favg[24]+0.3061862178478971*favg[6]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[19]+alphaDrag[0]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[3])+alphaDrag[9]*(0.273861278752583*favg[19]+0.1581138830084189*favg[3])+alphaDrag[6]*(0.3061862178478971*favg[18]+0.1767766952966368*favg[2])+alphaDrag[5]*(0.3061862178478971*favg[17]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[2]*favg[9]+0.1767766952966368*alphaDrag[1]*favg[8]-0.5*fjump[3]; 
  Ghat[4] = ((-1.684024198163434*nuVtSqSum[0]*fr[20])-1.684024198163434*nuVtSqSum[0]*fl[20]+1.325825214724776*nuVtSqSum[3]*fr[12]-1.325825214724776*nuVtSqSum[3]*fl[12]+1.325825214724776*nuVtSqSum[2]*fr[11]-1.325825214724776*nuVtSqSum[2]*fl[11]+1.325825214724776*nuVtSqSum[1]*fr[10]-1.325825214724776*nuVtSqSum[1]*fl[10]+1.325825214724776*nuVtSqSum[0]*fr[4]-1.325825214724776*nuVtSqSum[0]*fl[4])*rdv-0.8660254037844386*fjump[20]+alphaDrag[0]*(0.3061862178478971*favg[20]+0.1767766952966368*favg[4])+0.1767766952966368*alphaDrag[3]*favg[12]+0.1767766952966368*alphaDrag[2]*favg[11]+0.1767766952966368*alphaDrag[1]*favg[10]-0.5*fjump[4]; 
  Ghat[5] = ((-1.684024198163434*nuVtSqSum[0]*fr[21])-1.684024198163434*nuVtSqSum[0]*fl[21]+1.325825214724776*nuVtSqSum[3]*fr[15]-1.325825214724776*nuVtSqSum[3]*fl[15]+1.325825214724776*nuVtSqSum[2]*fr[14]-1.325825214724776*nuVtSqSum[2]*fl[14]+1.325825214724776*nuVtSqSum[1]*fr[13]-1.325825214724776*nuVtSqSum[1]*fl[13]+1.325825214724776*nuVtSqSum[0]*fr[5]-1.325825214724776*nuVtSqSum[0]*fl[5])*rdv-0.8660254037844386*fjump[21]+alphaDrag[0]*(0.3061862178478971*favg[21]+0.1767766952966368*favg[5])+0.1767766952966368*alphaDrag[3]*favg[15]+0.1767766952966368*alphaDrag[2]*favg[14]+0.1767766952966368*alphaDrag[1]*favg[13]-0.5*fjump[5]; 
  Ghat[7] = (0.9486832980505137*nuVtSqSum[4]*fr[27]-0.9486832980505137*nuVtSqSum[4]*fl[27]+1.185854122563142*nuVtSqSum[4]*fr[23]-1.185854122563142*nuVtSqSum[4]*fl[23]+1.185854122563142*nuVtSqSum[4]*fr[22]-1.185854122563142*nuVtSqSum[4]*fl[22]-1.684024198163434*nuVtSqSum[1]*fr[18]-1.684024198163434*nuVtSqSum[1]*fl[18]-1.684024198163434*nuVtSqSum[2]*fr[17]-1.684024198163434*nuVtSqSum[2]*fl[17]+1.325825214724776*nuVtSqSum[5]*fr[9]-1.325825214724776*nuVtSqSum[5]*fl[9]+1.185854122563142*fr[7]*nuVtSqSum[8]-1.185854122563142*fl[7]*nuVtSqSum[8]+1.325825214724776*nuVtSqSum[6]*fr[8]-1.325825214724776*nuVtSqSum[6]*fl[8]+1.185854122563142*fr[7]*nuVtSqSum[7]-1.185854122563142*fl[7]*nuVtSqSum[7]+1.325825214724776*nuVtSqSum[0]*fr[7]-1.325825214724776*nuVtSqSum[0]*fl[7]-1.684024198163434*nuVtSqSum[4]*fr[6]-1.684024198163434*nuVtSqSum[4]*fl[6]+1.325825214724776*fr[0]*nuVtSqSum[4]-1.325825214724776*fl[0]*nuVtSqSum[4]+1.325825214724776*fr[1]*nuVtSqSum[2]-1.325825214724776*fl[1]*nuVtSqSum[2]+1.325825214724776*nuVtSqSum[1]*fr[2]-1.325825214724776*nuVtSqSum[1]*fl[2])*rdv+alphaDrag[4]*(0.3952847075210473*favg[27]+0.1581138830084189*favg[23]+0.1581138830084189*favg[22]+0.3061862178478971*favg[6]+0.1767766952966368*favg[0])+alphaDrag[1]*(0.3061862178478971*favg[18]+0.1767766952966368*favg[2])+alphaDrag[2]*(0.3061862178478971*favg[17]+0.1767766952966368*favg[1])+0.1767766952966368*alphaDrag[5]*favg[9]+0.1581138830084189*favg[7]*alphaDrag[8]+0.1767766952966368*alphaDrag[6]*favg[8]+0.1581138830084189*favg[7]*alphaDrag[7]-0.5*fjump[7]+0.1767766952966368*alphaDrag[0]*favg[7]; 
  Ghat[8] = (0.9486832980505137*nuVtSqSum[5]*fr[27]-0.9486832980505137*nuVtSqSum[5]*fl[27]+1.185854122563142*nuVtSqSum[5]*fr[24]-1.185854122563142*nuVtSqSum[5]*fl[24]+1.185854122563142*nuVtSqSum[5]*fr[22]-1.185854122563142*nuVtSqSum[5]*fl[22]-1.684024198163434*nuVtSqSum[1]*fr[19]-1.684024198163434*nuVtSqSum[1]*fl[19]-1.684024198163434*nuVtSqSum[3]*fr[17]-1.684024198163434*nuVtSqSum[3]*fl[17]+1.185854122563142*fr[8]*nuVtSqSum[9]-1.185854122563142*fl[8]*nuVtSqSum[9]+1.325825214724776*nuVtSqSum[4]*fr[9]-1.325825214724776*nuVtSqSum[4]*fl[9]+1.185854122563142*nuVtSqSum[7]*fr[8]+1.325825214724776*nuVtSqSum[0]*fr[8]-1.185854122563142*nuVtSqSum[7]*fl[8]-1.325825214724776*nuVtSqSum[0]*fl[8]+1.325825214724776*nuVtSqSum[6]*fr[7]-1.325825214724776*nuVtSqSum[6]*fl[7]-1.684024198163434*nuVtSqSum[5]*fr[6]-1.684024198163434*nuVtSqSum[5]*fl[6]+1.325825214724776*fr[0]*nuVtSqSum[5]-1.325825214724776*fl[0]*nuVtSqSum[5]+1.325825214724776*fr[1]*nuVtSqSum[3]-1.325825214724776*fl[1]*nuVtSqSum[3]+1.325825214724776*nuVtSqSum[1]*fr[3]-1.325825214724776*nuVtSqSum[1]*fl[3])*rdv+alphaDrag[5]*(0.3952847075210473*favg[27]+0.1581138830084189*favg[24]+0.1581138830084189*favg[22]+0.3061862178478971*favg[6]+0.1767766952966368*favg[0])+alphaDrag[1]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[3])+alphaDrag[3]*(0.3061862178478971*favg[17]+0.1767766952966368*favg[1])+0.1581138830084189*favg[8]*alphaDrag[9]+0.1767766952966368*alphaDrag[4]*favg[9]-0.5*fjump[8]+0.1581138830084189*alphaDrag[7]*favg[8]+0.1767766952966368*alphaDrag[0]*favg[8]+0.1767766952966368*alphaDrag[6]*favg[7]; 
  Ghat[9] = (0.9486832980505137*nuVtSqSum[6]*fr[27]-0.9486832980505137*nuVtSqSum[6]*fl[27]+1.185854122563142*nuVtSqSum[6]*fr[24]-1.185854122563142*nuVtSqSum[6]*fl[24]+1.185854122563142*nuVtSqSum[6]*fr[23]-1.185854122563142*nuVtSqSum[6]*fl[23]-1.684024198163434*nuVtSqSum[2]*fr[19]-1.684024198163434*nuVtSqSum[2]*fl[19]-1.684024198163434*nuVtSqSum[3]*fr[18]-1.684024198163434*nuVtSqSum[3]*fl[18]+1.185854122563142*fr[9]*nuVtSqSum[9]-1.185854122563142*fl[9]*nuVtSqSum[9]+1.185854122563142*nuVtSqSum[8]*fr[9]+1.325825214724776*nuVtSqSum[0]*fr[9]-1.185854122563142*nuVtSqSum[8]*fl[9]-1.325825214724776*nuVtSqSum[0]*fl[9]+1.325825214724776*nuVtSqSum[4]*fr[8]-1.325825214724776*nuVtSqSum[4]*fl[8]+1.325825214724776*nuVtSqSum[5]*fr[7]-1.325825214724776*nuVtSqSum[5]*fl[7]-1.684024198163434*fr[6]*nuVtSqSum[6]-1.684024198163434*fl[6]*nuVtSqSum[6]+1.325825214724776*fr[0]*nuVtSqSum[6]-1.325825214724776*fl[0]*nuVtSqSum[6]+1.325825214724776*fr[2]*nuVtSqSum[3]-1.325825214724776*fl[2]*nuVtSqSum[3]+1.325825214724776*nuVtSqSum[2]*fr[3]-1.325825214724776*nuVtSqSum[2]*fl[3])*rdv+alphaDrag[6]*(0.3952847075210473*favg[27]+0.1581138830084189*favg[24]+0.1581138830084189*favg[23]+0.3061862178478971*favg[6]+0.1767766952966368*favg[0])+alphaDrag[2]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[3])+alphaDrag[3]*(0.3061862178478971*favg[18]+0.1767766952966368*favg[2])+0.1581138830084189*favg[9]*alphaDrag[9]-0.5*fjump[9]+0.1581138830084189*alphaDrag[8]*favg[9]+0.1767766952966368*alphaDrag[0]*favg[9]+0.1767766952966368*alphaDrag[4]*favg[8]+0.1767766952966368*alphaDrag[5]*favg[7]; 
  Ghat[10] = ((-1.684024198163434*nuVtSqSum[1]*fr[20])-1.684024198163434*nuVtSqSum[1]*fl[20]+1.325825214724776*nuVtSqSum[5]*fr[12]-1.325825214724776*nuVtSqSum[5]*fl[12]+1.325825214724776*nuVtSqSum[4]*fr[11]-1.325825214724776*nuVtSqSum[4]*fl[11]+1.185854122563142*nuVtSqSum[7]*fr[10]+1.325825214724776*nuVtSqSum[0]*fr[10]-1.185854122563142*nuVtSqSum[7]*fl[10]-1.325825214724776*nuVtSqSum[0]*fl[10]+1.325825214724776*nuVtSqSum[1]*fr[4]-1.325825214724776*nuVtSqSum[1]*fl[4])*rdv+alphaDrag[1]*(0.3061862178478971*favg[20]+0.1767766952966368*favg[4])+0.1767766952966368*alphaDrag[5]*favg[12]+0.1767766952966368*alphaDrag[4]*favg[11]-0.5*fjump[10]+0.1581138830084189*alphaDrag[7]*favg[10]+0.1767766952966368*alphaDrag[0]*favg[10]; 
  Ghat[11] = ((-1.684024198163434*nuVtSqSum[2]*fr[20])-1.684024198163434*nuVtSqSum[2]*fl[20]+1.325825214724776*nuVtSqSum[6]*fr[12]-1.325825214724776*nuVtSqSum[6]*fl[12]+1.185854122563142*nuVtSqSum[8]*fr[11]+1.325825214724776*nuVtSqSum[0]*fr[11]-1.185854122563142*nuVtSqSum[8]*fl[11]-1.325825214724776*nuVtSqSum[0]*fl[11]+1.325825214724776*nuVtSqSum[4]*fr[10]-1.325825214724776*nuVtSqSum[4]*fl[10]+1.325825214724776*nuVtSqSum[2]*fr[4]-1.325825214724776*nuVtSqSum[2]*fl[4])*rdv+alphaDrag[2]*(0.3061862178478971*favg[20]+0.1767766952966368*favg[4])+0.1767766952966368*alphaDrag[6]*favg[12]-0.5*fjump[11]+0.1581138830084189*alphaDrag[8]*favg[11]+0.1767766952966368*alphaDrag[0]*favg[11]+0.1767766952966368*alphaDrag[4]*favg[10]; 
  Ghat[12] = ((-1.684024198163434*nuVtSqSum[3]*fr[20])-1.684024198163434*nuVtSqSum[3]*fl[20]+1.185854122563142*nuVtSqSum[9]*fr[12]+1.325825214724776*nuVtSqSum[0]*fr[12]-1.185854122563142*nuVtSqSum[9]*fl[12]-1.325825214724776*nuVtSqSum[0]*fl[12]+1.325825214724776*nuVtSqSum[6]*fr[11]-1.325825214724776*nuVtSqSum[6]*fl[11]+1.325825214724776*nuVtSqSum[5]*fr[10]-1.325825214724776*nuVtSqSum[5]*fl[10]+1.325825214724776*nuVtSqSum[3]*fr[4]-1.325825214724776*nuVtSqSum[3]*fl[4])*rdv+alphaDrag[3]*(0.3061862178478971*favg[20]+0.1767766952966368*favg[4])-0.5*fjump[12]+0.1581138830084189*alphaDrag[9]*favg[12]+0.1767766952966368*alphaDrag[0]*favg[12]+0.1767766952966368*alphaDrag[6]*favg[11]+0.1767766952966368*alphaDrag[5]*favg[10]; 
  Ghat[13] = ((-1.684024198163434*nuVtSqSum[1]*fr[21])-1.684024198163434*nuVtSqSum[1]*fl[21]+1.325825214724776*nuVtSqSum[5]*fr[15]-1.325825214724776*nuVtSqSum[5]*fl[15]+1.325825214724776*nuVtSqSum[4]*fr[14]-1.325825214724776*nuVtSqSum[4]*fl[14]+1.185854122563142*nuVtSqSum[7]*fr[13]+1.325825214724776*nuVtSqSum[0]*fr[13]-1.185854122563142*nuVtSqSum[7]*fl[13]-1.325825214724776*nuVtSqSum[0]*fl[13]+1.325825214724776*nuVtSqSum[1]*fr[5]-1.325825214724776*nuVtSqSum[1]*fl[5])*rdv+alphaDrag[1]*(0.3061862178478971*favg[21]+0.1767766952966368*favg[5])+0.1767766952966368*alphaDrag[5]*favg[15]+0.1767766952966368*alphaDrag[4]*favg[14]-0.5*fjump[13]+0.1581138830084189*alphaDrag[7]*favg[13]+0.1767766952966368*alphaDrag[0]*favg[13]; 
  Ghat[14] = ((-1.684024198163434*nuVtSqSum[2]*fr[21])-1.684024198163434*nuVtSqSum[2]*fl[21]+1.325825214724776*nuVtSqSum[6]*fr[15]-1.325825214724776*nuVtSqSum[6]*fl[15]+1.185854122563142*nuVtSqSum[8]*fr[14]+1.325825214724776*nuVtSqSum[0]*fr[14]-1.185854122563142*nuVtSqSum[8]*fl[14]-1.325825214724776*nuVtSqSum[0]*fl[14]+1.325825214724776*nuVtSqSum[4]*fr[13]-1.325825214724776*nuVtSqSum[4]*fl[13]+1.325825214724776*nuVtSqSum[2]*fr[5]-1.325825214724776*nuVtSqSum[2]*fl[5])*rdv+alphaDrag[2]*(0.3061862178478971*favg[21]+0.1767766952966368*favg[5])+0.1767766952966368*alphaDrag[6]*favg[15]-0.5*fjump[14]+0.1581138830084189*alphaDrag[8]*favg[14]+0.1767766952966368*alphaDrag[0]*favg[14]+0.1767766952966368*alphaDrag[4]*favg[13]; 
  Ghat[15] = ((-1.684024198163434*nuVtSqSum[3]*fr[21])-1.684024198163434*nuVtSqSum[3]*fl[21]+1.185854122563142*nuVtSqSum[9]*fr[15]+1.325825214724776*nuVtSqSum[0]*fr[15]-1.185854122563142*nuVtSqSum[9]*fl[15]-1.325825214724776*nuVtSqSum[0]*fl[15]+1.325825214724776*nuVtSqSum[6]*fr[14]-1.325825214724776*nuVtSqSum[6]*fl[14]+1.325825214724776*nuVtSqSum[5]*fr[13]-1.325825214724776*nuVtSqSum[5]*fl[13]+1.325825214724776*nuVtSqSum[3]*fr[5]-1.325825214724776*nuVtSqSum[3]*fl[5])*rdv+alphaDrag[3]*(0.3061862178478971*favg[21]+0.1767766952966368*favg[5])-0.5*fjump[15]+0.1581138830084189*alphaDrag[9]*favg[15]+0.1767766952966368*alphaDrag[0]*favg[15]+0.1767766952966368*alphaDrag[6]*favg[14]+0.1767766952966368*alphaDrag[5]*favg[13]; 
  Ghat[16] = (1.325825214724776*nuVtSqSum[0]*fr[16]-1.325825214724776*nuVtSqSum[0]*fl[16])*rdv-0.5*fjump[16]+0.1767766952966368*alphaDrag[0]*favg[16]; 
  Ghat[22] = (0.9486832980505137*nuVtSqSum[7]*fr[27]-0.9486832980505137*nuVtSqSum[7]*fl[27]+0.8470386589736728*nuVtSqSum[7]*fr[22]+1.325825214724776*nuVtSqSum[0]*fr[22]-0.8470386589736728*nuVtSqSum[7]*fl[22]-1.325825214724776*nuVtSqSum[0]*fl[22]-1.506237033139206*nuVtSqSum[1]*fr[17]-1.506237033139206*nuVtSqSum[1]*fl[17]+1.185854122563142*nuVtSqSum[5]*fr[8]-1.185854122563142*nuVtSqSum[5]*fl[8]-1.684024198163434*fr[6]*nuVtSqSum[7]-1.684024198163434*fl[6]*nuVtSqSum[7]+1.325825214724776*fr[0]*nuVtSqSum[7]-1.325825214724776*fl[0]*nuVtSqSum[7]+1.185854122563142*nuVtSqSum[4]*fr[7]-1.185854122563142*nuVtSqSum[4]*fl[7]+1.185854122563142*fr[1]*nuVtSqSum[1]-1.185854122563142*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[7]*(0.3952847075210473*favg[27]+0.1129384878631564*favg[22]+0.3061862178478971*favg[6]+0.1767766952966368*favg[0])-0.5*fjump[22]+0.1767766952966368*alphaDrag[0]*favg[22]+alphaDrag[1]*(0.273861278752583*favg[17]+0.1581138830084189*favg[1])+0.1581138830084189*alphaDrag[5]*favg[8]+0.1581138830084189*alphaDrag[4]*favg[7]; 
  Ghat[23] = (0.9486832980505137*nuVtSqSum[8]*fr[27]-0.9486832980505137*nuVtSqSum[8]*fl[27]+0.8470386589736728*nuVtSqSum[8]*fr[23]+1.325825214724776*nuVtSqSum[0]*fr[23]-0.8470386589736728*nuVtSqSum[8]*fl[23]-1.325825214724776*nuVtSqSum[0]*fl[23]-1.506237033139206*nuVtSqSum[2]*fr[18]-1.506237033139206*nuVtSqSum[2]*fl[18]+1.185854122563142*nuVtSqSum[6]*fr[9]-1.185854122563142*nuVtSqSum[6]*fl[9]-1.684024198163434*fr[6]*nuVtSqSum[8]-1.684024198163434*fl[6]*nuVtSqSum[8]+1.325825214724776*fr[0]*nuVtSqSum[8]-1.325825214724776*fl[0]*nuVtSqSum[8]+1.185854122563142*nuVtSqSum[4]*fr[7]-1.185854122563142*nuVtSqSum[4]*fl[7]+1.185854122563142*fr[2]*nuVtSqSum[2]-1.185854122563142*fl[2]*nuVtSqSum[2])*rdv+alphaDrag[8]*(0.3952847075210473*favg[27]+0.1129384878631564*favg[23]+0.3061862178478971*favg[6]+0.1767766952966368*favg[0])-0.5*fjump[23]+0.1767766952966368*alphaDrag[0]*favg[23]+alphaDrag[2]*(0.273861278752583*favg[18]+0.1581138830084189*favg[2])+0.1581138830084189*alphaDrag[6]*favg[9]+0.1581138830084189*alphaDrag[4]*favg[7]; 
  Ghat[24] = (0.9486832980505137*nuVtSqSum[9]*fr[27]-0.9486832980505137*nuVtSqSum[9]*fl[27]+0.8470386589736728*nuVtSqSum[9]*fr[24]+1.325825214724776*nuVtSqSum[0]*fr[24]-0.8470386589736728*nuVtSqSum[9]*fl[24]-1.325825214724776*nuVtSqSum[0]*fl[24]-1.506237033139206*nuVtSqSum[3]*fr[19]-1.506237033139206*nuVtSqSum[3]*fl[19]-1.684024198163434*fr[6]*nuVtSqSum[9]-1.684024198163434*fl[6]*nuVtSqSum[9]+1.325825214724776*fr[0]*nuVtSqSum[9]-1.325825214724776*fl[0]*nuVtSqSum[9]+1.185854122563142*nuVtSqSum[6]*fr[9]-1.185854122563142*nuVtSqSum[6]*fl[9]+1.185854122563142*nuVtSqSum[5]*fr[8]-1.185854122563142*nuVtSqSum[5]*fl[8]+1.185854122563142*fr[3]*nuVtSqSum[3]-1.185854122563142*fl[3]*nuVtSqSum[3])*rdv+alphaDrag[9]*(0.3952847075210473*favg[27]+0.1129384878631564*favg[24]+0.3061862178478971*favg[6]+0.1767766952966368*favg[0])-0.5*fjump[24]+0.1767766952966368*alphaDrag[0]*favg[24]+alphaDrag[3]*(0.273861278752583*favg[19]+0.1581138830084189*favg[3])+0.1581138830084189*alphaDrag[6]*favg[9]+0.1581138830084189*alphaDrag[5]*favg[8]; 
  Ghat[25] = (1.325825214724776*nuVtSqSum[0]*fr[25]-1.325825214724776*nuVtSqSum[0]*fl[25])*rdv-0.5*fjump[25]+0.1767766952966368*alphaDrag[0]*favg[25]; 
  Ghat[26] = (1.325825214724776*nuVtSqSum[0]*fr[26]-1.325825214724776*nuVtSqSum[0]*fl[26])*rdv-0.5*fjump[26]+0.1767766952966368*alphaDrag[0]*favg[26]; 

  double incr1[28]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[0]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = -0.5*Ghat[9]; 
  incr1[10] = -0.5*Ghat[10]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = -0.5*Ghat[13]; 
  incr1[14] = -0.5*Ghat[14]; 
  incr1[15] = -0.5*Ghat[15]; 
  incr1[16] = -0.5*Ghat[16]; 
  incr1[17] = 0.8660254037844386*Ghat[1]; 
  incr1[18] = 0.8660254037844386*Ghat[2]; 
  incr1[19] = 0.8660254037844386*Ghat[3]; 
  incr1[20] = 0.8660254037844386*Ghat[4]; 
  incr1[21] = 0.8660254037844386*Ghat[5]; 
  incr1[22] = -0.5*Ghat[22]; 
  incr1[23] = -0.5*Ghat[23]; 
  incr1[24] = -0.5*Ghat[24]; 
  incr1[25] = -0.5*Ghat[25]; 
  incr1[26] = -0.5*Ghat[26]; 
  incr1[27] = -1.118033988749895*Ghat[0]; 

  double incr2[28]; 
  incr2[6] = nuVtSqSum[0]*(0.1497678868178187*(fr[27]+fl[27])-0.215446597392776*fr[6]+0.215446597392776*fl[6]+0.1530931089239486*(fr[0]+fl[0]))+0.1530931089239486*(nuVtSqSum[9]*(fr[24]+fl[24])+nuVtSqSum[8]*(fr[23]+fl[23])+nuVtSqSum[7]*(fr[22]+fl[22]))+nuVtSqSum[3]*((-0.215446597392776*fr[19])+0.215446597392776*fl[19]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[2]*((-0.215446597392776*fr[18])+0.215446597392776*fl[18]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.215446597392776*fr[17])+0.215446597392776*fl[17]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[6]*(fr[9]+fl[9])+nuVtSqSum[5]*(fr[8]+fl[8])+nuVtSqSum[4]*(fr[7]+fl[7])); 
  incr2[17] = nuVtSqSum[1]*(0.1497678868178187*(fr[27]+fl[27])+0.1369306393762915*(fr[22]+fl[22])-0.215446597392776*fr[6]+0.215446597392776*fl[6]+0.1530931089239486*(fr[0]+fl[0]))+nuVtSqSum[5]*((-0.215446597392776*fr[19])+0.215446597392776*fl[19]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[4]*((-0.215446597392776*fr[18])+0.215446597392776*fl[18]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[7]*((-0.1927012949165104*fr[17])+0.1927012949165104*fl[17]+0.1369306393762915*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.215446597392776*fr[17])+0.215446597392776*fl[17]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[3]*(fr[8]+fl[8])+nuVtSqSum[2]*(fr[7]+fl[7])); 
  incr2[18] = nuVtSqSum[2]*(0.1497678868178187*(fr[27]+fl[27])+0.1369306393762915*(fr[23]+fl[23])-0.215446597392776*fr[6]+0.215446597392776*fl[6]+0.1530931089239486*(fr[0]+fl[0]))+nuVtSqSum[6]*((-0.215446597392776*fr[19])+0.215446597392776*fl[19]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[8]*((-0.1927012949165104*fr[18])+0.1927012949165104*fl[18]+0.1369306393762915*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.215446597392776*fr[18])+0.215446597392776*fl[18]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[4]*((-0.215446597392776*fr[17])+0.215446597392776*fl[17]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[3]*(fr[9]+fl[9])+nuVtSqSum[1]*(fr[7]+fl[7])); 
  incr2[19] = nuVtSqSum[3]*(0.1497678868178187*(fr[27]+fl[27])+0.1369306393762915*(fr[24]+fl[24])-0.215446597392776*fr[6]+0.215446597392776*fl[6]+0.1530931089239486*(fr[0]+fl[0]))+nuVtSqSum[9]*((-0.1927012949165104*fr[19])+0.1927012949165104*fl[19]+0.1369306393762915*(fr[3]+fl[3]))+nuVtSqSum[0]*((-0.215446597392776*fr[19])+0.215446597392776*fl[19]+0.1530931089239486*(fr[3]+fl[3]))+nuVtSqSum[6]*((-0.215446597392776*fr[18])+0.215446597392776*fl[18]+0.1530931089239486*(fr[2]+fl[2]))+nuVtSqSum[5]*((-0.215446597392776*fr[17])+0.215446597392776*fl[17]+0.1530931089239486*(fr[1]+fl[1]))+0.1530931089239486*(nuVtSqSum[2]*(fr[9]+fl[9])+nuVtSqSum[1]*(fr[8]+fl[8])); 
  incr2[20] = nuVtSqSum[0]*((-0.215446597392776*fr[20])+0.215446597392776*fl[20]+0.1530931089239486*(fr[4]+fl[4]))+0.1530931089239486*(nuVtSqSum[3]*(fr[12]+fl[12])+nuVtSqSum[2]*(fr[11]+fl[11])+nuVtSqSum[1]*(fr[10]+fl[10])); 
  incr2[21] = nuVtSqSum[0]*((-0.215446597392776*fr[21])+0.215446597392776*fl[21]+0.1530931089239486*(fr[5]+fl[5]))+0.1530931089239486*(nuVtSqSum[3]*(fr[15]+fl[15])+nuVtSqSum[2]*(fr[14]+fl[14])+nuVtSqSum[1]*(fr[13]+fl[13])); 
  incr2[27] = nuVtSqSum[0]*((-0.5800485314420891*(fr[27]+fl[27]))+0.8344210836992756*fr[6]-0.8344210836992756*fl[6]-0.592927061281571*(fr[0]+fl[0]))-0.592927061281571*(nuVtSqSum[9]*(fr[24]+fl[24])+nuVtSqSum[8]*(fr[23]+fl[23])+nuVtSqSum[7]*(fr[22]+fl[22]))+nuVtSqSum[3]*(0.8344210836992756*fr[19]-0.8344210836992756*fl[19]-0.592927061281571*(fr[3]+fl[3]))+nuVtSqSum[2]*(0.8344210836992756*fr[18]-0.8344210836992756*fl[18]-0.592927061281571*(fr[2]+fl[2]))+nuVtSqSum[1]*(0.8344210836992756*fr[17]-0.8344210836992756*fl[17]-0.592927061281571*(fr[1]+fl[1]))-0.592927061281571*(nuVtSqSum[6]*(fr[9]+fl[9])+nuVtSqSum[5]*(fr[8]+fl[8])+nuVtSqSum[4]*(fr[7]+fl[7])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr1[9]*rdv2R; 
  outr[10] += incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr1[13]*rdv2R; 
  outr[14] += incr1[14]*rdv2R; 
  outr[15] += incr1[15]*rdv2R; 
  outr[16] += incr1[16]*rdv2R; 
  outr[17] += incr2[17]*rdvSq4R+incr1[17]*rdv2R; 
  outr[18] += incr2[18]*rdvSq4R+incr1[18]*rdv2R; 
  outr[19] += incr2[19]*rdvSq4R+incr1[19]*rdv2R; 
  outr[20] += incr2[20]*rdvSq4R+incr1[20]*rdv2R; 
  outr[21] += incr2[21]*rdvSq4R+incr1[21]*rdv2R; 
  outr[22] += incr1[22]*rdv2R; 
  outr[23] += incr1[23]*rdv2R; 
  outr[24] += incr1[24]*rdv2R; 
  outr[25] += incr1[25]*rdv2R; 
  outr[26] += incr1[26]*rdv2R; 
  outr[27] += incr2[27]*rdvSq4R+incr1[27]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += -1.0*incr1[9]*rdv2L; 
  outl[10] += -1.0*incr1[10]*rdv2L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += -1.0*incr1[13]*rdv2L; 
  outl[14] += -1.0*incr1[14]*rdv2L; 
  outl[15] += -1.0*incr1[15]*rdv2L; 
  outl[16] += -1.0*incr1[16]*rdv2L; 
  outl[17] += incr1[17]*rdv2L-1.0*incr2[17]*rdvSq4L; 
  outl[18] += incr1[18]*rdv2L-1.0*incr2[18]*rdvSq4L; 
  outl[19] += incr1[19]*rdv2L-1.0*incr2[19]*rdvSq4L; 
  outl[20] += incr1[20]*rdv2L-1.0*incr2[20]*rdvSq4L; 
  outl[21] += incr1[21]*rdv2L-1.0*incr2[21]*rdvSq4L; 
  outl[22] += -1.0*incr1[22]*rdv2L; 
  outl[23] += -1.0*incr1[23]*rdv2L; 
  outl[24] += -1.0*incr1[24]*rdv2L; 
  outl[25] += -1.0*incr1[25]*rdv2L; 
  outl[26] += -1.0*incr1[26]*rdv2L; 
  outl[27] += incr2[27]*rdvSq4L-1.0*incr1[27]*rdv2L; 

  return std::abs((0.3952847075210473*sumNuUz[9])/nuSum+(0.3952847075210473*sumNuUz[8])/nuSum+(0.3952847075210473*sumNuUz[7])/nuSum-(0.3535533905932737*sumNuUz[0])/nuSum+wl[5]); 
} 
