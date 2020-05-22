#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[9]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[15]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 

  double fjump[15]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(-1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(-1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(1*fr[14])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 

  double Ghat[15]; 
  for(unsigned short int i=0; i<15; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*nuVtSqSum[0]*fr[12]-1.897366596101028*nuVtSqSum[0]*fl[12]+2.651650429449552*nuVtSqSum[2]*fr[11]-2.651650429449552*nuVtSqSum[2]*fl[11]-3.368048396326869*nuVtSqSum[1]*fr[5]-3.368048396326869*nuVtSqSum[1]*fl[5]-3.368048396326869*nuVtSqSum[0]*fr[2]-3.368048396326869*nuVtSqSum[0]*fl[2]+2.651650429449552*fr[1]*nuVtSqSum[1]-2.651650429449552*fl[1]*nuVtSqSum[1]+2.651650429449552*fr[0]*nuVtSqSum[0]-2.651650429449552*fl[0]*nuVtSqSum[0])*rdv-1.118033988749895*fjump[12]+alphaDrag[0]*(0.7905694150420947*favg[12]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[2]*favg[11]+alphaDrag[1]*(0.6123724356957944*favg[5]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (1.897366596101028*nuVtSqSum[1]*fr[12]-1.897366596101028*nuVtSqSum[1]*fl[12]+2.371708245126284*nuVtSqSum[1]*fr[11]-2.371708245126284*nuVtSqSum[1]*fl[11]-3.012474066278413*nuVtSqSum[2]*fr[5]-3.368048396326869*nuVtSqSum[0]*fr[5]-3.012474066278413*nuVtSqSum[2]*fl[5]-3.368048396326869*nuVtSqSum[0]*fl[5]+2.371708245126284*fr[1]*nuVtSqSum[2]-2.371708245126284*fl[1]*nuVtSqSum[2]-3.368048396326869*nuVtSqSum[1]*fr[2]-3.368048396326869*nuVtSqSum[1]*fl[2]+2.651650429449552*fr[0]*nuVtSqSum[1]-2.651650429449552*fl[0]*nuVtSqSum[1]+2.651650429449552*nuVtSqSum[0]*fr[1]-2.651650429449552*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.7905694150420947*favg[12]+0.3162277660168379*favg[11]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[5]+alphaDrag[0]*(0.6123724356957944*favg[5]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.5477225575051661*favg[5]+0.3162277660168379*favg[1])-0.5*fjump[1]; 
  Ghat[3] = ((-3.368048396326869*nuVtSqSum[0]*fr[7])-3.368048396326869*nuVtSqSum[0]*fl[7]+2.651650429449552*nuVtSqSum[1]*fr[6]-2.651650429449552*nuVtSqSum[1]*fl[6]+2.651650429449552*nuVtSqSum[0]*fr[3]-2.651650429449552*nuVtSqSum[0]*fl[3])*rdv-0.8660254037844386*fjump[7]+alphaDrag[0]*(0.6123724356957944*favg[7]+0.3535533905932737*favg[3])+0.3535533905932737*alphaDrag[1]*favg[6]-0.5*fjump[3]; 
  Ghat[4] = ((-3.368048396326869*nuVtSqSum[0]*fr[9])-3.368048396326869*nuVtSqSum[0]*fl[9]+2.651650429449552*nuVtSqSum[1]*fr[8]-2.651650429449552*nuVtSqSum[1]*fl[8]+2.651650429449552*nuVtSqSum[0]*fr[4]-2.651650429449552*nuVtSqSum[0]*fl[4])*rdv-0.8660254037844386*fjump[9]+alphaDrag[0]*(0.6123724356957944*favg[9]+0.3535533905932737*favg[4])+0.3535533905932737*alphaDrag[1]*favg[8]-0.5*fjump[4]; 
  Ghat[6] = ((-3.368048396326869*nuVtSqSum[1]*fr[7])-3.368048396326869*nuVtSqSum[1]*fl[7]+2.371708245126284*nuVtSqSum[2]*fr[6]+2.651650429449552*nuVtSqSum[0]*fr[6]-2.371708245126284*nuVtSqSum[2]*fl[6]-2.651650429449552*nuVtSqSum[0]*fl[6]+2.651650429449552*nuVtSqSum[1]*fr[3]-2.651650429449552*nuVtSqSum[1]*fl[3])*rdv+alphaDrag[1]*(0.6123724356957944*favg[7]+0.3535533905932737*favg[3])-0.5*fjump[6]+0.3162277660168379*alphaDrag[2]*favg[6]+0.3535533905932737*alphaDrag[0]*favg[6]; 
  Ghat[8] = ((-3.368048396326869*nuVtSqSum[1]*fr[9])-3.368048396326869*nuVtSqSum[1]*fl[9]+2.371708245126284*nuVtSqSum[2]*fr[8]+2.651650429449552*nuVtSqSum[0]*fr[8]-2.371708245126284*nuVtSqSum[2]*fl[8]-2.651650429449552*nuVtSqSum[0]*fl[8]+2.651650429449552*nuVtSqSum[1]*fr[4]-2.651650429449552*nuVtSqSum[1]*fl[4])*rdv+alphaDrag[1]*(0.6123724356957944*favg[9]+0.3535533905932737*favg[4])-0.5*fjump[8]+0.3162277660168379*alphaDrag[2]*favg[8]+0.3535533905932737*alphaDrag[0]*favg[8]; 
  Ghat[10] = (2.651650429449552*nuVtSqSum[0]*fr[10]-2.651650429449552*nuVtSqSum[0]*fl[10])*rdv-0.5*fjump[10]+0.3535533905932737*alphaDrag[0]*favg[10]; 
  Ghat[11] = (1.897366596101028*nuVtSqSum[2]*fr[12]-1.897366596101028*nuVtSqSum[2]*fl[12]+1.694077317947346*nuVtSqSum[2]*fr[11]+2.651650429449552*nuVtSqSum[0]*fr[11]-1.694077317947346*nuVtSqSum[2]*fl[11]-2.651650429449552*nuVtSqSum[0]*fl[11]-3.012474066278413*nuVtSqSum[1]*fr[5]-3.012474066278413*nuVtSqSum[1]*fl[5]-3.368048396326869*fr[2]*nuVtSqSum[2]-3.368048396326869*fl[2]*nuVtSqSum[2]+2.651650429449552*fr[0]*nuVtSqSum[2]-2.651650429449552*fl[0]*nuVtSqSum[2]+2.371708245126284*fr[1]*nuVtSqSum[1]-2.371708245126284*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[2]*(0.7905694150420947*favg[12]+0.2258769757263128*favg[11]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[11]+0.3535533905932737*alphaDrag[0]*favg[11]+alphaDrag[1]*(0.5477225575051661*favg[5]+0.3162277660168379*favg[1]); 
  Ghat[13] = (2.651650429449552*nuVtSqSum[0]*fr[13]-2.651650429449552*nuVtSqSum[0]*fl[13])*rdv-0.5*fjump[13]+0.3535533905932737*alphaDrag[0]*favg[13]; 
  Ghat[14] = (2.651650429449552*nuVtSqSum[0]*fr[14]-2.651650429449552*nuVtSqSum[0]*fl[14])*rdv-0.5*fjump[14]+0.3535533905932737*alphaDrag[0]*favg[14]; 

  double incr1[15]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = -0.5*Ghat[6]; 
  incr1[7] = 0.8660254037844386*Ghat[3]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = 0.8660254037844386*Ghat[4]; 
  incr1[10] = -0.5*Ghat[10]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = -1.118033988749895*Ghat[0]; 
  incr1[13] = -0.5*Ghat[13]; 
  incr1[14] = -0.5*Ghat[14]; 

  double incr2[15]; 
  incr2[2] = nuVtSqSum[0]*(0.2995357736356374*(fr[12]+fl[12])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*nuVtSqSum[2]*(fr[11]+fl[11])+nuVtSqSum[1]*((-0.430893194785552*fr[5])+0.430893194785552*fl[5]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[5] = nuVtSqSum[1]*(0.2995357736356374*(fr[12]+fl[12])+0.273861278752583*(fr[11]+fl[11])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[2]*((-0.3854025898330209*fr[5])+0.3854025898330209*fl[5]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.430893194785552*fr[5])+0.430893194785552*fl[5]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[7] = nuVtSqSum[0]*((-0.430893194785552*fr[7])+0.430893194785552*fl[7]+0.3061862178478971*(fr[3]+fl[3]))+0.3061862178478971*nuVtSqSum[1]*(fr[6]+fl[6]); 
  incr2[9] = nuVtSqSum[0]*((-0.430893194785552*fr[9])+0.430893194785552*fl[9]+0.3061862178478971*(fr[4]+fl[4]))+0.3061862178478971*nuVtSqSum[1]*(fr[8]+fl[8]); 
  incr2[12] = nuVtSqSum[0]*((-1.160097062884178*(fr[12]+fl[12]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*nuVtSqSum[2]*(fr[11]+fl[11])+nuVtSqSum[1]*(1.668842167398551*fr[5]-1.668842167398551*fl[5]-1.185854122563142*(fr[1]+fl[1])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr2[9]*rdvSq4R+incr1[9]*rdv2R; 
  outr[10] += incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr2[12]*rdvSq4R+incr1[12]*rdv2R; 
  outr[13] += incr1[13]*rdv2R; 
  outr[14] += incr1[14]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += incr1[9]*rdv2L-1.0*incr2[9]*rdvSq4L; 
  outl[10] += -1.0*incr1[10]*rdv2L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += incr2[12]*rdvSq4L-1.0*incr1[12]*rdv2L; 
  outl[13] += -1.0*incr1[13]*rdv2L; 
  outl[14] += -1.0*incr1[14]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUx[2])/nuSum-(0.7071067811865475*sumNuUx[0])/nuSum+wl[1]); 
} 
double VmLBOconstNuSurf1x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[9]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUy = &nuUSum[3]; 

  double favg[15]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 

  double fjump[15]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(-1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(1*fr[14])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[2]*nuSum+0.7071067811865475*dxvl[2]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 

  double Ghat[15]; 
  for(unsigned short int i=0; i<15; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*nuVtSqSum[0]*fr[13]-1.897366596101028*nuVtSqSum[0]*fl[13]+2.651650429449552*nuVtSqSum[2]*fr[11]-2.651650429449552*nuVtSqSum[2]*fl[11]-3.368048396326869*nuVtSqSum[1]*fr[6]-3.368048396326869*nuVtSqSum[1]*fl[6]-3.368048396326869*nuVtSqSum[0]*fr[3]-3.368048396326869*nuVtSqSum[0]*fl[3]+2.651650429449552*fr[1]*nuVtSqSum[1]-2.651650429449552*fl[1]*nuVtSqSum[1]+2.651650429449552*fr[0]*nuVtSqSum[0]-2.651650429449552*fl[0]*nuVtSqSum[0])*rdv-1.118033988749895*fjump[13]+alphaDrag[0]*(0.7905694150420947*favg[13]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[2]*favg[11]+alphaDrag[1]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = (1.897366596101028*nuVtSqSum[1]*fr[13]-1.897366596101028*nuVtSqSum[1]*fl[13]+2.371708245126284*nuVtSqSum[1]*fr[11]-2.371708245126284*nuVtSqSum[1]*fl[11]-3.012474066278413*nuVtSqSum[2]*fr[6]-3.368048396326869*nuVtSqSum[0]*fr[6]-3.012474066278413*nuVtSqSum[2]*fl[6]-3.368048396326869*nuVtSqSum[0]*fl[6]-3.368048396326869*nuVtSqSum[1]*fr[3]-3.368048396326869*nuVtSqSum[1]*fl[3]+2.371708245126284*fr[1]*nuVtSqSum[2]-2.371708245126284*fl[1]*nuVtSqSum[2]+2.651650429449552*fr[0]*nuVtSqSum[1]-2.651650429449552*fl[0]*nuVtSqSum[1]+2.651650429449552*nuVtSqSum[0]*fr[1]-2.651650429449552*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.7905694150420947*favg[13]+0.3162277660168379*favg[11]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[6]+alphaDrag[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.5477225575051661*favg[6]+0.3162277660168379*favg[1])-0.5*fjump[1]; 
  Ghat[2] = ((-3.368048396326869*nuVtSqSum[0]*fr[7])-3.368048396326869*nuVtSqSum[0]*fl[7]+2.651650429449552*nuVtSqSum[1]*fr[5]-2.651650429449552*nuVtSqSum[1]*fl[5]+2.651650429449552*nuVtSqSum[0]*fr[2]-2.651650429449552*nuVtSqSum[0]*fl[2])*rdv-0.8660254037844386*fjump[7]+alphaDrag[0]*(0.6123724356957944*favg[7]+0.3535533905932737*favg[2])+0.3535533905932737*alphaDrag[1]*favg[5]-0.5*fjump[2]; 
  Ghat[4] = ((-3.368048396326869*nuVtSqSum[0]*fr[10])-3.368048396326869*nuVtSqSum[0]*fl[10]+2.651650429449552*nuVtSqSum[1]*fr[8]-2.651650429449552*nuVtSqSum[1]*fl[8]+2.651650429449552*nuVtSqSum[0]*fr[4]-2.651650429449552*nuVtSqSum[0]*fl[4])*rdv-0.8660254037844386*fjump[10]+alphaDrag[0]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[4])+0.3535533905932737*alphaDrag[1]*favg[8]-0.5*fjump[4]; 
  Ghat[5] = ((-3.368048396326869*nuVtSqSum[1]*fr[7])-3.368048396326869*nuVtSqSum[1]*fl[7]+2.371708245126284*nuVtSqSum[2]*fr[5]+2.651650429449552*nuVtSqSum[0]*fr[5]-2.371708245126284*nuVtSqSum[2]*fl[5]-2.651650429449552*nuVtSqSum[0]*fl[5]+2.651650429449552*nuVtSqSum[1]*fr[2]-2.651650429449552*nuVtSqSum[1]*fl[2])*rdv+alphaDrag[1]*(0.6123724356957944*favg[7]+0.3535533905932737*favg[2])-0.5*fjump[5]+0.3162277660168379*alphaDrag[2]*favg[5]+0.3535533905932737*alphaDrag[0]*favg[5]; 
  Ghat[8] = ((-3.368048396326869*nuVtSqSum[1]*fr[10])-3.368048396326869*nuVtSqSum[1]*fl[10]+2.371708245126284*nuVtSqSum[2]*fr[8]+2.651650429449552*nuVtSqSum[0]*fr[8]-2.371708245126284*nuVtSqSum[2]*fl[8]-2.651650429449552*nuVtSqSum[0]*fl[8]+2.651650429449552*nuVtSqSum[1]*fr[4]-2.651650429449552*nuVtSqSum[1]*fl[4])*rdv+alphaDrag[1]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[4])-0.5*fjump[8]+0.3162277660168379*alphaDrag[2]*favg[8]+0.3535533905932737*alphaDrag[0]*favg[8]; 
  Ghat[9] = (2.651650429449552*nuVtSqSum[0]*fr[9]-2.651650429449552*nuVtSqSum[0]*fl[9])*rdv-0.5*fjump[9]+0.3535533905932737*alphaDrag[0]*favg[9]; 
  Ghat[11] = (1.897366596101028*nuVtSqSum[2]*fr[13]-1.897366596101028*nuVtSqSum[2]*fl[13]+1.694077317947346*nuVtSqSum[2]*fr[11]+2.651650429449552*nuVtSqSum[0]*fr[11]-1.694077317947346*nuVtSqSum[2]*fl[11]-2.651650429449552*nuVtSqSum[0]*fl[11]-3.012474066278413*nuVtSqSum[1]*fr[6]-3.012474066278413*nuVtSqSum[1]*fl[6]-3.368048396326869*nuVtSqSum[2]*fr[3]-3.368048396326869*nuVtSqSum[2]*fl[3]+2.651650429449552*fr[0]*nuVtSqSum[2]-2.651650429449552*fl[0]*nuVtSqSum[2]+2.371708245126284*fr[1]*nuVtSqSum[1]-2.371708245126284*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[2]*(0.7905694150420947*favg[13]+0.2258769757263128*favg[11]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.5*fjump[11]+0.3535533905932737*alphaDrag[0]*favg[11]+alphaDrag[1]*(0.5477225575051661*favg[6]+0.3162277660168379*favg[1]); 
  Ghat[12] = (2.651650429449552*nuVtSqSum[0]*fr[12]-2.651650429449552*nuVtSqSum[0]*fl[12])*rdv-0.5*fjump[12]+0.3535533905932737*alphaDrag[0]*favg[12]; 
  Ghat[14] = (2.651650429449552*nuVtSqSum[0]*fr[14]-2.651650429449552*nuVtSqSum[0]*fl[14])*rdv-0.5*fjump[14]+0.3535533905932737*alphaDrag[0]*favg[14]; 

  double incr1[15]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[1]; 
  incr1[7] = 0.8660254037844386*Ghat[2]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = -0.5*Ghat[9]; 
  incr1[10] = 0.8660254037844386*Ghat[4]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = -1.118033988749895*Ghat[0]; 
  incr1[14] = -0.5*Ghat[14]; 

  double incr2[15]; 
  incr2[3] = nuVtSqSum[0]*(0.2995357736356374*(fr[13]+fl[13])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*nuVtSqSum[2]*(fr[11]+fl[11])+nuVtSqSum[1]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[6] = nuVtSqSum[1]*(0.2995357736356374*(fr[13]+fl[13])+0.273861278752583*(fr[11]+fl[11])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[2]*((-0.3854025898330209*fr[6])+0.3854025898330209*fl[6]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[7] = nuVtSqSum[0]*((-0.430893194785552*fr[7])+0.430893194785552*fl[7]+0.3061862178478971*(fr[2]+fl[2]))+0.3061862178478971*nuVtSqSum[1]*(fr[5]+fl[5]); 
  incr2[10] = nuVtSqSum[0]*((-0.430893194785552*fr[10])+0.430893194785552*fl[10]+0.3061862178478971*(fr[4]+fl[4]))+0.3061862178478971*nuVtSqSum[1]*(fr[8]+fl[8]); 
  incr2[13] = nuVtSqSum[0]*((-1.160097062884178*(fr[13]+fl[13]))+1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*nuVtSqSum[2]*(fr[11]+fl[11])+nuVtSqSum[1]*(1.668842167398551*fr[6]-1.668842167398551*fl[6]-1.185854122563142*(fr[1]+fl[1])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr1[9]*rdv2R; 
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr1[14]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += -1.0*incr1[9]*rdv2L; 
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += incr2[13]*rdvSq4L-1.0*incr1[13]*rdv2L; 
  outl[14] += -1.0*incr1[14]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUy[2])/nuSum-(0.7071067811865475*sumNuUy[0])/nuSum+wl[2]); 
} 
double VmLBOconstNuSurf1x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[9]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[3]; 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  const double *sumNuUz = &nuUSum[6]; 

  double favg[15]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = -1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 

  double fjump[15]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(-1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(-1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(1*fr[14])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[3]*nuSum+0.7071067811865475*dxvl[3]*nuSum-1.0*sumNuUz[0]; 
  alphaDrag[1] = -1.0*sumNuUz[1]; 
  alphaDrag[2] = -1.0*sumNuUz[2]; 

  double Ghat[15]; 
  for(unsigned short int i=0; i<15; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*nuVtSqSum[0]*fr[14]-1.897366596101028*nuVtSqSum[0]*fl[14]+2.651650429449552*nuVtSqSum[2]*fr[11]-2.651650429449552*nuVtSqSum[2]*fl[11]-3.368048396326869*nuVtSqSum[1]*fr[8]-3.368048396326869*nuVtSqSum[1]*fl[8]-3.368048396326869*nuVtSqSum[0]*fr[4]-3.368048396326869*nuVtSqSum[0]*fl[4]+2.651650429449552*fr[1]*nuVtSqSum[1]-2.651650429449552*fl[1]*nuVtSqSum[1]+2.651650429449552*fr[0]*nuVtSqSum[0]-2.651650429449552*fl[0]*nuVtSqSum[0])*rdv-1.118033988749895*fjump[14]+alphaDrag[0]*(0.7905694150420947*favg[14]+0.6123724356957944*favg[4]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[2]*favg[11]+alphaDrag[1]*(0.6123724356957944*favg[8]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[4]-0.5*fjump[0]; 
  Ghat[1] = (1.897366596101028*nuVtSqSum[1]*fr[14]-1.897366596101028*nuVtSqSum[1]*fl[14]+2.371708245126284*nuVtSqSum[1]*fr[11]-2.371708245126284*nuVtSqSum[1]*fl[11]-3.012474066278413*nuVtSqSum[2]*fr[8]-3.368048396326869*nuVtSqSum[0]*fr[8]-3.012474066278413*nuVtSqSum[2]*fl[8]-3.368048396326869*nuVtSqSum[0]*fl[8]-3.368048396326869*nuVtSqSum[1]*fr[4]-3.368048396326869*nuVtSqSum[1]*fl[4]+2.371708245126284*fr[1]*nuVtSqSum[2]-2.371708245126284*fl[1]*nuVtSqSum[2]+2.651650429449552*fr[0]*nuVtSqSum[1]-2.651650429449552*fl[0]*nuVtSqSum[1]+2.651650429449552*nuVtSqSum[0]*fr[1]-2.651650429449552*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.7905694150420947*favg[14]+0.3162277660168379*favg[11]+0.6123724356957944*favg[4]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[8]+alphaDrag[0]*(0.6123724356957944*favg[8]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.5477225575051661*favg[8]+0.3162277660168379*favg[1])-0.5*fjump[1]; 
  Ghat[2] = ((-3.368048396326869*nuVtSqSum[0]*fr[9])-3.368048396326869*nuVtSqSum[0]*fl[9]+2.651650429449552*nuVtSqSum[1]*fr[5]-2.651650429449552*nuVtSqSum[1]*fl[5]+2.651650429449552*nuVtSqSum[0]*fr[2]-2.651650429449552*nuVtSqSum[0]*fl[2])*rdv-0.8660254037844386*fjump[9]+alphaDrag[0]*(0.6123724356957944*favg[9]+0.3535533905932737*favg[2])+0.3535533905932737*alphaDrag[1]*favg[5]-0.5*fjump[2]; 
  Ghat[3] = ((-3.368048396326869*nuVtSqSum[0]*fr[10])-3.368048396326869*nuVtSqSum[0]*fl[10]+2.651650429449552*nuVtSqSum[1]*fr[6]-2.651650429449552*nuVtSqSum[1]*fl[6]+2.651650429449552*nuVtSqSum[0]*fr[3]-2.651650429449552*nuVtSqSum[0]*fl[3])*rdv-0.8660254037844386*fjump[10]+alphaDrag[0]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[3])+0.3535533905932737*alphaDrag[1]*favg[6]-0.5*fjump[3]; 
  Ghat[5] = ((-3.368048396326869*nuVtSqSum[1]*fr[9])-3.368048396326869*nuVtSqSum[1]*fl[9]+2.371708245126284*nuVtSqSum[2]*fr[5]+2.651650429449552*nuVtSqSum[0]*fr[5]-2.371708245126284*nuVtSqSum[2]*fl[5]-2.651650429449552*nuVtSqSum[0]*fl[5]+2.651650429449552*nuVtSqSum[1]*fr[2]-2.651650429449552*nuVtSqSum[1]*fl[2])*rdv+alphaDrag[1]*(0.6123724356957944*favg[9]+0.3535533905932737*favg[2])-0.5*fjump[5]+0.3162277660168379*alphaDrag[2]*favg[5]+0.3535533905932737*alphaDrag[0]*favg[5]; 
  Ghat[6] = ((-3.368048396326869*nuVtSqSum[1]*fr[10])-3.368048396326869*nuVtSqSum[1]*fl[10]+2.371708245126284*nuVtSqSum[2]*fr[6]+2.651650429449552*nuVtSqSum[0]*fr[6]-2.371708245126284*nuVtSqSum[2]*fl[6]-2.651650429449552*nuVtSqSum[0]*fl[6]+2.651650429449552*nuVtSqSum[1]*fr[3]-2.651650429449552*nuVtSqSum[1]*fl[3])*rdv+alphaDrag[1]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[3])-0.5*fjump[6]+0.3162277660168379*alphaDrag[2]*favg[6]+0.3535533905932737*alphaDrag[0]*favg[6]; 
  Ghat[7] = (2.651650429449552*nuVtSqSum[0]*fr[7]-2.651650429449552*nuVtSqSum[0]*fl[7])*rdv-0.5*fjump[7]+0.3535533905932737*alphaDrag[0]*favg[7]; 
  Ghat[11] = (1.897366596101028*nuVtSqSum[2]*fr[14]-1.897366596101028*nuVtSqSum[2]*fl[14]+1.694077317947346*nuVtSqSum[2]*fr[11]+2.651650429449552*nuVtSqSum[0]*fr[11]-1.694077317947346*nuVtSqSum[2]*fl[11]-2.651650429449552*nuVtSqSum[0]*fl[11]-3.012474066278413*nuVtSqSum[1]*fr[8]-3.012474066278413*nuVtSqSum[1]*fl[8]-3.368048396326869*nuVtSqSum[2]*fr[4]-3.368048396326869*nuVtSqSum[2]*fl[4]+2.651650429449552*fr[0]*nuVtSqSum[2]-2.651650429449552*fl[0]*nuVtSqSum[2]+2.371708245126284*fr[1]*nuVtSqSum[1]-2.371708245126284*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[2]*(0.7905694150420947*favg[14]+0.2258769757263128*favg[11]+0.6123724356957944*favg[4]+0.3535533905932737*favg[0])-0.5*fjump[11]+0.3535533905932737*alphaDrag[0]*favg[11]+alphaDrag[1]*(0.5477225575051661*favg[8]+0.3162277660168379*favg[1]); 
  Ghat[12] = (2.651650429449552*nuVtSqSum[0]*fr[12]-2.651650429449552*nuVtSqSum[0]*fl[12])*rdv-0.5*fjump[12]+0.3535533905932737*alphaDrag[0]*favg[12]; 
  Ghat[13] = (2.651650429449552*nuVtSqSum[0]*fr[13]-2.651650429449552*nuVtSqSum[0]*fl[13])*rdv-0.5*fjump[13]+0.3535533905932737*alphaDrag[0]*favg[13]; 

  double incr1[15]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[0]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = -0.5*Ghat[6]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = 0.8660254037844386*Ghat[1]; 
  incr1[9] = 0.8660254037844386*Ghat[2]; 
  incr1[10] = 0.8660254037844386*Ghat[3]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = -0.5*Ghat[13]; 
  incr1[14] = -1.118033988749895*Ghat[0]; 

  double incr2[15]; 
  incr2[4] = nuVtSqSum[0]*(0.2995357736356374*(fr[14]+fl[14])-0.430893194785552*fr[4]+0.430893194785552*fl[4]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*nuVtSqSum[2]*(fr[11]+fl[11])+nuVtSqSum[1]*((-0.430893194785552*fr[8])+0.430893194785552*fl[8]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[8] = nuVtSqSum[1]*(0.2995357736356374*(fr[14]+fl[14])+0.273861278752583*(fr[11]+fl[11])-0.430893194785552*fr[4]+0.430893194785552*fl[4]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[2]*((-0.3854025898330209*fr[8])+0.3854025898330209*fl[8]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.430893194785552*fr[8])+0.430893194785552*fl[8]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[9] = nuVtSqSum[0]*((-0.430893194785552*fr[9])+0.430893194785552*fl[9]+0.3061862178478971*(fr[2]+fl[2]))+0.3061862178478971*nuVtSqSum[1]*(fr[5]+fl[5]); 
  incr2[10] = nuVtSqSum[0]*((-0.430893194785552*fr[10])+0.430893194785552*fl[10]+0.3061862178478971*(fr[3]+fl[3]))+0.3061862178478971*nuVtSqSum[1]*(fr[6]+fl[6]); 
  incr2[14] = nuVtSqSum[0]*((-1.160097062884178*(fr[14]+fl[14]))+1.668842167398551*fr[4]-1.668842167398551*fl[4]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*nuVtSqSum[2]*(fr[11]+fl[11])+nuVtSqSum[1]*(1.668842167398551*fr[8]-1.668842167398551*fl[8]-1.185854122563142*(fr[1]+fl[1])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr2[8]*rdvSq4R+incr1[8]*rdv2R; 
  outr[9] += incr2[9]*rdvSq4R+incr1[9]*rdv2R; 
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += incr1[8]*rdv2L-1.0*incr2[8]*rdvSq4L; 
  outl[9] += incr1[9]*rdv2L-1.0*incr2[9]*rdvSq4L; 
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += -1.0*incr1[13]*rdv2L; 
  outl[14] += incr2[14]*rdvSq4L-1.0*incr1[14]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUz[2])/nuSum-(0.7071067811865475*sumNuUz[0])/nuSum+wl[3]); 
} 
