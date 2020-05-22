#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[10]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 

  double fjump[10]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(1*fr[9])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 

  double Ghat[10]; 
  for(unsigned short int i=0; i<10; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*nuVtSqSum[0]*fr[8]-1.897366596101028*nuVtSqSum[0]*fl[8]+2.651650429449552*nuVtSqSum[2]*fr[7]-2.651650429449552*nuVtSqSum[2]*fl[7]-3.368048396326869*nuVtSqSum[1]*fr[4]-3.368048396326869*nuVtSqSum[1]*fl[4]-3.368048396326869*nuVtSqSum[0]*fr[2]-3.368048396326869*nuVtSqSum[0]*fl[2]+2.651650429449552*fr[1]*nuVtSqSum[1]-2.651650429449552*fl[1]*nuVtSqSum[1]+2.651650429449552*fr[0]*nuVtSqSum[0]-2.651650429449552*fl[0]*nuVtSqSum[0])*rdv-1.118033988749895*fjump[8]+alphaDrag[0]*(0.7905694150420947*favg[8]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[2]*favg[7]+alphaDrag[1]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (1.897366596101028*nuVtSqSum[1]*fr[8]-1.897366596101028*nuVtSqSum[1]*fl[8]+2.371708245126284*nuVtSqSum[1]*fr[7]-2.371708245126284*nuVtSqSum[1]*fl[7]-3.012474066278413*nuVtSqSum[2]*fr[4]-3.368048396326869*nuVtSqSum[0]*fr[4]-3.012474066278413*nuVtSqSum[2]*fl[4]-3.368048396326869*nuVtSqSum[0]*fl[4]+2.371708245126284*fr[1]*nuVtSqSum[2]-2.371708245126284*fl[1]*nuVtSqSum[2]-3.368048396326869*nuVtSqSum[1]*fr[2]-3.368048396326869*nuVtSqSum[1]*fl[2]+2.651650429449552*fr[0]*nuVtSqSum[1]-2.651650429449552*fl[0]*nuVtSqSum[1]+2.651650429449552*nuVtSqSum[0]*fr[1]-2.651650429449552*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.7905694150420947*favg[8]+0.3162277660168379*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[4]+alphaDrag[0]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.5477225575051661*favg[4]+0.3162277660168379*favg[1])-0.5*fjump[1]; 
  Ghat[3] = ((-3.368048396326869*nuVtSqSum[0]*fr[6])-3.368048396326869*nuVtSqSum[0]*fl[6]+2.651650429449552*nuVtSqSum[1]*fr[5]-2.651650429449552*nuVtSqSum[1]*fl[5]+2.651650429449552*nuVtSqSum[0]*fr[3]-2.651650429449552*nuVtSqSum[0]*fl[3])*rdv-0.8660254037844386*fjump[6]+alphaDrag[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[3])+0.3535533905932737*alphaDrag[1]*favg[5]-0.5*fjump[3]; 
  Ghat[5] = ((-3.368048396326869*nuVtSqSum[1]*fr[6])-3.368048396326869*nuVtSqSum[1]*fl[6]+2.371708245126284*nuVtSqSum[2]*fr[5]+2.651650429449552*nuVtSqSum[0]*fr[5]-2.371708245126284*nuVtSqSum[2]*fl[5]-2.651650429449552*nuVtSqSum[0]*fl[5]+2.651650429449552*nuVtSqSum[1]*fr[3]-2.651650429449552*nuVtSqSum[1]*fl[3])*rdv+alphaDrag[1]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[5]+0.3162277660168379*alphaDrag[2]*favg[5]+0.3535533905932737*alphaDrag[0]*favg[5]; 
  Ghat[7] = (1.897366596101028*nuVtSqSum[2]*fr[8]-1.897366596101028*nuVtSqSum[2]*fl[8]+1.694077317947346*nuVtSqSum[2]*fr[7]+2.651650429449552*nuVtSqSum[0]*fr[7]-1.694077317947346*nuVtSqSum[2]*fl[7]-2.651650429449552*nuVtSqSum[0]*fl[7]-3.012474066278413*nuVtSqSum[1]*fr[4]-3.012474066278413*nuVtSqSum[1]*fl[4]-3.368048396326869*fr[2]*nuVtSqSum[2]-3.368048396326869*fl[2]*nuVtSqSum[2]+2.651650429449552*fr[0]*nuVtSqSum[2]-2.651650429449552*fl[0]*nuVtSqSum[2]+2.371708245126284*fr[1]*nuVtSqSum[1]-2.371708245126284*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[2]*(0.7905694150420947*favg[8]+0.2258769757263128*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[7]+0.3535533905932737*alphaDrag[0]*favg[7]+alphaDrag[1]*(0.5477225575051661*favg[4]+0.3162277660168379*favg[1]); 
  Ghat[9] = (2.651650429449552*nuVtSqSum[0]*fr[9]-2.651650429449552*nuVtSqSum[0]*fl[9])*rdv-0.5*fjump[9]+0.3535533905932737*alphaDrag[0]*favg[9]; 

  double incr1[10]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[1]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[3]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = -1.118033988749895*Ghat[0]; 
  incr1[9] = -0.5*Ghat[9]; 

  double incr2[10]; 
  incr2[2] = nuVtSqSum[0]*(0.2995357736356374*(fr[8]+fl[8])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*nuVtSqSum[2]*(fr[7]+fl[7])+nuVtSqSum[1]*((-0.430893194785552*fr[4])+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[4] = nuVtSqSum[1]*(0.2995357736356374*(fr[8]+fl[8])+0.273861278752583*(fr[7]+fl[7])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[2]*((-0.3854025898330209*fr[4])+0.3854025898330209*fl[4]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.430893194785552*fr[4])+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[6] = nuVtSqSum[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3]))+0.3061862178478971*nuVtSqSum[1]*(fr[5]+fl[5]); 
  incr2[8] = nuVtSqSum[0]*((-1.160097062884178*(fr[8]+fl[8]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*nuVtSqSum[2]*(fr[7]+fl[7])+nuVtSqSum[1]*(1.668842167398551*fr[4]-1.668842167398551*fl[4]-1.185854122563142*(fr[1]+fl[1])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr2[8]*rdvSq4R+incr1[8]*rdv2R; 
  outr[9] += incr1[9]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += incr2[8]*rdvSq4L-1.0*incr1[8]*rdv2L; 
  outl[9] += -1.0*incr1[9]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUx[2])/nuSum-(0.7071067811865475*sumNuUx[0])/nuSum+wl[1]); 
} 
double VmLBOconstNuSurf1x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUy = &nuUSum[3]; 

  double favg[10]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 

  double fjump[10]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(1*fr[9])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[2]*nuSum+0.7071067811865475*dxvl[2]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 

  double Ghat[10]; 
  for(unsigned short int i=0; i<10; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*nuVtSqSum[0]*fr[9]-1.897366596101028*nuVtSqSum[0]*fl[9]+2.651650429449552*nuVtSqSum[2]*fr[7]-2.651650429449552*nuVtSqSum[2]*fl[7]-3.368048396326869*nuVtSqSum[1]*fr[5]-3.368048396326869*nuVtSqSum[1]*fl[5]-3.368048396326869*nuVtSqSum[0]*fr[3]-3.368048396326869*nuVtSqSum[0]*fl[3]+2.651650429449552*fr[1]*nuVtSqSum[1]-2.651650429449552*fl[1]*nuVtSqSum[1]+2.651650429449552*fr[0]*nuVtSqSum[0]-2.651650429449552*fl[0]*nuVtSqSum[0])*rdv-1.118033988749895*fjump[9]+alphaDrag[0]*(0.7905694150420947*favg[9]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[2]*favg[7]+alphaDrag[1]*(0.6123724356957944*favg[5]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = (1.897366596101028*nuVtSqSum[1]*fr[9]-1.897366596101028*nuVtSqSum[1]*fl[9]+2.371708245126284*nuVtSqSum[1]*fr[7]-2.371708245126284*nuVtSqSum[1]*fl[7]-3.012474066278413*nuVtSqSum[2]*fr[5]-3.368048396326869*nuVtSqSum[0]*fr[5]-3.012474066278413*nuVtSqSum[2]*fl[5]-3.368048396326869*nuVtSqSum[0]*fl[5]-3.368048396326869*nuVtSqSum[1]*fr[3]-3.368048396326869*nuVtSqSum[1]*fl[3]+2.371708245126284*fr[1]*nuVtSqSum[2]-2.371708245126284*fl[1]*nuVtSqSum[2]+2.651650429449552*fr[0]*nuVtSqSum[1]-2.651650429449552*fl[0]*nuVtSqSum[1]+2.651650429449552*nuVtSqSum[0]*fr[1]-2.651650429449552*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.7905694150420947*favg[9]+0.3162277660168379*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[5]+alphaDrag[0]*(0.6123724356957944*favg[5]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.5477225575051661*favg[5]+0.3162277660168379*favg[1])-0.5*fjump[1]; 
  Ghat[2] = ((-3.368048396326869*nuVtSqSum[0]*fr[6])-3.368048396326869*nuVtSqSum[0]*fl[6]+2.651650429449552*nuVtSqSum[1]*fr[4]-2.651650429449552*nuVtSqSum[1]*fl[4]+2.651650429449552*nuVtSqSum[0]*fr[2]-2.651650429449552*nuVtSqSum[0]*fl[2])*rdv-0.8660254037844386*fjump[6]+alphaDrag[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[2])+0.3535533905932737*alphaDrag[1]*favg[4]-0.5*fjump[2]; 
  Ghat[4] = ((-3.368048396326869*nuVtSqSum[1]*fr[6])-3.368048396326869*nuVtSqSum[1]*fl[6]+2.371708245126284*nuVtSqSum[2]*fr[4]+2.651650429449552*nuVtSqSum[0]*fr[4]-2.371708245126284*nuVtSqSum[2]*fl[4]-2.651650429449552*nuVtSqSum[0]*fl[4]+2.651650429449552*nuVtSqSum[1]*fr[2]-2.651650429449552*nuVtSqSum[1]*fl[2])*rdv+alphaDrag[1]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-0.5*fjump[4]+0.3162277660168379*alphaDrag[2]*favg[4]+0.3535533905932737*alphaDrag[0]*favg[4]; 
  Ghat[7] = (1.897366596101028*nuVtSqSum[2]*fr[9]-1.897366596101028*nuVtSqSum[2]*fl[9]+1.694077317947346*nuVtSqSum[2]*fr[7]+2.651650429449552*nuVtSqSum[0]*fr[7]-1.694077317947346*nuVtSqSum[2]*fl[7]-2.651650429449552*nuVtSqSum[0]*fl[7]-3.012474066278413*nuVtSqSum[1]*fr[5]-3.012474066278413*nuVtSqSum[1]*fl[5]-3.368048396326869*nuVtSqSum[2]*fr[3]-3.368048396326869*nuVtSqSum[2]*fl[3]+2.651650429449552*fr[0]*nuVtSqSum[2]-2.651650429449552*fl[0]*nuVtSqSum[2]+2.371708245126284*fr[1]*nuVtSqSum[1]-2.371708245126284*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[2]*(0.7905694150420947*favg[9]+0.2258769757263128*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.5*fjump[7]+0.3535533905932737*alphaDrag[0]*favg[7]+alphaDrag[1]*(0.5477225575051661*favg[5]+0.3162277660168379*favg[1]); 
  Ghat[8] = (2.651650429449552*nuVtSqSum[0]*fr[8]-2.651650429449552*nuVtSqSum[0]*fl[8])*rdv-0.5*fjump[8]+0.3535533905932737*alphaDrag[0]*favg[8]; 

  double incr1[10]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = 0.8660254037844386*Ghat[2]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = -1.118033988749895*Ghat[0]; 

  double incr2[10]; 
  incr2[3] = nuVtSqSum[0]*(0.2995357736356374*(fr[9]+fl[9])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*nuVtSqSum[2]*(fr[7]+fl[7])+nuVtSqSum[1]*((-0.430893194785552*fr[5])+0.430893194785552*fl[5]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[5] = nuVtSqSum[1]*(0.2995357736356374*(fr[9]+fl[9])+0.273861278752583*(fr[7]+fl[7])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[2]*((-0.3854025898330209*fr[5])+0.3854025898330209*fl[5]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.430893194785552*fr[5])+0.430893194785552*fl[5]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[6] = nuVtSqSum[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[2]+fl[2]))+0.3061862178478971*nuVtSqSum[1]*(fr[4]+fl[4]); 
  incr2[9] = nuVtSqSum[0]*((-1.160097062884178*(fr[9]+fl[9]))+1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*nuVtSqSum[2]*(fr[7]+fl[7])+nuVtSqSum[1]*(1.668842167398551*fr[5]-1.668842167398551*fl[5]-1.185854122563142*(fr[1]+fl[1])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr2[9]*rdvSq4R+incr1[9]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += incr2[9]*rdvSq4L-1.0*incr1[9]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUy[2])/nuSum-(0.7071067811865475*sumNuUy[0])/nuSum+wl[2]); 
} 
