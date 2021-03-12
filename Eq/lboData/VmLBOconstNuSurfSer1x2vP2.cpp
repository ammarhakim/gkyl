#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x2vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[20]; 
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
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = -1*fr[16]+fl[16]; 
  favg[17] = -1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = -1*fr[19]+fl[19]; 

  double fjump[20]; 
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
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = nuSum*vMuMidMax*(fl[16]-(-1*fr[16])); 
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(-1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(-1*fr[19])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 

  double Ghat[20]; 
  for(unsigned short int i=0; i<20; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (0.9486832980505138*nuVtSqSum[1]*fr[12]-0.9486832980505138*nuVtSqSum[1]*fl[12]-1.684024198163434*nuVtSqSum[2]*fr[11]-1.684024198163434*nuVtSqSum[2]*fl[11]+0.9486832980505137*nuVtSqSum[0]*fr[8]-0.9486832980505137*nuVtSqSum[0]*fl[8]+1.325825214724776*nuVtSqSum[2]*fr[7]-1.325825214724776*nuVtSqSum[2]*fl[7]-1.684024198163434*nuVtSqSum[1]*fr[4]-1.684024198163434*nuVtSqSum[1]*fl[4]-1.684024198163434*nuVtSqSum[0]*fr[2]-1.684024198163434*nuVtSqSum[0]*fl[2]+1.325825214724776*fr[1]*nuVtSqSum[1]-1.325825214724776*fl[1]*nuVtSqSum[1]+1.325825214724776*fr[0]*nuVtSqSum[0]-1.325825214724776*fl[0]*nuVtSqSum[0])*rdv2L+alphaDrag[1]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.6123724356957944*favg[11]+0.3535533905932737*favg[7])-1.118033988749895*fjump[8]+alphaDrag[0]*(0.7905694150420947*favg[8]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (0.848528137423857*nuVtSqSum[2]*fr[12]+0.9486832980505138*nuVtSqSum[0]*fr[12]-0.848528137423857*nuVtSqSum[2]*fl[12]-0.9486832980505138*nuVtSqSum[0]*fl[12]-1.506237033139207*nuVtSqSum[1]*fr[11]-1.506237033139207*nuVtSqSum[1]*fl[11]+0.9486832980505137*nuVtSqSum[1]*fr[8]-0.9486832980505137*nuVtSqSum[1]*fl[8]+1.185854122563142*nuVtSqSum[1]*fr[7]-1.185854122563142*nuVtSqSum[1]*fl[7]-1.506237033139206*nuVtSqSum[2]*fr[4]-1.684024198163434*nuVtSqSum[0]*fr[4]-1.506237033139206*nuVtSqSum[2]*fl[4]-1.684024198163434*nuVtSqSum[0]*fl[4]+1.185854122563142*fr[1]*nuVtSqSum[2]-1.185854122563142*fl[1]*nuVtSqSum[2]-1.684024198163434*nuVtSqSum[1]*fr[2]-1.684024198163434*nuVtSqSum[1]*fl[2]+1.325825214724776*fr[0]*nuVtSqSum[1]-1.325825214724776*fl[0]*nuVtSqSum[1]+1.325825214724776*nuVtSqSum[0]*fr[1]-1.325825214724776*nuVtSqSum[0]*fl[1])*rdv2L-1.118033988749895*fjump[12]+alphaDrag[0]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])+alphaDrag[1]*(0.5477225575051661*favg[11]+0.7905694150420947*favg[8]+0.3162277660168379*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[4]-0.5*fjump[1]; 
  Ghat[3] = (0.9486832980505137*nuVtSqSum[1]*fr[18]-0.9486832980505137*nuVtSqSum[1]*fl[18]-1.684024198163434*nuVtSqSum[2]*fr[17]-1.684024198163434*nuVtSqSum[2]*fl[17]+0.9486832980505138*nuVtSqSum[0]*fr[14]-0.9486832980505138*nuVtSqSum[0]*fl[14]+1.325825214724776*nuVtSqSum[2]*fr[13]-1.325825214724776*nuVtSqSum[2]*fl[13]-1.684024198163434*nuVtSqSum[1]*fr[10]-1.684024198163434*nuVtSqSum[1]*fl[10]-1.684024198163434*nuVtSqSum[0]*fr[6]-1.684024198163434*nuVtSqSum[0]*fl[6]+1.325825214724776*nuVtSqSum[1]*fr[5]-1.325825214724776*nuVtSqSum[1]*fl[5]+1.325825214724776*nuVtSqSum[0]*fr[3]-1.325825214724776*nuVtSqSum[0]*fl[3])*rdv2L+alphaDrag[1]*(0.7905694150420947*favg[18]+0.6123724356957944*favg[10]+0.3535533905932737*favg[5])+alphaDrag[2]*(0.6123724356957944*favg[17]+0.3535533905932737*favg[13])-1.118033988749895*fjump[14]+alphaDrag[0]*(0.7905694150420948*favg[14]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.8660254037844386*fjump[6]-0.5*fjump[3]; 
  Ghat[5] = (0.848528137423857*nuVtSqSum[2]*fr[18]+0.9486832980505137*nuVtSqSum[0]*fr[18]-0.848528137423857*nuVtSqSum[2]*fl[18]-0.9486832980505137*nuVtSqSum[0]*fl[18]-1.506237033139206*nuVtSqSum[1]*fr[17]-1.506237033139206*nuVtSqSum[1]*fl[17]+0.9486832980505138*nuVtSqSum[1]*fr[14]-0.9486832980505138*nuVtSqSum[1]*fl[14]+1.185854122563142*nuVtSqSum[1]*fr[13]-1.185854122563142*nuVtSqSum[1]*fl[13]-1.506237033139206*nuVtSqSum[2]*fr[10]-1.684024198163434*nuVtSqSum[0]*fr[10]-1.506237033139206*nuVtSqSum[2]*fl[10]-1.684024198163434*nuVtSqSum[0]*fl[10]-1.684024198163434*nuVtSqSum[1]*fr[6]-1.684024198163434*nuVtSqSum[1]*fl[6]+1.185854122563142*nuVtSqSum[2]*fr[5]+1.325825214724776*nuVtSqSum[0]*fr[5]-1.185854122563142*nuVtSqSum[2]*fl[5]-1.325825214724776*nuVtSqSum[0]*fl[5]+1.325825214724776*nuVtSqSum[1]*fr[3]-1.325825214724776*nuVtSqSum[1]*fl[3])*rdv2L-1.118033988749895*fjump[18]+alphaDrag[0]*(0.7905694150420947*favg[18]+0.6123724356957944*favg[10]+0.3535533905932737*favg[5])+alphaDrag[2]*(0.7071067811865475*favg[18]+0.5477225575051661*favg[10]+0.3162277660168379*favg[5])+alphaDrag[1]*(0.5477225575051661*favg[17]+0.7905694150420948*favg[14]+0.3162277660168379*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.8660254037844386*fjump[10]-0.5*fjump[5]; 
  Ghat[7] = (0.848528137423857*nuVtSqSum[1]*fr[12]-0.848528137423857*nuVtSqSum[1]*fl[12]-1.075883595099433*nuVtSqSum[2]*fr[11]-1.684024198163434*nuVtSqSum[0]*fr[11]-1.075883595099433*nuVtSqSum[2]*fl[11]-1.684024198163434*nuVtSqSum[0]*fl[11]+0.9486832980505137*nuVtSqSum[2]*fr[8]-0.9486832980505137*nuVtSqSum[2]*fl[8]+0.8470386589736728*nuVtSqSum[2]*fr[7]+1.325825214724776*nuVtSqSum[0]*fr[7]-0.8470386589736728*nuVtSqSum[2]*fl[7]-1.325825214724776*nuVtSqSum[0]*fl[7]-1.506237033139206*nuVtSqSum[1]*fr[4]-1.506237033139206*nuVtSqSum[1]*fl[4]-1.684024198163434*fr[2]*nuVtSqSum[2]-1.684024198163434*fl[2]*nuVtSqSum[2]+1.325825214724776*fr[0]*nuVtSqSum[2]-1.325825214724776*fl[0]*nuVtSqSum[2]+1.185854122563142*fr[1]*nuVtSqSum[1]-1.185854122563142*fl[1]*nuVtSqSum[1])*rdv2L+alphaDrag[1]*(0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])-0.8660254037844387*fjump[11]+alphaDrag[0]*(0.6123724356957944*favg[11]+0.3535533905932737*favg[7])+alphaDrag[2]*(0.3912303982179757*favg[11]+0.7905694150420947*favg[8]+0.2258769757263128*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[7]; 
  Ghat[9] = ((-1.684024198163434*nuVtSqSum[1]*fr[19])-1.684024198163434*nuVtSqSum[1]*fl[19]-1.684024198163434*nuVtSqSum[0]*fr[16]-1.684024198163434*nuVtSqSum[0]*fl[16]+1.325825214724776*nuVtSqSum[1]*fr[15]-1.325825214724776*nuVtSqSum[1]*fl[15]+1.325825214724776*nuVtSqSum[0]*fr[9]-1.325825214724776*nuVtSqSum[0]*fl[9])*rdv2L+alphaDrag[1]*(0.6123724356957944*favg[19]+0.3535533905932737*favg[15])-0.8660254037844387*fjump[16]+alphaDrag[0]*(0.6123724356957944*favg[16]+0.3535533905932737*favg[9])-0.5*fjump[9]; 
  Ghat[13] = (0.848528137423857*nuVtSqSum[1]*fr[18]-0.848528137423857*nuVtSqSum[1]*fl[18]-1.075883595099433*nuVtSqSum[2]*fr[17]-1.684024198163434*nuVtSqSum[0]*fr[17]-1.075883595099433*nuVtSqSum[2]*fl[17]-1.684024198163434*nuVtSqSum[0]*fl[17]+0.9486832980505137*nuVtSqSum[2]*fr[14]-0.9486832980505137*nuVtSqSum[2]*fl[14]+0.8470386589736728*nuVtSqSum[2]*fr[13]+1.325825214724776*nuVtSqSum[0]*fr[13]-0.8470386589736728*nuVtSqSum[2]*fl[13]-1.325825214724776*nuVtSqSum[0]*fl[13]-1.506237033139207*nuVtSqSum[1]*fr[10]-1.506237033139207*nuVtSqSum[1]*fl[10]-1.684024198163434*nuVtSqSum[2]*fr[6]-1.684024198163434*nuVtSqSum[2]*fl[6]+1.185854122563142*nuVtSqSum[1]*fr[5]-1.185854122563142*nuVtSqSum[1]*fl[5]+1.325825214724776*nuVtSqSum[2]*fr[3]-1.325825214724776*nuVtSqSum[2]*fl[3])*rdv2L+alphaDrag[1]*(0.7071067811865475*favg[18]+0.5477225575051661*favg[10]+0.3162277660168379*favg[5])-0.8660254037844387*fjump[17]+alphaDrag[0]*(0.6123724356957944*favg[17]+0.3535533905932737*favg[13])+alphaDrag[2]*(0.3912303982179757*favg[17]+0.7905694150420947*favg[14]+0.2258769757263128*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[13]; 
  Ghat[15] = ((-1.506237033139207*nuVtSqSum[2]*fr[19])-1.684024198163434*nuVtSqSum[0]*fr[19]-1.506237033139207*nuVtSqSum[2]*fl[19]-1.684024198163434*nuVtSqSum[0]*fl[19]-1.684024198163434*nuVtSqSum[1]*fr[16]-1.684024198163434*nuVtSqSum[1]*fl[16]+1.185854122563142*nuVtSqSum[2]*fr[15]+1.325825214724776*nuVtSqSum[0]*fr[15]-1.185854122563142*nuVtSqSum[2]*fl[15]-1.325825214724776*nuVtSqSum[0]*fl[15]+1.325825214724776*nuVtSqSum[1]*fr[9]-1.325825214724776*nuVtSqSum[1]*fl[9])*rdv2L-0.8660254037844387*fjump[19]+alphaDrag[0]*(0.6123724356957944*favg[19]+0.3535533905932737*favg[15])+alphaDrag[2]*(0.5477225575051661*favg[19]+0.3162277660168379*favg[15])+alphaDrag[1]*(0.6123724356957944*favg[16]+0.3535533905932737*favg[9])-0.5*fjump[15]; 

  double incr1[20]; 
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
  incr1[10] = 0.8660254037844386*Ghat[5]; 
  incr1[11] = 0.8660254037844387*Ghat[7]; 
  incr1[12] = -1.118033988749895*Ghat[1]; 
  incr1[13] = -0.5*Ghat[13]; 
  incr1[14] = -1.118033988749895*Ghat[3]; 
  incr1[15] = -0.5*Ghat[15]; 
  incr1[16] = 0.8660254037844387*Ghat[9]; 
  incr1[17] = 0.8660254037844387*Ghat[13]; 
  incr1[18] = -1.118033988749895*Ghat[5]; 
  incr1[19] = 0.8660254037844387*Ghat[15]; 

  double incr2[20]; 
  incr2[2] = nuVtSqSum[1]*(0.2995357736356374*(fr[12]+fl[12])-0.430893194785552*fr[4]+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.430893194785552*fr[11])+0.430893194785552*fl[11]+0.3061862178478971*(fr[7]+fl[7]))+nuVtSqSum[0]*(0.2995357736356374*(fr[8]+fl[8])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[4] = nuVtSqSum[0]*(0.2995357736356374*(fr[12]+fl[12])-0.430893194785552*fr[4]+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[2]*(0.2679129406169099*(fr[12]+fl[12])-0.3854025898330209*fr[4]+0.3854025898330209*fl[4]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.3854025898330209*fr[11])+0.3854025898330209*fl[11]+0.2995357736356374*(fr[8]+fl[8])+0.273861278752583*(fr[7]+fl[7])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[6] = nuVtSqSum[1]*(0.2995357736356374*(fr[18]+fl[18])-0.430893194785552*fr[10]+0.430893194785552*fl[10]+0.3061862178478971*(fr[5]+fl[5]))+nuVtSqSum[2]*((-0.430893194785552*fr[17])+0.430893194785552*fl[17]+0.3061862178478971*(fr[13]+fl[13]))+nuVtSqSum[0]*(0.2995357736356374*(fr[14]+fl[14])-0.430893194785552*fr[6]+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 
  incr2[8] = nuVtSqSum[1]*((-1.160097062884178*(fr[12]+fl[12]))+1.668842167398551*fr[4]-1.668842167398551*fl[4]-1.185854122563142*(fr[1]+fl[1]))+nuVtSqSum[2]*(1.668842167398552*fr[11]-1.668842167398552*fl[11]-1.185854122563142*(fr[7]+fl[7]))+nuVtSqSum[0]*((-1.160097062884178*(fr[8]+fl[8]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0])); 
  incr2[10] = nuVtSqSum[0]*(0.2995357736356374*(fr[18]+fl[18])-0.430893194785552*fr[10]+0.430893194785552*fl[10]+0.3061862178478971*(fr[5]+fl[5]))+nuVtSqSum[2]*(0.2679129406169099*(fr[18]+fl[18])-0.3854025898330209*fr[10]+0.3854025898330209*fl[10]+0.273861278752583*(fr[5]+fl[5]))+nuVtSqSum[1]*((-0.3854025898330209*fr[17])+0.3854025898330209*fl[17]+0.2995357736356374*(fr[14]+fl[14])+0.273861278752583*(fr[13]+fl[13])-0.430893194785552*fr[6]+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 
  incr2[11] = nuVtSqSum[1]*(0.2679129406169099*(fr[12]+fl[12])-0.3854025898330209*fr[4]+0.3854025898330209*fl[4]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.2752875641664436*fr[11])+0.2752875641664436*fl[11]+0.2995357736356374*(fr[8]+fl[8])+0.1956151991089878*(fr[7]+fl[7])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[0]*((-0.430893194785552*fr[11])+0.430893194785552*fl[11]+0.3061862178478971*(fr[7]+fl[7])); 
  incr2[12] = nuVtSqSum[2]*((-1.037622357242749*(fr[12]+fl[12]))+1.492657812008498*fr[4]-1.492657812008498*fl[4]-1.060660171779821*(fr[1]+fl[1]))+nuVtSqSum[0]*((-1.160097062884178*(fr[12]+fl[12]))+1.668842167398552*fr[4]-1.668842167398552*fl[4]-1.185854122563142*(fr[1]+fl[1]))+nuVtSqSum[1]*(1.492657812008498*fr[11]-1.492657812008498*fl[11]-1.160097062884178*(fr[8]+fl[8])-1.060660171779821*(fr[7]+fl[7])+1.668842167398552*fr[2]-1.668842167398552*fl[2]-1.185854122563142*(fr[0]+fl[0])); 
  incr2[14] = nuVtSqSum[1]*((-1.160097062884178*(fr[18]+fl[18]))+1.668842167398552*fr[10]-1.668842167398552*fl[10]-1.185854122563142*(fr[5]+fl[5]))+nuVtSqSum[2]*(1.668842167398552*fr[17]-1.668842167398552*fl[17]-1.185854122563142*(fr[13]+fl[13]))+nuVtSqSum[0]*((-1.160097062884178*(fr[14]+fl[14]))+1.668842167398552*fr[6]-1.668842167398552*fl[6]-1.185854122563142*(fr[3]+fl[3])); 
  incr2[16] = nuVtSqSum[1]*((-0.430893194785552*fr[19])+0.430893194785552*fl[19]+0.3061862178478971*(fr[15]+fl[15]))+nuVtSqSum[0]*((-0.430893194785552*fr[16])+0.430893194785552*fl[16]+0.3061862178478971*(fr[9]+fl[9])); 
  incr2[17] = nuVtSqSum[1]*(0.2679129406169099*(fr[18]+fl[18])-0.3854025898330209*fr[10]+0.3854025898330209*fl[10]+0.273861278752583*(fr[5]+fl[5]))+nuVtSqSum[2]*((-0.2752875641664436*fr[17])+0.2752875641664436*fl[17]+0.2995357736356374*(fr[14]+fl[14])+0.1956151991089878*(fr[13]+fl[13])-0.430893194785552*fr[6]+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3]))+nuVtSqSum[0]*((-0.430893194785552*fr[17])+0.430893194785552*fl[17]+0.3061862178478971*(fr[13]+fl[13])); 
  incr2[18] = nuVtSqSum[2]*((-1.037622357242749*(fr[18]+fl[18]))+1.492657812008498*fr[10]-1.492657812008498*fl[10]-1.060660171779821*(fr[5]+fl[5]))+nuVtSqSum[0]*((-1.160097062884178*(fr[18]+fl[18]))+1.668842167398551*fr[10]-1.668842167398551*fl[10]-1.185854122563142*(fr[5]+fl[5]))+nuVtSqSum[1]*(1.492657812008498*fr[17]-1.492657812008498*fl[17]-1.160097062884178*(fr[14]+fl[14])-1.060660171779821*(fr[13]+fl[13])+1.668842167398551*fr[6]-1.668842167398551*fl[6]-1.185854122563142*(fr[3]+fl[3])); 
  incr2[19] = nuVtSqSum[2]*((-0.3854025898330209*fr[19])+0.3854025898330209*fl[19]+0.273861278752583*(fr[15]+fl[15]))+nuVtSqSum[0]*((-0.430893194785552*fr[19])+0.430893194785552*fl[19]+0.3061862178478971*(fr[15]+fl[15]))+nuVtSqSum[1]*((-0.430893194785552*fr[16])+0.430893194785552*fl[16]+0.3061862178478971*(fr[9]+fl[9])); 

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
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr2[11]*rdvSq4R+incr1[11]*rdv2R; 
  outr[12] += incr2[12]*rdvSq4R+incr1[12]*rdv2R; 
  outr[13] += incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr1[15]*rdv2R; 
  outr[16] += incr2[16]*rdvSq4R+incr1[16]*rdv2R; 
  outr[17] += incr2[17]*rdvSq4R+incr1[17]*rdv2R; 
  outr[18] += incr2[18]*rdvSq4R+incr1[18]*rdv2R; 
  outr[19] += incr2[19]*rdvSq4R+incr1[19]*rdv2R; 

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
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += incr1[11]*rdv2L-1.0*incr2[11]*rdvSq4L; 
  outl[12] += incr2[12]*rdvSq4L-1.0*incr1[12]*rdv2L; 
  outl[13] += -1.0*incr1[13]*rdv2L; 
  outl[14] += incr2[14]*rdvSq4L-1.0*incr1[14]*rdv2L; 
  outl[15] += -1.0*incr1[15]*rdv2L; 
  outl[16] += incr1[16]*rdv2L-1.0*incr2[16]*rdvSq4L; 
  outl[17] += incr1[17]*rdv2L-1.0*incr2[17]*rdvSq4L; 
  outl[18] += incr2[18]*rdvSq4L-1.0*incr1[18]*rdv2L; 
  outl[19] += incr1[19]*rdv2L-1.0*incr2[19]*rdvSq4L; 

  return std::abs((0.7905694150420947*sumNuUx[2])/nuSum-(0.7071067811865475*sumNuUx[0])/nuSum+wl[1]); 
} 
double VmLBOconstNuSurf1x2vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUy = &nuUSum[3]; 

  double favg[20]; 
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
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = -1*fr[17]+fl[17]; 
  favg[18] = -1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 

  double fjump[20]; 
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
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(-1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = nuSum*vMuMidMax*(fl[16]-(1*fr[16])); 
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(-1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(-1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(1*fr[19])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[2]*nuSum+0.7071067811865475*dxvl[2]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 

  double Ghat[20]; 
  for(unsigned short int i=0; i<20; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (0.9486832980505138*nuVtSqSum[1]*fr[15]-0.9486832980505138*nuVtSqSum[1]*fl[15]-1.684024198163434*nuVtSqSum[2]*fr[13]-1.684024198163434*nuVtSqSum[2]*fl[13]+0.9486832980505137*nuVtSqSum[0]*fr[9]-0.9486832980505137*nuVtSqSum[0]*fl[9]+1.325825214724776*nuVtSqSum[2]*fr[7]-1.325825214724776*nuVtSqSum[2]*fl[7]-1.684024198163434*nuVtSqSum[1]*fr[5]-1.684024198163434*nuVtSqSum[1]*fl[5]-1.684024198163434*nuVtSqSum[0]*fr[3]-1.684024198163434*nuVtSqSum[0]*fl[3]+1.325825214724776*fr[1]*nuVtSqSum[1]-1.325825214724776*fl[1]*nuVtSqSum[1]+1.325825214724776*fr[0]*nuVtSqSum[0]-1.325825214724776*fl[0]*nuVtSqSum[0])*rdv2L+alphaDrag[1]*(0.7905694150420948*favg[15]+0.6123724356957944*favg[5]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.6123724356957944*favg[13]+0.3535533905932737*favg[7])-1.118033988749895*fjump[9]+alphaDrag[0]*(0.7905694150420947*favg[9]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = (0.848528137423857*nuVtSqSum[2]*fr[15]+0.9486832980505138*nuVtSqSum[0]*fr[15]-0.848528137423857*nuVtSqSum[2]*fl[15]-0.9486832980505138*nuVtSqSum[0]*fl[15]-1.506237033139207*nuVtSqSum[1]*fr[13]-1.506237033139207*nuVtSqSum[1]*fl[13]+0.9486832980505137*nuVtSqSum[1]*fr[9]-0.9486832980505137*nuVtSqSum[1]*fl[9]+1.185854122563142*nuVtSqSum[1]*fr[7]-1.185854122563142*nuVtSqSum[1]*fl[7]-1.506237033139206*nuVtSqSum[2]*fr[5]-1.684024198163434*nuVtSqSum[0]*fr[5]-1.506237033139206*nuVtSqSum[2]*fl[5]-1.684024198163434*nuVtSqSum[0]*fl[5]-1.684024198163434*nuVtSqSum[1]*fr[3]-1.684024198163434*nuVtSqSum[1]*fl[3]+1.185854122563142*fr[1]*nuVtSqSum[2]-1.185854122563142*fl[1]*nuVtSqSum[2]+1.325825214724776*fr[0]*nuVtSqSum[1]-1.325825214724776*fl[0]*nuVtSqSum[1]+1.325825214724776*nuVtSqSum[0]*fr[1]-1.325825214724776*nuVtSqSum[0]*fl[1])*rdv2L-1.118033988749895*fjump[15]+alphaDrag[0]*(0.7905694150420948*favg[15]+0.6123724356957944*favg[5]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.7071067811865475*favg[15]+0.5477225575051661*favg[5]+0.3162277660168379*favg[1])+alphaDrag[1]*(0.5477225575051661*favg[13]+0.7905694150420947*favg[9]+0.3162277660168379*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[5]-0.5*fjump[1]; 
  Ghat[2] = (0.9486832980505137*nuVtSqSum[1]*fr[19]-0.9486832980505137*nuVtSqSum[1]*fl[19]-1.684024198163434*nuVtSqSum[2]*fr[17]-1.684024198163434*nuVtSqSum[2]*fl[17]+0.9486832980505138*nuVtSqSum[0]*fr[16]-0.9486832980505138*nuVtSqSum[0]*fl[16]+1.325825214724776*nuVtSqSum[2]*fr[11]-1.325825214724776*nuVtSqSum[2]*fl[11]-1.684024198163434*nuVtSqSum[1]*fr[10]-1.684024198163434*nuVtSqSum[1]*fl[10]-1.684024198163434*nuVtSqSum[0]*fr[6]-1.684024198163434*nuVtSqSum[0]*fl[6]+1.325825214724776*nuVtSqSum[1]*fr[4]-1.325825214724776*nuVtSqSum[1]*fl[4]+1.325825214724776*nuVtSqSum[0]*fr[2]-1.325825214724776*nuVtSqSum[0]*fl[2])*rdv2L+alphaDrag[1]*(0.7905694150420947*favg[19]+0.6123724356957944*favg[10]+0.3535533905932737*favg[4])+alphaDrag[2]*(0.6123724356957944*favg[17]+0.3535533905932737*favg[11])-1.118033988749895*fjump[16]+alphaDrag[0]*(0.7905694150420948*favg[16]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-0.8660254037844386*fjump[6]-0.5*fjump[2]; 
  Ghat[4] = (0.848528137423857*nuVtSqSum[2]*fr[19]+0.9486832980505137*nuVtSqSum[0]*fr[19]-0.848528137423857*nuVtSqSum[2]*fl[19]-0.9486832980505137*nuVtSqSum[0]*fl[19]-1.506237033139206*nuVtSqSum[1]*fr[17]-1.506237033139206*nuVtSqSum[1]*fl[17]+0.9486832980505138*nuVtSqSum[1]*fr[16]-0.9486832980505138*nuVtSqSum[1]*fl[16]+1.185854122563142*nuVtSqSum[1]*fr[11]-1.185854122563142*nuVtSqSum[1]*fl[11]-1.506237033139206*nuVtSqSum[2]*fr[10]-1.684024198163434*nuVtSqSum[0]*fr[10]-1.506237033139206*nuVtSqSum[2]*fl[10]-1.684024198163434*nuVtSqSum[0]*fl[10]-1.684024198163434*nuVtSqSum[1]*fr[6]-1.684024198163434*nuVtSqSum[1]*fl[6]+1.185854122563142*nuVtSqSum[2]*fr[4]+1.325825214724776*nuVtSqSum[0]*fr[4]-1.185854122563142*nuVtSqSum[2]*fl[4]-1.325825214724776*nuVtSqSum[0]*fl[4]+1.325825214724776*nuVtSqSum[1]*fr[2]-1.325825214724776*nuVtSqSum[1]*fl[2])*rdv2L-1.118033988749895*fjump[19]+alphaDrag[0]*(0.7905694150420947*favg[19]+0.6123724356957944*favg[10]+0.3535533905932737*favg[4])+alphaDrag[2]*(0.7071067811865475*favg[19]+0.5477225575051661*favg[10]+0.3162277660168379*favg[4])+alphaDrag[1]*(0.5477225575051661*favg[17]+0.7905694150420948*favg[16]+0.3162277660168379*favg[11]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-0.8660254037844386*fjump[10]-0.5*fjump[4]; 
  Ghat[7] = (0.848528137423857*nuVtSqSum[1]*fr[15]-0.848528137423857*nuVtSqSum[1]*fl[15]-1.075883595099433*nuVtSqSum[2]*fr[13]-1.684024198163434*nuVtSqSum[0]*fr[13]-1.075883595099433*nuVtSqSum[2]*fl[13]-1.684024198163434*nuVtSqSum[0]*fl[13]+0.9486832980505137*nuVtSqSum[2]*fr[9]-0.9486832980505137*nuVtSqSum[2]*fl[9]+0.8470386589736728*nuVtSqSum[2]*fr[7]+1.325825214724776*nuVtSqSum[0]*fr[7]-0.8470386589736728*nuVtSqSum[2]*fl[7]-1.325825214724776*nuVtSqSum[0]*fl[7]-1.506237033139206*nuVtSqSum[1]*fr[5]-1.506237033139206*nuVtSqSum[1]*fl[5]-1.684024198163434*nuVtSqSum[2]*fr[3]-1.684024198163434*nuVtSqSum[2]*fl[3]+1.325825214724776*fr[0]*nuVtSqSum[2]-1.325825214724776*fl[0]*nuVtSqSum[2]+1.185854122563142*fr[1]*nuVtSqSum[1]-1.185854122563142*fl[1]*nuVtSqSum[1])*rdv2L+alphaDrag[1]*(0.7071067811865475*favg[15]+0.5477225575051661*favg[5]+0.3162277660168379*favg[1])-0.8660254037844387*fjump[13]+alphaDrag[0]*(0.6123724356957944*favg[13]+0.3535533905932737*favg[7])+alphaDrag[2]*(0.3912303982179757*favg[13]+0.7905694150420947*favg[9]+0.2258769757263128*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.5*fjump[7]; 
  Ghat[8] = ((-1.684024198163434*nuVtSqSum[1]*fr[18])-1.684024198163434*nuVtSqSum[1]*fl[18]-1.684024198163434*nuVtSqSum[0]*fr[14]-1.684024198163434*nuVtSqSum[0]*fl[14]+1.325825214724776*nuVtSqSum[1]*fr[12]-1.325825214724776*nuVtSqSum[1]*fl[12]+1.325825214724776*nuVtSqSum[0]*fr[8]-1.325825214724776*nuVtSqSum[0]*fl[8])*rdv2L+alphaDrag[1]*(0.6123724356957944*favg[18]+0.3535533905932737*favg[12])-0.8660254037844387*fjump[14]+alphaDrag[0]*(0.6123724356957944*favg[14]+0.3535533905932737*favg[8])-0.5*fjump[8]; 
  Ghat[11] = (0.848528137423857*nuVtSqSum[1]*fr[19]-0.848528137423857*nuVtSqSum[1]*fl[19]-1.075883595099433*nuVtSqSum[2]*fr[17]-1.684024198163434*nuVtSqSum[0]*fr[17]-1.075883595099433*nuVtSqSum[2]*fl[17]-1.684024198163434*nuVtSqSum[0]*fl[17]+0.9486832980505137*nuVtSqSum[2]*fr[16]-0.9486832980505137*nuVtSqSum[2]*fl[16]+0.8470386589736728*nuVtSqSum[2]*fr[11]+1.325825214724776*nuVtSqSum[0]*fr[11]-0.8470386589736728*nuVtSqSum[2]*fl[11]-1.325825214724776*nuVtSqSum[0]*fl[11]-1.506237033139207*nuVtSqSum[1]*fr[10]-1.506237033139207*nuVtSqSum[1]*fl[10]-1.684024198163434*nuVtSqSum[2]*fr[6]-1.684024198163434*nuVtSqSum[2]*fl[6]+1.185854122563142*nuVtSqSum[1]*fr[4]-1.185854122563142*nuVtSqSum[1]*fl[4]+1.325825214724776*fr[2]*nuVtSqSum[2]-1.325825214724776*fl[2]*nuVtSqSum[2])*rdv2L+alphaDrag[1]*(0.7071067811865475*favg[19]+0.5477225575051661*favg[10]+0.3162277660168379*favg[4])-0.8660254037844387*fjump[17]+alphaDrag[0]*(0.6123724356957944*favg[17]+0.3535533905932737*favg[11])+alphaDrag[2]*(0.3912303982179757*favg[17]+0.7905694150420947*favg[16]+0.2258769757263128*favg[11]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-0.5*fjump[11]; 
  Ghat[12] = ((-1.506237033139207*nuVtSqSum[2]*fr[18])-1.684024198163434*nuVtSqSum[0]*fr[18]-1.506237033139207*nuVtSqSum[2]*fl[18]-1.684024198163434*nuVtSqSum[0]*fl[18]-1.684024198163434*nuVtSqSum[1]*fr[14]-1.684024198163434*nuVtSqSum[1]*fl[14]+1.185854122563142*nuVtSqSum[2]*fr[12]+1.325825214724776*nuVtSqSum[0]*fr[12]-1.185854122563142*nuVtSqSum[2]*fl[12]-1.325825214724776*nuVtSqSum[0]*fl[12]+1.325825214724776*nuVtSqSum[1]*fr[8]-1.325825214724776*nuVtSqSum[1]*fl[8])*rdv2L-0.8660254037844387*fjump[18]+alphaDrag[0]*(0.6123724356957944*favg[18]+0.3535533905932737*favg[12])+alphaDrag[2]*(0.5477225575051661*favg[18]+0.3162277660168379*favg[12])+alphaDrag[1]*(0.6123724356957944*favg[14]+0.3535533905932737*favg[8])-0.5*fjump[12]; 

  double incr1[20]; 
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
  incr1[10] = 0.8660254037844386*Ghat[4]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = 0.8660254037844387*Ghat[7]; 
  incr1[14] = 0.8660254037844387*Ghat[8]; 
  incr1[15] = -1.118033988749895*Ghat[1]; 
  incr1[16] = -1.118033988749895*Ghat[2]; 
  incr1[17] = 0.8660254037844387*Ghat[11]; 
  incr1[18] = 0.8660254037844387*Ghat[12]; 
  incr1[19] = -1.118033988749895*Ghat[4]; 

  double incr2[20]; 
  incr2[3] = nuVtSqSum[1]*(0.2995357736356374*(fr[15]+fl[15])-0.430893194785552*fr[5]+0.430893194785552*fl[5]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.430893194785552*fr[13])+0.430893194785552*fl[13]+0.3061862178478971*(fr[7]+fl[7]))+nuVtSqSum[0]*(0.2995357736356374*(fr[9]+fl[9])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[5] = nuVtSqSum[0]*(0.2995357736356374*(fr[15]+fl[15])-0.430893194785552*fr[5]+0.430893194785552*fl[5]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[2]*(0.2679129406169099*(fr[15]+fl[15])-0.3854025898330209*fr[5]+0.3854025898330209*fl[5]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.3854025898330209*fr[13])+0.3854025898330209*fl[13]+0.2995357736356374*(fr[9]+fl[9])+0.273861278752583*(fr[7]+fl[7])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[6] = nuVtSqSum[1]*(0.2995357736356374*(fr[19]+fl[19])-0.430893194785552*fr[10]+0.430893194785552*fl[10]+0.3061862178478971*(fr[4]+fl[4]))+nuVtSqSum[2]*((-0.430893194785552*fr[17])+0.430893194785552*fl[17]+0.3061862178478971*(fr[11]+fl[11]))+nuVtSqSum[0]*(0.2995357736356374*(fr[16]+fl[16])-0.430893194785552*fr[6]+0.430893194785552*fl[6]+0.3061862178478971*(fr[2]+fl[2])); 
  incr2[9] = nuVtSqSum[1]*((-1.160097062884178*(fr[15]+fl[15]))+1.668842167398551*fr[5]-1.668842167398551*fl[5]-1.185854122563142*(fr[1]+fl[1]))+nuVtSqSum[2]*(1.668842167398552*fr[13]-1.668842167398552*fl[13]-1.185854122563142*(fr[7]+fl[7]))+nuVtSqSum[0]*((-1.160097062884178*(fr[9]+fl[9]))+1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[0]+fl[0])); 
  incr2[10] = nuVtSqSum[0]*(0.2995357736356374*(fr[19]+fl[19])-0.430893194785552*fr[10]+0.430893194785552*fl[10]+0.3061862178478971*(fr[4]+fl[4]))+nuVtSqSum[2]*(0.2679129406169099*(fr[19]+fl[19])-0.3854025898330209*fr[10]+0.3854025898330209*fl[10]+0.273861278752583*(fr[4]+fl[4]))+nuVtSqSum[1]*((-0.3854025898330209*fr[17])+0.3854025898330209*fl[17]+0.2995357736356374*(fr[16]+fl[16])+0.273861278752583*(fr[11]+fl[11])-0.430893194785552*fr[6]+0.430893194785552*fl[6]+0.3061862178478971*(fr[2]+fl[2])); 
  incr2[13] = nuVtSqSum[1]*(0.2679129406169099*(fr[15]+fl[15])-0.3854025898330209*fr[5]+0.3854025898330209*fl[5]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.2752875641664436*fr[13])+0.2752875641664436*fl[13]+0.2995357736356374*(fr[9]+fl[9])+0.1956151991089878*(fr[7]+fl[7])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[0]*((-0.430893194785552*fr[13])+0.430893194785552*fl[13]+0.3061862178478971*(fr[7]+fl[7])); 
  incr2[14] = nuVtSqSum[1]*((-0.430893194785552*fr[18])+0.430893194785552*fl[18]+0.3061862178478971*(fr[12]+fl[12]))+nuVtSqSum[0]*((-0.430893194785552*fr[14])+0.430893194785552*fl[14]+0.3061862178478971*(fr[8]+fl[8])); 
  incr2[15] = nuVtSqSum[2]*((-1.037622357242749*(fr[15]+fl[15]))+1.492657812008498*fr[5]-1.492657812008498*fl[5]-1.060660171779821*(fr[1]+fl[1]))+nuVtSqSum[0]*((-1.160097062884178*(fr[15]+fl[15]))+1.668842167398552*fr[5]-1.668842167398552*fl[5]-1.185854122563142*(fr[1]+fl[1]))+nuVtSqSum[1]*(1.492657812008498*fr[13]-1.492657812008498*fl[13]-1.160097062884178*(fr[9]+fl[9])-1.060660171779821*(fr[7]+fl[7])+1.668842167398552*fr[3]-1.668842167398552*fl[3]-1.185854122563142*(fr[0]+fl[0])); 
  incr2[16] = nuVtSqSum[1]*((-1.160097062884178*(fr[19]+fl[19]))+1.668842167398552*fr[10]-1.668842167398552*fl[10]-1.185854122563142*(fr[4]+fl[4]))+nuVtSqSum[2]*(1.668842167398552*fr[17]-1.668842167398552*fl[17]-1.185854122563142*(fr[11]+fl[11]))+nuVtSqSum[0]*((-1.160097062884178*(fr[16]+fl[16]))+1.668842167398552*fr[6]-1.668842167398552*fl[6]-1.185854122563142*(fr[2]+fl[2])); 
  incr2[17] = nuVtSqSum[1]*(0.2679129406169099*(fr[19]+fl[19])-0.3854025898330209*fr[10]+0.3854025898330209*fl[10]+0.273861278752583*(fr[4]+fl[4]))+nuVtSqSum[2]*((-0.2752875641664436*fr[17])+0.2752875641664436*fl[17]+0.2995357736356374*(fr[16]+fl[16])+0.1956151991089878*(fr[11]+fl[11])-0.430893194785552*fr[6]+0.430893194785552*fl[6]+0.3061862178478971*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.430893194785552*fr[17])+0.430893194785552*fl[17]+0.3061862178478971*(fr[11]+fl[11])); 
  incr2[18] = nuVtSqSum[2]*((-0.3854025898330209*fr[18])+0.3854025898330209*fl[18]+0.273861278752583*(fr[12]+fl[12]))+nuVtSqSum[0]*((-0.430893194785552*fr[18])+0.430893194785552*fl[18]+0.3061862178478971*(fr[12]+fl[12]))+nuVtSqSum[1]*((-0.430893194785552*fr[14])+0.430893194785552*fl[14]+0.3061862178478971*(fr[8]+fl[8])); 
  incr2[19] = nuVtSqSum[2]*((-1.037622357242749*(fr[19]+fl[19]))+1.492657812008498*fr[10]-1.492657812008498*fl[10]-1.060660171779821*(fr[4]+fl[4]))+nuVtSqSum[0]*((-1.160097062884178*(fr[19]+fl[19]))+1.668842167398551*fr[10]-1.668842167398551*fl[10]-1.185854122563142*(fr[4]+fl[4]))+nuVtSqSum[1]*(1.492657812008498*fr[17]-1.492657812008498*fl[17]-1.160097062884178*(fr[16]+fl[16])-1.060660171779821*(fr[11]+fl[11])+1.668842167398551*fr[6]-1.668842167398551*fl[6]-1.185854122563142*(fr[2]+fl[2])); 

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
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 
  outr[16] += incr2[16]*rdvSq4R+incr1[16]*rdv2R; 
  outr[17] += incr2[17]*rdvSq4R+incr1[17]*rdv2R; 
  outr[18] += incr2[18]*rdvSq4R+incr1[18]*rdv2R; 
  outr[19] += incr2[19]*rdvSq4R+incr1[19]*rdv2R; 

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
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += incr1[13]*rdv2L-1.0*incr2[13]*rdvSq4L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += incr2[15]*rdvSq4L-1.0*incr1[15]*rdv2L; 
  outl[16] += incr2[16]*rdvSq4L-1.0*incr1[16]*rdv2L; 
  outl[17] += incr1[17]*rdv2L-1.0*incr2[17]*rdvSq4L; 
  outl[18] += incr1[18]*rdv2L-1.0*incr2[18]*rdvSq4L; 
  outl[19] += incr2[19]*rdvSq4L-1.0*incr1[19]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUy[2])/nuSum-(0.7071067811865475*sumNuUy[0])/nuSum+wl[2]); 
} 
