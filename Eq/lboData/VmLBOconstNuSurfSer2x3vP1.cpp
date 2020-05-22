#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf2x3vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:          Cell-center coordinates. 
  // dxv[5]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[12]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[32]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = -1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = -1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = -1*fr[18]+fl[18]; 
  favg[19] = -1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  favg[21] = -1*fr[21]+fl[21]; 
  favg[22] = -1*fr[22]+fl[22]; 
  favg[23] = 1*fr[23]+fl[23]; 
  favg[24] = 1*fr[24]+fl[24]; 
  favg[25] = -1*fr[25]+fl[25]; 
  favg[26] = -1*fr[26]+fl[26]; 
  favg[27] = -1*fr[27]+fl[27]; 
  favg[28] = 1*fr[28]+fl[28]; 
  favg[29] = -1*fr[29]+fl[29]; 
  favg[30] = -1*fr[30]+fl[30]; 
  favg[31] = -1*fr[31]+fl[31]; 

  double fjump[32]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(-1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(-1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = nuSum*vMuMidMax*(fl[16]-(-1*fr[16])); 
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(-1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(-1*fr[19])); 
  fjump[20] = nuSum*vMuMidMax*(fl[20]-(1*fr[20])); 
  fjump[21] = nuSum*vMuMidMax*(fl[21]-(-1*fr[21])); 
  fjump[22] = nuSum*vMuMidMax*(fl[22]-(-1*fr[22])); 
  fjump[23] = nuSum*vMuMidMax*(fl[23]-(1*fr[23])); 
  fjump[24] = nuSum*vMuMidMax*(fl[24]-(1*fr[24])); 
  fjump[25] = nuSum*vMuMidMax*(fl[25]-(-1*fr[25])); 
  fjump[26] = nuSum*vMuMidMax*(fl[26]-(-1*fr[26])); 
  fjump[27] = nuSum*vMuMidMax*(fl[27]-(-1*fr[27])); 
  fjump[28] = nuSum*vMuMidMax*(fl[28]-(1*fr[28])); 
  fjump[29] = nuSum*vMuMidMax*(fl[29]-(-1*fr[29])); 
  fjump[30] = nuSum*vMuMidMax*(fl[30]-(-1*fr[30])); 
  fjump[31] = nuSum*vMuMidMax*(fl[31]-(-1*fr[31])); 

  double alphaDrag[4]; 
  alphaDrag[0] = 2.0*wl[2]*nuSum+dxvl[2]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 
  alphaDrag[3] = -1.0*sumNuUx[3]; 

  double Ghat[32]; 
  for(unsigned short int i=0; i<32; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.082531754730548*nuVtSqSum[3]*fr[16])-1.082531754730548*nuVtSqSum[3]*fl[16]-1.082531754730548*nuVtSqSum[2]*fr[8]-1.082531754730548*nuVtSqSum[2]*fl[8]-1.082531754730548*nuVtSqSum[1]*fr[7]-1.082531754730548*nuVtSqSum[1]*fl[7]+1.125*nuVtSqSum[3]*fr[6]-1.125*nuVtSqSum[3]*fl[6]-1.082531754730548*nuVtSqSum[0]*fr[3]-1.082531754730548*nuVtSqSum[0]*fl[3]+1.125*fr[2]*nuVtSqSum[2]-1.125*fl[2]*nuVtSqSum[2]+1.125*fr[1]*nuVtSqSum[1]-1.125*fl[1]*nuVtSqSum[1]+1.125*fr[0]*nuVtSqSum[0]-1.125*fl[0]*nuVtSqSum[0])*rdv+alphaDrag[3]*(0.4330127018922193*favg[16]+0.25*favg[6])+alphaDrag[2]*(0.4330127018922193*favg[8]+0.25*favg[2])+alphaDrag[1]*(0.4330127018922193*favg[7]+0.25*favg[1])-0.8660254037844386*fjump[3]+alphaDrag[0]*(0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[0]; 
  Ghat[1] = ((-1.082531754730548*nuVtSqSum[2]*fr[16])-1.082531754730548*nuVtSqSum[2]*fl[16]-1.082531754730548*nuVtSqSum[3]*fr[8]-1.082531754730548*nuVtSqSum[3]*fl[8]-1.082531754730548*nuVtSqSum[0]*fr[7]-1.082531754730548*nuVtSqSum[0]*fl[7]+1.125*nuVtSqSum[2]*fr[6]-1.125*nuVtSqSum[2]*fl[6]+1.125*fr[2]*nuVtSqSum[3]-1.125*fl[2]*nuVtSqSum[3]-1.082531754730548*nuVtSqSum[1]*fr[3]-1.082531754730548*nuVtSqSum[1]*fl[3]+1.125*fr[0]*nuVtSqSum[1]-1.125*fl[0]*nuVtSqSum[1]+1.125*nuVtSqSum[0]*fr[1]-1.125*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[2]*(0.4330127018922193*favg[16]+0.25*favg[6])+alphaDrag[3]*(0.4330127018922193*favg[8]+0.25*favg[2])-0.8660254037844386*fjump[7]+alphaDrag[0]*(0.4330127018922193*favg[7]+0.25*favg[1])+alphaDrag[1]*(0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[1]; 
  Ghat[2] = ((-1.082531754730548*nuVtSqSum[1]*fr[16])-1.082531754730548*nuVtSqSum[1]*fl[16]-1.082531754730548*nuVtSqSum[0]*fr[8]-1.082531754730548*nuVtSqSum[0]*fl[8]-1.082531754730548*nuVtSqSum[3]*fr[7]-1.082531754730548*nuVtSqSum[3]*fl[7]+1.125*nuVtSqSum[1]*fr[6]-1.125*nuVtSqSum[1]*fl[6]+1.125*fr[1]*nuVtSqSum[3]-1.125*fl[1]*nuVtSqSum[3]-1.082531754730548*nuVtSqSum[2]*fr[3]-1.082531754730548*nuVtSqSum[2]*fl[3]+1.125*fr[0]*nuVtSqSum[2]-1.125*fl[0]*nuVtSqSum[2]+1.125*nuVtSqSum[0]*fr[2]-1.125*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[1]*(0.4330127018922193*favg[16]+0.25*favg[6])-0.8660254037844386*fjump[8]+alphaDrag[0]*(0.4330127018922193*favg[8]+0.25*favg[2])+alphaDrag[3]*(0.4330127018922193*favg[7]+0.25*favg[1])+alphaDrag[2]*(0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[2]; 
  Ghat[4] = ((-1.082531754730548*nuVtSqSum[3]*fr[26])-1.082531754730548*nuVtSqSum[3]*fl[26]-1.082531754730548*nuVtSqSum[2]*fr[19]-1.082531754730548*nuVtSqSum[2]*fl[19]-1.082531754730548*nuVtSqSum[1]*fr[18]-1.082531754730548*nuVtSqSum[1]*fl[18]+1.125*nuVtSqSum[3]*fr[17]-1.125*nuVtSqSum[3]*fl[17]-1.082531754730548*nuVtSqSum[0]*fr[11]-1.082531754730548*nuVtSqSum[0]*fl[11]+1.125*nuVtSqSum[2]*fr[10]-1.125*nuVtSqSum[2]*fl[10]+1.125*nuVtSqSum[1]*fr[9]-1.125*nuVtSqSum[1]*fl[9]+1.125*nuVtSqSum[0]*fr[4]-1.125*nuVtSqSum[0]*fl[4])*rdv+alphaDrag[3]*(0.4330127018922193*favg[26]+0.25*favg[17])+alphaDrag[2]*(0.4330127018922193*favg[19]+0.25*favg[10])+alphaDrag[1]*(0.4330127018922193*favg[18]+0.25*favg[9])-0.8660254037844386*fjump[11]+alphaDrag[0]*(0.4330127018922193*favg[11]+0.25*favg[4])-0.5*fjump[4]; 
  Ghat[5] = ((-1.082531754730548*nuVtSqSum[3]*fr[27])-1.082531754730548*nuVtSqSum[3]*fl[27]-1.082531754730548*nuVtSqSum[2]*fr[22]-1.082531754730548*nuVtSqSum[2]*fl[22]-1.082531754730548*nuVtSqSum[1]*fr[21]-1.082531754730548*nuVtSqSum[1]*fl[21]+1.125*nuVtSqSum[3]*fr[20]-1.125*nuVtSqSum[3]*fl[20]-1.082531754730548*nuVtSqSum[0]*fr[14]-1.082531754730548*nuVtSqSum[0]*fl[14]+1.125*nuVtSqSum[2]*fr[13]-1.125*nuVtSqSum[2]*fl[13]+1.125*nuVtSqSum[1]*fr[12]-1.125*nuVtSqSum[1]*fl[12]+1.125*nuVtSqSum[0]*fr[5]-1.125*nuVtSqSum[0]*fl[5])*rdv+alphaDrag[3]*(0.4330127018922193*favg[27]+0.25*favg[20])+alphaDrag[2]*(0.4330127018922193*favg[22]+0.25*favg[13])+alphaDrag[1]*(0.4330127018922193*favg[21]+0.25*favg[12])-0.8660254037844386*fjump[14]+alphaDrag[0]*(0.4330127018922193*favg[14]+0.25*favg[5])-0.5*fjump[5]; 
  Ghat[6] = ((-1.082531754730548*nuVtSqSum[0]*fr[16])-1.082531754730548*nuVtSqSum[0]*fl[16]-1.082531754730548*nuVtSqSum[1]*fr[8]-1.082531754730548*nuVtSqSum[1]*fl[8]-1.082531754730548*nuVtSqSum[2]*fr[7]-1.082531754730548*nuVtSqSum[2]*fl[7]+1.125*nuVtSqSum[0]*fr[6]-1.125*nuVtSqSum[0]*fl[6]-1.082531754730548*fr[3]*nuVtSqSum[3]-1.082531754730548*fl[3]*nuVtSqSum[3]+1.125*fr[0]*nuVtSqSum[3]-1.125*fl[0]*nuVtSqSum[3]+1.125*fr[1]*nuVtSqSum[2]-1.125*fl[1]*nuVtSqSum[2]+1.125*nuVtSqSum[1]*fr[2]-1.125*nuVtSqSum[1]*fl[2])*rdv-0.8660254037844386*fjump[16]+alphaDrag[0]*(0.4330127018922193*favg[16]+0.25*favg[6])+alphaDrag[1]*(0.4330127018922193*favg[8]+0.25*favg[2])+alphaDrag[2]*(0.4330127018922193*favg[7]+0.25*favg[1])-0.5*fjump[6]+alphaDrag[3]*(0.4330127018922193*favg[3]+0.25*favg[0]); 
  Ghat[9] = ((-1.082531754730548*nuVtSqSum[2]*fr[26])-1.082531754730548*nuVtSqSum[2]*fl[26]-1.082531754730548*nuVtSqSum[3]*fr[19]-1.082531754730548*nuVtSqSum[3]*fl[19]-1.082531754730548*nuVtSqSum[0]*fr[18]-1.082531754730548*nuVtSqSum[0]*fl[18]+1.125*nuVtSqSum[2]*fr[17]-1.125*nuVtSqSum[2]*fl[17]-1.082531754730548*nuVtSqSum[1]*fr[11]-1.082531754730548*nuVtSqSum[1]*fl[11]+1.125*nuVtSqSum[3]*fr[10]-1.125*nuVtSqSum[3]*fl[10]+1.125*nuVtSqSum[0]*fr[9]-1.125*nuVtSqSum[0]*fl[9]+1.125*nuVtSqSum[1]*fr[4]-1.125*nuVtSqSum[1]*fl[4])*rdv+alphaDrag[2]*(0.4330127018922193*favg[26]+0.25*favg[17])+alphaDrag[3]*(0.4330127018922193*favg[19]+0.25*favg[10])-0.8660254037844386*fjump[18]+alphaDrag[0]*(0.4330127018922193*favg[18]+0.25*favg[9])+alphaDrag[1]*(0.4330127018922193*favg[11]+0.25*favg[4])-0.5*fjump[9]; 
  Ghat[10] = ((-1.082531754730548*nuVtSqSum[1]*fr[26])-1.082531754730548*nuVtSqSum[1]*fl[26]-1.082531754730548*nuVtSqSum[0]*fr[19]-1.082531754730548*nuVtSqSum[0]*fl[19]-1.082531754730548*nuVtSqSum[3]*fr[18]-1.082531754730548*nuVtSqSum[3]*fl[18]+1.125*nuVtSqSum[1]*fr[17]-1.125*nuVtSqSum[1]*fl[17]-1.082531754730548*nuVtSqSum[2]*fr[11]-1.082531754730548*nuVtSqSum[2]*fl[11]+1.125*nuVtSqSum[0]*fr[10]-1.125*nuVtSqSum[0]*fl[10]+1.125*nuVtSqSum[3]*fr[9]-1.125*nuVtSqSum[3]*fl[9]+1.125*nuVtSqSum[2]*fr[4]-1.125*nuVtSqSum[2]*fl[4])*rdv+alphaDrag[1]*(0.4330127018922193*favg[26]+0.25*favg[17])-0.8660254037844386*fjump[19]+alphaDrag[0]*(0.4330127018922193*favg[19]+0.25*favg[10])+alphaDrag[3]*(0.4330127018922193*favg[18]+0.25*favg[9])+alphaDrag[2]*(0.4330127018922193*favg[11]+0.25*favg[4])-0.5*fjump[10]; 
  Ghat[12] = ((-1.082531754730548*nuVtSqSum[2]*fr[27])-1.082531754730548*nuVtSqSum[2]*fl[27]-1.082531754730548*nuVtSqSum[3]*fr[22]-1.082531754730548*nuVtSqSum[3]*fl[22]-1.082531754730548*nuVtSqSum[0]*fr[21]-1.082531754730548*nuVtSqSum[0]*fl[21]+1.125*nuVtSqSum[2]*fr[20]-1.125*nuVtSqSum[2]*fl[20]-1.082531754730548*nuVtSqSum[1]*fr[14]-1.082531754730548*nuVtSqSum[1]*fl[14]+1.125*nuVtSqSum[3]*fr[13]-1.125*nuVtSqSum[3]*fl[13]+1.125*nuVtSqSum[0]*fr[12]-1.125*nuVtSqSum[0]*fl[12]+1.125*nuVtSqSum[1]*fr[5]-1.125*nuVtSqSum[1]*fl[5])*rdv+alphaDrag[2]*(0.4330127018922193*favg[27]+0.25*favg[20])+alphaDrag[3]*(0.4330127018922193*favg[22]+0.25*favg[13])-0.8660254037844386*fjump[21]+alphaDrag[0]*(0.4330127018922193*favg[21]+0.25*favg[12])+alphaDrag[1]*(0.4330127018922193*favg[14]+0.25*favg[5])-0.5*fjump[12]; 
  Ghat[13] = ((-1.082531754730548*nuVtSqSum[1]*fr[27])-1.082531754730548*nuVtSqSum[1]*fl[27]-1.082531754730548*nuVtSqSum[0]*fr[22]-1.082531754730548*nuVtSqSum[0]*fl[22]-1.082531754730548*nuVtSqSum[3]*fr[21]-1.082531754730548*nuVtSqSum[3]*fl[21]+1.125*nuVtSqSum[1]*fr[20]-1.125*nuVtSqSum[1]*fl[20]-1.082531754730548*nuVtSqSum[2]*fr[14]-1.082531754730548*nuVtSqSum[2]*fl[14]+1.125*nuVtSqSum[0]*fr[13]-1.125*nuVtSqSum[0]*fl[13]+1.125*nuVtSqSum[3]*fr[12]-1.125*nuVtSqSum[3]*fl[12]+1.125*nuVtSqSum[2]*fr[5]-1.125*nuVtSqSum[2]*fl[5])*rdv+alphaDrag[1]*(0.4330127018922193*favg[27]+0.25*favg[20])-0.8660254037844386*fjump[22]+alphaDrag[0]*(0.4330127018922193*favg[22]+0.25*favg[13])+alphaDrag[3]*(0.4330127018922193*favg[21]+0.25*favg[12])+alphaDrag[2]*(0.4330127018922193*favg[14]+0.25*favg[5])-0.5*fjump[13]; 
  Ghat[15] = ((-1.082531754730548*nuVtSqSum[3]*fr[31])-1.082531754730548*nuVtSqSum[3]*fl[31]-1.082531754730548*nuVtSqSum[2]*fr[30]-1.082531754730548*nuVtSqSum[2]*fl[30]-1.082531754730548*nuVtSqSum[1]*fr[29]-1.082531754730548*nuVtSqSum[1]*fl[29]+1.125*nuVtSqSum[3]*fr[28]-1.125*nuVtSqSum[3]*fl[28]-1.082531754730548*nuVtSqSum[0]*fr[25]-1.082531754730548*nuVtSqSum[0]*fl[25]+1.125*nuVtSqSum[2]*fr[24]-1.125*nuVtSqSum[2]*fl[24]+1.125*nuVtSqSum[1]*fr[23]-1.125*nuVtSqSum[1]*fl[23]+1.125*nuVtSqSum[0]*fr[15]-1.125*nuVtSqSum[0]*fl[15])*rdv+alphaDrag[3]*(0.4330127018922193*favg[31]+0.25*favg[28])+alphaDrag[2]*(0.4330127018922193*favg[30]+0.25*favg[24])+alphaDrag[1]*(0.4330127018922193*favg[29]+0.25*favg[23])-0.8660254037844386*fjump[25]+alphaDrag[0]*(0.4330127018922193*favg[25]+0.25*favg[15])-0.5*fjump[15]; 
  Ghat[17] = ((-1.082531754730548*nuVtSqSum[0]*fr[26])-1.082531754730548*nuVtSqSum[0]*fl[26]-1.082531754730548*nuVtSqSum[1]*fr[19]-1.082531754730548*nuVtSqSum[1]*fl[19]-1.082531754730548*nuVtSqSum[2]*fr[18]-1.082531754730548*nuVtSqSum[2]*fl[18]+1.125*nuVtSqSum[0]*fr[17]-1.125*nuVtSqSum[0]*fl[17]-1.082531754730548*nuVtSqSum[3]*fr[11]-1.082531754730548*nuVtSqSum[3]*fl[11]+1.125*nuVtSqSum[1]*fr[10]-1.125*nuVtSqSum[1]*fl[10]+1.125*nuVtSqSum[2]*fr[9]-1.125*nuVtSqSum[2]*fl[9]+1.125*nuVtSqSum[3]*fr[4]-1.125*nuVtSqSum[3]*fl[4])*rdv-0.8660254037844386*fjump[26]+alphaDrag[0]*(0.4330127018922193*favg[26]+0.25*favg[17])+alphaDrag[1]*(0.4330127018922193*favg[19]+0.25*favg[10])+alphaDrag[2]*(0.4330127018922193*favg[18]+0.25*favg[9])-0.5*fjump[17]+alphaDrag[3]*(0.4330127018922193*favg[11]+0.25*favg[4]); 
  Ghat[20] = ((-1.082531754730548*nuVtSqSum[0]*fr[27])-1.082531754730548*nuVtSqSum[0]*fl[27]-1.082531754730548*nuVtSqSum[1]*fr[22]-1.082531754730548*nuVtSqSum[1]*fl[22]-1.082531754730548*nuVtSqSum[2]*fr[21]-1.082531754730548*nuVtSqSum[2]*fl[21]+1.125*nuVtSqSum[0]*fr[20]-1.125*nuVtSqSum[0]*fl[20]-1.082531754730548*nuVtSqSum[3]*fr[14]-1.082531754730548*nuVtSqSum[3]*fl[14]+1.125*nuVtSqSum[1]*fr[13]-1.125*nuVtSqSum[1]*fl[13]+1.125*nuVtSqSum[2]*fr[12]-1.125*nuVtSqSum[2]*fl[12]+1.125*nuVtSqSum[3]*fr[5]-1.125*nuVtSqSum[3]*fl[5])*rdv-0.8660254037844386*fjump[27]+alphaDrag[0]*(0.4330127018922193*favg[27]+0.25*favg[20])+alphaDrag[1]*(0.4330127018922193*favg[22]+0.25*favg[13])+alphaDrag[2]*(0.4330127018922193*favg[21]+0.25*favg[12])-0.5*fjump[20]+alphaDrag[3]*(0.4330127018922193*favg[14]+0.25*favg[5]); 
  Ghat[23] = ((-1.082531754730548*nuVtSqSum[2]*fr[31])-1.082531754730548*nuVtSqSum[2]*fl[31]-1.082531754730548*nuVtSqSum[3]*fr[30]-1.082531754730548*nuVtSqSum[3]*fl[30]-1.082531754730548*nuVtSqSum[0]*fr[29]-1.082531754730548*nuVtSqSum[0]*fl[29]+1.125*nuVtSqSum[2]*fr[28]-1.125*nuVtSqSum[2]*fl[28]-1.082531754730548*nuVtSqSum[1]*fr[25]-1.082531754730548*nuVtSqSum[1]*fl[25]+1.125*nuVtSqSum[3]*fr[24]-1.125*nuVtSqSum[3]*fl[24]+1.125*nuVtSqSum[0]*fr[23]-1.125*nuVtSqSum[0]*fl[23]+1.125*nuVtSqSum[1]*fr[15]-1.125*nuVtSqSum[1]*fl[15])*rdv+alphaDrag[2]*(0.4330127018922193*favg[31]+0.25*favg[28])+alphaDrag[3]*(0.4330127018922193*favg[30]+0.25*favg[24])-0.8660254037844386*fjump[29]+alphaDrag[0]*(0.4330127018922193*favg[29]+0.25*favg[23])+alphaDrag[1]*(0.4330127018922193*favg[25]+0.25*favg[15])-0.5*fjump[23]; 
  Ghat[24] = ((-1.082531754730548*nuVtSqSum[1]*fr[31])-1.082531754730548*nuVtSqSum[1]*fl[31]-1.082531754730548*nuVtSqSum[0]*fr[30]-1.082531754730548*nuVtSqSum[0]*fl[30]-1.082531754730548*nuVtSqSum[3]*fr[29]-1.082531754730548*nuVtSqSum[3]*fl[29]+1.125*nuVtSqSum[1]*fr[28]-1.125*nuVtSqSum[1]*fl[28]-1.082531754730548*nuVtSqSum[2]*fr[25]-1.082531754730548*nuVtSqSum[2]*fl[25]+1.125*nuVtSqSum[0]*fr[24]-1.125*nuVtSqSum[0]*fl[24]+1.125*nuVtSqSum[3]*fr[23]-1.125*nuVtSqSum[3]*fl[23]+1.125*nuVtSqSum[2]*fr[15]-1.125*nuVtSqSum[2]*fl[15])*rdv+alphaDrag[1]*(0.4330127018922193*favg[31]+0.25*favg[28])-0.8660254037844386*fjump[30]+alphaDrag[0]*(0.4330127018922193*favg[30]+0.25*favg[24])+alphaDrag[3]*(0.4330127018922193*favg[29]+0.25*favg[23])+alphaDrag[2]*(0.4330127018922193*favg[25]+0.25*favg[15])-0.5*fjump[24]; 
  Ghat[28] = ((-1.082531754730548*nuVtSqSum[0]*fr[31])-1.082531754730548*nuVtSqSum[0]*fl[31]-1.082531754730548*nuVtSqSum[1]*fr[30]-1.082531754730548*nuVtSqSum[1]*fl[30]-1.082531754730548*nuVtSqSum[2]*fr[29]-1.082531754730548*nuVtSqSum[2]*fl[29]+1.125*nuVtSqSum[0]*fr[28]-1.125*nuVtSqSum[0]*fl[28]-1.082531754730548*nuVtSqSum[3]*fr[25]-1.082531754730548*nuVtSqSum[3]*fl[25]+1.125*nuVtSqSum[1]*fr[24]-1.125*nuVtSqSum[1]*fl[24]+1.125*nuVtSqSum[2]*fr[23]-1.125*nuVtSqSum[2]*fl[23]+1.125*nuVtSqSum[3]*fr[15]-1.125*nuVtSqSum[3]*fl[15])*rdv-0.8660254037844386*fjump[31]+alphaDrag[0]*(0.4330127018922193*favg[31]+0.25*favg[28])+alphaDrag[1]*(0.4330127018922193*favg[30]+0.25*favg[24])+alphaDrag[2]*(0.4330127018922193*favg[29]+0.25*favg[23])-0.5*fjump[28]+alphaDrag[3]*(0.4330127018922193*favg[25]+0.25*favg[15]); 

  double incr1[32]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = -0.5*Ghat[6]; 
  incr1[7] = 0.8660254037844386*Ghat[1]; 
  incr1[8] = 0.8660254037844386*Ghat[2]; 
  incr1[9] = -0.5*Ghat[9]; 
  incr1[10] = -0.5*Ghat[10]; 
  incr1[11] = 0.8660254037844386*Ghat[4]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = -0.5*Ghat[13]; 
  incr1[14] = 0.8660254037844386*Ghat[5]; 
  incr1[15] = -0.5*Ghat[15]; 
  incr1[16] = 0.8660254037844386*Ghat[6]; 
  incr1[17] = -0.5*Ghat[17]; 
  incr1[18] = 0.8660254037844386*Ghat[9]; 
  incr1[19] = 0.8660254037844386*Ghat[10]; 
  incr1[20] = -0.5*Ghat[20]; 
  incr1[21] = 0.8660254037844386*Ghat[12]; 
  incr1[22] = 0.8660254037844386*Ghat[13]; 
  incr1[23] = -0.5*Ghat[23]; 
  incr1[24] = -0.5*Ghat[24]; 
  incr1[25] = 0.8660254037844386*Ghat[15]; 
  incr1[26] = 0.8660254037844386*Ghat[17]; 
  incr1[27] = 0.8660254037844386*Ghat[20]; 
  incr1[28] = -0.5*Ghat[28]; 
  incr1[29] = 0.8660254037844386*Ghat[23]; 
  incr1[30] = 0.8660254037844386*Ghat[24]; 
  incr1[31] = 0.8660254037844386*Ghat[28]; 

  double incr2[32]; 
  incr2[3] = nuVtSqSum[3]*((-0.25*fr[16])+0.25*fl[16]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[2]*((-0.25*fr[8])+0.25*fl[8]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[7] = nuVtSqSum[2]*((-0.25*fr[16])+0.25*fl[16]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[3]*((-0.25*fr[8])+0.25*fl[8]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[8] = nuVtSqSum[1]*((-0.25*fr[16])+0.25*fl[16]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[0]*((-0.25*fr[8])+0.25*fl[8]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[3]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[11] = nuVtSqSum[3]*((-0.25*fr[26])+0.25*fl[26]+0.2165063509461096*(fr[17]+fl[17]))+nuVtSqSum[2]*((-0.25*fr[19])+0.25*fl[19]+0.2165063509461096*(fr[10]+fl[10]))+nuVtSqSum[1]*((-0.25*fr[18])+0.25*fl[18]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[0]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[14] = nuVtSqSum[3]*((-0.25*fr[27])+0.25*fl[27]+0.2165063509461096*(fr[20]+fl[20]))+nuVtSqSum[2]*((-0.25*fr[22])+0.25*fl[22]+0.2165063509461096*(fr[13]+fl[13]))+nuVtSqSum[1]*((-0.25*fr[21])+0.25*fl[21]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[0]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[5]+fl[5])); 
  incr2[16] = nuVtSqSum[0]*((-0.25*fr[16])+0.25*fl[16]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[1]*((-0.25*fr[8])+0.25*fl[8]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[2]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[1]+fl[1]))+((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0]))*nuVtSqSum[3]; 
  incr2[18] = nuVtSqSum[2]*((-0.25*fr[26])+0.25*fl[26]+0.2165063509461096*(fr[17]+fl[17]))+nuVtSqSum[3]*((-0.25*fr[19])+0.25*fl[19]+0.2165063509461096*(fr[10]+fl[10]))+nuVtSqSum[0]*((-0.25*fr[18])+0.25*fl[18]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[1]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[19] = nuVtSqSum[1]*((-0.25*fr[26])+0.25*fl[26]+0.2165063509461096*(fr[17]+fl[17]))+nuVtSqSum[0]*((-0.25*fr[19])+0.25*fl[19]+0.2165063509461096*(fr[10]+fl[10]))+nuVtSqSum[3]*((-0.25*fr[18])+0.25*fl[18]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[2]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[21] = nuVtSqSum[2]*((-0.25*fr[27])+0.25*fl[27]+0.2165063509461096*(fr[20]+fl[20]))+nuVtSqSum[3]*((-0.25*fr[22])+0.25*fl[22]+0.2165063509461096*(fr[13]+fl[13]))+nuVtSqSum[0]*((-0.25*fr[21])+0.25*fl[21]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[1]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[5]+fl[5])); 
  incr2[22] = nuVtSqSum[1]*((-0.25*fr[27])+0.25*fl[27]+0.2165063509461096*(fr[20]+fl[20]))+nuVtSqSum[0]*((-0.25*fr[22])+0.25*fl[22]+0.2165063509461096*(fr[13]+fl[13]))+nuVtSqSum[3]*((-0.25*fr[21])+0.25*fl[21]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[2]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[5]+fl[5])); 
  incr2[25] = nuVtSqSum[3]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[28]+fl[28]))+nuVtSqSum[2]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[24]+fl[24]))+nuVtSqSum[1]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[23]+fl[23]))+nuVtSqSum[0]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[15]+fl[15])); 
  incr2[26] = nuVtSqSum[0]*((-0.25*fr[26])+0.25*fl[26]+0.2165063509461096*(fr[17]+fl[17]))+nuVtSqSum[1]*((-0.25*fr[19])+0.25*fl[19]+0.2165063509461096*(fr[10]+fl[10]))+nuVtSqSum[2]*((-0.25*fr[18])+0.25*fl[18]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[3]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[27] = nuVtSqSum[0]*((-0.25*fr[27])+0.25*fl[27]+0.2165063509461096*(fr[20]+fl[20]))+nuVtSqSum[1]*((-0.25*fr[22])+0.25*fl[22]+0.2165063509461096*(fr[13]+fl[13]))+nuVtSqSum[2]*((-0.25*fr[21])+0.25*fl[21]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[3]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[5]+fl[5])); 
  incr2[29] = nuVtSqSum[2]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[28]+fl[28]))+nuVtSqSum[3]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[24]+fl[24]))+nuVtSqSum[0]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[23]+fl[23]))+nuVtSqSum[1]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[15]+fl[15])); 
  incr2[30] = nuVtSqSum[1]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[28]+fl[28]))+nuVtSqSum[0]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[24]+fl[24]))+nuVtSqSum[3]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[23]+fl[23]))+nuVtSqSum[2]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[15]+fl[15])); 
  incr2[31] = nuVtSqSum[0]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[28]+fl[28]))+nuVtSqSum[1]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[24]+fl[24]))+nuVtSqSum[2]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[23]+fl[23]))+nuVtSqSum[3]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[15]+fl[15])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr2[8]*rdvSq4R+incr1[8]*rdv2R; 
  outr[9] += incr1[9]*rdv2R; 
  outr[10] += incr1[10]*rdv2R; 
  outr[11] += incr2[11]*rdvSq4R+incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr1[15]*rdv2R; 
  outr[16] += incr2[16]*rdvSq4R+incr1[16]*rdv2R; 
  outr[17] += incr1[17]*rdv2R; 
  outr[18] += incr2[18]*rdvSq4R+incr1[18]*rdv2R; 
  outr[19] += incr2[19]*rdvSq4R+incr1[19]*rdv2R; 
  outr[20] += incr1[20]*rdv2R; 
  outr[21] += incr2[21]*rdvSq4R+incr1[21]*rdv2R; 
  outr[22] += incr2[22]*rdvSq4R+incr1[22]*rdv2R; 
  outr[23] += incr1[23]*rdv2R; 
  outr[24] += incr1[24]*rdv2R; 
  outr[25] += incr2[25]*rdvSq4R+incr1[25]*rdv2R; 
  outr[26] += incr2[26]*rdvSq4R+incr1[26]*rdv2R; 
  outr[27] += incr2[27]*rdvSq4R+incr1[27]*rdv2R; 
  outr[28] += incr1[28]*rdv2R; 
  outr[29] += incr2[29]*rdvSq4R+incr1[29]*rdv2R; 
  outr[30] += incr2[30]*rdvSq4R+incr1[30]*rdv2R; 
  outr[31] += incr2[31]*rdvSq4R+incr1[31]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 
  outl[8] += incr1[8]*rdv2L-1.0*incr2[8]*rdvSq4L; 
  outl[9] += -1.0*incr1[9]*rdv2L; 
  outl[10] += -1.0*incr1[10]*rdv2L; 
  outl[11] += incr1[11]*rdv2L-1.0*incr2[11]*rdvSq4L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += -1.0*incr1[13]*rdv2L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += -1.0*incr1[15]*rdv2L; 
  outl[16] += incr1[16]*rdv2L-1.0*incr2[16]*rdvSq4L; 
  outl[17] += -1.0*incr1[17]*rdv2L; 
  outl[18] += incr1[18]*rdv2L-1.0*incr2[18]*rdvSq4L; 
  outl[19] += incr1[19]*rdv2L-1.0*incr2[19]*rdvSq4L; 
  outl[20] += -1.0*incr1[20]*rdv2L; 
  outl[21] += incr1[21]*rdv2L-1.0*incr2[21]*rdvSq4L; 
  outl[22] += incr1[22]*rdv2L-1.0*incr2[22]*rdvSq4L; 
  outl[23] += -1.0*incr1[23]*rdv2L; 
  outl[24] += -1.0*incr1[24]*rdv2L; 
  outl[25] += incr1[25]*rdv2L-1.0*incr2[25]*rdvSq4L; 
  outl[26] += incr1[26]*rdv2L-1.0*incr2[26]*rdvSq4L; 
  outl[27] += incr1[27]*rdv2L-1.0*incr2[27]*rdvSq4L; 
  outl[28] += -1.0*incr1[28]*rdv2L; 
  outl[29] += incr1[29]*rdv2L-1.0*incr2[29]*rdvSq4L; 
  outl[30] += incr1[30]*rdv2L-1.0*incr2[30]*rdvSq4L; 
  outl[31] += incr1[31]*rdv2L-1.0*incr2[31]*rdvSq4L; 

  return std::abs(wl[2]-(0.5*sumNuUx[0])/nuSum); 
} 
double VmLBOconstNuSurf2x3vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:          Cell-center coordinates. 
  // dxv[5]:        Cell spacing. 
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

  const double *sumNuUy = &nuUSum[4]; 

  double favg[32]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = -1*fr[17]+fl[17]; 
  favg[18] = -1*fr[18]+fl[18]; 
  favg[19] = -1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  favg[21] = 1*fr[21]+fl[21]; 
  favg[22] = 1*fr[22]+fl[22]; 
  favg[23] = -1*fr[23]+fl[23]; 
  favg[24] = -1*fr[24]+fl[24]; 
  favg[25] = -1*fr[25]+fl[25]; 
  favg[26] = -1*fr[26]+fl[26]; 
  favg[27] = 1*fr[27]+fl[27]; 
  favg[28] = -1*fr[28]+fl[28]; 
  favg[29] = -1*fr[29]+fl[29]; 
  favg[30] = -1*fr[30]+fl[30]; 
  favg[31] = -1*fr[31]+fl[31]; 

  double fjump[32]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(-1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(-1*fr[15])); 
  fjump[16] = nuSum*vMuMidMax*(fl[16]-(1*fr[16])); 
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(-1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(-1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(-1*fr[19])); 
  fjump[20] = nuSum*vMuMidMax*(fl[20]-(1*fr[20])); 
  fjump[21] = nuSum*vMuMidMax*(fl[21]-(1*fr[21])); 
  fjump[22] = nuSum*vMuMidMax*(fl[22]-(1*fr[22])); 
  fjump[23] = nuSum*vMuMidMax*(fl[23]-(-1*fr[23])); 
  fjump[24] = nuSum*vMuMidMax*(fl[24]-(-1*fr[24])); 
  fjump[25] = nuSum*vMuMidMax*(fl[25]-(-1*fr[25])); 
  fjump[26] = nuSum*vMuMidMax*(fl[26]-(-1*fr[26])); 
  fjump[27] = nuSum*vMuMidMax*(fl[27]-(1*fr[27])); 
  fjump[28] = nuSum*vMuMidMax*(fl[28]-(-1*fr[28])); 
  fjump[29] = nuSum*vMuMidMax*(fl[29]-(-1*fr[29])); 
  fjump[30] = nuSum*vMuMidMax*(fl[30]-(-1*fr[30])); 
  fjump[31] = nuSum*vMuMidMax*(fl[31]-(-1*fr[31])); 

  double alphaDrag[4]; 
  alphaDrag[0] = 2.0*wl[3]*nuSum+dxvl[3]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 
  alphaDrag[3] = -1.0*sumNuUy[3]; 

  double Ghat[32]; 
  for(unsigned short int i=0; i<32; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.082531754730548*nuVtSqSum[3]*fr[17])-1.082531754730548*nuVtSqSum[3]*fl[17]-1.082531754730548*nuVtSqSum[2]*fr[10]-1.082531754730548*nuVtSqSum[2]*fl[10]-1.082531754730548*nuVtSqSum[1]*fr[9]-1.082531754730548*nuVtSqSum[1]*fl[9]+1.125*nuVtSqSum[3]*fr[6]-1.125*nuVtSqSum[3]*fl[6]-1.082531754730548*nuVtSqSum[0]*fr[4]-1.082531754730548*nuVtSqSum[0]*fl[4]+1.125*fr[2]*nuVtSqSum[2]-1.125*fl[2]*nuVtSqSum[2]+1.125*fr[1]*nuVtSqSum[1]-1.125*fl[1]*nuVtSqSum[1]+1.125*fr[0]*nuVtSqSum[0]-1.125*fl[0]*nuVtSqSum[0])*rdv+alphaDrag[3]*(0.4330127018922193*favg[17]+0.25*favg[6])+alphaDrag[2]*(0.4330127018922193*favg[10]+0.25*favg[2])+alphaDrag[1]*(0.4330127018922193*favg[9]+0.25*favg[1])-0.8660254037844386*fjump[4]+alphaDrag[0]*(0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[0]; 
  Ghat[1] = ((-1.082531754730548*nuVtSqSum[2]*fr[17])-1.082531754730548*nuVtSqSum[2]*fl[17]-1.082531754730548*nuVtSqSum[3]*fr[10]-1.082531754730548*nuVtSqSum[3]*fl[10]-1.082531754730548*nuVtSqSum[0]*fr[9]-1.082531754730548*nuVtSqSum[0]*fl[9]+1.125*nuVtSqSum[2]*fr[6]-1.125*nuVtSqSum[2]*fl[6]-1.082531754730548*nuVtSqSum[1]*fr[4]-1.082531754730548*nuVtSqSum[1]*fl[4]+1.125*fr[2]*nuVtSqSum[3]-1.125*fl[2]*nuVtSqSum[3]+1.125*fr[0]*nuVtSqSum[1]-1.125*fl[0]*nuVtSqSum[1]+1.125*nuVtSqSum[0]*fr[1]-1.125*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[2]*(0.4330127018922193*favg[17]+0.25*favg[6])+alphaDrag[3]*(0.4330127018922193*favg[10]+0.25*favg[2])-0.8660254037844386*fjump[9]+alphaDrag[0]*(0.4330127018922193*favg[9]+0.25*favg[1])+alphaDrag[1]*(0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[1]; 
  Ghat[2] = ((-1.082531754730548*nuVtSqSum[1]*fr[17])-1.082531754730548*nuVtSqSum[1]*fl[17]-1.082531754730548*nuVtSqSum[0]*fr[10]-1.082531754730548*nuVtSqSum[0]*fl[10]-1.082531754730548*nuVtSqSum[3]*fr[9]-1.082531754730548*nuVtSqSum[3]*fl[9]+1.125*nuVtSqSum[1]*fr[6]-1.125*nuVtSqSum[1]*fl[6]-1.082531754730548*nuVtSqSum[2]*fr[4]-1.082531754730548*nuVtSqSum[2]*fl[4]+1.125*fr[1]*nuVtSqSum[3]-1.125*fl[1]*nuVtSqSum[3]+1.125*fr[0]*nuVtSqSum[2]-1.125*fl[0]*nuVtSqSum[2]+1.125*nuVtSqSum[0]*fr[2]-1.125*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[1]*(0.4330127018922193*favg[17]+0.25*favg[6])-0.8660254037844386*fjump[10]+alphaDrag[0]*(0.4330127018922193*favg[10]+0.25*favg[2])+alphaDrag[3]*(0.4330127018922193*favg[9]+0.25*favg[1])+alphaDrag[2]*(0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[2]; 
  Ghat[3] = ((-1.082531754730548*nuVtSqSum[3]*fr[26])-1.082531754730548*nuVtSqSum[3]*fl[26]-1.082531754730548*nuVtSqSum[2]*fr[19]-1.082531754730548*nuVtSqSum[2]*fl[19]-1.082531754730548*nuVtSqSum[1]*fr[18]-1.082531754730548*nuVtSqSum[1]*fl[18]+1.125*nuVtSqSum[3]*fr[16]-1.125*nuVtSqSum[3]*fl[16]-1.082531754730548*nuVtSqSum[0]*fr[11]-1.082531754730548*nuVtSqSum[0]*fl[11]+1.125*nuVtSqSum[2]*fr[8]-1.125*nuVtSqSum[2]*fl[8]+1.125*nuVtSqSum[1]*fr[7]-1.125*nuVtSqSum[1]*fl[7]+1.125*nuVtSqSum[0]*fr[3]-1.125*nuVtSqSum[0]*fl[3])*rdv+alphaDrag[3]*(0.4330127018922193*favg[26]+0.25*favg[16])+alphaDrag[2]*(0.4330127018922193*favg[19]+0.25*favg[8])+alphaDrag[1]*(0.4330127018922193*favg[18]+0.25*favg[7])-0.8660254037844386*fjump[11]+alphaDrag[0]*(0.4330127018922193*favg[11]+0.25*favg[3])-0.5*fjump[3]; 
  Ghat[5] = ((-1.082531754730548*nuVtSqSum[3]*fr[28])-1.082531754730548*nuVtSqSum[3]*fl[28]-1.082531754730548*nuVtSqSum[2]*fr[24]-1.082531754730548*nuVtSqSum[2]*fl[24]-1.082531754730548*nuVtSqSum[1]*fr[23]-1.082531754730548*nuVtSqSum[1]*fl[23]+1.125*nuVtSqSum[3]*fr[20]-1.125*nuVtSqSum[3]*fl[20]-1.082531754730548*nuVtSqSum[0]*fr[15]-1.082531754730548*nuVtSqSum[0]*fl[15]+1.125*nuVtSqSum[2]*fr[13]-1.125*nuVtSqSum[2]*fl[13]+1.125*nuVtSqSum[1]*fr[12]-1.125*nuVtSqSum[1]*fl[12]+1.125*nuVtSqSum[0]*fr[5]-1.125*nuVtSqSum[0]*fl[5])*rdv+alphaDrag[3]*(0.4330127018922193*favg[28]+0.25*favg[20])+alphaDrag[2]*(0.4330127018922193*favg[24]+0.25*favg[13])+alphaDrag[1]*(0.4330127018922193*favg[23]+0.25*favg[12])-0.8660254037844386*fjump[15]+alphaDrag[0]*(0.4330127018922193*favg[15]+0.25*favg[5])-0.5*fjump[5]; 
  Ghat[6] = ((-1.082531754730548*nuVtSqSum[0]*fr[17])-1.082531754730548*nuVtSqSum[0]*fl[17]-1.082531754730548*nuVtSqSum[1]*fr[10]-1.082531754730548*nuVtSqSum[1]*fl[10]-1.082531754730548*nuVtSqSum[2]*fr[9]-1.082531754730548*nuVtSqSum[2]*fl[9]+1.125*nuVtSqSum[0]*fr[6]-1.125*nuVtSqSum[0]*fl[6]-1.082531754730548*nuVtSqSum[3]*fr[4]-1.082531754730548*nuVtSqSum[3]*fl[4]+1.125*fr[0]*nuVtSqSum[3]-1.125*fl[0]*nuVtSqSum[3]+1.125*fr[1]*nuVtSqSum[2]-1.125*fl[1]*nuVtSqSum[2]+1.125*nuVtSqSum[1]*fr[2]-1.125*nuVtSqSum[1]*fl[2])*rdv-0.8660254037844386*fjump[17]+alphaDrag[0]*(0.4330127018922193*favg[17]+0.25*favg[6])+alphaDrag[1]*(0.4330127018922193*favg[10]+0.25*favg[2])+alphaDrag[2]*(0.4330127018922193*favg[9]+0.25*favg[1])-0.5*fjump[6]+alphaDrag[3]*(0.4330127018922193*favg[4]+0.25*favg[0]); 
  Ghat[7] = ((-1.082531754730548*nuVtSqSum[2]*fr[26])-1.082531754730548*nuVtSqSum[2]*fl[26]-1.082531754730548*nuVtSqSum[3]*fr[19]-1.082531754730548*nuVtSqSum[3]*fl[19]-1.082531754730548*nuVtSqSum[0]*fr[18]-1.082531754730548*nuVtSqSum[0]*fl[18]+1.125*nuVtSqSum[2]*fr[16]-1.125*nuVtSqSum[2]*fl[16]-1.082531754730548*nuVtSqSum[1]*fr[11]-1.082531754730548*nuVtSqSum[1]*fl[11]+1.125*nuVtSqSum[3]*fr[8]-1.125*nuVtSqSum[3]*fl[8]+1.125*nuVtSqSum[0]*fr[7]-1.125*nuVtSqSum[0]*fl[7]+1.125*nuVtSqSum[1]*fr[3]-1.125*nuVtSqSum[1]*fl[3])*rdv+alphaDrag[2]*(0.4330127018922193*favg[26]+0.25*favg[16])+alphaDrag[3]*(0.4330127018922193*favg[19]+0.25*favg[8])-0.8660254037844386*fjump[18]+alphaDrag[0]*(0.4330127018922193*favg[18]+0.25*favg[7])+alphaDrag[1]*(0.4330127018922193*favg[11]+0.25*favg[3])-0.5*fjump[7]; 
  Ghat[8] = ((-1.082531754730548*nuVtSqSum[1]*fr[26])-1.082531754730548*nuVtSqSum[1]*fl[26]-1.082531754730548*nuVtSqSum[0]*fr[19]-1.082531754730548*nuVtSqSum[0]*fl[19]-1.082531754730548*nuVtSqSum[3]*fr[18]-1.082531754730548*nuVtSqSum[3]*fl[18]+1.125*nuVtSqSum[1]*fr[16]-1.125*nuVtSqSum[1]*fl[16]-1.082531754730548*nuVtSqSum[2]*fr[11]-1.082531754730548*nuVtSqSum[2]*fl[11]+1.125*nuVtSqSum[0]*fr[8]-1.125*nuVtSqSum[0]*fl[8]+1.125*nuVtSqSum[3]*fr[7]-1.125*nuVtSqSum[3]*fl[7]+1.125*nuVtSqSum[2]*fr[3]-1.125*nuVtSqSum[2]*fl[3])*rdv+alphaDrag[1]*(0.4330127018922193*favg[26]+0.25*favg[16])-0.8660254037844386*fjump[19]+alphaDrag[0]*(0.4330127018922193*favg[19]+0.25*favg[8])+alphaDrag[3]*(0.4330127018922193*favg[18]+0.25*favg[7])+alphaDrag[2]*(0.4330127018922193*favg[11]+0.25*favg[3])-0.5*fjump[8]; 
  Ghat[12] = ((-1.082531754730548*nuVtSqSum[2]*fr[28])-1.082531754730548*nuVtSqSum[2]*fl[28]-1.082531754730548*nuVtSqSum[3]*fr[24]-1.082531754730548*nuVtSqSum[3]*fl[24]-1.082531754730548*nuVtSqSum[0]*fr[23]-1.082531754730548*nuVtSqSum[0]*fl[23]+1.125*nuVtSqSum[2]*fr[20]-1.125*nuVtSqSum[2]*fl[20]-1.082531754730548*nuVtSqSum[1]*fr[15]-1.082531754730548*nuVtSqSum[1]*fl[15]+1.125*nuVtSqSum[3]*fr[13]-1.125*nuVtSqSum[3]*fl[13]+1.125*nuVtSqSum[0]*fr[12]-1.125*nuVtSqSum[0]*fl[12]+1.125*nuVtSqSum[1]*fr[5]-1.125*nuVtSqSum[1]*fl[5])*rdv+alphaDrag[2]*(0.4330127018922193*favg[28]+0.25*favg[20])+alphaDrag[3]*(0.4330127018922193*favg[24]+0.25*favg[13])-0.8660254037844386*fjump[23]+alphaDrag[0]*(0.4330127018922193*favg[23]+0.25*favg[12])+alphaDrag[1]*(0.4330127018922193*favg[15]+0.25*favg[5])-0.5*fjump[12]; 
  Ghat[13] = ((-1.082531754730548*nuVtSqSum[1]*fr[28])-1.082531754730548*nuVtSqSum[1]*fl[28]-1.082531754730548*nuVtSqSum[0]*fr[24]-1.082531754730548*nuVtSqSum[0]*fl[24]-1.082531754730548*nuVtSqSum[3]*fr[23]-1.082531754730548*nuVtSqSum[3]*fl[23]+1.125*nuVtSqSum[1]*fr[20]-1.125*nuVtSqSum[1]*fl[20]-1.082531754730548*nuVtSqSum[2]*fr[15]-1.082531754730548*nuVtSqSum[2]*fl[15]+1.125*nuVtSqSum[0]*fr[13]-1.125*nuVtSqSum[0]*fl[13]+1.125*nuVtSqSum[3]*fr[12]-1.125*nuVtSqSum[3]*fl[12]+1.125*nuVtSqSum[2]*fr[5]-1.125*nuVtSqSum[2]*fl[5])*rdv+alphaDrag[1]*(0.4330127018922193*favg[28]+0.25*favg[20])-0.8660254037844386*fjump[24]+alphaDrag[0]*(0.4330127018922193*favg[24]+0.25*favg[13])+alphaDrag[3]*(0.4330127018922193*favg[23]+0.25*favg[12])+alphaDrag[2]*(0.4330127018922193*favg[15]+0.25*favg[5])-0.5*fjump[13]; 
  Ghat[14] = ((-1.082531754730548*nuVtSqSum[3]*fr[31])-1.082531754730548*nuVtSqSum[3]*fl[31]-1.082531754730548*nuVtSqSum[2]*fr[30]-1.082531754730548*nuVtSqSum[2]*fl[30]-1.082531754730548*nuVtSqSum[1]*fr[29]-1.082531754730548*nuVtSqSum[1]*fl[29]+1.125*nuVtSqSum[3]*fr[27]-1.125*nuVtSqSum[3]*fl[27]-1.082531754730548*nuVtSqSum[0]*fr[25]-1.082531754730548*nuVtSqSum[0]*fl[25]+1.125*nuVtSqSum[2]*fr[22]-1.125*nuVtSqSum[2]*fl[22]+1.125*nuVtSqSum[1]*fr[21]-1.125*nuVtSqSum[1]*fl[21]+1.125*nuVtSqSum[0]*fr[14]-1.125*nuVtSqSum[0]*fl[14])*rdv+alphaDrag[3]*(0.4330127018922193*favg[31]+0.25*favg[27])+alphaDrag[2]*(0.4330127018922193*favg[30]+0.25*favg[22])+alphaDrag[1]*(0.4330127018922193*favg[29]+0.25*favg[21])-0.8660254037844386*fjump[25]+alphaDrag[0]*(0.4330127018922193*favg[25]+0.25*favg[14])-0.5*fjump[14]; 
  Ghat[16] = ((-1.082531754730548*nuVtSqSum[0]*fr[26])-1.082531754730548*nuVtSqSum[0]*fl[26]-1.082531754730548*nuVtSqSum[1]*fr[19]-1.082531754730548*nuVtSqSum[1]*fl[19]-1.082531754730548*nuVtSqSum[2]*fr[18]-1.082531754730548*nuVtSqSum[2]*fl[18]+1.125*nuVtSqSum[0]*fr[16]-1.125*nuVtSqSum[0]*fl[16]-1.082531754730548*nuVtSqSum[3]*fr[11]-1.082531754730548*nuVtSqSum[3]*fl[11]+1.125*nuVtSqSum[1]*fr[8]-1.125*nuVtSqSum[1]*fl[8]+1.125*nuVtSqSum[2]*fr[7]-1.125*nuVtSqSum[2]*fl[7]+1.125*fr[3]*nuVtSqSum[3]-1.125*fl[3]*nuVtSqSum[3])*rdv-0.8660254037844386*fjump[26]+alphaDrag[0]*(0.4330127018922193*favg[26]+0.25*favg[16])+alphaDrag[1]*(0.4330127018922193*favg[19]+0.25*favg[8])+alphaDrag[2]*(0.4330127018922193*favg[18]+0.25*favg[7])-0.5*fjump[16]+alphaDrag[3]*(0.4330127018922193*favg[11]+0.25*favg[3]); 
  Ghat[20] = ((-1.082531754730548*nuVtSqSum[0]*fr[28])-1.082531754730548*nuVtSqSum[0]*fl[28]-1.082531754730548*nuVtSqSum[1]*fr[24]-1.082531754730548*nuVtSqSum[1]*fl[24]-1.082531754730548*nuVtSqSum[2]*fr[23]-1.082531754730548*nuVtSqSum[2]*fl[23]+1.125*nuVtSqSum[0]*fr[20]-1.125*nuVtSqSum[0]*fl[20]-1.082531754730548*nuVtSqSum[3]*fr[15]-1.082531754730548*nuVtSqSum[3]*fl[15]+1.125*nuVtSqSum[1]*fr[13]-1.125*nuVtSqSum[1]*fl[13]+1.125*nuVtSqSum[2]*fr[12]-1.125*nuVtSqSum[2]*fl[12]+1.125*nuVtSqSum[3]*fr[5]-1.125*nuVtSqSum[3]*fl[5])*rdv-0.8660254037844386*fjump[28]+alphaDrag[0]*(0.4330127018922193*favg[28]+0.25*favg[20])+alphaDrag[1]*(0.4330127018922193*favg[24]+0.25*favg[13])+alphaDrag[2]*(0.4330127018922193*favg[23]+0.25*favg[12])-0.5*fjump[20]+alphaDrag[3]*(0.4330127018922193*favg[15]+0.25*favg[5]); 
  Ghat[21] = ((-1.082531754730548*nuVtSqSum[2]*fr[31])-1.082531754730548*nuVtSqSum[2]*fl[31]-1.082531754730548*nuVtSqSum[3]*fr[30]-1.082531754730548*nuVtSqSum[3]*fl[30]-1.082531754730548*nuVtSqSum[0]*fr[29]-1.082531754730548*nuVtSqSum[0]*fl[29]+1.125*nuVtSqSum[2]*fr[27]-1.125*nuVtSqSum[2]*fl[27]-1.082531754730548*nuVtSqSum[1]*fr[25]-1.082531754730548*nuVtSqSum[1]*fl[25]+1.125*nuVtSqSum[3]*fr[22]-1.125*nuVtSqSum[3]*fl[22]+1.125*nuVtSqSum[0]*fr[21]-1.125*nuVtSqSum[0]*fl[21]+1.125*nuVtSqSum[1]*fr[14]-1.125*nuVtSqSum[1]*fl[14])*rdv+alphaDrag[2]*(0.4330127018922193*favg[31]+0.25*favg[27])+alphaDrag[3]*(0.4330127018922193*favg[30]+0.25*favg[22])-0.8660254037844386*fjump[29]+alphaDrag[0]*(0.4330127018922193*favg[29]+0.25*favg[21])+alphaDrag[1]*(0.4330127018922193*favg[25]+0.25*favg[14])-0.5*fjump[21]; 
  Ghat[22] = ((-1.082531754730548*nuVtSqSum[1]*fr[31])-1.082531754730548*nuVtSqSum[1]*fl[31]-1.082531754730548*nuVtSqSum[0]*fr[30]-1.082531754730548*nuVtSqSum[0]*fl[30]-1.082531754730548*nuVtSqSum[3]*fr[29]-1.082531754730548*nuVtSqSum[3]*fl[29]+1.125*nuVtSqSum[1]*fr[27]-1.125*nuVtSqSum[1]*fl[27]-1.082531754730548*nuVtSqSum[2]*fr[25]-1.082531754730548*nuVtSqSum[2]*fl[25]+1.125*nuVtSqSum[0]*fr[22]-1.125*nuVtSqSum[0]*fl[22]+1.125*nuVtSqSum[3]*fr[21]-1.125*nuVtSqSum[3]*fl[21]+1.125*nuVtSqSum[2]*fr[14]-1.125*nuVtSqSum[2]*fl[14])*rdv+alphaDrag[1]*(0.4330127018922193*favg[31]+0.25*favg[27])-0.8660254037844386*fjump[30]+alphaDrag[0]*(0.4330127018922193*favg[30]+0.25*favg[22])+alphaDrag[3]*(0.4330127018922193*favg[29]+0.25*favg[21])+alphaDrag[2]*(0.4330127018922193*favg[25]+0.25*favg[14])-0.5*fjump[22]; 
  Ghat[27] = ((-1.082531754730548*nuVtSqSum[0]*fr[31])-1.082531754730548*nuVtSqSum[0]*fl[31]-1.082531754730548*nuVtSqSum[1]*fr[30]-1.082531754730548*nuVtSqSum[1]*fl[30]-1.082531754730548*nuVtSqSum[2]*fr[29]-1.082531754730548*nuVtSqSum[2]*fl[29]+1.125*nuVtSqSum[0]*fr[27]-1.125*nuVtSqSum[0]*fl[27]-1.082531754730548*nuVtSqSum[3]*fr[25]-1.082531754730548*nuVtSqSum[3]*fl[25]+1.125*nuVtSqSum[1]*fr[22]-1.125*nuVtSqSum[1]*fl[22]+1.125*nuVtSqSum[2]*fr[21]-1.125*nuVtSqSum[2]*fl[21]+1.125*nuVtSqSum[3]*fr[14]-1.125*nuVtSqSum[3]*fl[14])*rdv-0.8660254037844386*fjump[31]+alphaDrag[0]*(0.4330127018922193*favg[31]+0.25*favg[27])+alphaDrag[1]*(0.4330127018922193*favg[30]+0.25*favg[22])+alphaDrag[2]*(0.4330127018922193*favg[29]+0.25*favg[21])-0.5*fjump[27]+alphaDrag[3]*(0.4330127018922193*favg[25]+0.25*favg[14]); 

  double incr1[32]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[0]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = -0.5*Ghat[6]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = 0.8660254037844386*Ghat[1]; 
  incr1[10] = 0.8660254037844386*Ghat[2]; 
  incr1[11] = 0.8660254037844386*Ghat[3]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = -0.5*Ghat[13]; 
  incr1[14] = -0.5*Ghat[14]; 
  incr1[15] = 0.8660254037844386*Ghat[5]; 
  incr1[16] = -0.5*Ghat[16]; 
  incr1[17] = 0.8660254037844386*Ghat[6]; 
  incr1[18] = 0.8660254037844386*Ghat[7]; 
  incr1[19] = 0.8660254037844386*Ghat[8]; 
  incr1[20] = -0.5*Ghat[20]; 
  incr1[21] = -0.5*Ghat[21]; 
  incr1[22] = -0.5*Ghat[22]; 
  incr1[23] = 0.8660254037844386*Ghat[12]; 
  incr1[24] = 0.8660254037844386*Ghat[13]; 
  incr1[25] = 0.8660254037844386*Ghat[14]; 
  incr1[26] = 0.8660254037844386*Ghat[16]; 
  incr1[27] = -0.5*Ghat[27]; 
  incr1[28] = 0.8660254037844386*Ghat[20]; 
  incr1[29] = 0.8660254037844386*Ghat[21]; 
  incr1[30] = 0.8660254037844386*Ghat[22]; 
  incr1[31] = 0.8660254037844386*Ghat[27]; 

  double incr2[32]; 
  incr2[4] = nuVtSqSum[3]*((-0.25*fr[17])+0.25*fl[17]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[2]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.25*fr[9])+0.25*fl[9]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.25*fr[4])+0.25*fl[4]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[9] = nuVtSqSum[2]*((-0.25*fr[17])+0.25*fl[17]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[3]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.25*fr[9])+0.25*fl[9]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.25*fr[4])+0.25*fl[4]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[10] = nuVtSqSum[1]*((-0.25*fr[17])+0.25*fl[17]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[0]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[3]*((-0.25*fr[9])+0.25*fl[9]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.25*fr[4])+0.25*fl[4]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[11] = nuVtSqSum[3]*((-0.25*fr[26])+0.25*fl[26]+0.2165063509461096*(fr[16]+fl[16]))+nuVtSqSum[2]*((-0.25*fr[19])+0.25*fl[19]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[1]*((-0.25*fr[18])+0.25*fl[18]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[0]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[3]+fl[3])); 
  incr2[15] = nuVtSqSum[3]*((-0.25*fr[28])+0.25*fl[28]+0.2165063509461096*(fr[20]+fl[20]))+nuVtSqSum[2]*((-0.25*fr[24])+0.25*fl[24]+0.2165063509461096*(fr[13]+fl[13]))+nuVtSqSum[1]*((-0.25*fr[23])+0.25*fl[23]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[0]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[5]+fl[5])); 
  incr2[17] = nuVtSqSum[0]*((-0.25*fr[17])+0.25*fl[17]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[1]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[2]*((-0.25*fr[9])+0.25*fl[9]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[3]*((-0.25*fr[4])+0.25*fl[4]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[18] = nuVtSqSum[2]*((-0.25*fr[26])+0.25*fl[26]+0.2165063509461096*(fr[16]+fl[16]))+nuVtSqSum[3]*((-0.25*fr[19])+0.25*fl[19]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[0]*((-0.25*fr[18])+0.25*fl[18]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[1]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[3]+fl[3])); 
  incr2[19] = nuVtSqSum[1]*((-0.25*fr[26])+0.25*fl[26]+0.2165063509461096*(fr[16]+fl[16]))+nuVtSqSum[0]*((-0.25*fr[19])+0.25*fl[19]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[3]*((-0.25*fr[18])+0.25*fl[18]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[2]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[3]+fl[3])); 
  incr2[23] = nuVtSqSum[2]*((-0.25*fr[28])+0.25*fl[28]+0.2165063509461096*(fr[20]+fl[20]))+nuVtSqSum[3]*((-0.25*fr[24])+0.25*fl[24]+0.2165063509461096*(fr[13]+fl[13]))+nuVtSqSum[0]*((-0.25*fr[23])+0.25*fl[23]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[1]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[5]+fl[5])); 
  incr2[24] = nuVtSqSum[1]*((-0.25*fr[28])+0.25*fl[28]+0.2165063509461096*(fr[20]+fl[20]))+nuVtSqSum[0]*((-0.25*fr[24])+0.25*fl[24]+0.2165063509461096*(fr[13]+fl[13]))+nuVtSqSum[3]*((-0.25*fr[23])+0.25*fl[23]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[2]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[5]+fl[5])); 
  incr2[25] = nuVtSqSum[3]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[27]+fl[27]))+nuVtSqSum[2]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[22]+fl[22]))+nuVtSqSum[1]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[21]+fl[21]))+nuVtSqSum[0]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[14]+fl[14])); 
  incr2[26] = nuVtSqSum[0]*((-0.25*fr[26])+0.25*fl[26]+0.2165063509461096*(fr[16]+fl[16]))+nuVtSqSum[1]*((-0.25*fr[19])+0.25*fl[19]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[2]*((-0.25*fr[18])+0.25*fl[18]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[3]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[3]+fl[3])); 
  incr2[28] = nuVtSqSum[0]*((-0.25*fr[28])+0.25*fl[28]+0.2165063509461096*(fr[20]+fl[20]))+nuVtSqSum[1]*((-0.25*fr[24])+0.25*fl[24]+0.2165063509461096*(fr[13]+fl[13]))+nuVtSqSum[2]*((-0.25*fr[23])+0.25*fl[23]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[3]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[5]+fl[5])); 
  incr2[29] = nuVtSqSum[2]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[27]+fl[27]))+nuVtSqSum[3]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[22]+fl[22]))+nuVtSqSum[0]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[21]+fl[21]))+nuVtSqSum[1]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[14]+fl[14])); 
  incr2[30] = nuVtSqSum[1]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[27]+fl[27]))+nuVtSqSum[0]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[22]+fl[22]))+nuVtSqSum[3]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[21]+fl[21]))+nuVtSqSum[2]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[14]+fl[14])); 
  incr2[31] = nuVtSqSum[0]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[27]+fl[27]))+nuVtSqSum[1]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[22]+fl[22]))+nuVtSqSum[2]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[21]+fl[21]))+nuVtSqSum[3]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[14]+fl[14])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr2[9]*rdvSq4R+incr1[9]*rdv2R; 
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr2[11]*rdvSq4R+incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr1[13]*rdv2R; 
  outr[14] += incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 
  outr[16] += incr1[16]*rdv2R; 
  outr[17] += incr2[17]*rdvSq4R+incr1[17]*rdv2R; 
  outr[18] += incr2[18]*rdvSq4R+incr1[18]*rdv2R; 
  outr[19] += incr2[19]*rdvSq4R+incr1[19]*rdv2R; 
  outr[20] += incr1[20]*rdv2R; 
  outr[21] += incr1[21]*rdv2R; 
  outr[22] += incr1[22]*rdv2R; 
  outr[23] += incr2[23]*rdvSq4R+incr1[23]*rdv2R; 
  outr[24] += incr2[24]*rdvSq4R+incr1[24]*rdv2R; 
  outr[25] += incr2[25]*rdvSq4R+incr1[25]*rdv2R; 
  outr[26] += incr2[26]*rdvSq4R+incr1[26]*rdv2R; 
  outr[27] += incr1[27]*rdv2R; 
  outr[28] += incr2[28]*rdvSq4R+incr1[28]*rdv2R; 
  outr[29] += incr2[29]*rdvSq4R+incr1[29]*rdv2R; 
  outr[30] += incr2[30]*rdvSq4R+incr1[30]*rdv2R; 
  outr[31] += incr2[31]*rdvSq4R+incr1[31]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += incr1[9]*rdv2L-1.0*incr2[9]*rdvSq4L; 
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += incr1[11]*rdv2L-1.0*incr2[11]*rdvSq4L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += -1.0*incr1[13]*rdv2L; 
  outl[14] += -1.0*incr1[14]*rdv2L; 
  outl[15] += incr1[15]*rdv2L-1.0*incr2[15]*rdvSq4L; 
  outl[16] += -1.0*incr1[16]*rdv2L; 
  outl[17] += incr1[17]*rdv2L-1.0*incr2[17]*rdvSq4L; 
  outl[18] += incr1[18]*rdv2L-1.0*incr2[18]*rdvSq4L; 
  outl[19] += incr1[19]*rdv2L-1.0*incr2[19]*rdvSq4L; 
  outl[20] += -1.0*incr1[20]*rdv2L; 
  outl[21] += -1.0*incr1[21]*rdv2L; 
  outl[22] += -1.0*incr1[22]*rdv2L; 
  outl[23] += incr1[23]*rdv2L-1.0*incr2[23]*rdvSq4L; 
  outl[24] += incr1[24]*rdv2L-1.0*incr2[24]*rdvSq4L; 
  outl[25] += incr1[25]*rdv2L-1.0*incr2[25]*rdvSq4L; 
  outl[26] += incr1[26]*rdv2L-1.0*incr2[26]*rdvSq4L; 
  outl[27] += -1.0*incr1[27]*rdv2L; 
  outl[28] += incr1[28]*rdv2L-1.0*incr2[28]*rdvSq4L; 
  outl[29] += incr1[29]*rdv2L-1.0*incr2[29]*rdvSq4L; 
  outl[30] += incr1[30]*rdv2L-1.0*incr2[30]*rdvSq4L; 
  outl[31] += incr1[31]*rdv2L-1.0*incr2[31]*rdvSq4L; 

  return std::abs(wl[3]-(0.5*sumNuUy[0])/nuSum); 
} 
double VmLBOconstNuSurf2x3vSer_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:          Cell-center coordinates. 
  // dxv[5]:        Cell spacing. 
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

  const double *sumNuUz = &nuUSum[8]; 

  double favg[32]; 
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
  favg[12] = -1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = -1*fr[20]+fl[20]; 
  favg[21] = -1*fr[21]+fl[21]; 
  favg[22] = -1*fr[22]+fl[22]; 
  favg[23] = -1*fr[23]+fl[23]; 
  favg[24] = -1*fr[24]+fl[24]; 
  favg[25] = -1*fr[25]+fl[25]; 
  favg[26] = 1*fr[26]+fl[26]; 
  favg[27] = -1*fr[27]+fl[27]; 
  favg[28] = -1*fr[28]+fl[28]; 
  favg[29] = -1*fr[29]+fl[29]; 
  favg[30] = -1*fr[30]+fl[30]; 
  favg[31] = -1*fr[31]+fl[31]; 

  double fjump[32]; 
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
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(-1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(-1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(-1*fr[15])); 
  fjump[16] = nuSum*vMuMidMax*(fl[16]-(1*fr[16])); 
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(1*fr[19])); 
  fjump[20] = nuSum*vMuMidMax*(fl[20]-(-1*fr[20])); 
  fjump[21] = nuSum*vMuMidMax*(fl[21]-(-1*fr[21])); 
  fjump[22] = nuSum*vMuMidMax*(fl[22]-(-1*fr[22])); 
  fjump[23] = nuSum*vMuMidMax*(fl[23]-(-1*fr[23])); 
  fjump[24] = nuSum*vMuMidMax*(fl[24]-(-1*fr[24])); 
  fjump[25] = nuSum*vMuMidMax*(fl[25]-(-1*fr[25])); 
  fjump[26] = nuSum*vMuMidMax*(fl[26]-(1*fr[26])); 
  fjump[27] = nuSum*vMuMidMax*(fl[27]-(-1*fr[27])); 
  fjump[28] = nuSum*vMuMidMax*(fl[28]-(-1*fr[28])); 
  fjump[29] = nuSum*vMuMidMax*(fl[29]-(-1*fr[29])); 
  fjump[30] = nuSum*vMuMidMax*(fl[30]-(-1*fr[30])); 
  fjump[31] = nuSum*vMuMidMax*(fl[31]-(-1*fr[31])); 

  double alphaDrag[4]; 
  alphaDrag[0] = 2.0*wl[4]*nuSum+dxvl[4]*nuSum-1.0*sumNuUz[0]; 
  alphaDrag[1] = -1.0*sumNuUz[1]; 
  alphaDrag[2] = -1.0*sumNuUz[2]; 
  alphaDrag[3] = -1.0*sumNuUz[3]; 

  double Ghat[32]; 
  for(unsigned short int i=0; i<32; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.082531754730548*nuVtSqSum[3]*fr[20])-1.082531754730548*nuVtSqSum[3]*fl[20]-1.082531754730548*nuVtSqSum[2]*fr[13]-1.082531754730548*nuVtSqSum[2]*fl[13]-1.082531754730548*nuVtSqSum[1]*fr[12]-1.082531754730548*nuVtSqSum[1]*fl[12]+1.125*nuVtSqSum[3]*fr[6]-1.125*nuVtSqSum[3]*fl[6]-1.082531754730548*nuVtSqSum[0]*fr[5]-1.082531754730548*nuVtSqSum[0]*fl[5]+1.125*fr[2]*nuVtSqSum[2]-1.125*fl[2]*nuVtSqSum[2]+1.125*fr[1]*nuVtSqSum[1]-1.125*fl[1]*nuVtSqSum[1]+1.125*fr[0]*nuVtSqSum[0]-1.125*fl[0]*nuVtSqSum[0])*rdv+alphaDrag[3]*(0.4330127018922193*favg[20]+0.25*favg[6])+alphaDrag[2]*(0.4330127018922193*favg[13]+0.25*favg[2])+alphaDrag[1]*(0.4330127018922193*favg[12]+0.25*favg[1])-0.8660254037844386*fjump[5]+alphaDrag[0]*(0.4330127018922193*favg[5]+0.25*favg[0])-0.5*fjump[0]; 
  Ghat[1] = ((-1.082531754730548*nuVtSqSum[2]*fr[20])-1.082531754730548*nuVtSqSum[2]*fl[20]-1.082531754730548*nuVtSqSum[3]*fr[13]-1.082531754730548*nuVtSqSum[3]*fl[13]-1.082531754730548*nuVtSqSum[0]*fr[12]-1.082531754730548*nuVtSqSum[0]*fl[12]+1.125*nuVtSqSum[2]*fr[6]-1.125*nuVtSqSum[2]*fl[6]-1.082531754730548*nuVtSqSum[1]*fr[5]-1.082531754730548*nuVtSqSum[1]*fl[5]+1.125*fr[2]*nuVtSqSum[3]-1.125*fl[2]*nuVtSqSum[3]+1.125*fr[0]*nuVtSqSum[1]-1.125*fl[0]*nuVtSqSum[1]+1.125*nuVtSqSum[0]*fr[1]-1.125*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[2]*(0.4330127018922193*favg[20]+0.25*favg[6])+alphaDrag[3]*(0.4330127018922193*favg[13]+0.25*favg[2])-0.8660254037844386*fjump[12]+alphaDrag[0]*(0.4330127018922193*favg[12]+0.25*favg[1])+alphaDrag[1]*(0.4330127018922193*favg[5]+0.25*favg[0])-0.5*fjump[1]; 
  Ghat[2] = ((-1.082531754730548*nuVtSqSum[1]*fr[20])-1.082531754730548*nuVtSqSum[1]*fl[20]-1.082531754730548*nuVtSqSum[0]*fr[13]-1.082531754730548*nuVtSqSum[0]*fl[13]-1.082531754730548*nuVtSqSum[3]*fr[12]-1.082531754730548*nuVtSqSum[3]*fl[12]+1.125*nuVtSqSum[1]*fr[6]-1.125*nuVtSqSum[1]*fl[6]-1.082531754730548*nuVtSqSum[2]*fr[5]-1.082531754730548*nuVtSqSum[2]*fl[5]+1.125*fr[1]*nuVtSqSum[3]-1.125*fl[1]*nuVtSqSum[3]+1.125*fr[0]*nuVtSqSum[2]-1.125*fl[0]*nuVtSqSum[2]+1.125*nuVtSqSum[0]*fr[2]-1.125*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[1]*(0.4330127018922193*favg[20]+0.25*favg[6])-0.8660254037844386*fjump[13]+alphaDrag[0]*(0.4330127018922193*favg[13]+0.25*favg[2])+alphaDrag[3]*(0.4330127018922193*favg[12]+0.25*favg[1])+alphaDrag[2]*(0.4330127018922193*favg[5]+0.25*favg[0])-0.5*fjump[2]; 
  Ghat[3] = ((-1.082531754730548*nuVtSqSum[3]*fr[27])-1.082531754730548*nuVtSqSum[3]*fl[27]-1.082531754730548*nuVtSqSum[2]*fr[22]-1.082531754730548*nuVtSqSum[2]*fl[22]-1.082531754730548*nuVtSqSum[1]*fr[21]-1.082531754730548*nuVtSqSum[1]*fl[21]+1.125*nuVtSqSum[3]*fr[16]-1.125*nuVtSqSum[3]*fl[16]-1.082531754730548*nuVtSqSum[0]*fr[14]-1.082531754730548*nuVtSqSum[0]*fl[14]+1.125*nuVtSqSum[2]*fr[8]-1.125*nuVtSqSum[2]*fl[8]+1.125*nuVtSqSum[1]*fr[7]-1.125*nuVtSqSum[1]*fl[7]+1.125*nuVtSqSum[0]*fr[3]-1.125*nuVtSqSum[0]*fl[3])*rdv+alphaDrag[3]*(0.4330127018922193*favg[27]+0.25*favg[16])+alphaDrag[2]*(0.4330127018922193*favg[22]+0.25*favg[8])+alphaDrag[1]*(0.4330127018922193*favg[21]+0.25*favg[7])-0.8660254037844386*fjump[14]+alphaDrag[0]*(0.4330127018922193*favg[14]+0.25*favg[3])-0.5*fjump[3]; 
  Ghat[4] = ((-1.082531754730548*nuVtSqSum[3]*fr[28])-1.082531754730548*nuVtSqSum[3]*fl[28]-1.082531754730548*nuVtSqSum[2]*fr[24]-1.082531754730548*nuVtSqSum[2]*fl[24]-1.082531754730548*nuVtSqSum[1]*fr[23]-1.082531754730548*nuVtSqSum[1]*fl[23]+1.125*nuVtSqSum[3]*fr[17]-1.125*nuVtSqSum[3]*fl[17]-1.082531754730548*nuVtSqSum[0]*fr[15]-1.082531754730548*nuVtSqSum[0]*fl[15]+1.125*nuVtSqSum[2]*fr[10]-1.125*nuVtSqSum[2]*fl[10]+1.125*nuVtSqSum[1]*fr[9]-1.125*nuVtSqSum[1]*fl[9]+1.125*nuVtSqSum[0]*fr[4]-1.125*nuVtSqSum[0]*fl[4])*rdv+alphaDrag[3]*(0.4330127018922193*favg[28]+0.25*favg[17])+alphaDrag[2]*(0.4330127018922193*favg[24]+0.25*favg[10])+alphaDrag[1]*(0.4330127018922193*favg[23]+0.25*favg[9])-0.8660254037844386*fjump[15]+alphaDrag[0]*(0.4330127018922193*favg[15]+0.25*favg[4])-0.5*fjump[4]; 
  Ghat[6] = ((-1.082531754730548*nuVtSqSum[0]*fr[20])-1.082531754730548*nuVtSqSum[0]*fl[20]-1.082531754730548*nuVtSqSum[1]*fr[13]-1.082531754730548*nuVtSqSum[1]*fl[13]-1.082531754730548*nuVtSqSum[2]*fr[12]-1.082531754730548*nuVtSqSum[2]*fl[12]+1.125*nuVtSqSum[0]*fr[6]-1.125*nuVtSqSum[0]*fl[6]-1.082531754730548*nuVtSqSum[3]*fr[5]-1.082531754730548*nuVtSqSum[3]*fl[5]+1.125*fr[0]*nuVtSqSum[3]-1.125*fl[0]*nuVtSqSum[3]+1.125*fr[1]*nuVtSqSum[2]-1.125*fl[1]*nuVtSqSum[2]+1.125*nuVtSqSum[1]*fr[2]-1.125*nuVtSqSum[1]*fl[2])*rdv-0.8660254037844386*fjump[20]+alphaDrag[0]*(0.4330127018922193*favg[20]+0.25*favg[6])+alphaDrag[1]*(0.4330127018922193*favg[13]+0.25*favg[2])+alphaDrag[2]*(0.4330127018922193*favg[12]+0.25*favg[1])-0.5*fjump[6]+alphaDrag[3]*(0.4330127018922193*favg[5]+0.25*favg[0]); 
  Ghat[7] = ((-1.082531754730548*nuVtSqSum[2]*fr[27])-1.082531754730548*nuVtSqSum[2]*fl[27]-1.082531754730548*nuVtSqSum[3]*fr[22]-1.082531754730548*nuVtSqSum[3]*fl[22]-1.082531754730548*nuVtSqSum[0]*fr[21]-1.082531754730548*nuVtSqSum[0]*fl[21]+1.125*nuVtSqSum[2]*fr[16]-1.125*nuVtSqSum[2]*fl[16]-1.082531754730548*nuVtSqSum[1]*fr[14]-1.082531754730548*nuVtSqSum[1]*fl[14]+1.125*nuVtSqSum[3]*fr[8]-1.125*nuVtSqSum[3]*fl[8]+1.125*nuVtSqSum[0]*fr[7]-1.125*nuVtSqSum[0]*fl[7]+1.125*nuVtSqSum[1]*fr[3]-1.125*nuVtSqSum[1]*fl[3])*rdv+alphaDrag[2]*(0.4330127018922193*favg[27]+0.25*favg[16])+alphaDrag[3]*(0.4330127018922193*favg[22]+0.25*favg[8])-0.8660254037844386*fjump[21]+alphaDrag[0]*(0.4330127018922193*favg[21]+0.25*favg[7])+alphaDrag[1]*(0.4330127018922193*favg[14]+0.25*favg[3])-0.5*fjump[7]; 
  Ghat[8] = ((-1.082531754730548*nuVtSqSum[1]*fr[27])-1.082531754730548*nuVtSqSum[1]*fl[27]-1.082531754730548*nuVtSqSum[0]*fr[22]-1.082531754730548*nuVtSqSum[0]*fl[22]-1.082531754730548*nuVtSqSum[3]*fr[21]-1.082531754730548*nuVtSqSum[3]*fl[21]+1.125*nuVtSqSum[1]*fr[16]-1.125*nuVtSqSum[1]*fl[16]-1.082531754730548*nuVtSqSum[2]*fr[14]-1.082531754730548*nuVtSqSum[2]*fl[14]+1.125*nuVtSqSum[0]*fr[8]-1.125*nuVtSqSum[0]*fl[8]+1.125*nuVtSqSum[3]*fr[7]-1.125*nuVtSqSum[3]*fl[7]+1.125*nuVtSqSum[2]*fr[3]-1.125*nuVtSqSum[2]*fl[3])*rdv+alphaDrag[1]*(0.4330127018922193*favg[27]+0.25*favg[16])-0.8660254037844386*fjump[22]+alphaDrag[0]*(0.4330127018922193*favg[22]+0.25*favg[8])+alphaDrag[3]*(0.4330127018922193*favg[21]+0.25*favg[7])+alphaDrag[2]*(0.4330127018922193*favg[14]+0.25*favg[3])-0.5*fjump[8]; 
  Ghat[9] = ((-1.082531754730548*nuVtSqSum[2]*fr[28])-1.082531754730548*nuVtSqSum[2]*fl[28]-1.082531754730548*nuVtSqSum[3]*fr[24]-1.082531754730548*nuVtSqSum[3]*fl[24]-1.082531754730548*nuVtSqSum[0]*fr[23]-1.082531754730548*nuVtSqSum[0]*fl[23]+1.125*nuVtSqSum[2]*fr[17]-1.125*nuVtSqSum[2]*fl[17]-1.082531754730548*nuVtSqSum[1]*fr[15]-1.082531754730548*nuVtSqSum[1]*fl[15]+1.125*nuVtSqSum[3]*fr[10]-1.125*nuVtSqSum[3]*fl[10]+1.125*nuVtSqSum[0]*fr[9]-1.125*nuVtSqSum[0]*fl[9]+1.125*nuVtSqSum[1]*fr[4]-1.125*nuVtSqSum[1]*fl[4])*rdv+alphaDrag[2]*(0.4330127018922193*favg[28]+0.25*favg[17])+alphaDrag[3]*(0.4330127018922193*favg[24]+0.25*favg[10])-0.8660254037844386*fjump[23]+alphaDrag[0]*(0.4330127018922193*favg[23]+0.25*favg[9])+alphaDrag[1]*(0.4330127018922193*favg[15]+0.25*favg[4])-0.5*fjump[9]; 
  Ghat[10] = ((-1.082531754730548*nuVtSqSum[1]*fr[28])-1.082531754730548*nuVtSqSum[1]*fl[28]-1.082531754730548*nuVtSqSum[0]*fr[24]-1.082531754730548*nuVtSqSum[0]*fl[24]-1.082531754730548*nuVtSqSum[3]*fr[23]-1.082531754730548*nuVtSqSum[3]*fl[23]+1.125*nuVtSqSum[1]*fr[17]-1.125*nuVtSqSum[1]*fl[17]-1.082531754730548*nuVtSqSum[2]*fr[15]-1.082531754730548*nuVtSqSum[2]*fl[15]+1.125*nuVtSqSum[0]*fr[10]-1.125*nuVtSqSum[0]*fl[10]+1.125*nuVtSqSum[3]*fr[9]-1.125*nuVtSqSum[3]*fl[9]+1.125*nuVtSqSum[2]*fr[4]-1.125*nuVtSqSum[2]*fl[4])*rdv+alphaDrag[1]*(0.4330127018922193*favg[28]+0.25*favg[17])-0.8660254037844386*fjump[24]+alphaDrag[0]*(0.4330127018922193*favg[24]+0.25*favg[10])+alphaDrag[3]*(0.4330127018922193*favg[23]+0.25*favg[9])+alphaDrag[2]*(0.4330127018922193*favg[15]+0.25*favg[4])-0.5*fjump[10]; 
  Ghat[11] = ((-1.082531754730548*nuVtSqSum[3]*fr[31])-1.082531754730548*nuVtSqSum[3]*fl[31]-1.082531754730548*nuVtSqSum[2]*fr[30]-1.082531754730548*nuVtSqSum[2]*fl[30]-1.082531754730548*nuVtSqSum[1]*fr[29]-1.082531754730548*nuVtSqSum[1]*fl[29]+1.125*nuVtSqSum[3]*fr[26]-1.125*nuVtSqSum[3]*fl[26]-1.082531754730548*nuVtSqSum[0]*fr[25]-1.082531754730548*nuVtSqSum[0]*fl[25]+1.125*nuVtSqSum[2]*fr[19]-1.125*nuVtSqSum[2]*fl[19]+1.125*nuVtSqSum[1]*fr[18]-1.125*nuVtSqSum[1]*fl[18]+1.125*nuVtSqSum[0]*fr[11]-1.125*nuVtSqSum[0]*fl[11])*rdv+alphaDrag[3]*(0.4330127018922193*favg[31]+0.25*favg[26])+alphaDrag[2]*(0.4330127018922193*favg[30]+0.25*favg[19])+alphaDrag[1]*(0.4330127018922193*favg[29]+0.25*favg[18])-0.8660254037844386*fjump[25]+alphaDrag[0]*(0.4330127018922193*favg[25]+0.25*favg[11])-0.5*fjump[11]; 
  Ghat[16] = ((-1.082531754730548*nuVtSqSum[0]*fr[27])-1.082531754730548*nuVtSqSum[0]*fl[27]-1.082531754730548*nuVtSqSum[1]*fr[22]-1.082531754730548*nuVtSqSum[1]*fl[22]-1.082531754730548*nuVtSqSum[2]*fr[21]-1.082531754730548*nuVtSqSum[2]*fl[21]+1.125*nuVtSqSum[0]*fr[16]-1.125*nuVtSqSum[0]*fl[16]-1.082531754730548*nuVtSqSum[3]*fr[14]-1.082531754730548*nuVtSqSum[3]*fl[14]+1.125*nuVtSqSum[1]*fr[8]-1.125*nuVtSqSum[1]*fl[8]+1.125*nuVtSqSum[2]*fr[7]-1.125*nuVtSqSum[2]*fl[7]+1.125*fr[3]*nuVtSqSum[3]-1.125*fl[3]*nuVtSqSum[3])*rdv-0.8660254037844386*fjump[27]+alphaDrag[0]*(0.4330127018922193*favg[27]+0.25*favg[16])+alphaDrag[1]*(0.4330127018922193*favg[22]+0.25*favg[8])+alphaDrag[2]*(0.4330127018922193*favg[21]+0.25*favg[7])-0.5*fjump[16]+alphaDrag[3]*(0.4330127018922193*favg[14]+0.25*favg[3]); 
  Ghat[17] = ((-1.082531754730548*nuVtSqSum[0]*fr[28])-1.082531754730548*nuVtSqSum[0]*fl[28]-1.082531754730548*nuVtSqSum[1]*fr[24]-1.082531754730548*nuVtSqSum[1]*fl[24]-1.082531754730548*nuVtSqSum[2]*fr[23]-1.082531754730548*nuVtSqSum[2]*fl[23]+1.125*nuVtSqSum[0]*fr[17]-1.125*nuVtSqSum[0]*fl[17]-1.082531754730548*nuVtSqSum[3]*fr[15]-1.082531754730548*nuVtSqSum[3]*fl[15]+1.125*nuVtSqSum[1]*fr[10]-1.125*nuVtSqSum[1]*fl[10]+1.125*nuVtSqSum[2]*fr[9]-1.125*nuVtSqSum[2]*fl[9]+1.125*nuVtSqSum[3]*fr[4]-1.125*nuVtSqSum[3]*fl[4])*rdv-0.8660254037844386*fjump[28]+alphaDrag[0]*(0.4330127018922193*favg[28]+0.25*favg[17])+alphaDrag[1]*(0.4330127018922193*favg[24]+0.25*favg[10])+alphaDrag[2]*(0.4330127018922193*favg[23]+0.25*favg[9])-0.5*fjump[17]+alphaDrag[3]*(0.4330127018922193*favg[15]+0.25*favg[4]); 
  Ghat[18] = ((-1.082531754730548*nuVtSqSum[2]*fr[31])-1.082531754730548*nuVtSqSum[2]*fl[31]-1.082531754730548*nuVtSqSum[3]*fr[30]-1.082531754730548*nuVtSqSum[3]*fl[30]-1.082531754730548*nuVtSqSum[0]*fr[29]-1.082531754730548*nuVtSqSum[0]*fl[29]+1.125*nuVtSqSum[2]*fr[26]-1.125*nuVtSqSum[2]*fl[26]-1.082531754730548*nuVtSqSum[1]*fr[25]-1.082531754730548*nuVtSqSum[1]*fl[25]+1.125*nuVtSqSum[3]*fr[19]-1.125*nuVtSqSum[3]*fl[19]+1.125*nuVtSqSum[0]*fr[18]-1.125*nuVtSqSum[0]*fl[18]+1.125*nuVtSqSum[1]*fr[11]-1.125*nuVtSqSum[1]*fl[11])*rdv+alphaDrag[2]*(0.4330127018922193*favg[31]+0.25*favg[26])+alphaDrag[3]*(0.4330127018922193*favg[30]+0.25*favg[19])-0.8660254037844386*fjump[29]+alphaDrag[0]*(0.4330127018922193*favg[29]+0.25*favg[18])+alphaDrag[1]*(0.4330127018922193*favg[25]+0.25*favg[11])-0.5*fjump[18]; 
  Ghat[19] = ((-1.082531754730548*nuVtSqSum[1]*fr[31])-1.082531754730548*nuVtSqSum[1]*fl[31]-1.082531754730548*nuVtSqSum[0]*fr[30]-1.082531754730548*nuVtSqSum[0]*fl[30]-1.082531754730548*nuVtSqSum[3]*fr[29]-1.082531754730548*nuVtSqSum[3]*fl[29]+1.125*nuVtSqSum[1]*fr[26]-1.125*nuVtSqSum[1]*fl[26]-1.082531754730548*nuVtSqSum[2]*fr[25]-1.082531754730548*nuVtSqSum[2]*fl[25]+1.125*nuVtSqSum[0]*fr[19]-1.125*nuVtSqSum[0]*fl[19]+1.125*nuVtSqSum[3]*fr[18]-1.125*nuVtSqSum[3]*fl[18]+1.125*nuVtSqSum[2]*fr[11]-1.125*nuVtSqSum[2]*fl[11])*rdv+alphaDrag[1]*(0.4330127018922193*favg[31]+0.25*favg[26])-0.8660254037844386*fjump[30]+alphaDrag[0]*(0.4330127018922193*favg[30]+0.25*favg[19])+alphaDrag[3]*(0.4330127018922193*favg[29]+0.25*favg[18])+alphaDrag[2]*(0.4330127018922193*favg[25]+0.25*favg[11])-0.5*fjump[19]; 
  Ghat[26] = ((-1.082531754730548*nuVtSqSum[0]*fr[31])-1.082531754730548*nuVtSqSum[0]*fl[31]-1.082531754730548*nuVtSqSum[1]*fr[30]-1.082531754730548*nuVtSqSum[1]*fl[30]-1.082531754730548*nuVtSqSum[2]*fr[29]-1.082531754730548*nuVtSqSum[2]*fl[29]+1.125*nuVtSqSum[0]*fr[26]-1.125*nuVtSqSum[0]*fl[26]-1.082531754730548*nuVtSqSum[3]*fr[25]-1.082531754730548*nuVtSqSum[3]*fl[25]+1.125*nuVtSqSum[1]*fr[19]-1.125*nuVtSqSum[1]*fl[19]+1.125*nuVtSqSum[2]*fr[18]-1.125*nuVtSqSum[2]*fl[18]+1.125*nuVtSqSum[3]*fr[11]-1.125*nuVtSqSum[3]*fl[11])*rdv-0.8660254037844386*fjump[31]+alphaDrag[0]*(0.4330127018922193*favg[31]+0.25*favg[26])+alphaDrag[1]*(0.4330127018922193*favg[30]+0.25*favg[19])+alphaDrag[2]*(0.4330127018922193*favg[29]+0.25*favg[18])-0.5*fjump[26]+alphaDrag[3]*(0.4330127018922193*favg[25]+0.25*favg[11]); 

  double incr1[32]; 
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
  incr1[12] = 0.8660254037844386*Ghat[1]; 
  incr1[13] = 0.8660254037844386*Ghat[2]; 
  incr1[14] = 0.8660254037844386*Ghat[3]; 
  incr1[15] = 0.8660254037844386*Ghat[4]; 
  incr1[16] = -0.5*Ghat[16]; 
  incr1[17] = -0.5*Ghat[17]; 
  incr1[18] = -0.5*Ghat[18]; 
  incr1[19] = -0.5*Ghat[19]; 
  incr1[20] = 0.8660254037844386*Ghat[6]; 
  incr1[21] = 0.8660254037844386*Ghat[7]; 
  incr1[22] = 0.8660254037844386*Ghat[8]; 
  incr1[23] = 0.8660254037844386*Ghat[9]; 
  incr1[24] = 0.8660254037844386*Ghat[10]; 
  incr1[25] = 0.8660254037844386*Ghat[11]; 
  incr1[26] = -0.5*Ghat[26]; 
  incr1[27] = 0.8660254037844386*Ghat[16]; 
  incr1[28] = 0.8660254037844386*Ghat[17]; 
  incr1[29] = 0.8660254037844386*Ghat[18]; 
  incr1[30] = 0.8660254037844386*Ghat[19]; 
  incr1[31] = 0.8660254037844386*Ghat[26]; 

  double incr2[32]; 
  incr2[5] = nuVtSqSum[3]*((-0.25*fr[20])+0.25*fl[20]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[2]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.25*fr[12])+0.25*fl[12]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.25*fr[5])+0.25*fl[5]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[12] = nuVtSqSum[2]*((-0.25*fr[20])+0.25*fl[20]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[3]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.25*fr[12])+0.25*fl[12]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.25*fr[5])+0.25*fl[5]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[13] = nuVtSqSum[1]*((-0.25*fr[20])+0.25*fl[20]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[0]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[3]*((-0.25*fr[12])+0.25*fl[12]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.25*fr[5])+0.25*fl[5]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[14] = nuVtSqSum[3]*((-0.25*fr[27])+0.25*fl[27]+0.2165063509461096*(fr[16]+fl[16]))+nuVtSqSum[2]*((-0.25*fr[22])+0.25*fl[22]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[1]*((-0.25*fr[21])+0.25*fl[21]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[0]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[3]+fl[3])); 
  incr2[15] = nuVtSqSum[3]*((-0.25*fr[28])+0.25*fl[28]+0.2165063509461096*(fr[17]+fl[17]))+nuVtSqSum[2]*((-0.25*fr[24])+0.25*fl[24]+0.2165063509461096*(fr[10]+fl[10]))+nuVtSqSum[1]*((-0.25*fr[23])+0.25*fl[23]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[0]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[20] = nuVtSqSum[0]*((-0.25*fr[20])+0.25*fl[20]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[1]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[2]*((-0.25*fr[12])+0.25*fl[12]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[3]*((-0.25*fr[5])+0.25*fl[5]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[21] = nuVtSqSum[2]*((-0.25*fr[27])+0.25*fl[27]+0.2165063509461096*(fr[16]+fl[16]))+nuVtSqSum[3]*((-0.25*fr[22])+0.25*fl[22]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[0]*((-0.25*fr[21])+0.25*fl[21]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[1]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[3]+fl[3])); 
  incr2[22] = nuVtSqSum[1]*((-0.25*fr[27])+0.25*fl[27]+0.2165063509461096*(fr[16]+fl[16]))+nuVtSqSum[0]*((-0.25*fr[22])+0.25*fl[22]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[3]*((-0.25*fr[21])+0.25*fl[21]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[2]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[3]+fl[3])); 
  incr2[23] = nuVtSqSum[2]*((-0.25*fr[28])+0.25*fl[28]+0.2165063509461096*(fr[17]+fl[17]))+nuVtSqSum[3]*((-0.25*fr[24])+0.25*fl[24]+0.2165063509461096*(fr[10]+fl[10]))+nuVtSqSum[0]*((-0.25*fr[23])+0.25*fl[23]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[1]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[24] = nuVtSqSum[1]*((-0.25*fr[28])+0.25*fl[28]+0.2165063509461096*(fr[17]+fl[17]))+nuVtSqSum[0]*((-0.25*fr[24])+0.25*fl[24]+0.2165063509461096*(fr[10]+fl[10]))+nuVtSqSum[3]*((-0.25*fr[23])+0.25*fl[23]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[2]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[25] = nuVtSqSum[3]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[26]+fl[26]))+nuVtSqSum[2]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[19]+fl[19]))+nuVtSqSum[1]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[18]+fl[18]))+nuVtSqSum[0]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[11]+fl[11])); 
  incr2[27] = nuVtSqSum[0]*((-0.25*fr[27])+0.25*fl[27]+0.2165063509461096*(fr[16]+fl[16]))+nuVtSqSum[1]*((-0.25*fr[22])+0.25*fl[22]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[2]*((-0.25*fr[21])+0.25*fl[21]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[3]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[3]+fl[3])); 
  incr2[28] = nuVtSqSum[0]*((-0.25*fr[28])+0.25*fl[28]+0.2165063509461096*(fr[17]+fl[17]))+nuVtSqSum[1]*((-0.25*fr[24])+0.25*fl[24]+0.2165063509461096*(fr[10]+fl[10]))+nuVtSqSum[2]*((-0.25*fr[23])+0.25*fl[23]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[3]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[29] = nuVtSqSum[2]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[26]+fl[26]))+nuVtSqSum[3]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[19]+fl[19]))+nuVtSqSum[0]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[18]+fl[18]))+nuVtSqSum[1]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[11]+fl[11])); 
  incr2[30] = nuVtSqSum[1]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[26]+fl[26]))+nuVtSqSum[0]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[19]+fl[19]))+nuVtSqSum[3]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[18]+fl[18]))+nuVtSqSum[2]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[11]+fl[11])); 
  incr2[31] = nuVtSqSum[0]*((-0.25*fr[31])+0.25*fl[31]+0.2165063509461096*(fr[26]+fl[26]))+nuVtSqSum[1]*((-0.25*fr[30])+0.25*fl[30]+0.2165063509461096*(fr[19]+fl[19]))+nuVtSqSum[2]*((-0.25*fr[29])+0.25*fl[29]+0.2165063509461096*(fr[18]+fl[18]))+nuVtSqSum[3]*((-0.25*fr[25])+0.25*fl[25]+0.2165063509461096*(fr[11]+fl[11])); 

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
  outr[12] += incr2[12]*rdvSq4R+incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 
  outr[16] += incr1[16]*rdv2R; 
  outr[17] += incr1[17]*rdv2R; 
  outr[18] += incr1[18]*rdv2R; 
  outr[19] += incr1[19]*rdv2R; 
  outr[20] += incr2[20]*rdvSq4R+incr1[20]*rdv2R; 
  outr[21] += incr2[21]*rdvSq4R+incr1[21]*rdv2R; 
  outr[22] += incr2[22]*rdvSq4R+incr1[22]*rdv2R; 
  outr[23] += incr2[23]*rdvSq4R+incr1[23]*rdv2R; 
  outr[24] += incr2[24]*rdvSq4R+incr1[24]*rdv2R; 
  outr[25] += incr2[25]*rdvSq4R+incr1[25]*rdv2R; 
  outr[26] += incr1[26]*rdv2R; 
  outr[27] += incr2[27]*rdvSq4R+incr1[27]*rdv2R; 
  outr[28] += incr2[28]*rdvSq4R+incr1[28]*rdv2R; 
  outr[29] += incr2[29]*rdvSq4R+incr1[29]*rdv2R; 
  outr[30] += incr2[30]*rdvSq4R+incr1[30]*rdv2R; 
  outr[31] += incr2[31]*rdvSq4R+incr1[31]*rdv2R; 

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
  outl[12] += incr1[12]*rdv2L-1.0*incr2[12]*rdvSq4L; 
  outl[13] += incr1[13]*rdv2L-1.0*incr2[13]*rdvSq4L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += incr1[15]*rdv2L-1.0*incr2[15]*rdvSq4L; 
  outl[16] += -1.0*incr1[16]*rdv2L; 
  outl[17] += -1.0*incr1[17]*rdv2L; 
  outl[18] += -1.0*incr1[18]*rdv2L; 
  outl[19] += -1.0*incr1[19]*rdv2L; 
  outl[20] += incr1[20]*rdv2L-1.0*incr2[20]*rdvSq4L; 
  outl[21] += incr1[21]*rdv2L-1.0*incr2[21]*rdvSq4L; 
  outl[22] += incr1[22]*rdv2L-1.0*incr2[22]*rdvSq4L; 
  outl[23] += incr1[23]*rdv2L-1.0*incr2[23]*rdvSq4L; 
  outl[24] += incr1[24]*rdv2L-1.0*incr2[24]*rdvSq4L; 
  outl[25] += incr1[25]*rdv2L-1.0*incr2[25]*rdvSq4L; 
  outl[26] += -1.0*incr1[26]*rdv2L; 
  outl[27] += incr1[27]*rdv2L-1.0*incr2[27]*rdvSq4L; 
  outl[28] += incr1[28]*rdv2L-1.0*incr2[28]*rdvSq4L; 
  outl[29] += incr1[29]*rdv2L-1.0*incr2[29]*rdvSq4L; 
  outl[30] += incr1[30]*rdv2L-1.0*incr2[30]*rdvSq4L; 
  outl[31] += incr1[31]*rdv2L-1.0*incr2[31]*rdvSq4L; 

  return std::abs(wl[4]-(0.5*sumNuUz[0])/nuSum); 
} 
