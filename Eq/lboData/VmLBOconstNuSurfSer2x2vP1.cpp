#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf2x2vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[8]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[16]; 
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
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 

  double fjump[16]; 
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
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(-1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(-1*fr[15])); 

  double alphaDrag[4]; 
  alphaDrag[0] = 2.0*wl[2]*nuSum+dxvl[2]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 
  alphaDrag[3] = -1.0*sumNuUx[3]; 

  double Ghat[16]; 
  for(unsigned short int i=0; i<16; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.082531754730548*nuVtSqSum[3]*fr[11])-1.082531754730548*nuVtSqSum[3]*fl[11]-1.082531754730548*nuVtSqSum[2]*fr[7]-1.082531754730548*nuVtSqSum[2]*fl[7]-1.082531754730548*nuVtSqSum[1]*fr[6]-1.082531754730548*nuVtSqSum[1]*fl[6]+1.125*nuVtSqSum[3]*fr[5]-1.125*nuVtSqSum[3]*fl[5]-1.082531754730548*nuVtSqSum[0]*fr[3]-1.082531754730548*nuVtSqSum[0]*fl[3]+1.125*fr[2]*nuVtSqSum[2]-1.125*fl[2]*nuVtSqSum[2]+1.125*fr[1]*nuVtSqSum[1]-1.125*fl[1]*nuVtSqSum[1]+1.125*fr[0]*nuVtSqSum[0]-1.125*fl[0]*nuVtSqSum[0])*rdv+alphaDrag[3]*(0.4330127018922193*favg[11]+0.25*favg[5])+alphaDrag[2]*(0.4330127018922193*favg[7]+0.25*favg[2])+alphaDrag[1]*(0.4330127018922193*favg[6]+0.25*favg[1])-0.8660254037844386*fjump[3]+alphaDrag[0]*(0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[0]; 
  Ghat[1] = ((-1.082531754730548*nuVtSqSum[2]*fr[11])-1.082531754730548*nuVtSqSum[2]*fl[11]-1.082531754730548*nuVtSqSum[3]*fr[7]-1.082531754730548*nuVtSqSum[3]*fl[7]-1.082531754730548*nuVtSqSum[0]*fr[6]-1.082531754730548*nuVtSqSum[0]*fl[6]+1.125*nuVtSqSum[2]*fr[5]-1.125*nuVtSqSum[2]*fl[5]+1.125*fr[2]*nuVtSqSum[3]-1.125*fl[2]*nuVtSqSum[3]-1.082531754730548*nuVtSqSum[1]*fr[3]-1.082531754730548*nuVtSqSum[1]*fl[3]+1.125*fr[0]*nuVtSqSum[1]-1.125*fl[0]*nuVtSqSum[1]+1.125*nuVtSqSum[0]*fr[1]-1.125*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[2]*(0.4330127018922193*favg[11]+0.25*favg[5])+alphaDrag[3]*(0.4330127018922193*favg[7]+0.25*favg[2])-0.8660254037844386*fjump[6]+alphaDrag[0]*(0.4330127018922193*favg[6]+0.25*favg[1])+alphaDrag[1]*(0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[1]; 
  Ghat[2] = ((-1.082531754730548*nuVtSqSum[1]*fr[11])-1.082531754730548*nuVtSqSum[1]*fl[11]-1.082531754730548*nuVtSqSum[0]*fr[7]-1.082531754730548*nuVtSqSum[0]*fl[7]-1.082531754730548*nuVtSqSum[3]*fr[6]-1.082531754730548*nuVtSqSum[3]*fl[6]+1.125*nuVtSqSum[1]*fr[5]-1.125*nuVtSqSum[1]*fl[5]+1.125*fr[1]*nuVtSqSum[3]-1.125*fl[1]*nuVtSqSum[3]-1.082531754730548*nuVtSqSum[2]*fr[3]-1.082531754730548*nuVtSqSum[2]*fl[3]+1.125*fr[0]*nuVtSqSum[2]-1.125*fl[0]*nuVtSqSum[2]+1.125*nuVtSqSum[0]*fr[2]-1.125*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[1]*(0.4330127018922193*favg[11]+0.25*favg[5])-0.8660254037844386*fjump[7]+alphaDrag[0]*(0.4330127018922193*favg[7]+0.25*favg[2])+alphaDrag[3]*(0.4330127018922193*favg[6]+0.25*favg[1])+alphaDrag[2]*(0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[2]; 
  Ghat[4] = ((-1.082531754730548*nuVtSqSum[3]*fr[15])-1.082531754730548*nuVtSqSum[3]*fl[15]-1.082531754730548*nuVtSqSum[2]*fr[14]-1.082531754730548*nuVtSqSum[2]*fl[14]-1.082531754730548*nuVtSqSum[1]*fr[13]-1.082531754730548*nuVtSqSum[1]*fl[13]+1.125*nuVtSqSum[3]*fr[12]-1.125*nuVtSqSum[3]*fl[12]-1.082531754730548*nuVtSqSum[0]*fr[10]-1.082531754730548*nuVtSqSum[0]*fl[10]+1.125*nuVtSqSum[2]*fr[9]-1.125*nuVtSqSum[2]*fl[9]+1.125*nuVtSqSum[1]*fr[8]-1.125*nuVtSqSum[1]*fl[8]+1.125*nuVtSqSum[0]*fr[4]-1.125*nuVtSqSum[0]*fl[4])*rdv+alphaDrag[3]*(0.4330127018922193*favg[15]+0.25*favg[12])+alphaDrag[2]*(0.4330127018922193*favg[14]+0.25*favg[9])+alphaDrag[1]*(0.4330127018922193*favg[13]+0.25*favg[8])-0.8660254037844386*fjump[10]+alphaDrag[0]*(0.4330127018922193*favg[10]+0.25*favg[4])-0.5*fjump[4]; 
  Ghat[5] = ((-1.082531754730548*nuVtSqSum[0]*fr[11])-1.082531754730548*nuVtSqSum[0]*fl[11]-1.082531754730548*nuVtSqSum[1]*fr[7]-1.082531754730548*nuVtSqSum[1]*fl[7]-1.082531754730548*nuVtSqSum[2]*fr[6]-1.082531754730548*nuVtSqSum[2]*fl[6]+1.125*nuVtSqSum[0]*fr[5]-1.125*nuVtSqSum[0]*fl[5]-1.082531754730548*fr[3]*nuVtSqSum[3]-1.082531754730548*fl[3]*nuVtSqSum[3]+1.125*fr[0]*nuVtSqSum[3]-1.125*fl[0]*nuVtSqSum[3]+1.125*fr[1]*nuVtSqSum[2]-1.125*fl[1]*nuVtSqSum[2]+1.125*nuVtSqSum[1]*fr[2]-1.125*nuVtSqSum[1]*fl[2])*rdv-0.8660254037844386*fjump[11]+alphaDrag[0]*(0.4330127018922193*favg[11]+0.25*favg[5])+alphaDrag[1]*(0.4330127018922193*favg[7]+0.25*favg[2])+alphaDrag[2]*(0.4330127018922193*favg[6]+0.25*favg[1])-0.5*fjump[5]+alphaDrag[3]*(0.4330127018922193*favg[3]+0.25*favg[0]); 
  Ghat[8] = ((-1.082531754730548*nuVtSqSum[2]*fr[15])-1.082531754730548*nuVtSqSum[2]*fl[15]-1.082531754730548*nuVtSqSum[3]*fr[14]-1.082531754730548*nuVtSqSum[3]*fl[14]-1.082531754730548*nuVtSqSum[0]*fr[13]-1.082531754730548*nuVtSqSum[0]*fl[13]+1.125*nuVtSqSum[2]*fr[12]-1.125*nuVtSqSum[2]*fl[12]-1.082531754730548*nuVtSqSum[1]*fr[10]-1.082531754730548*nuVtSqSum[1]*fl[10]+1.125*nuVtSqSum[3]*fr[9]-1.125*nuVtSqSum[3]*fl[9]+1.125*nuVtSqSum[0]*fr[8]-1.125*nuVtSqSum[0]*fl[8]+1.125*nuVtSqSum[1]*fr[4]-1.125*nuVtSqSum[1]*fl[4])*rdv+alphaDrag[2]*(0.4330127018922193*favg[15]+0.25*favg[12])+alphaDrag[3]*(0.4330127018922193*favg[14]+0.25*favg[9])-0.8660254037844386*fjump[13]+alphaDrag[0]*(0.4330127018922193*favg[13]+0.25*favg[8])+alphaDrag[1]*(0.4330127018922193*favg[10]+0.25*favg[4])-0.5*fjump[8]; 
  Ghat[9] = ((-1.082531754730548*nuVtSqSum[1]*fr[15])-1.082531754730548*nuVtSqSum[1]*fl[15]-1.082531754730548*nuVtSqSum[0]*fr[14]-1.082531754730548*nuVtSqSum[0]*fl[14]-1.082531754730548*nuVtSqSum[3]*fr[13]-1.082531754730548*nuVtSqSum[3]*fl[13]+1.125*nuVtSqSum[1]*fr[12]-1.125*nuVtSqSum[1]*fl[12]-1.082531754730548*nuVtSqSum[2]*fr[10]-1.082531754730548*nuVtSqSum[2]*fl[10]+1.125*nuVtSqSum[0]*fr[9]-1.125*nuVtSqSum[0]*fl[9]+1.125*nuVtSqSum[3]*fr[8]-1.125*nuVtSqSum[3]*fl[8]+1.125*nuVtSqSum[2]*fr[4]-1.125*nuVtSqSum[2]*fl[4])*rdv+alphaDrag[1]*(0.4330127018922193*favg[15]+0.25*favg[12])-0.8660254037844386*fjump[14]+alphaDrag[0]*(0.4330127018922193*favg[14]+0.25*favg[9])+alphaDrag[3]*(0.4330127018922193*favg[13]+0.25*favg[8])+alphaDrag[2]*(0.4330127018922193*favg[10]+0.25*favg[4])-0.5*fjump[9]; 
  Ghat[12] = ((-1.082531754730548*nuVtSqSum[0]*fr[15])-1.082531754730548*nuVtSqSum[0]*fl[15]-1.082531754730548*nuVtSqSum[1]*fr[14]-1.082531754730548*nuVtSqSum[1]*fl[14]-1.082531754730548*nuVtSqSum[2]*fr[13]-1.082531754730548*nuVtSqSum[2]*fl[13]+1.125*nuVtSqSum[0]*fr[12]-1.125*nuVtSqSum[0]*fl[12]-1.082531754730548*nuVtSqSum[3]*fr[10]-1.082531754730548*nuVtSqSum[3]*fl[10]+1.125*nuVtSqSum[1]*fr[9]-1.125*nuVtSqSum[1]*fl[9]+1.125*nuVtSqSum[2]*fr[8]-1.125*nuVtSqSum[2]*fl[8]+1.125*nuVtSqSum[3]*fr[4]-1.125*nuVtSqSum[3]*fl[4])*rdv-0.8660254037844386*fjump[15]+alphaDrag[0]*(0.4330127018922193*favg[15]+0.25*favg[12])+alphaDrag[1]*(0.4330127018922193*favg[14]+0.25*favg[9])+alphaDrag[2]*(0.4330127018922193*favg[13]+0.25*favg[8])-0.5*fjump[12]+alphaDrag[3]*(0.4330127018922193*favg[10]+0.25*favg[4]); 

  double incr1[16]; 
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
  incr1[11] = 0.8660254037844386*Ghat[5]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = 0.8660254037844386*Ghat[8]; 
  incr1[14] = 0.8660254037844386*Ghat[9]; 
  incr1[15] = 0.8660254037844386*Ghat[12]; 

  double incr2[16]; 
  incr2[3] = nuVtSqSum[3]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[2]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[6] = nuVtSqSum[2]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[3]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[7] = nuVtSqSum[1]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[0]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[3]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[10] = nuVtSqSum[3]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[2]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[1]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[0]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[11] = nuVtSqSum[0]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[1]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[2]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0]))*nuVtSqSum[3]; 
  incr2[13] = nuVtSqSum[2]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[3]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[0]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[1]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[14] = nuVtSqSum[1]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[0]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[3]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[2]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[15] = nuVtSqSum[0]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[1]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[2]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[3]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 

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
  outr[11] += incr2[11]*rdvSq4R+incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 

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
  outl[11] += incr1[11]*rdv2L-1.0*incr2[11]*rdvSq4L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += incr1[13]*rdv2L-1.0*incr2[13]*rdvSq4L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += incr1[15]*rdv2L-1.0*incr2[15]*rdvSq4L; 

  return std::abs(wl[2]-(0.5*sumNuUx[0])/nuSum); 
} 
double VmLBOconstNuSurf2x2vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[8]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[3]; 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  const double *sumNuUy = &nuUSum[4]; 

  double favg[16]; 
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
  favg[12] = -1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 

  double fjump[16]; 
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
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(-1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(-1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(-1*fr[15])); 

  double alphaDrag[4]; 
  alphaDrag[0] = 2.0*wl[3]*nuSum+dxvl[3]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 
  alphaDrag[3] = -1.0*sumNuUy[3]; 

  double Ghat[16]; 
  for(unsigned short int i=0; i<16; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.082531754730548*nuVtSqSum[3]*fr[12])-1.082531754730548*nuVtSqSum[3]*fl[12]-1.082531754730548*nuVtSqSum[2]*fr[9]-1.082531754730548*nuVtSqSum[2]*fl[9]-1.082531754730548*nuVtSqSum[1]*fr[8]-1.082531754730548*nuVtSqSum[1]*fl[8]+1.125*nuVtSqSum[3]*fr[5]-1.125*nuVtSqSum[3]*fl[5]-1.082531754730548*nuVtSqSum[0]*fr[4]-1.082531754730548*nuVtSqSum[0]*fl[4]+1.125*fr[2]*nuVtSqSum[2]-1.125*fl[2]*nuVtSqSum[2]+1.125*fr[1]*nuVtSqSum[1]-1.125*fl[1]*nuVtSqSum[1]+1.125*fr[0]*nuVtSqSum[0]-1.125*fl[0]*nuVtSqSum[0])*rdv+alphaDrag[3]*(0.4330127018922193*favg[12]+0.25*favg[5])+alphaDrag[2]*(0.4330127018922193*favg[9]+0.25*favg[2])+alphaDrag[1]*(0.4330127018922193*favg[8]+0.25*favg[1])-0.8660254037844386*fjump[4]+alphaDrag[0]*(0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[0]; 
  Ghat[1] = ((-1.082531754730548*nuVtSqSum[2]*fr[12])-1.082531754730548*nuVtSqSum[2]*fl[12]-1.082531754730548*nuVtSqSum[3]*fr[9]-1.082531754730548*nuVtSqSum[3]*fl[9]-1.082531754730548*nuVtSqSum[0]*fr[8]-1.082531754730548*nuVtSqSum[0]*fl[8]+1.125*nuVtSqSum[2]*fr[5]-1.125*nuVtSqSum[2]*fl[5]-1.082531754730548*nuVtSqSum[1]*fr[4]-1.082531754730548*nuVtSqSum[1]*fl[4]+1.125*fr[2]*nuVtSqSum[3]-1.125*fl[2]*nuVtSqSum[3]+1.125*fr[0]*nuVtSqSum[1]-1.125*fl[0]*nuVtSqSum[1]+1.125*nuVtSqSum[0]*fr[1]-1.125*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[2]*(0.4330127018922193*favg[12]+0.25*favg[5])+alphaDrag[3]*(0.4330127018922193*favg[9]+0.25*favg[2])-0.8660254037844386*fjump[8]+alphaDrag[0]*(0.4330127018922193*favg[8]+0.25*favg[1])+alphaDrag[1]*(0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[1]; 
  Ghat[2] = ((-1.082531754730548*nuVtSqSum[1]*fr[12])-1.082531754730548*nuVtSqSum[1]*fl[12]-1.082531754730548*nuVtSqSum[0]*fr[9]-1.082531754730548*nuVtSqSum[0]*fl[9]-1.082531754730548*nuVtSqSum[3]*fr[8]-1.082531754730548*nuVtSqSum[3]*fl[8]+1.125*nuVtSqSum[1]*fr[5]-1.125*nuVtSqSum[1]*fl[5]-1.082531754730548*nuVtSqSum[2]*fr[4]-1.082531754730548*nuVtSqSum[2]*fl[4]+1.125*fr[1]*nuVtSqSum[3]-1.125*fl[1]*nuVtSqSum[3]+1.125*fr[0]*nuVtSqSum[2]-1.125*fl[0]*nuVtSqSum[2]+1.125*nuVtSqSum[0]*fr[2]-1.125*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[1]*(0.4330127018922193*favg[12]+0.25*favg[5])-0.8660254037844386*fjump[9]+alphaDrag[0]*(0.4330127018922193*favg[9]+0.25*favg[2])+alphaDrag[3]*(0.4330127018922193*favg[8]+0.25*favg[1])+alphaDrag[2]*(0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[2]; 
  Ghat[3] = ((-1.082531754730548*nuVtSqSum[3]*fr[15])-1.082531754730548*nuVtSqSum[3]*fl[15]-1.082531754730548*nuVtSqSum[2]*fr[14]-1.082531754730548*nuVtSqSum[2]*fl[14]-1.082531754730548*nuVtSqSum[1]*fr[13]-1.082531754730548*nuVtSqSum[1]*fl[13]+1.125*nuVtSqSum[3]*fr[11]-1.125*nuVtSqSum[3]*fl[11]-1.082531754730548*nuVtSqSum[0]*fr[10]-1.082531754730548*nuVtSqSum[0]*fl[10]+1.125*nuVtSqSum[2]*fr[7]-1.125*nuVtSqSum[2]*fl[7]+1.125*nuVtSqSum[1]*fr[6]-1.125*nuVtSqSum[1]*fl[6]+1.125*nuVtSqSum[0]*fr[3]-1.125*nuVtSqSum[0]*fl[3])*rdv+alphaDrag[3]*(0.4330127018922193*favg[15]+0.25*favg[11])+alphaDrag[2]*(0.4330127018922193*favg[14]+0.25*favg[7])+alphaDrag[1]*(0.4330127018922193*favg[13]+0.25*favg[6])-0.8660254037844386*fjump[10]+alphaDrag[0]*(0.4330127018922193*favg[10]+0.25*favg[3])-0.5*fjump[3]; 
  Ghat[5] = ((-1.082531754730548*nuVtSqSum[0]*fr[12])-1.082531754730548*nuVtSqSum[0]*fl[12]-1.082531754730548*nuVtSqSum[1]*fr[9]-1.082531754730548*nuVtSqSum[1]*fl[9]-1.082531754730548*nuVtSqSum[2]*fr[8]-1.082531754730548*nuVtSqSum[2]*fl[8]+1.125*nuVtSqSum[0]*fr[5]-1.125*nuVtSqSum[0]*fl[5]-1.082531754730548*nuVtSqSum[3]*fr[4]-1.082531754730548*nuVtSqSum[3]*fl[4]+1.125*fr[0]*nuVtSqSum[3]-1.125*fl[0]*nuVtSqSum[3]+1.125*fr[1]*nuVtSqSum[2]-1.125*fl[1]*nuVtSqSum[2]+1.125*nuVtSqSum[1]*fr[2]-1.125*nuVtSqSum[1]*fl[2])*rdv-0.8660254037844386*fjump[12]+alphaDrag[0]*(0.4330127018922193*favg[12]+0.25*favg[5])+alphaDrag[1]*(0.4330127018922193*favg[9]+0.25*favg[2])+alphaDrag[2]*(0.4330127018922193*favg[8]+0.25*favg[1])-0.5*fjump[5]+alphaDrag[3]*(0.4330127018922193*favg[4]+0.25*favg[0]); 
  Ghat[6] = ((-1.082531754730548*nuVtSqSum[2]*fr[15])-1.082531754730548*nuVtSqSum[2]*fl[15]-1.082531754730548*nuVtSqSum[3]*fr[14]-1.082531754730548*nuVtSqSum[3]*fl[14]-1.082531754730548*nuVtSqSum[0]*fr[13]-1.082531754730548*nuVtSqSum[0]*fl[13]+1.125*nuVtSqSum[2]*fr[11]-1.125*nuVtSqSum[2]*fl[11]-1.082531754730548*nuVtSqSum[1]*fr[10]-1.082531754730548*nuVtSqSum[1]*fl[10]+1.125*nuVtSqSum[3]*fr[7]-1.125*nuVtSqSum[3]*fl[7]+1.125*nuVtSqSum[0]*fr[6]-1.125*nuVtSqSum[0]*fl[6]+1.125*nuVtSqSum[1]*fr[3]-1.125*nuVtSqSum[1]*fl[3])*rdv+alphaDrag[2]*(0.4330127018922193*favg[15]+0.25*favg[11])+alphaDrag[3]*(0.4330127018922193*favg[14]+0.25*favg[7])-0.8660254037844386*fjump[13]+alphaDrag[0]*(0.4330127018922193*favg[13]+0.25*favg[6])+alphaDrag[1]*(0.4330127018922193*favg[10]+0.25*favg[3])-0.5*fjump[6]; 
  Ghat[7] = ((-1.082531754730548*nuVtSqSum[1]*fr[15])-1.082531754730548*nuVtSqSum[1]*fl[15]-1.082531754730548*nuVtSqSum[0]*fr[14]-1.082531754730548*nuVtSqSum[0]*fl[14]-1.082531754730548*nuVtSqSum[3]*fr[13]-1.082531754730548*nuVtSqSum[3]*fl[13]+1.125*nuVtSqSum[1]*fr[11]-1.125*nuVtSqSum[1]*fl[11]-1.082531754730548*nuVtSqSum[2]*fr[10]-1.082531754730548*nuVtSqSum[2]*fl[10]+1.125*nuVtSqSum[0]*fr[7]-1.125*nuVtSqSum[0]*fl[7]+1.125*nuVtSqSum[3]*fr[6]-1.125*nuVtSqSum[3]*fl[6]+1.125*nuVtSqSum[2]*fr[3]-1.125*nuVtSqSum[2]*fl[3])*rdv+alphaDrag[1]*(0.4330127018922193*favg[15]+0.25*favg[11])-0.8660254037844386*fjump[14]+alphaDrag[0]*(0.4330127018922193*favg[14]+0.25*favg[7])+alphaDrag[3]*(0.4330127018922193*favg[13]+0.25*favg[6])+alphaDrag[2]*(0.4330127018922193*favg[10]+0.25*favg[3])-0.5*fjump[7]; 
  Ghat[11] = ((-1.082531754730548*nuVtSqSum[0]*fr[15])-1.082531754730548*nuVtSqSum[0]*fl[15]-1.082531754730548*nuVtSqSum[1]*fr[14]-1.082531754730548*nuVtSqSum[1]*fl[14]-1.082531754730548*nuVtSqSum[2]*fr[13]-1.082531754730548*nuVtSqSum[2]*fl[13]+1.125*nuVtSqSum[0]*fr[11]-1.125*nuVtSqSum[0]*fl[11]-1.082531754730548*nuVtSqSum[3]*fr[10]-1.082531754730548*nuVtSqSum[3]*fl[10]+1.125*nuVtSqSum[1]*fr[7]-1.125*nuVtSqSum[1]*fl[7]+1.125*nuVtSqSum[2]*fr[6]-1.125*nuVtSqSum[2]*fl[6]+1.125*fr[3]*nuVtSqSum[3]-1.125*fl[3]*nuVtSqSum[3])*rdv-0.8660254037844386*fjump[15]+alphaDrag[0]*(0.4330127018922193*favg[15]+0.25*favg[11])+alphaDrag[1]*(0.4330127018922193*favg[14]+0.25*favg[7])+alphaDrag[2]*(0.4330127018922193*favg[13]+0.25*favg[6])-0.5*fjump[11]+alphaDrag[3]*(0.4330127018922193*favg[10]+0.25*favg[3]); 

  double incr1[16]; 
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
  incr1[12] = 0.8660254037844386*Ghat[5]; 
  incr1[13] = 0.8660254037844386*Ghat[6]; 
  incr1[14] = 0.8660254037844386*Ghat[7]; 
  incr1[15] = 0.8660254037844386*Ghat[11]; 

  double incr2[16]; 
  incr2[4] = nuVtSqSum[3]*((-0.25*fr[12])+0.25*fl[12]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[2]*((-0.25*fr[9])+0.25*fl[9]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.25*fr[8])+0.25*fl[8]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.25*fr[4])+0.25*fl[4]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[8] = nuVtSqSum[2]*((-0.25*fr[12])+0.25*fl[12]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[3]*((-0.25*fr[9])+0.25*fl[9]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.25*fr[8])+0.25*fl[8]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.25*fr[4])+0.25*fl[4]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[9] = nuVtSqSum[1]*((-0.25*fr[12])+0.25*fl[12]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[0]*((-0.25*fr[9])+0.25*fl[9]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[3]*((-0.25*fr[8])+0.25*fl[8]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.25*fr[4])+0.25*fl[4]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[10] = nuVtSqSum[3]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[11]+fl[11]))+nuVtSqSum[2]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[1]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[0]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[3]+fl[3])); 
  incr2[12] = nuVtSqSum[0]*((-0.25*fr[12])+0.25*fl[12]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[1]*((-0.25*fr[9])+0.25*fl[9]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[2]*((-0.25*fr[8])+0.25*fl[8]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[3]*((-0.25*fr[4])+0.25*fl[4]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[13] = nuVtSqSum[2]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[11]+fl[11]))+nuVtSqSum[3]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[0]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[1]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[3]+fl[3])); 
  incr2[14] = nuVtSqSum[1]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[11]+fl[11]))+nuVtSqSum[0]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[3]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[2]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[3]+fl[3])); 
  incr2[15] = nuVtSqSum[0]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[11]+fl[11]))+nuVtSqSum[1]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[7]+fl[7]))+nuVtSqSum[2]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[6]+fl[6]))+nuVtSqSum[3]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[3]+fl[3])); 

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
  outr[12] += incr2[12]*rdvSq4R+incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 

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
  outl[12] += incr1[12]*rdv2L-1.0*incr2[12]*rdvSq4L; 
  outl[13] += incr1[13]*rdv2L-1.0*incr2[13]*rdvSq4L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += incr1[15]*rdv2L-1.0*incr2[15]*rdvSq4L; 

  return std::abs(wl[3]-(0.5*sumNuUy[0])/nuSum); 
} 
