#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf2x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:          Cell-center coordinates. 
  // dxv[5]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[18]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[6]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[21]; 
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
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 

  double fjump[21]; 
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
  fjump[16] = nuSum*vMuMidMax*(fl[16]-(1*fr[16])); 
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(1*fr[19])); 
  fjump[20] = nuSum*vMuMidMax*(fl[20]-(1*fr[20])); 

  double alphaDrag[6]; 
  alphaDrag[0] = 2.0*wl[2]*nuSum+dxvl[2]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 
  alphaDrag[3] = -1.0*sumNuUx[3]; 
  alphaDrag[4] = -1.0*sumNuUx[4]; 
  alphaDrag[5] = -1.0*sumNuUx[5]; 

  double Ghat[21]; 
  for(unsigned short int i=0; i<21; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.341640786499874*nuVtSqSum[0]*fr[18]-1.341640786499874*nuVtSqSum[0]*fl[18]+1.875*nuVtSqSum[5]*fr[17]-1.875*nuVtSqSum[5]*fl[17]+1.875*nuVtSqSum[4]*fr[16]-1.875*nuVtSqSum[4]*fl[16]-2.381569860407206*nuVtSqSum[2]*fr[8]-2.381569860407206*nuVtSqSum[2]*fl[8]-2.381569860407206*nuVtSqSum[1]*fr[7]-2.381569860407206*nuVtSqSum[1]*fl[7]+1.875*nuVtSqSum[3]*fr[6]-1.875*nuVtSqSum[3]*fl[6]-2.381569860407206*nuVtSqSum[0]*fr[3]-2.381569860407206*nuVtSqSum[0]*fl[3]+1.875*fr[2]*nuVtSqSum[2]-1.875*fl[2]*nuVtSqSum[2]+1.875*fr[1]*nuVtSqSum[1]-1.875*fl[1]*nuVtSqSum[1]+1.875*fr[0]*nuVtSqSum[0]-1.875*fl[0]*nuVtSqSum[0])*rdv-1.118033988749895*fjump[18]+alphaDrag[0]*(0.5590169943749475*favg[18]+0.4330127018922193*favg[3]+0.25*favg[0])+0.25*alphaDrag[5]*favg[17]+0.25*alphaDrag[4]*favg[16]+alphaDrag[2]*(0.4330127018922193*favg[8]+0.25*favg[2])+alphaDrag[1]*(0.4330127018922193*favg[7]+0.25*favg[1])+0.25*alphaDrag[3]*favg[6]-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = (1.341640786499874*nuVtSqSum[1]*fr[18]-1.341640786499874*nuVtSqSum[1]*fl[18]+1.677050983124842*nuVtSqSum[1]*fr[16]-1.677050983124842*nuVtSqSum[1]*fl[16]-2.381569860407206*nuVtSqSum[3]*fr[8]-2.381569860407206*nuVtSqSum[3]*fl[8]-2.130140840414079*nuVtSqSum[4]*fr[7]-2.381569860407206*nuVtSqSum[0]*fr[7]-2.130140840414079*nuVtSqSum[4]*fl[7]-2.381569860407206*nuVtSqSum[0]*fl[7]+1.875*nuVtSqSum[2]*fr[6]-1.875*nuVtSqSum[2]*fl[6]+1.677050983124842*fr[1]*nuVtSqSum[4]-1.677050983124842*fl[1]*nuVtSqSum[4]+1.875*fr[2]*nuVtSqSum[3]-1.875*fl[2]*nuVtSqSum[3]-2.381569860407206*nuVtSqSum[1]*fr[3]-2.381569860407206*nuVtSqSum[1]*fl[3]+1.875*fr[0]*nuVtSqSum[1]-1.875*fl[0]*nuVtSqSum[1]+1.875*nuVtSqSum[0]*fr[1]-1.875*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.5590169943749475*favg[18]+0.223606797749979*favg[16]+0.4330127018922193*favg[3]+0.25*favg[0])+alphaDrag[3]*(0.4330127018922193*favg[8]+0.25*favg[2])-0.8660254037844386*fjump[7]+alphaDrag[0]*(0.4330127018922193*favg[7]+0.25*favg[1])+alphaDrag[4]*(0.3872983346207416*favg[7]+0.223606797749979*favg[1])+0.25*alphaDrag[2]*favg[6]-0.5*fjump[1]; 
  Ghat[2] = (1.341640786499874*nuVtSqSum[2]*fr[18]-1.341640786499874*nuVtSqSum[2]*fl[18]+1.677050983124842*nuVtSqSum[2]*fr[17]-1.677050983124842*nuVtSqSum[2]*fl[17]-2.130140840414079*nuVtSqSum[5]*fr[8]-2.381569860407206*nuVtSqSum[0]*fr[8]-2.130140840414079*nuVtSqSum[5]*fl[8]-2.381569860407206*nuVtSqSum[0]*fl[8]-2.381569860407206*nuVtSqSum[3]*fr[7]-2.381569860407206*nuVtSqSum[3]*fl[7]+1.875*nuVtSqSum[1]*fr[6]-1.875*nuVtSqSum[1]*fl[6]+1.677050983124842*fr[2]*nuVtSqSum[5]-1.677050983124842*fl[2]*nuVtSqSum[5]+1.875*fr[1]*nuVtSqSum[3]-1.875*fl[1]*nuVtSqSum[3]-2.381569860407206*nuVtSqSum[2]*fr[3]-2.381569860407206*nuVtSqSum[2]*fl[3]+1.875*fr[0]*nuVtSqSum[2]-1.875*fl[0]*nuVtSqSum[2]+1.875*nuVtSqSum[0]*fr[2]-1.875*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[2]*(0.5590169943749475*favg[18]+0.223606797749979*favg[17]+0.4330127018922193*favg[3]+0.25*favg[0])-0.8660254037844386*fjump[8]+alphaDrag[0]*(0.4330127018922193*favg[8]+0.25*favg[2])+alphaDrag[5]*(0.3872983346207416*favg[8]+0.223606797749979*favg[2])+alphaDrag[3]*(0.4330127018922193*favg[7]+0.25*favg[1])+0.25*alphaDrag[1]*favg[6]-0.5*fjump[2]; 
  Ghat[4] = ((-2.381569860407206*nuVtSqSum[0]*fr[11])-2.381569860407206*nuVtSqSum[0]*fl[11]+1.875*nuVtSqSum[2]*fr[10]-1.875*nuVtSqSum[2]*fl[10]+1.875*nuVtSqSum[1]*fr[9]-1.875*nuVtSqSum[1]*fl[9]+1.875*nuVtSqSum[0]*fr[4]-1.875*nuVtSqSum[0]*fl[4])*rdv-0.8660254037844386*fjump[11]+alphaDrag[0]*(0.4330127018922193*favg[11]+0.25*favg[4])+0.25*alphaDrag[2]*favg[10]+0.25*alphaDrag[1]*favg[9]-0.5*fjump[4]; 
  Ghat[5] = ((-2.381569860407206*nuVtSqSum[0]*fr[14])-2.381569860407206*nuVtSqSum[0]*fl[14]+1.875*nuVtSqSum[2]*fr[13]-1.875*nuVtSqSum[2]*fl[13]+1.875*nuVtSqSum[1]*fr[12]-1.875*nuVtSqSum[1]*fl[12]+1.875*nuVtSqSum[0]*fr[5]-1.875*nuVtSqSum[0]*fl[5])*rdv-0.8660254037844386*fjump[14]+alphaDrag[0]*(0.4330127018922193*favg[14]+0.25*favg[5])+0.25*alphaDrag[2]*favg[13]+0.25*alphaDrag[1]*favg[12]-0.5*fjump[5]; 
  Ghat[6] = (1.341640786499874*nuVtSqSum[3]*fr[18]-1.341640786499874*nuVtSqSum[3]*fl[18]+1.677050983124842*nuVtSqSum[3]*fr[17]-1.677050983124842*nuVtSqSum[3]*fl[17]+1.677050983124842*nuVtSqSum[3]*fr[16]-1.677050983124842*nuVtSqSum[3]*fl[16]-2.381569860407206*nuVtSqSum[1]*fr[8]-2.381569860407206*nuVtSqSum[1]*fl[8]-2.381569860407206*nuVtSqSum[2]*fr[7]-2.381569860407206*nuVtSqSum[2]*fl[7]+1.677050983124842*nuVtSqSum[5]*fr[6]+1.677050983124842*nuVtSqSum[4]*fr[6]+1.875*nuVtSqSum[0]*fr[6]-1.677050983124842*nuVtSqSum[5]*fl[6]-1.677050983124842*nuVtSqSum[4]*fl[6]-1.875*nuVtSqSum[0]*fl[6]-2.381569860407206*fr[3]*nuVtSqSum[3]-2.381569860407206*fl[3]*nuVtSqSum[3]+1.875*fr[0]*nuVtSqSum[3]-1.875*fl[0]*nuVtSqSum[3]+1.875*fr[1]*nuVtSqSum[2]-1.875*fl[1]*nuVtSqSum[2]+1.875*nuVtSqSum[1]*fr[2]-1.875*nuVtSqSum[1]*fl[2])*rdv+alphaDrag[3]*(0.5590169943749475*favg[18]+0.223606797749979*favg[17]+0.223606797749979*favg[16]+0.4330127018922193*favg[3]+0.25*favg[0])+alphaDrag[1]*(0.4330127018922193*favg[8]+0.25*favg[2])+alphaDrag[2]*(0.4330127018922193*favg[7]+0.25*favg[1])-0.5*fjump[6]+0.223606797749979*alphaDrag[5]*favg[6]+0.223606797749979*alphaDrag[4]*favg[6]+0.25*alphaDrag[0]*favg[6]; 
  Ghat[9] = ((-2.381569860407206*nuVtSqSum[1]*fr[11])-2.381569860407206*nuVtSqSum[1]*fl[11]+1.875*nuVtSqSum[3]*fr[10]-1.875*nuVtSqSum[3]*fl[10]+1.677050983124842*nuVtSqSum[4]*fr[9]+1.875*nuVtSqSum[0]*fr[9]-1.677050983124842*nuVtSqSum[4]*fl[9]-1.875*nuVtSqSum[0]*fl[9]+1.875*nuVtSqSum[1]*fr[4]-1.875*nuVtSqSum[1]*fl[4])*rdv+alphaDrag[1]*(0.4330127018922193*favg[11]+0.25*favg[4])+0.25*alphaDrag[3]*favg[10]-0.5*fjump[9]+0.223606797749979*alphaDrag[4]*favg[9]+0.25*alphaDrag[0]*favg[9]; 
  Ghat[10] = ((-2.381569860407206*nuVtSqSum[2]*fr[11])-2.381569860407206*nuVtSqSum[2]*fl[11]+1.677050983124842*nuVtSqSum[5]*fr[10]+1.875*nuVtSqSum[0]*fr[10]-1.677050983124842*nuVtSqSum[5]*fl[10]-1.875*nuVtSqSum[0]*fl[10]+1.875*nuVtSqSum[3]*fr[9]-1.875*nuVtSqSum[3]*fl[9]+1.875*nuVtSqSum[2]*fr[4]-1.875*nuVtSqSum[2]*fl[4])*rdv+alphaDrag[2]*(0.4330127018922193*favg[11]+0.25*favg[4])-0.5*fjump[10]+0.223606797749979*alphaDrag[5]*favg[10]+0.25*alphaDrag[0]*favg[10]+0.25*alphaDrag[3]*favg[9]; 
  Ghat[12] = ((-2.381569860407206*nuVtSqSum[1]*fr[14])-2.381569860407206*nuVtSqSum[1]*fl[14]+1.875*nuVtSqSum[3]*fr[13]-1.875*nuVtSqSum[3]*fl[13]+1.677050983124842*nuVtSqSum[4]*fr[12]+1.875*nuVtSqSum[0]*fr[12]-1.677050983124842*nuVtSqSum[4]*fl[12]-1.875*nuVtSqSum[0]*fl[12]+1.875*nuVtSqSum[1]*fr[5]-1.875*nuVtSqSum[1]*fl[5])*rdv+alphaDrag[1]*(0.4330127018922193*favg[14]+0.25*favg[5])+0.25*alphaDrag[3]*favg[13]-0.5*fjump[12]+0.223606797749979*alphaDrag[4]*favg[12]+0.25*alphaDrag[0]*favg[12]; 
  Ghat[13] = ((-2.381569860407206*nuVtSqSum[2]*fr[14])-2.381569860407206*nuVtSqSum[2]*fl[14]+1.677050983124842*nuVtSqSum[5]*fr[13]+1.875*nuVtSqSum[0]*fr[13]-1.677050983124842*nuVtSqSum[5]*fl[13]-1.875*nuVtSqSum[0]*fl[13]+1.875*nuVtSqSum[3]*fr[12]-1.875*nuVtSqSum[3]*fl[12]+1.875*nuVtSqSum[2]*fr[5]-1.875*nuVtSqSum[2]*fl[5])*rdv+alphaDrag[2]*(0.4330127018922193*favg[14]+0.25*favg[5])-0.5*fjump[13]+0.223606797749979*alphaDrag[5]*favg[13]+0.25*alphaDrag[0]*favg[13]+0.25*alphaDrag[3]*favg[12]; 
  Ghat[15] = (1.875*nuVtSqSum[0]*fr[15]-1.875*nuVtSqSum[0]*fl[15])*rdv-0.5*fjump[15]+0.25*alphaDrag[0]*favg[15]; 
  Ghat[16] = (1.341640786499874*nuVtSqSum[4]*fr[18]-1.341640786499874*nuVtSqSum[4]*fl[18]+1.197893559374888*nuVtSqSum[4]*fr[16]+1.875*nuVtSqSum[0]*fr[16]-1.197893559374888*nuVtSqSum[4]*fl[16]-1.875*nuVtSqSum[0]*fl[16]-2.130140840414079*nuVtSqSum[1]*fr[7]-2.130140840414079*nuVtSqSum[1]*fl[7]+1.677050983124842*nuVtSqSum[3]*fr[6]-1.677050983124842*nuVtSqSum[3]*fl[6]-2.381569860407206*fr[3]*nuVtSqSum[4]-2.381569860407206*fl[3]*nuVtSqSum[4]+1.875*fr[0]*nuVtSqSum[4]-1.875*fl[0]*nuVtSqSum[4]+1.677050983124842*fr[1]*nuVtSqSum[1]-1.677050983124842*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[4]*(0.5590169943749475*favg[18]+0.159719141249985*favg[16]+0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[16]+0.25*alphaDrag[0]*favg[16]+alphaDrag[1]*(0.3872983346207416*favg[7]+0.223606797749979*favg[1])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[17] = (1.341640786499874*nuVtSqSum[5]*fr[18]-1.341640786499874*nuVtSqSum[5]*fl[18]+1.197893559374888*nuVtSqSum[5]*fr[17]+1.875*nuVtSqSum[0]*fr[17]-1.197893559374888*nuVtSqSum[5]*fl[17]-1.875*nuVtSqSum[0]*fl[17]-2.130140840414079*nuVtSqSum[2]*fr[8]-2.130140840414079*nuVtSqSum[2]*fl[8]+1.677050983124842*nuVtSqSum[3]*fr[6]-1.677050983124842*nuVtSqSum[3]*fl[6]-2.381569860407206*fr[3]*nuVtSqSum[5]-2.381569860407206*fl[3]*nuVtSqSum[5]+1.875*fr[0]*nuVtSqSum[5]-1.875*fl[0]*nuVtSqSum[5]+1.677050983124842*fr[2]*nuVtSqSum[2]-1.677050983124842*fl[2]*nuVtSqSum[2])*rdv+alphaDrag[5]*(0.5590169943749475*favg[18]+0.159719141249985*favg[17]+0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[17]+0.25*alphaDrag[0]*favg[17]+alphaDrag[2]*(0.3872983346207416*favg[8]+0.223606797749979*favg[2])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[19] = (1.875*nuVtSqSum[0]*fr[19]-1.875*nuVtSqSum[0]*fl[19])*rdv-0.5*fjump[19]+0.25*alphaDrag[0]*favg[19]; 
  Ghat[20] = (1.875*nuVtSqSum[0]*fr[20]-1.875*nuVtSqSum[0]*fl[20])*rdv-0.5*fjump[20]+0.25*alphaDrag[0]*favg[20]; 

  double incr1[21]; 
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
  incr1[16] = -0.5*Ghat[16]; 
  incr1[17] = -0.5*Ghat[17]; 
  incr1[18] = -1.118033988749895*Ghat[0]; 
  incr1[19] = -0.5*Ghat[19]; 
  incr1[20] = -0.5*Ghat[20]; 

  double incr2[21]; 
  incr2[3] = nuVtSqSum[0]*(0.2118037767457181*(fr[18]+fl[18])-0.3046875*fr[3]+0.3046875*fl[3]+0.2165063509461096*(fr[0]+fl[0]))+0.2165063509461096*(nuVtSqSum[5]*(fr[17]+fl[17])+nuVtSqSum[4]*(fr[16]+fl[16]))+nuVtSqSum[2]*((-0.3046875*fr[8])+0.3046875*fl[8]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.3046875*fr[7])+0.3046875*fl[7]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*nuVtSqSum[3]*(fr[6]+fl[6]); 
  incr2[7] = nuVtSqSum[1]*(0.2118037767457181*(fr[18]+fl[18])+0.1936491673103708*(fr[16]+fl[16])-0.3046875*fr[3]+0.3046875*fl[3]+0.2165063509461096*(fr[0]+fl[0]))+nuVtSqSum[3]*((-0.3046875*fr[8])+0.3046875*fl[8]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[4]*((-0.2725207847577868*fr[7])+0.2725207847577868*fl[7]+0.1936491673103708*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.3046875*fr[7])+0.3046875*fl[7]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*nuVtSqSum[2]*(fr[6]+fl[6]); 
  incr2[8] = nuVtSqSum[2]*(0.2118037767457181*(fr[18]+fl[18])+0.1936491673103708*(fr[17]+fl[17])-0.3046875*fr[3]+0.3046875*fl[3]+0.2165063509461096*(fr[0]+fl[0]))+nuVtSqSum[5]*((-0.2725207847577868*fr[8])+0.2725207847577868*fl[8]+0.1936491673103708*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.3046875*fr[8])+0.3046875*fl[8]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[3]*((-0.3046875*fr[7])+0.3046875*fl[7]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*nuVtSqSum[1]*(fr[6]+fl[6]); 
  incr2[11] = nuVtSqSum[0]*((-0.3046875*fr[11])+0.3046875*fl[11]+0.2165063509461096*(fr[4]+fl[4]))+0.2165063509461096*(nuVtSqSum[2]*(fr[10]+fl[10])+nuVtSqSum[1]*(fr[9]+fl[9])); 
  incr2[14] = nuVtSqSum[0]*((-0.3046875*fr[14])+0.3046875*fl[14]+0.2165063509461096*(fr[5]+fl[5]))+0.2165063509461096*(nuVtSqSum[2]*(fr[13]+fl[13])+nuVtSqSum[1]*(fr[12]+fl[12])); 
  incr2[18] = nuVtSqSum[0]*((-0.8203125*(fr[18]+fl[18]))+1.180049613297572*fr[3]-1.180049613297572*fl[3]-0.8385254915624212*(fr[0]+fl[0]))-0.8385254915624212*(nuVtSqSum[5]*(fr[17]+fl[17])+nuVtSqSum[4]*(fr[16]+fl[16]))+nuVtSqSum[2]*(1.180049613297572*fr[8]-1.180049613297572*fl[8]-0.8385254915624212*(fr[2]+fl[2]))+nuVtSqSum[1]*(1.180049613297572*fr[7]-1.180049613297572*fl[7]-0.8385254915624212*(fr[1]+fl[1]))-0.8385254915624212*nuVtSqSum[3]*(fr[6]+fl[6]); 

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
  outr[16] += incr1[16]*rdv2R; 
  outr[17] += incr1[17]*rdv2R; 
  outr[18] += incr2[18]*rdvSq4R+incr1[18]*rdv2R; 
  outr[19] += incr1[19]*rdv2R; 
  outr[20] += incr1[20]*rdv2R; 

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
  outl[16] += -1.0*incr1[16]*rdv2L; 
  outl[17] += -1.0*incr1[17]*rdv2L; 
  outl[18] += incr2[18]*rdvSq4L-1.0*incr1[18]*rdv2L; 
  outl[19] += -1.0*incr1[19]*rdv2L; 
  outl[20] += -1.0*incr1[20]*rdv2L; 

  return std::abs((0.5590169943749475*sumNuUx[5])/nuSum+(0.5590169943749475*sumNuUx[4])/nuSum-(0.5*sumNuUx[0])/nuSum+wl[2]); 
} 
double VmLBOconstNuSurf2x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:          Cell-center coordinates. 
  // dxv[5]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[18]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[6]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[3]; 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  const double *sumNuUy = &nuUSum[6]; 

  double favg[21]; 
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
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 

  double fjump[21]; 
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
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(1*fr[19])); 
  fjump[20] = nuSum*vMuMidMax*(fl[20]-(1*fr[20])); 

  double alphaDrag[6]; 
  alphaDrag[0] = 2.0*wl[3]*nuSum+dxvl[3]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 
  alphaDrag[3] = -1.0*sumNuUy[3]; 
  alphaDrag[4] = -1.0*sumNuUy[4]; 
  alphaDrag[5] = -1.0*sumNuUy[5]; 

  double Ghat[21]; 
  for(unsigned short int i=0; i<21; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.341640786499874*nuVtSqSum[0]*fr[19]-1.341640786499874*nuVtSqSum[0]*fl[19]+1.875*nuVtSqSum[5]*fr[17]-1.875*nuVtSqSum[5]*fl[17]+1.875*nuVtSqSum[4]*fr[16]-1.875*nuVtSqSum[4]*fl[16]-2.381569860407206*nuVtSqSum[2]*fr[10]-2.381569860407206*nuVtSqSum[2]*fl[10]-2.381569860407206*nuVtSqSum[1]*fr[9]-2.381569860407206*nuVtSqSum[1]*fl[9]+1.875*nuVtSqSum[3]*fr[6]-1.875*nuVtSqSum[3]*fl[6]-2.381569860407206*nuVtSqSum[0]*fr[4]-2.381569860407206*nuVtSqSum[0]*fl[4]+1.875*fr[2]*nuVtSqSum[2]-1.875*fl[2]*nuVtSqSum[2]+1.875*fr[1]*nuVtSqSum[1]-1.875*fl[1]*nuVtSqSum[1]+1.875*fr[0]*nuVtSqSum[0]-1.875*fl[0]*nuVtSqSum[0])*rdv-1.118033988749895*fjump[19]+alphaDrag[0]*(0.5590169943749475*favg[19]+0.4330127018922193*favg[4]+0.25*favg[0])+0.25*alphaDrag[5]*favg[17]+0.25*alphaDrag[4]*favg[16]+alphaDrag[2]*(0.4330127018922193*favg[10]+0.25*favg[2])+alphaDrag[1]*(0.4330127018922193*favg[9]+0.25*favg[1])+0.25*alphaDrag[3]*favg[6]-0.8660254037844386*fjump[4]-0.5*fjump[0]; 
  Ghat[1] = (1.341640786499874*nuVtSqSum[1]*fr[19]-1.341640786499874*nuVtSqSum[1]*fl[19]+1.677050983124842*nuVtSqSum[1]*fr[16]-1.677050983124842*nuVtSqSum[1]*fl[16]-2.381569860407206*nuVtSqSum[3]*fr[10]-2.381569860407206*nuVtSqSum[3]*fl[10]-2.130140840414079*nuVtSqSum[4]*fr[9]-2.381569860407206*nuVtSqSum[0]*fr[9]-2.130140840414079*nuVtSqSum[4]*fl[9]-2.381569860407206*nuVtSqSum[0]*fl[9]+1.875*nuVtSqSum[2]*fr[6]-1.875*nuVtSqSum[2]*fl[6]+1.677050983124842*fr[1]*nuVtSqSum[4]-1.677050983124842*fl[1]*nuVtSqSum[4]-2.381569860407206*nuVtSqSum[1]*fr[4]-2.381569860407206*nuVtSqSum[1]*fl[4]+1.875*fr[2]*nuVtSqSum[3]-1.875*fl[2]*nuVtSqSum[3]+1.875*fr[0]*nuVtSqSum[1]-1.875*fl[0]*nuVtSqSum[1]+1.875*nuVtSqSum[0]*fr[1]-1.875*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.5590169943749475*favg[19]+0.223606797749979*favg[16]+0.4330127018922193*favg[4]+0.25*favg[0])+alphaDrag[3]*(0.4330127018922193*favg[10]+0.25*favg[2])-0.8660254037844386*fjump[9]+alphaDrag[0]*(0.4330127018922193*favg[9]+0.25*favg[1])+alphaDrag[4]*(0.3872983346207416*favg[9]+0.223606797749979*favg[1])+0.25*alphaDrag[2]*favg[6]-0.5*fjump[1]; 
  Ghat[2] = (1.341640786499874*nuVtSqSum[2]*fr[19]-1.341640786499874*nuVtSqSum[2]*fl[19]+1.677050983124842*nuVtSqSum[2]*fr[17]-1.677050983124842*nuVtSqSum[2]*fl[17]-2.130140840414079*nuVtSqSum[5]*fr[10]-2.381569860407206*nuVtSqSum[0]*fr[10]-2.130140840414079*nuVtSqSum[5]*fl[10]-2.381569860407206*nuVtSqSum[0]*fl[10]-2.381569860407206*nuVtSqSum[3]*fr[9]-2.381569860407206*nuVtSqSum[3]*fl[9]+1.875*nuVtSqSum[1]*fr[6]-1.875*nuVtSqSum[1]*fl[6]+1.677050983124842*fr[2]*nuVtSqSum[5]-1.677050983124842*fl[2]*nuVtSqSum[5]-2.381569860407206*nuVtSqSum[2]*fr[4]-2.381569860407206*nuVtSqSum[2]*fl[4]+1.875*fr[1]*nuVtSqSum[3]-1.875*fl[1]*nuVtSqSum[3]+1.875*fr[0]*nuVtSqSum[2]-1.875*fl[0]*nuVtSqSum[2]+1.875*nuVtSqSum[0]*fr[2]-1.875*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[2]*(0.5590169943749475*favg[19]+0.223606797749979*favg[17]+0.4330127018922193*favg[4]+0.25*favg[0])-0.8660254037844386*fjump[10]+alphaDrag[0]*(0.4330127018922193*favg[10]+0.25*favg[2])+alphaDrag[5]*(0.3872983346207416*favg[10]+0.223606797749979*favg[2])+alphaDrag[3]*(0.4330127018922193*favg[9]+0.25*favg[1])+0.25*alphaDrag[1]*favg[6]-0.5*fjump[2]; 
  Ghat[3] = ((-2.381569860407206*nuVtSqSum[0]*fr[11])-2.381569860407206*nuVtSqSum[0]*fl[11]+1.875*nuVtSqSum[2]*fr[8]-1.875*nuVtSqSum[2]*fl[8]+1.875*nuVtSqSum[1]*fr[7]-1.875*nuVtSqSum[1]*fl[7]+1.875*nuVtSqSum[0]*fr[3]-1.875*nuVtSqSum[0]*fl[3])*rdv-0.8660254037844386*fjump[11]+alphaDrag[0]*(0.4330127018922193*favg[11]+0.25*favg[3])+0.25*alphaDrag[2]*favg[8]+0.25*alphaDrag[1]*favg[7]-0.5*fjump[3]; 
  Ghat[5] = ((-2.381569860407206*nuVtSqSum[0]*fr[15])-2.381569860407206*nuVtSqSum[0]*fl[15]+1.875*nuVtSqSum[2]*fr[13]-1.875*nuVtSqSum[2]*fl[13]+1.875*nuVtSqSum[1]*fr[12]-1.875*nuVtSqSum[1]*fl[12]+1.875*nuVtSqSum[0]*fr[5]-1.875*nuVtSqSum[0]*fl[5])*rdv-0.8660254037844386*fjump[15]+alphaDrag[0]*(0.4330127018922193*favg[15]+0.25*favg[5])+0.25*alphaDrag[2]*favg[13]+0.25*alphaDrag[1]*favg[12]-0.5*fjump[5]; 
  Ghat[6] = (1.341640786499874*nuVtSqSum[3]*fr[19]-1.341640786499874*nuVtSqSum[3]*fl[19]+1.677050983124842*nuVtSqSum[3]*fr[17]-1.677050983124842*nuVtSqSum[3]*fl[17]+1.677050983124842*nuVtSqSum[3]*fr[16]-1.677050983124842*nuVtSqSum[3]*fl[16]-2.381569860407206*nuVtSqSum[1]*fr[10]-2.381569860407206*nuVtSqSum[1]*fl[10]-2.381569860407206*nuVtSqSum[2]*fr[9]-2.381569860407206*nuVtSqSum[2]*fl[9]+1.677050983124842*nuVtSqSum[5]*fr[6]+1.677050983124842*nuVtSqSum[4]*fr[6]+1.875*nuVtSqSum[0]*fr[6]-1.677050983124842*nuVtSqSum[5]*fl[6]-1.677050983124842*nuVtSqSum[4]*fl[6]-1.875*nuVtSqSum[0]*fl[6]-2.381569860407206*nuVtSqSum[3]*fr[4]-2.381569860407206*nuVtSqSum[3]*fl[4]+1.875*fr[0]*nuVtSqSum[3]-1.875*fl[0]*nuVtSqSum[3]+1.875*fr[1]*nuVtSqSum[2]-1.875*fl[1]*nuVtSqSum[2]+1.875*nuVtSqSum[1]*fr[2]-1.875*nuVtSqSum[1]*fl[2])*rdv+alphaDrag[3]*(0.5590169943749475*favg[19]+0.223606797749979*favg[17]+0.223606797749979*favg[16]+0.4330127018922193*favg[4]+0.25*favg[0])+alphaDrag[1]*(0.4330127018922193*favg[10]+0.25*favg[2])+alphaDrag[2]*(0.4330127018922193*favg[9]+0.25*favg[1])-0.5*fjump[6]+0.223606797749979*alphaDrag[5]*favg[6]+0.223606797749979*alphaDrag[4]*favg[6]+0.25*alphaDrag[0]*favg[6]; 
  Ghat[7] = ((-2.381569860407206*nuVtSqSum[1]*fr[11])-2.381569860407206*nuVtSqSum[1]*fl[11]+1.875*nuVtSqSum[3]*fr[8]-1.875*nuVtSqSum[3]*fl[8]+1.677050983124842*nuVtSqSum[4]*fr[7]+1.875*nuVtSqSum[0]*fr[7]-1.677050983124842*nuVtSqSum[4]*fl[7]-1.875*nuVtSqSum[0]*fl[7]+1.875*nuVtSqSum[1]*fr[3]-1.875*nuVtSqSum[1]*fl[3])*rdv+alphaDrag[1]*(0.4330127018922193*favg[11]+0.25*favg[3])+0.25*alphaDrag[3]*favg[8]-0.5*fjump[7]+0.223606797749979*alphaDrag[4]*favg[7]+0.25*alphaDrag[0]*favg[7]; 
  Ghat[8] = ((-2.381569860407206*nuVtSqSum[2]*fr[11])-2.381569860407206*nuVtSqSum[2]*fl[11]+1.677050983124842*nuVtSqSum[5]*fr[8]+1.875*nuVtSqSum[0]*fr[8]-1.677050983124842*nuVtSqSum[5]*fl[8]-1.875*nuVtSqSum[0]*fl[8]+1.875*nuVtSqSum[3]*fr[7]-1.875*nuVtSqSum[3]*fl[7]+1.875*nuVtSqSum[2]*fr[3]-1.875*nuVtSqSum[2]*fl[3])*rdv+alphaDrag[2]*(0.4330127018922193*favg[11]+0.25*favg[3])-0.5*fjump[8]+0.223606797749979*alphaDrag[5]*favg[8]+0.25*alphaDrag[0]*favg[8]+0.25*alphaDrag[3]*favg[7]; 
  Ghat[12] = ((-2.381569860407206*nuVtSqSum[1]*fr[15])-2.381569860407206*nuVtSqSum[1]*fl[15]+1.875*nuVtSqSum[3]*fr[13]-1.875*nuVtSqSum[3]*fl[13]+1.677050983124842*nuVtSqSum[4]*fr[12]+1.875*nuVtSqSum[0]*fr[12]-1.677050983124842*nuVtSqSum[4]*fl[12]-1.875*nuVtSqSum[0]*fl[12]+1.875*nuVtSqSum[1]*fr[5]-1.875*nuVtSqSum[1]*fl[5])*rdv+alphaDrag[1]*(0.4330127018922193*favg[15]+0.25*favg[5])+0.25*alphaDrag[3]*favg[13]-0.5*fjump[12]+0.223606797749979*alphaDrag[4]*favg[12]+0.25*alphaDrag[0]*favg[12]; 
  Ghat[13] = ((-2.381569860407206*nuVtSqSum[2]*fr[15])-2.381569860407206*nuVtSqSum[2]*fl[15]+1.677050983124842*nuVtSqSum[5]*fr[13]+1.875*nuVtSqSum[0]*fr[13]-1.677050983124842*nuVtSqSum[5]*fl[13]-1.875*nuVtSqSum[0]*fl[13]+1.875*nuVtSqSum[3]*fr[12]-1.875*nuVtSqSum[3]*fl[12]+1.875*nuVtSqSum[2]*fr[5]-1.875*nuVtSqSum[2]*fl[5])*rdv+alphaDrag[2]*(0.4330127018922193*favg[15]+0.25*favg[5])-0.5*fjump[13]+0.223606797749979*alphaDrag[5]*favg[13]+0.25*alphaDrag[0]*favg[13]+0.25*alphaDrag[3]*favg[12]; 
  Ghat[14] = (1.875*nuVtSqSum[0]*fr[14]-1.875*nuVtSqSum[0]*fl[14])*rdv-0.5*fjump[14]+0.25*alphaDrag[0]*favg[14]; 
  Ghat[16] = (1.341640786499874*nuVtSqSum[4]*fr[19]-1.341640786499874*nuVtSqSum[4]*fl[19]+1.197893559374888*nuVtSqSum[4]*fr[16]+1.875*nuVtSqSum[0]*fr[16]-1.197893559374888*nuVtSqSum[4]*fl[16]-1.875*nuVtSqSum[0]*fl[16]-2.130140840414079*nuVtSqSum[1]*fr[9]-2.130140840414079*nuVtSqSum[1]*fl[9]+1.677050983124842*nuVtSqSum[3]*fr[6]-1.677050983124842*nuVtSqSum[3]*fl[6]-2.381569860407206*fr[4]*nuVtSqSum[4]-2.381569860407206*fl[4]*nuVtSqSum[4]+1.875*fr[0]*nuVtSqSum[4]-1.875*fl[0]*nuVtSqSum[4]+1.677050983124842*fr[1]*nuVtSqSum[1]-1.677050983124842*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[4]*(0.5590169943749475*favg[19]+0.159719141249985*favg[16]+0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[16]+0.25*alphaDrag[0]*favg[16]+alphaDrag[1]*(0.3872983346207416*favg[9]+0.223606797749979*favg[1])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[17] = (1.341640786499874*nuVtSqSum[5]*fr[19]-1.341640786499874*nuVtSqSum[5]*fl[19]+1.197893559374888*nuVtSqSum[5]*fr[17]+1.875*nuVtSqSum[0]*fr[17]-1.197893559374888*nuVtSqSum[5]*fl[17]-1.875*nuVtSqSum[0]*fl[17]-2.130140840414079*nuVtSqSum[2]*fr[10]-2.130140840414079*nuVtSqSum[2]*fl[10]+1.677050983124842*nuVtSqSum[3]*fr[6]-1.677050983124842*nuVtSqSum[3]*fl[6]-2.381569860407206*fr[4]*nuVtSqSum[5]-2.381569860407206*fl[4]*nuVtSqSum[5]+1.875*fr[0]*nuVtSqSum[5]-1.875*fl[0]*nuVtSqSum[5]+1.677050983124842*fr[2]*nuVtSqSum[2]-1.677050983124842*fl[2]*nuVtSqSum[2])*rdv+alphaDrag[5]*(0.5590169943749475*favg[19]+0.159719141249985*favg[17]+0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[17]+0.25*alphaDrag[0]*favg[17]+alphaDrag[2]*(0.3872983346207416*favg[10]+0.223606797749979*favg[2])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[18] = (1.875*nuVtSqSum[0]*fr[18]-1.875*nuVtSqSum[0]*fl[18])*rdv-0.5*fjump[18]+0.25*alphaDrag[0]*favg[18]; 
  Ghat[20] = (1.875*nuVtSqSum[0]*fr[20]-1.875*nuVtSqSum[0]*fl[20])*rdv-0.5*fjump[20]+0.25*alphaDrag[0]*favg[20]; 

  double incr1[21]; 
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
  incr1[17] = -0.5*Ghat[17]; 
  incr1[18] = -0.5*Ghat[18]; 
  incr1[19] = -1.118033988749895*Ghat[0]; 
  incr1[20] = -0.5*Ghat[20]; 

  double incr2[21]; 
  incr2[4] = nuVtSqSum[0]*(0.2118037767457181*(fr[19]+fl[19])-0.3046875*fr[4]+0.3046875*fl[4]+0.2165063509461096*(fr[0]+fl[0]))+0.2165063509461096*(nuVtSqSum[5]*(fr[17]+fl[17])+nuVtSqSum[4]*(fr[16]+fl[16]))+nuVtSqSum[2]*((-0.3046875*fr[10])+0.3046875*fl[10]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.3046875*fr[9])+0.3046875*fl[9]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*nuVtSqSum[3]*(fr[6]+fl[6]); 
  incr2[9] = nuVtSqSum[1]*(0.2118037767457181*(fr[19]+fl[19])+0.1936491673103708*(fr[16]+fl[16])-0.3046875*fr[4]+0.3046875*fl[4]+0.2165063509461096*(fr[0]+fl[0]))+nuVtSqSum[3]*((-0.3046875*fr[10])+0.3046875*fl[10]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[4]*((-0.2725207847577868*fr[9])+0.2725207847577868*fl[9]+0.1936491673103708*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.3046875*fr[9])+0.3046875*fl[9]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*nuVtSqSum[2]*(fr[6]+fl[6]); 
  incr2[10] = nuVtSqSum[2]*(0.2118037767457181*(fr[19]+fl[19])+0.1936491673103708*(fr[17]+fl[17])-0.3046875*fr[4]+0.3046875*fl[4]+0.2165063509461096*(fr[0]+fl[0]))+nuVtSqSum[5]*((-0.2725207847577868*fr[10])+0.2725207847577868*fl[10]+0.1936491673103708*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.3046875*fr[10])+0.3046875*fl[10]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[3]*((-0.3046875*fr[9])+0.3046875*fl[9]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*nuVtSqSum[1]*(fr[6]+fl[6]); 
  incr2[11] = nuVtSqSum[0]*((-0.3046875*fr[11])+0.3046875*fl[11]+0.2165063509461096*(fr[3]+fl[3]))+0.2165063509461096*(nuVtSqSum[2]*(fr[8]+fl[8])+nuVtSqSum[1]*(fr[7]+fl[7])); 
  incr2[15] = nuVtSqSum[0]*((-0.3046875*fr[15])+0.3046875*fl[15]+0.2165063509461096*(fr[5]+fl[5]))+0.2165063509461096*(nuVtSqSum[2]*(fr[13]+fl[13])+nuVtSqSum[1]*(fr[12]+fl[12])); 
  incr2[19] = nuVtSqSum[0]*((-0.8203125*(fr[19]+fl[19]))+1.180049613297572*fr[4]-1.180049613297572*fl[4]-0.8385254915624212*(fr[0]+fl[0]))-0.8385254915624212*(nuVtSqSum[5]*(fr[17]+fl[17])+nuVtSqSum[4]*(fr[16]+fl[16]))+nuVtSqSum[2]*(1.180049613297572*fr[10]-1.180049613297572*fl[10]-0.8385254915624212*(fr[2]+fl[2]))+nuVtSqSum[1]*(1.180049613297572*fr[9]-1.180049613297572*fl[9]-0.8385254915624212*(fr[1]+fl[1]))-0.8385254915624212*nuVtSqSum[3]*(fr[6]+fl[6]); 

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
  outr[17] += incr1[17]*rdv2R; 
  outr[18] += incr1[18]*rdv2R; 
  outr[19] += incr2[19]*rdvSq4R+incr1[19]*rdv2R; 
  outr[20] += incr1[20]*rdv2R; 

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
  outl[17] += -1.0*incr1[17]*rdv2L; 
  outl[18] += -1.0*incr1[18]*rdv2L; 
  outl[19] += incr2[19]*rdvSq4L-1.0*incr1[19]*rdv2L; 
  outl[20] += -1.0*incr1[20]*rdv2L; 

  return std::abs((0.5590169943749475*sumNuUy[5])/nuSum+(0.5590169943749475*sumNuUy[4])/nuSum-(0.5*sumNuUy[0])/nuSum+wl[3]); 
} 
double VmLBOconstNuSurf2x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:          Cell-center coordinates. 
  // dxv[5]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[18]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[6]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[4]; 
  double rdv2L = 2.0/dxvl[4]; 
  double rdv2R = 2.0/dxvr[4]; 
  double rdvSq4L = 4.0/(dxvl[4]*dxvl[4]); 
  double rdvSq4R = 4.0/(dxvr[4]*dxvr[4]); 

  const double *sumNuUz = &nuUSum[12]; 

  double favg[21]; 
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
  favg[20] = 1*fr[20]+fl[20]; 

  double fjump[21]; 
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
  fjump[20] = nuSum*vMuMidMax*(fl[20]-(1*fr[20])); 

  double alphaDrag[6]; 
  alphaDrag[0] = 2.0*wl[4]*nuSum+dxvl[4]*nuSum-1.0*sumNuUz[0]; 
  alphaDrag[1] = -1.0*sumNuUz[1]; 
  alphaDrag[2] = -1.0*sumNuUz[2]; 
  alphaDrag[3] = -1.0*sumNuUz[3]; 
  alphaDrag[4] = -1.0*sumNuUz[4]; 
  alphaDrag[5] = -1.0*sumNuUz[5]; 

  double Ghat[21]; 
  for(unsigned short int i=0; i<21; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.341640786499874*nuVtSqSum[0]*fr[20]-1.341640786499874*nuVtSqSum[0]*fl[20]+1.875*nuVtSqSum[5]*fr[17]-1.875*nuVtSqSum[5]*fl[17]+1.875*nuVtSqSum[4]*fr[16]-1.875*nuVtSqSum[4]*fl[16]-2.381569860407206*nuVtSqSum[2]*fr[13]-2.381569860407206*nuVtSqSum[2]*fl[13]-2.381569860407206*nuVtSqSum[1]*fr[12]-2.381569860407206*nuVtSqSum[1]*fl[12]+1.875*nuVtSqSum[3]*fr[6]-1.875*nuVtSqSum[3]*fl[6]-2.381569860407206*nuVtSqSum[0]*fr[5]-2.381569860407206*nuVtSqSum[0]*fl[5]+1.875*fr[2]*nuVtSqSum[2]-1.875*fl[2]*nuVtSqSum[2]+1.875*fr[1]*nuVtSqSum[1]-1.875*fl[1]*nuVtSqSum[1]+1.875*fr[0]*nuVtSqSum[0]-1.875*fl[0]*nuVtSqSum[0])*rdv-1.118033988749895*fjump[20]+alphaDrag[0]*(0.5590169943749475*favg[20]+0.4330127018922193*favg[5]+0.25*favg[0])+0.25*alphaDrag[5]*favg[17]+0.25*alphaDrag[4]*favg[16]+alphaDrag[2]*(0.4330127018922193*favg[13]+0.25*favg[2])+alphaDrag[1]*(0.4330127018922193*favg[12]+0.25*favg[1])+0.25*alphaDrag[3]*favg[6]-0.8660254037844386*fjump[5]-0.5*fjump[0]; 
  Ghat[1] = (1.341640786499874*nuVtSqSum[1]*fr[20]-1.341640786499874*nuVtSqSum[1]*fl[20]+1.677050983124842*nuVtSqSum[1]*fr[16]-1.677050983124842*nuVtSqSum[1]*fl[16]-2.381569860407206*nuVtSqSum[3]*fr[13]-2.381569860407206*nuVtSqSum[3]*fl[13]-2.130140840414079*nuVtSqSum[4]*fr[12]-2.381569860407206*nuVtSqSum[0]*fr[12]-2.130140840414079*nuVtSqSum[4]*fl[12]-2.381569860407206*nuVtSqSum[0]*fl[12]+1.875*nuVtSqSum[2]*fr[6]-1.875*nuVtSqSum[2]*fl[6]-2.381569860407206*nuVtSqSum[1]*fr[5]-2.381569860407206*nuVtSqSum[1]*fl[5]+1.677050983124842*fr[1]*nuVtSqSum[4]-1.677050983124842*fl[1]*nuVtSqSum[4]+1.875*fr[2]*nuVtSqSum[3]-1.875*fl[2]*nuVtSqSum[3]+1.875*fr[0]*nuVtSqSum[1]-1.875*fl[0]*nuVtSqSum[1]+1.875*nuVtSqSum[0]*fr[1]-1.875*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.5590169943749475*favg[20]+0.223606797749979*favg[16]+0.4330127018922193*favg[5]+0.25*favg[0])+alphaDrag[3]*(0.4330127018922193*favg[13]+0.25*favg[2])-0.8660254037844386*fjump[12]+alphaDrag[0]*(0.4330127018922193*favg[12]+0.25*favg[1])+alphaDrag[4]*(0.3872983346207416*favg[12]+0.223606797749979*favg[1])+0.25*alphaDrag[2]*favg[6]-0.5*fjump[1]; 
  Ghat[2] = (1.341640786499874*nuVtSqSum[2]*fr[20]-1.341640786499874*nuVtSqSum[2]*fl[20]+1.677050983124842*nuVtSqSum[2]*fr[17]-1.677050983124842*nuVtSqSum[2]*fl[17]-2.130140840414079*nuVtSqSum[5]*fr[13]-2.381569860407206*nuVtSqSum[0]*fr[13]-2.130140840414079*nuVtSqSum[5]*fl[13]-2.381569860407206*nuVtSqSum[0]*fl[13]-2.381569860407206*nuVtSqSum[3]*fr[12]-2.381569860407206*nuVtSqSum[3]*fl[12]+1.875*nuVtSqSum[1]*fr[6]-1.875*nuVtSqSum[1]*fl[6]+1.677050983124842*fr[2]*nuVtSqSum[5]-1.677050983124842*fl[2]*nuVtSqSum[5]-2.381569860407206*nuVtSqSum[2]*fr[5]-2.381569860407206*nuVtSqSum[2]*fl[5]+1.875*fr[1]*nuVtSqSum[3]-1.875*fl[1]*nuVtSqSum[3]+1.875*fr[0]*nuVtSqSum[2]-1.875*fl[0]*nuVtSqSum[2]+1.875*nuVtSqSum[0]*fr[2]-1.875*nuVtSqSum[0]*fl[2])*rdv+alphaDrag[2]*(0.5590169943749475*favg[20]+0.223606797749979*favg[17]+0.4330127018922193*favg[5]+0.25*favg[0])-0.8660254037844386*fjump[13]+alphaDrag[0]*(0.4330127018922193*favg[13]+0.25*favg[2])+alphaDrag[5]*(0.3872983346207416*favg[13]+0.223606797749979*favg[2])+alphaDrag[3]*(0.4330127018922193*favg[12]+0.25*favg[1])+0.25*alphaDrag[1]*favg[6]-0.5*fjump[2]; 
  Ghat[3] = ((-2.381569860407206*nuVtSqSum[0]*fr[14])-2.381569860407206*nuVtSqSum[0]*fl[14]+1.875*nuVtSqSum[2]*fr[8]-1.875*nuVtSqSum[2]*fl[8]+1.875*nuVtSqSum[1]*fr[7]-1.875*nuVtSqSum[1]*fl[7]+1.875*nuVtSqSum[0]*fr[3]-1.875*nuVtSqSum[0]*fl[3])*rdv-0.8660254037844386*fjump[14]+alphaDrag[0]*(0.4330127018922193*favg[14]+0.25*favg[3])+0.25*alphaDrag[2]*favg[8]+0.25*alphaDrag[1]*favg[7]-0.5*fjump[3]; 
  Ghat[4] = ((-2.381569860407206*nuVtSqSum[0]*fr[15])-2.381569860407206*nuVtSqSum[0]*fl[15]+1.875*nuVtSqSum[2]*fr[10]-1.875*nuVtSqSum[2]*fl[10]+1.875*nuVtSqSum[1]*fr[9]-1.875*nuVtSqSum[1]*fl[9]+1.875*nuVtSqSum[0]*fr[4]-1.875*nuVtSqSum[0]*fl[4])*rdv-0.8660254037844386*fjump[15]+alphaDrag[0]*(0.4330127018922193*favg[15]+0.25*favg[4])+0.25*alphaDrag[2]*favg[10]+0.25*alphaDrag[1]*favg[9]-0.5*fjump[4]; 
  Ghat[6] = (1.341640786499874*nuVtSqSum[3]*fr[20]-1.341640786499874*nuVtSqSum[3]*fl[20]+1.677050983124842*nuVtSqSum[3]*fr[17]-1.677050983124842*nuVtSqSum[3]*fl[17]+1.677050983124842*nuVtSqSum[3]*fr[16]-1.677050983124842*nuVtSqSum[3]*fl[16]-2.381569860407206*nuVtSqSum[1]*fr[13]-2.381569860407206*nuVtSqSum[1]*fl[13]-2.381569860407206*nuVtSqSum[2]*fr[12]-2.381569860407206*nuVtSqSum[2]*fl[12]+1.677050983124842*nuVtSqSum[5]*fr[6]+1.677050983124842*nuVtSqSum[4]*fr[6]+1.875*nuVtSqSum[0]*fr[6]-1.677050983124842*nuVtSqSum[5]*fl[6]-1.677050983124842*nuVtSqSum[4]*fl[6]-1.875*nuVtSqSum[0]*fl[6]-2.381569860407206*nuVtSqSum[3]*fr[5]-2.381569860407206*nuVtSqSum[3]*fl[5]+1.875*fr[0]*nuVtSqSum[3]-1.875*fl[0]*nuVtSqSum[3]+1.875*fr[1]*nuVtSqSum[2]-1.875*fl[1]*nuVtSqSum[2]+1.875*nuVtSqSum[1]*fr[2]-1.875*nuVtSqSum[1]*fl[2])*rdv+alphaDrag[3]*(0.5590169943749475*favg[20]+0.223606797749979*favg[17]+0.223606797749979*favg[16]+0.4330127018922193*favg[5]+0.25*favg[0])+alphaDrag[1]*(0.4330127018922193*favg[13]+0.25*favg[2])+alphaDrag[2]*(0.4330127018922193*favg[12]+0.25*favg[1])-0.5*fjump[6]+0.223606797749979*alphaDrag[5]*favg[6]+0.223606797749979*alphaDrag[4]*favg[6]+0.25*alphaDrag[0]*favg[6]; 
  Ghat[7] = ((-2.381569860407206*nuVtSqSum[1]*fr[14])-2.381569860407206*nuVtSqSum[1]*fl[14]+1.875*nuVtSqSum[3]*fr[8]-1.875*nuVtSqSum[3]*fl[8]+1.677050983124842*nuVtSqSum[4]*fr[7]+1.875*nuVtSqSum[0]*fr[7]-1.677050983124842*nuVtSqSum[4]*fl[7]-1.875*nuVtSqSum[0]*fl[7]+1.875*nuVtSqSum[1]*fr[3]-1.875*nuVtSqSum[1]*fl[3])*rdv+alphaDrag[1]*(0.4330127018922193*favg[14]+0.25*favg[3])+0.25*alphaDrag[3]*favg[8]-0.5*fjump[7]+0.223606797749979*alphaDrag[4]*favg[7]+0.25*alphaDrag[0]*favg[7]; 
  Ghat[8] = ((-2.381569860407206*nuVtSqSum[2]*fr[14])-2.381569860407206*nuVtSqSum[2]*fl[14]+1.677050983124842*nuVtSqSum[5]*fr[8]+1.875*nuVtSqSum[0]*fr[8]-1.677050983124842*nuVtSqSum[5]*fl[8]-1.875*nuVtSqSum[0]*fl[8]+1.875*nuVtSqSum[3]*fr[7]-1.875*nuVtSqSum[3]*fl[7]+1.875*nuVtSqSum[2]*fr[3]-1.875*nuVtSqSum[2]*fl[3])*rdv+alphaDrag[2]*(0.4330127018922193*favg[14]+0.25*favg[3])-0.5*fjump[8]+0.223606797749979*alphaDrag[5]*favg[8]+0.25*alphaDrag[0]*favg[8]+0.25*alphaDrag[3]*favg[7]; 
  Ghat[9] = ((-2.381569860407206*nuVtSqSum[1]*fr[15])-2.381569860407206*nuVtSqSum[1]*fl[15]+1.875*nuVtSqSum[3]*fr[10]-1.875*nuVtSqSum[3]*fl[10]+1.677050983124842*nuVtSqSum[4]*fr[9]+1.875*nuVtSqSum[0]*fr[9]-1.677050983124842*nuVtSqSum[4]*fl[9]-1.875*nuVtSqSum[0]*fl[9]+1.875*nuVtSqSum[1]*fr[4]-1.875*nuVtSqSum[1]*fl[4])*rdv+alphaDrag[1]*(0.4330127018922193*favg[15]+0.25*favg[4])+0.25*alphaDrag[3]*favg[10]-0.5*fjump[9]+0.223606797749979*alphaDrag[4]*favg[9]+0.25*alphaDrag[0]*favg[9]; 
  Ghat[10] = ((-2.381569860407206*nuVtSqSum[2]*fr[15])-2.381569860407206*nuVtSqSum[2]*fl[15]+1.677050983124842*nuVtSqSum[5]*fr[10]+1.875*nuVtSqSum[0]*fr[10]-1.677050983124842*nuVtSqSum[5]*fl[10]-1.875*nuVtSqSum[0]*fl[10]+1.875*nuVtSqSum[3]*fr[9]-1.875*nuVtSqSum[3]*fl[9]+1.875*nuVtSqSum[2]*fr[4]-1.875*nuVtSqSum[2]*fl[4])*rdv+alphaDrag[2]*(0.4330127018922193*favg[15]+0.25*favg[4])-0.5*fjump[10]+0.223606797749979*alphaDrag[5]*favg[10]+0.25*alphaDrag[0]*favg[10]+0.25*alphaDrag[3]*favg[9]; 
  Ghat[11] = (1.875*nuVtSqSum[0]*fr[11]-1.875*nuVtSqSum[0]*fl[11])*rdv-0.5*fjump[11]+0.25*alphaDrag[0]*favg[11]; 
  Ghat[16] = (1.341640786499874*nuVtSqSum[4]*fr[20]-1.341640786499874*nuVtSqSum[4]*fl[20]+1.197893559374888*nuVtSqSum[4]*fr[16]+1.875*nuVtSqSum[0]*fr[16]-1.197893559374888*nuVtSqSum[4]*fl[16]-1.875*nuVtSqSum[0]*fl[16]-2.130140840414079*nuVtSqSum[1]*fr[12]-2.130140840414079*nuVtSqSum[1]*fl[12]+1.677050983124842*nuVtSqSum[3]*fr[6]-1.677050983124842*nuVtSqSum[3]*fl[6]-2.381569860407206*nuVtSqSum[4]*fr[5]-2.381569860407206*nuVtSqSum[4]*fl[5]+1.875*fr[0]*nuVtSqSum[4]-1.875*fl[0]*nuVtSqSum[4]+1.677050983124842*fr[1]*nuVtSqSum[1]-1.677050983124842*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[4]*(0.5590169943749475*favg[20]+0.159719141249985*favg[16]+0.4330127018922193*favg[5]+0.25*favg[0])-0.5*fjump[16]+0.25*alphaDrag[0]*favg[16]+alphaDrag[1]*(0.3872983346207416*favg[12]+0.223606797749979*favg[1])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[17] = (1.341640786499874*nuVtSqSum[5]*fr[20]-1.341640786499874*nuVtSqSum[5]*fl[20]+1.197893559374888*nuVtSqSum[5]*fr[17]+1.875*nuVtSqSum[0]*fr[17]-1.197893559374888*nuVtSqSum[5]*fl[17]-1.875*nuVtSqSum[0]*fl[17]-2.130140840414079*nuVtSqSum[2]*fr[13]-2.130140840414079*nuVtSqSum[2]*fl[13]+1.677050983124842*nuVtSqSum[3]*fr[6]-1.677050983124842*nuVtSqSum[3]*fl[6]-2.381569860407206*fr[5]*nuVtSqSum[5]-2.381569860407206*fl[5]*nuVtSqSum[5]+1.875*fr[0]*nuVtSqSum[5]-1.875*fl[0]*nuVtSqSum[5]+1.677050983124842*fr[2]*nuVtSqSum[2]-1.677050983124842*fl[2]*nuVtSqSum[2])*rdv+alphaDrag[5]*(0.5590169943749475*favg[20]+0.159719141249985*favg[17]+0.4330127018922193*favg[5]+0.25*favg[0])-0.5*fjump[17]+0.25*alphaDrag[0]*favg[17]+alphaDrag[2]*(0.3872983346207416*favg[13]+0.223606797749979*favg[2])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[18] = (1.875*nuVtSqSum[0]*fr[18]-1.875*nuVtSqSum[0]*fl[18])*rdv-0.5*fjump[18]+0.25*alphaDrag[0]*favg[18]; 
  Ghat[19] = (1.875*nuVtSqSum[0]*fr[19]-1.875*nuVtSqSum[0]*fl[19])*rdv-0.5*fjump[19]+0.25*alphaDrag[0]*favg[19]; 

  double incr1[21]; 
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
  incr1[20] = -1.118033988749895*Ghat[0]; 

  double incr2[21]; 
  incr2[5] = nuVtSqSum[0]*(0.2118037767457181*(fr[20]+fl[20])-0.3046875*fr[5]+0.3046875*fl[5]+0.2165063509461096*(fr[0]+fl[0]))+0.2165063509461096*(nuVtSqSum[5]*(fr[17]+fl[17])+nuVtSqSum[4]*(fr[16]+fl[16]))+nuVtSqSum[2]*((-0.3046875*fr[13])+0.3046875*fl[13]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.3046875*fr[12])+0.3046875*fl[12]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*nuVtSqSum[3]*(fr[6]+fl[6]); 
  incr2[12] = nuVtSqSum[1]*(0.2118037767457181*(fr[20]+fl[20])+0.1936491673103708*(fr[16]+fl[16])-0.3046875*fr[5]+0.3046875*fl[5]+0.2165063509461096*(fr[0]+fl[0]))+nuVtSqSum[3]*((-0.3046875*fr[13])+0.3046875*fl[13]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[4]*((-0.2725207847577868*fr[12])+0.2725207847577868*fl[12]+0.1936491673103708*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.3046875*fr[12])+0.3046875*fl[12]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*nuVtSqSum[2]*(fr[6]+fl[6]); 
  incr2[13] = nuVtSqSum[2]*(0.2118037767457181*(fr[20]+fl[20])+0.1936491673103708*(fr[17]+fl[17])-0.3046875*fr[5]+0.3046875*fl[5]+0.2165063509461096*(fr[0]+fl[0]))+nuVtSqSum[5]*((-0.2725207847577868*fr[13])+0.2725207847577868*fl[13]+0.1936491673103708*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.3046875*fr[13])+0.3046875*fl[13]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[3]*((-0.3046875*fr[12])+0.3046875*fl[12]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*nuVtSqSum[1]*(fr[6]+fl[6]); 
  incr2[14] = nuVtSqSum[0]*((-0.3046875*fr[14])+0.3046875*fl[14]+0.2165063509461096*(fr[3]+fl[3]))+0.2165063509461096*(nuVtSqSum[2]*(fr[8]+fl[8])+nuVtSqSum[1]*(fr[7]+fl[7])); 
  incr2[15] = nuVtSqSum[0]*((-0.3046875*fr[15])+0.3046875*fl[15]+0.2165063509461096*(fr[4]+fl[4]))+0.2165063509461096*(nuVtSqSum[2]*(fr[10]+fl[10])+nuVtSqSum[1]*(fr[9]+fl[9])); 
  incr2[20] = nuVtSqSum[0]*((-0.8203125*(fr[20]+fl[20]))+1.180049613297572*fr[5]-1.180049613297572*fl[5]-0.8385254915624212*(fr[0]+fl[0]))-0.8385254915624212*(nuVtSqSum[5]*(fr[17]+fl[17])+nuVtSqSum[4]*(fr[16]+fl[16]))+nuVtSqSum[2]*(1.180049613297572*fr[13]-1.180049613297572*fl[13]-0.8385254915624212*(fr[2]+fl[2]))+nuVtSqSum[1]*(1.180049613297572*fr[12]-1.180049613297572*fl[12]-0.8385254915624212*(fr[1]+fl[1]))-0.8385254915624212*nuVtSqSum[3]*(fr[6]+fl[6]); 

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
  outl[20] += incr2[20]*rdvSq4L-1.0*incr1[20]*rdv2L; 

  return std::abs((0.5590169943749475*sumNuUz[5])/nuSum+(0.5590169943749475*sumNuUz[4])/nuSum-(0.5*sumNuUz[0])/nuSum+wl[4]); 
} 
