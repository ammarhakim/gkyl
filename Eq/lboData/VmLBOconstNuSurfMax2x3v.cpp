#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf2x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates. 
  // dxv[5]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[9]:       bulk velocity (in 3 directions). 
  // vtSq[3]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2nuL = 2.0*nu/dxvl[2]; 
  double rdv2nuR = 2.0*nu/dxvr[2]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[2]*dxvl[2]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[2]*dxvr[2]); 

  const double *Ux = &u[0]; 

  double favg[6]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 

  double fjump[6]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(1*fr[5])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 2.0*wl[2]+dxvl[2]-1.0*Ux[0]; 
  alphaDrag[1] = -1.0*Ux[1]; 
  alphaDrag[2] = -1.0*Ux[2]; 

  double Ghat[6]; 
  for(unsigned short int i=0; i<6; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.082531754730548*vtSq[0]*fr[3])-1.082531754730548*vtSq[0]*fl[3]+1.125*fr[2]*vtSq[2]-1.125*fl[2]*vtSq[2]+1.125*fr[1]*vtSq[1]-1.125*fl[1]*vtSq[1]+1.125*fr[0]*vtSq[0]-1.125*fl[0]*vtSq[0])*rdv-0.8660254037844386*fjump[3]+alphaDrag[0]*(0.4330127018922193*favg[3]+0.25*favg[0])+0.25*alphaDrag[2]*favg[2]+0.25*alphaDrag[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = ((-1.082531754730548*vtSq[1]*fr[3])-1.082531754730548*vtSq[1]*fl[3]+1.125*fr[0]*vtSq[1]-1.125*fl[0]*vtSq[1]+1.125*vtSq[0]*fr[1]-1.125*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[1]+0.25*alphaDrag[0]*favg[1]; 
  Ghat[2] = ((-1.082531754730548*vtSq[2]*fr[3])-1.082531754730548*vtSq[2]*fl[3]+1.125*fr[0]*vtSq[2]-1.125*fl[0]*vtSq[2]+1.125*vtSq[0]*fr[2]-1.125*vtSq[0]*fl[2])*rdv+alphaDrag[2]*(0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[2]+0.25*alphaDrag[0]*favg[2]; 
  Ghat[4] = (1.125*vtSq[0]*fr[4]-1.125*vtSq[0]*fl[4])*rdv-0.5*fjump[4]+0.25*alphaDrag[0]*favg[4]; 
  Ghat[5] = (1.125*vtSq[0]*fr[5]-1.125*vtSq[0]*fl[5])*rdv-0.5*fjump[5]+0.25*alphaDrag[0]*favg[5]; 

  double incr1[6]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -0.5*Ghat[5]; 

  double incr2[6]; 
  incr2[3] = vtSq[0]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0]))+0.2165063509461096*((fr[2]+fl[2])*vtSq[2]+(fr[1]+fl[1])*vtSq[1]); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr1[2]*rdv2nuR; 
  outr[3] += incr2[3]*rdvSq4nuR+incr1[3]*rdv2nuR; 
  outr[4] += incr1[4]*rdv2nuR; 
  outr[5] += incr1[5]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += -1.0*incr1[2]*rdv2nuL; 
  outl[3] += incr1[3]*rdv2nuL-1.0*incr2[3]*rdvSq4nuL; 
  outl[4] += -1.0*incr1[4]*rdv2nuL; 
  outl[5] += -1.0*incr1[5]*rdv2nuL; 

  return std::abs(wl[2]-0.5*Ux[0]); 
} 
double VmLBOconstNuSurf2x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates. 
  // dxv[5]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[18]:       bulk velocity (in 3 directions). 
  // vtSq[6]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2nuL = 2.0*nu/dxvl[2]; 
  double rdv2nuR = 2.0*nu/dxvr[2]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[2]*dxvl[2]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[2]*dxvr[2]); 

  const double *Ux = &u[0]; 

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
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = vMuMidMax*(fl[7]-(-1*fr[7])); 
  fjump[8] = vMuMidMax*(fl[8]-(-1*fr[8])); 
  fjump[9] = vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = vMuMidMax*(fl[10]-(1*fr[10])); 
  fjump[11] = vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = vMuMidMax*(fl[16]-(1*fr[16])); 
  fjump[17] = vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = vMuMidMax*(fl[19]-(1*fr[19])); 
  fjump[20] = vMuMidMax*(fl[20]-(1*fr[20])); 

  double alphaDrag[6]; 
  alphaDrag[0] = 2.0*wl[2]+dxvl[2]-1.0*Ux[0]; 
  alphaDrag[1] = -1.0*Ux[1]; 
  alphaDrag[2] = -1.0*Ux[2]; 
  alphaDrag[3] = -1.0*Ux[3]; 
  alphaDrag[4] = -1.0*Ux[4]; 
  alphaDrag[5] = -1.0*Ux[5]; 

  double Ghat[21]; 
  for(unsigned short int i=0; i<21; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.341640786499874*vtSq[0]*fr[18]-1.341640786499874*vtSq[0]*fl[18]+1.875*vtSq[5]*fr[17]-1.875*vtSq[5]*fl[17]+1.875*vtSq[4]*fr[16]-1.875*vtSq[4]*fl[16]-2.381569860407206*vtSq[2]*fr[8]-2.381569860407206*vtSq[2]*fl[8]-2.381569860407206*vtSq[1]*fr[7]-2.381569860407206*vtSq[1]*fl[7]+1.875*vtSq[3]*fr[6]-1.875*vtSq[3]*fl[6]-2.381569860407206*vtSq[0]*fr[3]-2.381569860407206*vtSq[0]*fl[3]+1.875*fr[2]*vtSq[2]-1.875*fl[2]*vtSq[2]+1.875*fr[1]*vtSq[1]-1.875*fl[1]*vtSq[1]+1.875*fr[0]*vtSq[0]-1.875*fl[0]*vtSq[0])*rdv-1.118033988749895*fjump[18]+alphaDrag[0]*(0.5590169943749475*favg[18]+0.4330127018922193*favg[3]+0.25*favg[0])+0.25*alphaDrag[5]*favg[17]+0.25*alphaDrag[4]*favg[16]+alphaDrag[2]*(0.4330127018922193*favg[8]+0.25*favg[2])+alphaDrag[1]*(0.4330127018922193*favg[7]+0.25*favg[1])+0.25*alphaDrag[3]*favg[6]-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = (1.341640786499874*vtSq[1]*fr[18]-1.341640786499874*vtSq[1]*fl[18]+1.677050983124842*vtSq[1]*fr[16]-1.677050983124842*vtSq[1]*fl[16]-2.381569860407206*vtSq[3]*fr[8]-2.381569860407206*vtSq[3]*fl[8]-2.130140840414079*vtSq[4]*fr[7]-2.381569860407206*vtSq[0]*fr[7]-2.130140840414079*vtSq[4]*fl[7]-2.381569860407206*vtSq[0]*fl[7]+1.875*vtSq[2]*fr[6]-1.875*vtSq[2]*fl[6]+1.677050983124842*fr[1]*vtSq[4]-1.677050983124842*fl[1]*vtSq[4]+1.875*fr[2]*vtSq[3]-1.875*fl[2]*vtSq[3]-2.381569860407206*vtSq[1]*fr[3]-2.381569860407206*vtSq[1]*fl[3]+1.875*fr[0]*vtSq[1]-1.875*fl[0]*vtSq[1]+1.875*vtSq[0]*fr[1]-1.875*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.5590169943749475*favg[18]+0.223606797749979*favg[16]+0.4330127018922193*favg[3]+0.25*favg[0])+alphaDrag[3]*(0.4330127018922193*favg[8]+0.25*favg[2])-0.8660254037844386*fjump[7]+alphaDrag[0]*(0.4330127018922193*favg[7]+0.25*favg[1])+alphaDrag[4]*(0.3872983346207416*favg[7]+0.223606797749979*favg[1])+0.25*alphaDrag[2]*favg[6]-0.5*fjump[1]; 
  Ghat[2] = (1.341640786499874*vtSq[2]*fr[18]-1.341640786499874*vtSq[2]*fl[18]+1.677050983124842*vtSq[2]*fr[17]-1.677050983124842*vtSq[2]*fl[17]-2.130140840414079*vtSq[5]*fr[8]-2.381569860407206*vtSq[0]*fr[8]-2.130140840414079*vtSq[5]*fl[8]-2.381569860407206*vtSq[0]*fl[8]-2.381569860407206*vtSq[3]*fr[7]-2.381569860407206*vtSq[3]*fl[7]+1.875*vtSq[1]*fr[6]-1.875*vtSq[1]*fl[6]+1.677050983124842*fr[2]*vtSq[5]-1.677050983124842*fl[2]*vtSq[5]+1.875*fr[1]*vtSq[3]-1.875*fl[1]*vtSq[3]-2.381569860407206*vtSq[2]*fr[3]-2.381569860407206*vtSq[2]*fl[3]+1.875*fr[0]*vtSq[2]-1.875*fl[0]*vtSq[2]+1.875*vtSq[0]*fr[2]-1.875*vtSq[0]*fl[2])*rdv+alphaDrag[2]*(0.5590169943749475*favg[18]+0.223606797749979*favg[17]+0.4330127018922193*favg[3]+0.25*favg[0])-0.8660254037844386*fjump[8]+alphaDrag[0]*(0.4330127018922193*favg[8]+0.25*favg[2])+alphaDrag[5]*(0.3872983346207416*favg[8]+0.223606797749979*favg[2])+alphaDrag[3]*(0.4330127018922193*favg[7]+0.25*favg[1])+0.25*alphaDrag[1]*favg[6]-0.5*fjump[2]; 
  Ghat[4] = ((-2.381569860407206*vtSq[0]*fr[11])-2.381569860407206*vtSq[0]*fl[11]+1.875*vtSq[2]*fr[10]-1.875*vtSq[2]*fl[10]+1.875*vtSq[1]*fr[9]-1.875*vtSq[1]*fl[9]+1.875*vtSq[0]*fr[4]-1.875*vtSq[0]*fl[4])*rdv-0.8660254037844386*fjump[11]+alphaDrag[0]*(0.4330127018922193*favg[11]+0.25*favg[4])+0.25*alphaDrag[2]*favg[10]+0.25*alphaDrag[1]*favg[9]-0.5*fjump[4]; 
  Ghat[5] = ((-2.381569860407206*vtSq[0]*fr[14])-2.381569860407206*vtSq[0]*fl[14]+1.875*vtSq[2]*fr[13]-1.875*vtSq[2]*fl[13]+1.875*vtSq[1]*fr[12]-1.875*vtSq[1]*fl[12]+1.875*vtSq[0]*fr[5]-1.875*vtSq[0]*fl[5])*rdv-0.8660254037844386*fjump[14]+alphaDrag[0]*(0.4330127018922193*favg[14]+0.25*favg[5])+0.25*alphaDrag[2]*favg[13]+0.25*alphaDrag[1]*favg[12]-0.5*fjump[5]; 
  Ghat[6] = (1.341640786499874*vtSq[3]*fr[18]-1.341640786499874*vtSq[3]*fl[18]+1.677050983124842*vtSq[3]*fr[17]-1.677050983124842*vtSq[3]*fl[17]+1.677050983124842*vtSq[3]*fr[16]-1.677050983124842*vtSq[3]*fl[16]-2.381569860407206*vtSq[1]*fr[8]-2.381569860407206*vtSq[1]*fl[8]-2.381569860407206*vtSq[2]*fr[7]-2.381569860407206*vtSq[2]*fl[7]+1.677050983124842*vtSq[5]*fr[6]+1.677050983124842*vtSq[4]*fr[6]+1.875*vtSq[0]*fr[6]-1.677050983124842*vtSq[5]*fl[6]-1.677050983124842*vtSq[4]*fl[6]-1.875*vtSq[0]*fl[6]-2.381569860407206*fr[3]*vtSq[3]-2.381569860407206*fl[3]*vtSq[3]+1.875*fr[0]*vtSq[3]-1.875*fl[0]*vtSq[3]+1.875*fr[1]*vtSq[2]-1.875*fl[1]*vtSq[2]+1.875*vtSq[1]*fr[2]-1.875*vtSq[1]*fl[2])*rdv+alphaDrag[3]*(0.5590169943749475*favg[18]+0.223606797749979*favg[17]+0.223606797749979*favg[16]+0.4330127018922193*favg[3]+0.25*favg[0])+alphaDrag[1]*(0.4330127018922193*favg[8]+0.25*favg[2])+alphaDrag[2]*(0.4330127018922193*favg[7]+0.25*favg[1])-0.5*fjump[6]+0.223606797749979*alphaDrag[5]*favg[6]+0.223606797749979*alphaDrag[4]*favg[6]+0.25*alphaDrag[0]*favg[6]; 
  Ghat[9] = ((-2.381569860407206*vtSq[1]*fr[11])-2.381569860407206*vtSq[1]*fl[11]+1.875*vtSq[3]*fr[10]-1.875*vtSq[3]*fl[10]+1.677050983124842*vtSq[4]*fr[9]+1.875*vtSq[0]*fr[9]-1.677050983124842*vtSq[4]*fl[9]-1.875*vtSq[0]*fl[9]+1.875*vtSq[1]*fr[4]-1.875*vtSq[1]*fl[4])*rdv+alphaDrag[1]*(0.4330127018922193*favg[11]+0.25*favg[4])+0.25*alphaDrag[3]*favg[10]-0.5*fjump[9]+0.223606797749979*alphaDrag[4]*favg[9]+0.25*alphaDrag[0]*favg[9]; 
  Ghat[10] = ((-2.381569860407206*vtSq[2]*fr[11])-2.381569860407206*vtSq[2]*fl[11]+1.677050983124842*vtSq[5]*fr[10]+1.875*vtSq[0]*fr[10]-1.677050983124842*vtSq[5]*fl[10]-1.875*vtSq[0]*fl[10]+1.875*vtSq[3]*fr[9]-1.875*vtSq[3]*fl[9]+1.875*vtSq[2]*fr[4]-1.875*vtSq[2]*fl[4])*rdv+alphaDrag[2]*(0.4330127018922193*favg[11]+0.25*favg[4])-0.5*fjump[10]+0.223606797749979*alphaDrag[5]*favg[10]+0.25*alphaDrag[0]*favg[10]+0.25*alphaDrag[3]*favg[9]; 
  Ghat[12] = ((-2.381569860407206*vtSq[1]*fr[14])-2.381569860407206*vtSq[1]*fl[14]+1.875*vtSq[3]*fr[13]-1.875*vtSq[3]*fl[13]+1.677050983124842*vtSq[4]*fr[12]+1.875*vtSq[0]*fr[12]-1.677050983124842*vtSq[4]*fl[12]-1.875*vtSq[0]*fl[12]+1.875*vtSq[1]*fr[5]-1.875*vtSq[1]*fl[5])*rdv+alphaDrag[1]*(0.4330127018922193*favg[14]+0.25*favg[5])+0.25*alphaDrag[3]*favg[13]-0.5*fjump[12]+0.223606797749979*alphaDrag[4]*favg[12]+0.25*alphaDrag[0]*favg[12]; 
  Ghat[13] = ((-2.381569860407206*vtSq[2]*fr[14])-2.381569860407206*vtSq[2]*fl[14]+1.677050983124842*vtSq[5]*fr[13]+1.875*vtSq[0]*fr[13]-1.677050983124842*vtSq[5]*fl[13]-1.875*vtSq[0]*fl[13]+1.875*vtSq[3]*fr[12]-1.875*vtSq[3]*fl[12]+1.875*vtSq[2]*fr[5]-1.875*vtSq[2]*fl[5])*rdv+alphaDrag[2]*(0.4330127018922193*favg[14]+0.25*favg[5])-0.5*fjump[13]+0.223606797749979*alphaDrag[5]*favg[13]+0.25*alphaDrag[0]*favg[13]+0.25*alphaDrag[3]*favg[12]; 
  Ghat[15] = (1.875*vtSq[0]*fr[15]-1.875*vtSq[0]*fl[15])*rdv-0.5*fjump[15]+0.25*alphaDrag[0]*favg[15]; 
  Ghat[16] = (1.341640786499874*vtSq[4]*fr[18]-1.341640786499874*vtSq[4]*fl[18]+1.197893559374888*vtSq[4]*fr[16]+1.875*vtSq[0]*fr[16]-1.197893559374888*vtSq[4]*fl[16]-1.875*vtSq[0]*fl[16]-2.130140840414079*vtSq[1]*fr[7]-2.130140840414079*vtSq[1]*fl[7]+1.677050983124842*vtSq[3]*fr[6]-1.677050983124842*vtSq[3]*fl[6]-2.381569860407206*fr[3]*vtSq[4]-2.381569860407206*fl[3]*vtSq[4]+1.875*fr[0]*vtSq[4]-1.875*fl[0]*vtSq[4]+1.677050983124842*fr[1]*vtSq[1]-1.677050983124842*fl[1]*vtSq[1])*rdv+alphaDrag[4]*(0.5590169943749475*favg[18]+0.159719141249985*favg[16]+0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[16]+0.25*alphaDrag[0]*favg[16]+alphaDrag[1]*(0.3872983346207416*favg[7]+0.223606797749979*favg[1])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[17] = (1.341640786499874*vtSq[5]*fr[18]-1.341640786499874*vtSq[5]*fl[18]+1.197893559374888*vtSq[5]*fr[17]+1.875*vtSq[0]*fr[17]-1.197893559374888*vtSq[5]*fl[17]-1.875*vtSq[0]*fl[17]-2.130140840414079*vtSq[2]*fr[8]-2.130140840414079*vtSq[2]*fl[8]+1.677050983124842*vtSq[3]*fr[6]-1.677050983124842*vtSq[3]*fl[6]-2.381569860407206*fr[3]*vtSq[5]-2.381569860407206*fl[3]*vtSq[5]+1.875*fr[0]*vtSq[5]-1.875*fl[0]*vtSq[5]+1.677050983124842*fr[2]*vtSq[2]-1.677050983124842*fl[2]*vtSq[2])*rdv+alphaDrag[5]*(0.5590169943749475*favg[18]+0.159719141249985*favg[17]+0.4330127018922193*favg[3]+0.25*favg[0])-0.5*fjump[17]+0.25*alphaDrag[0]*favg[17]+alphaDrag[2]*(0.3872983346207416*favg[8]+0.223606797749979*favg[2])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[19] = (1.875*vtSq[0]*fr[19]-1.875*vtSq[0]*fl[19])*rdv-0.5*fjump[19]+0.25*alphaDrag[0]*favg[19]; 
  Ghat[20] = (1.875*vtSq[0]*fr[20]-1.875*vtSq[0]*fl[20])*rdv-0.5*fjump[20]+0.25*alphaDrag[0]*favg[20]; 

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
  incr2[3] = vtSq[0]*(0.2118037767457181*(fr[18]+fl[18])-0.3046875*fr[3]+0.3046875*fl[3]+0.2165063509461096*(fr[0]+fl[0]))+0.2165063509461096*(vtSq[5]*(fr[17]+fl[17])+vtSq[4]*(fr[16]+fl[16]))+vtSq[2]*((-0.3046875*fr[8])+0.3046875*fl[8]+0.2165063509461096*(fr[2]+fl[2]))+vtSq[1]*((-0.3046875*fr[7])+0.3046875*fl[7]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*vtSq[3]*(fr[6]+fl[6]); 
  incr2[7] = vtSq[1]*(0.2118037767457181*(fr[18]+fl[18])+0.1936491673103708*(fr[16]+fl[16])-0.3046875*fr[3]+0.3046875*fl[3]+0.2165063509461096*(fr[0]+fl[0]))+vtSq[3]*((-0.3046875*fr[8])+0.3046875*fl[8]+0.2165063509461096*(fr[2]+fl[2]))+vtSq[4]*((-0.2725207847577868*fr[7])+0.2725207847577868*fl[7]+0.1936491673103708*(fr[1]+fl[1]))+vtSq[0]*((-0.3046875*fr[7])+0.3046875*fl[7]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*vtSq[2]*(fr[6]+fl[6]); 
  incr2[8] = vtSq[2]*(0.2118037767457181*(fr[18]+fl[18])+0.1936491673103708*(fr[17]+fl[17])-0.3046875*fr[3]+0.3046875*fl[3]+0.2165063509461096*(fr[0]+fl[0]))+vtSq[5]*((-0.2725207847577868*fr[8])+0.2725207847577868*fl[8]+0.1936491673103708*(fr[2]+fl[2]))+vtSq[0]*((-0.3046875*fr[8])+0.3046875*fl[8]+0.2165063509461096*(fr[2]+fl[2]))+vtSq[3]*((-0.3046875*fr[7])+0.3046875*fl[7]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*vtSq[1]*(fr[6]+fl[6]); 
  incr2[11] = vtSq[0]*((-0.3046875*fr[11])+0.3046875*fl[11]+0.2165063509461096*(fr[4]+fl[4]))+0.2165063509461096*(vtSq[2]*(fr[10]+fl[10])+vtSq[1]*(fr[9]+fl[9])); 
  incr2[14] = vtSq[0]*((-0.3046875*fr[14])+0.3046875*fl[14]+0.2165063509461096*(fr[5]+fl[5]))+0.2165063509461096*(vtSq[2]*(fr[13]+fl[13])+vtSq[1]*(fr[12]+fl[12])); 
  incr2[18] = vtSq[0]*((-0.8203125*(fr[18]+fl[18]))+1.180049613297572*fr[3]-1.180049613297572*fl[3]-0.8385254915624212*(fr[0]+fl[0]))-0.8385254915624212*(vtSq[5]*(fr[17]+fl[17])+vtSq[4]*(fr[16]+fl[16]))+vtSq[2]*(1.180049613297572*fr[8]-1.180049613297572*fl[8]-0.8385254915624212*(fr[2]+fl[2]))+vtSq[1]*(1.180049613297572*fr[7]-1.180049613297572*fl[7]-0.8385254915624212*(fr[1]+fl[1]))-0.8385254915624212*vtSq[3]*(fr[6]+fl[6]); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr1[2]*rdv2nuR; 
  outr[3] += incr2[3]*rdvSq4nuR+incr1[3]*rdv2nuR; 
  outr[4] += incr1[4]*rdv2nuR; 
  outr[5] += incr1[5]*rdv2nuR; 
  outr[6] += incr1[6]*rdv2nuR; 
  outr[7] += incr2[7]*rdvSq4nuR+incr1[7]*rdv2nuR; 
  outr[8] += incr2[8]*rdvSq4nuR+incr1[8]*rdv2nuR; 
  outr[9] += incr1[9]*rdv2nuR; 
  outr[10] += incr1[10]*rdv2nuR; 
  outr[11] += incr2[11]*rdvSq4nuR+incr1[11]*rdv2nuR; 
  outr[12] += incr1[12]*rdv2nuR; 
  outr[13] += incr1[13]*rdv2nuR; 
  outr[14] += incr2[14]*rdvSq4nuR+incr1[14]*rdv2nuR; 
  outr[15] += incr1[15]*rdv2nuR; 
  outr[16] += incr1[16]*rdv2nuR; 
  outr[17] += incr1[17]*rdv2nuR; 
  outr[18] += incr2[18]*rdvSq4nuR+incr1[18]*rdv2nuR; 
  outr[19] += incr1[19]*rdv2nuR; 
  outr[20] += incr1[20]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += -1.0*incr1[2]*rdv2nuL; 
  outl[3] += incr1[3]*rdv2nuL-1.0*incr2[3]*rdvSq4nuL; 
  outl[4] += -1.0*incr1[4]*rdv2nuL; 
  outl[5] += -1.0*incr1[5]*rdv2nuL; 
  outl[6] += -1.0*incr1[6]*rdv2nuL; 
  outl[7] += incr1[7]*rdv2nuL-1.0*incr2[7]*rdvSq4nuL; 
  outl[8] += incr1[8]*rdv2nuL-1.0*incr2[8]*rdvSq4nuL; 
  outl[9] += -1.0*incr1[9]*rdv2nuL; 
  outl[10] += -1.0*incr1[10]*rdv2nuL; 
  outl[11] += incr1[11]*rdv2nuL-1.0*incr2[11]*rdvSq4nuL; 
  outl[12] += -1.0*incr1[12]*rdv2nuL; 
  outl[13] += -1.0*incr1[13]*rdv2nuL; 
  outl[14] += incr1[14]*rdv2nuL-1.0*incr2[14]*rdvSq4nuL; 
  outl[15] += -1.0*incr1[15]*rdv2nuL; 
  outl[16] += -1.0*incr1[16]*rdv2nuL; 
  outl[17] += -1.0*incr1[17]*rdv2nuL; 
  outl[18] += incr2[18]*rdvSq4nuL-1.0*incr1[18]*rdv2nuL; 
  outl[19] += -1.0*incr1[19]*rdv2nuL; 
  outl[20] += -1.0*incr1[20]*rdv2nuL; 

  return std::abs(0.5590169943749475*Ux[5]+0.5590169943749475*Ux[4]+wl[2]-0.5*Ux[0]); 
} 
double VmLBOconstNuSurf2x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates. 
  // dxv[5]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[9]:       bulk velocity (in 3 directions). 
  // vtSq[3]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[3]; 
  double rdv2nuL = 2.0*nu/dxvl[3]; 
  double rdv2nuR = 2.0*nu/dxvr[3]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[3]*dxvl[3]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[3]*dxvr[3]); 

  const double *Uy = &u[3]; 

  double favg[6]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 

  double fjump[6]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(1*fr[5])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 2.0*wl[3]+dxvl[3]-1.0*Uy[0]; 
  alphaDrag[1] = -1.0*Uy[1]; 
  alphaDrag[2] = -1.0*Uy[2]; 

  double Ghat[6]; 
  for(unsigned short int i=0; i<6; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.082531754730548*vtSq[0]*fr[4])-1.082531754730548*vtSq[0]*fl[4]+1.125*fr[2]*vtSq[2]-1.125*fl[2]*vtSq[2]+1.125*fr[1]*vtSq[1]-1.125*fl[1]*vtSq[1]+1.125*fr[0]*vtSq[0]-1.125*fl[0]*vtSq[0])*rdv-0.8660254037844386*fjump[4]+alphaDrag[0]*(0.4330127018922193*favg[4]+0.25*favg[0])+0.25*alphaDrag[2]*favg[2]+0.25*alphaDrag[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = ((-1.082531754730548*vtSq[1]*fr[4])-1.082531754730548*vtSq[1]*fl[4]+1.125*fr[0]*vtSq[1]-1.125*fl[0]*vtSq[1]+1.125*vtSq[0]*fr[1]-1.125*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[1]+0.25*alphaDrag[0]*favg[1]; 
  Ghat[2] = ((-1.082531754730548*vtSq[2]*fr[4])-1.082531754730548*vtSq[2]*fl[4]+1.125*fr[0]*vtSq[2]-1.125*fl[0]*vtSq[2]+1.125*vtSq[0]*fr[2]-1.125*vtSq[0]*fl[2])*rdv+alphaDrag[2]*(0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[2]+0.25*alphaDrag[0]*favg[2]; 
  Ghat[3] = (1.125*vtSq[0]*fr[3]-1.125*vtSq[0]*fl[3])*rdv-0.5*fjump[3]+0.25*alphaDrag[0]*favg[3]; 
  Ghat[5] = (1.125*vtSq[0]*fr[5]-1.125*vtSq[0]*fl[5])*rdv-0.5*fjump[5]+0.25*alphaDrag[0]*favg[5]; 

  double incr1[6]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[0]; 
  incr1[5] = -0.5*Ghat[5]; 

  double incr2[6]; 
  incr2[4] = vtSq[0]*((-0.25*fr[4])+0.25*fl[4]+0.2165063509461096*(fr[0]+fl[0]))+0.2165063509461096*((fr[2]+fl[2])*vtSq[2]+(fr[1]+fl[1])*vtSq[1]); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr1[2]*rdv2nuR; 
  outr[3] += incr1[3]*rdv2nuR; 
  outr[4] += incr2[4]*rdvSq4nuR+incr1[4]*rdv2nuR; 
  outr[5] += incr1[5]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += -1.0*incr1[2]*rdv2nuL; 
  outl[3] += -1.0*incr1[3]*rdv2nuL; 
  outl[4] += incr1[4]*rdv2nuL-1.0*incr2[4]*rdvSq4nuL; 
  outl[5] += -1.0*incr1[5]*rdv2nuL; 

  return std::abs(wl[3]-0.5*Uy[0]); 
} 
double VmLBOconstNuSurf2x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates. 
  // dxv[5]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[18]:       bulk velocity (in 3 directions). 
  // vtSq[6]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[3]; 
  double rdv2nuL = 2.0*nu/dxvl[3]; 
  double rdv2nuR = 2.0*nu/dxvr[3]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[3]*dxvl[3]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[3]*dxvr[3]); 

  const double *Uy = &u[6]; 

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
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = vMuMidMax*(fl[9]-(-1*fr[9])); 
  fjump[10] = vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = vMuMidMax*(fl[14]-(1*fr[14])); 
  fjump[15] = vMuMidMax*(fl[15]-(-1*fr[15])); 
  fjump[16] = vMuMidMax*(fl[16]-(1*fr[16])); 
  fjump[17] = vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = vMuMidMax*(fl[19]-(1*fr[19])); 
  fjump[20] = vMuMidMax*(fl[20]-(1*fr[20])); 

  double alphaDrag[6]; 
  alphaDrag[0] = 2.0*wl[3]+dxvl[3]-1.0*Uy[0]; 
  alphaDrag[1] = -1.0*Uy[1]; 
  alphaDrag[2] = -1.0*Uy[2]; 
  alphaDrag[3] = -1.0*Uy[3]; 
  alphaDrag[4] = -1.0*Uy[4]; 
  alphaDrag[5] = -1.0*Uy[5]; 

  double Ghat[21]; 
  for(unsigned short int i=0; i<21; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.341640786499874*vtSq[0]*fr[19]-1.341640786499874*vtSq[0]*fl[19]+1.875*vtSq[5]*fr[17]-1.875*vtSq[5]*fl[17]+1.875*vtSq[4]*fr[16]-1.875*vtSq[4]*fl[16]-2.381569860407206*vtSq[2]*fr[10]-2.381569860407206*vtSq[2]*fl[10]-2.381569860407206*vtSq[1]*fr[9]-2.381569860407206*vtSq[1]*fl[9]+1.875*vtSq[3]*fr[6]-1.875*vtSq[3]*fl[6]-2.381569860407206*vtSq[0]*fr[4]-2.381569860407206*vtSq[0]*fl[4]+1.875*fr[2]*vtSq[2]-1.875*fl[2]*vtSq[2]+1.875*fr[1]*vtSq[1]-1.875*fl[1]*vtSq[1]+1.875*fr[0]*vtSq[0]-1.875*fl[0]*vtSq[0])*rdv-1.118033988749895*fjump[19]+alphaDrag[0]*(0.5590169943749475*favg[19]+0.4330127018922193*favg[4]+0.25*favg[0])+0.25*alphaDrag[5]*favg[17]+0.25*alphaDrag[4]*favg[16]+alphaDrag[2]*(0.4330127018922193*favg[10]+0.25*favg[2])+alphaDrag[1]*(0.4330127018922193*favg[9]+0.25*favg[1])+0.25*alphaDrag[3]*favg[6]-0.8660254037844386*fjump[4]-0.5*fjump[0]; 
  Ghat[1] = (1.341640786499874*vtSq[1]*fr[19]-1.341640786499874*vtSq[1]*fl[19]+1.677050983124842*vtSq[1]*fr[16]-1.677050983124842*vtSq[1]*fl[16]-2.381569860407206*vtSq[3]*fr[10]-2.381569860407206*vtSq[3]*fl[10]-2.130140840414079*vtSq[4]*fr[9]-2.381569860407206*vtSq[0]*fr[9]-2.130140840414079*vtSq[4]*fl[9]-2.381569860407206*vtSq[0]*fl[9]+1.875*vtSq[2]*fr[6]-1.875*vtSq[2]*fl[6]+1.677050983124842*fr[1]*vtSq[4]-1.677050983124842*fl[1]*vtSq[4]-2.381569860407206*vtSq[1]*fr[4]-2.381569860407206*vtSq[1]*fl[4]+1.875*fr[2]*vtSq[3]-1.875*fl[2]*vtSq[3]+1.875*fr[0]*vtSq[1]-1.875*fl[0]*vtSq[1]+1.875*vtSq[0]*fr[1]-1.875*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.5590169943749475*favg[19]+0.223606797749979*favg[16]+0.4330127018922193*favg[4]+0.25*favg[0])+alphaDrag[3]*(0.4330127018922193*favg[10]+0.25*favg[2])-0.8660254037844386*fjump[9]+alphaDrag[0]*(0.4330127018922193*favg[9]+0.25*favg[1])+alphaDrag[4]*(0.3872983346207416*favg[9]+0.223606797749979*favg[1])+0.25*alphaDrag[2]*favg[6]-0.5*fjump[1]; 
  Ghat[2] = (1.341640786499874*vtSq[2]*fr[19]-1.341640786499874*vtSq[2]*fl[19]+1.677050983124842*vtSq[2]*fr[17]-1.677050983124842*vtSq[2]*fl[17]-2.130140840414079*vtSq[5]*fr[10]-2.381569860407206*vtSq[0]*fr[10]-2.130140840414079*vtSq[5]*fl[10]-2.381569860407206*vtSq[0]*fl[10]-2.381569860407206*vtSq[3]*fr[9]-2.381569860407206*vtSq[3]*fl[9]+1.875*vtSq[1]*fr[6]-1.875*vtSq[1]*fl[6]+1.677050983124842*fr[2]*vtSq[5]-1.677050983124842*fl[2]*vtSq[5]-2.381569860407206*vtSq[2]*fr[4]-2.381569860407206*vtSq[2]*fl[4]+1.875*fr[1]*vtSq[3]-1.875*fl[1]*vtSq[3]+1.875*fr[0]*vtSq[2]-1.875*fl[0]*vtSq[2]+1.875*vtSq[0]*fr[2]-1.875*vtSq[0]*fl[2])*rdv+alphaDrag[2]*(0.5590169943749475*favg[19]+0.223606797749979*favg[17]+0.4330127018922193*favg[4]+0.25*favg[0])-0.8660254037844386*fjump[10]+alphaDrag[0]*(0.4330127018922193*favg[10]+0.25*favg[2])+alphaDrag[5]*(0.3872983346207416*favg[10]+0.223606797749979*favg[2])+alphaDrag[3]*(0.4330127018922193*favg[9]+0.25*favg[1])+0.25*alphaDrag[1]*favg[6]-0.5*fjump[2]; 
  Ghat[3] = ((-2.381569860407206*vtSq[0]*fr[11])-2.381569860407206*vtSq[0]*fl[11]+1.875*vtSq[2]*fr[8]-1.875*vtSq[2]*fl[8]+1.875*vtSq[1]*fr[7]-1.875*vtSq[1]*fl[7]+1.875*vtSq[0]*fr[3]-1.875*vtSq[0]*fl[3])*rdv-0.8660254037844386*fjump[11]+alphaDrag[0]*(0.4330127018922193*favg[11]+0.25*favg[3])+0.25*alphaDrag[2]*favg[8]+0.25*alphaDrag[1]*favg[7]-0.5*fjump[3]; 
  Ghat[5] = ((-2.381569860407206*vtSq[0]*fr[15])-2.381569860407206*vtSq[0]*fl[15]+1.875*vtSq[2]*fr[13]-1.875*vtSq[2]*fl[13]+1.875*vtSq[1]*fr[12]-1.875*vtSq[1]*fl[12]+1.875*vtSq[0]*fr[5]-1.875*vtSq[0]*fl[5])*rdv-0.8660254037844386*fjump[15]+alphaDrag[0]*(0.4330127018922193*favg[15]+0.25*favg[5])+0.25*alphaDrag[2]*favg[13]+0.25*alphaDrag[1]*favg[12]-0.5*fjump[5]; 
  Ghat[6] = (1.341640786499874*vtSq[3]*fr[19]-1.341640786499874*vtSq[3]*fl[19]+1.677050983124842*vtSq[3]*fr[17]-1.677050983124842*vtSq[3]*fl[17]+1.677050983124842*vtSq[3]*fr[16]-1.677050983124842*vtSq[3]*fl[16]-2.381569860407206*vtSq[1]*fr[10]-2.381569860407206*vtSq[1]*fl[10]-2.381569860407206*vtSq[2]*fr[9]-2.381569860407206*vtSq[2]*fl[9]+1.677050983124842*vtSq[5]*fr[6]+1.677050983124842*vtSq[4]*fr[6]+1.875*vtSq[0]*fr[6]-1.677050983124842*vtSq[5]*fl[6]-1.677050983124842*vtSq[4]*fl[6]-1.875*vtSq[0]*fl[6]-2.381569860407206*vtSq[3]*fr[4]-2.381569860407206*vtSq[3]*fl[4]+1.875*fr[0]*vtSq[3]-1.875*fl[0]*vtSq[3]+1.875*fr[1]*vtSq[2]-1.875*fl[1]*vtSq[2]+1.875*vtSq[1]*fr[2]-1.875*vtSq[1]*fl[2])*rdv+alphaDrag[3]*(0.5590169943749475*favg[19]+0.223606797749979*favg[17]+0.223606797749979*favg[16]+0.4330127018922193*favg[4]+0.25*favg[0])+alphaDrag[1]*(0.4330127018922193*favg[10]+0.25*favg[2])+alphaDrag[2]*(0.4330127018922193*favg[9]+0.25*favg[1])-0.5*fjump[6]+0.223606797749979*alphaDrag[5]*favg[6]+0.223606797749979*alphaDrag[4]*favg[6]+0.25*alphaDrag[0]*favg[6]; 
  Ghat[7] = ((-2.381569860407206*vtSq[1]*fr[11])-2.381569860407206*vtSq[1]*fl[11]+1.875*vtSq[3]*fr[8]-1.875*vtSq[3]*fl[8]+1.677050983124842*vtSq[4]*fr[7]+1.875*vtSq[0]*fr[7]-1.677050983124842*vtSq[4]*fl[7]-1.875*vtSq[0]*fl[7]+1.875*vtSq[1]*fr[3]-1.875*vtSq[1]*fl[3])*rdv+alphaDrag[1]*(0.4330127018922193*favg[11]+0.25*favg[3])+0.25*alphaDrag[3]*favg[8]-0.5*fjump[7]+0.223606797749979*alphaDrag[4]*favg[7]+0.25*alphaDrag[0]*favg[7]; 
  Ghat[8] = ((-2.381569860407206*vtSq[2]*fr[11])-2.381569860407206*vtSq[2]*fl[11]+1.677050983124842*vtSq[5]*fr[8]+1.875*vtSq[0]*fr[8]-1.677050983124842*vtSq[5]*fl[8]-1.875*vtSq[0]*fl[8]+1.875*vtSq[3]*fr[7]-1.875*vtSq[3]*fl[7]+1.875*vtSq[2]*fr[3]-1.875*vtSq[2]*fl[3])*rdv+alphaDrag[2]*(0.4330127018922193*favg[11]+0.25*favg[3])-0.5*fjump[8]+0.223606797749979*alphaDrag[5]*favg[8]+0.25*alphaDrag[0]*favg[8]+0.25*alphaDrag[3]*favg[7]; 
  Ghat[12] = ((-2.381569860407206*vtSq[1]*fr[15])-2.381569860407206*vtSq[1]*fl[15]+1.875*vtSq[3]*fr[13]-1.875*vtSq[3]*fl[13]+1.677050983124842*vtSq[4]*fr[12]+1.875*vtSq[0]*fr[12]-1.677050983124842*vtSq[4]*fl[12]-1.875*vtSq[0]*fl[12]+1.875*vtSq[1]*fr[5]-1.875*vtSq[1]*fl[5])*rdv+alphaDrag[1]*(0.4330127018922193*favg[15]+0.25*favg[5])+0.25*alphaDrag[3]*favg[13]-0.5*fjump[12]+0.223606797749979*alphaDrag[4]*favg[12]+0.25*alphaDrag[0]*favg[12]; 
  Ghat[13] = ((-2.381569860407206*vtSq[2]*fr[15])-2.381569860407206*vtSq[2]*fl[15]+1.677050983124842*vtSq[5]*fr[13]+1.875*vtSq[0]*fr[13]-1.677050983124842*vtSq[5]*fl[13]-1.875*vtSq[0]*fl[13]+1.875*vtSq[3]*fr[12]-1.875*vtSq[3]*fl[12]+1.875*vtSq[2]*fr[5]-1.875*vtSq[2]*fl[5])*rdv+alphaDrag[2]*(0.4330127018922193*favg[15]+0.25*favg[5])-0.5*fjump[13]+0.223606797749979*alphaDrag[5]*favg[13]+0.25*alphaDrag[0]*favg[13]+0.25*alphaDrag[3]*favg[12]; 
  Ghat[14] = (1.875*vtSq[0]*fr[14]-1.875*vtSq[0]*fl[14])*rdv-0.5*fjump[14]+0.25*alphaDrag[0]*favg[14]; 
  Ghat[16] = (1.341640786499874*vtSq[4]*fr[19]-1.341640786499874*vtSq[4]*fl[19]+1.197893559374888*vtSq[4]*fr[16]+1.875*vtSq[0]*fr[16]-1.197893559374888*vtSq[4]*fl[16]-1.875*vtSq[0]*fl[16]-2.130140840414079*vtSq[1]*fr[9]-2.130140840414079*vtSq[1]*fl[9]+1.677050983124842*vtSq[3]*fr[6]-1.677050983124842*vtSq[3]*fl[6]-2.381569860407206*fr[4]*vtSq[4]-2.381569860407206*fl[4]*vtSq[4]+1.875*fr[0]*vtSq[4]-1.875*fl[0]*vtSq[4]+1.677050983124842*fr[1]*vtSq[1]-1.677050983124842*fl[1]*vtSq[1])*rdv+alphaDrag[4]*(0.5590169943749475*favg[19]+0.159719141249985*favg[16]+0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[16]+0.25*alphaDrag[0]*favg[16]+alphaDrag[1]*(0.3872983346207416*favg[9]+0.223606797749979*favg[1])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[17] = (1.341640786499874*vtSq[5]*fr[19]-1.341640786499874*vtSq[5]*fl[19]+1.197893559374888*vtSq[5]*fr[17]+1.875*vtSq[0]*fr[17]-1.197893559374888*vtSq[5]*fl[17]-1.875*vtSq[0]*fl[17]-2.130140840414079*vtSq[2]*fr[10]-2.130140840414079*vtSq[2]*fl[10]+1.677050983124842*vtSq[3]*fr[6]-1.677050983124842*vtSq[3]*fl[6]-2.381569860407206*fr[4]*vtSq[5]-2.381569860407206*fl[4]*vtSq[5]+1.875*fr[0]*vtSq[5]-1.875*fl[0]*vtSq[5]+1.677050983124842*fr[2]*vtSq[2]-1.677050983124842*fl[2]*vtSq[2])*rdv+alphaDrag[5]*(0.5590169943749475*favg[19]+0.159719141249985*favg[17]+0.4330127018922193*favg[4]+0.25*favg[0])-0.5*fjump[17]+0.25*alphaDrag[0]*favg[17]+alphaDrag[2]*(0.3872983346207416*favg[10]+0.223606797749979*favg[2])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[18] = (1.875*vtSq[0]*fr[18]-1.875*vtSq[0]*fl[18])*rdv-0.5*fjump[18]+0.25*alphaDrag[0]*favg[18]; 
  Ghat[20] = (1.875*vtSq[0]*fr[20]-1.875*vtSq[0]*fl[20])*rdv-0.5*fjump[20]+0.25*alphaDrag[0]*favg[20]; 

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
  incr2[4] = vtSq[0]*(0.2118037767457181*(fr[19]+fl[19])-0.3046875*fr[4]+0.3046875*fl[4]+0.2165063509461096*(fr[0]+fl[0]))+0.2165063509461096*(vtSq[5]*(fr[17]+fl[17])+vtSq[4]*(fr[16]+fl[16]))+vtSq[2]*((-0.3046875*fr[10])+0.3046875*fl[10]+0.2165063509461096*(fr[2]+fl[2]))+vtSq[1]*((-0.3046875*fr[9])+0.3046875*fl[9]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*vtSq[3]*(fr[6]+fl[6]); 
  incr2[9] = vtSq[1]*(0.2118037767457181*(fr[19]+fl[19])+0.1936491673103708*(fr[16]+fl[16])-0.3046875*fr[4]+0.3046875*fl[4]+0.2165063509461096*(fr[0]+fl[0]))+vtSq[3]*((-0.3046875*fr[10])+0.3046875*fl[10]+0.2165063509461096*(fr[2]+fl[2]))+vtSq[4]*((-0.2725207847577868*fr[9])+0.2725207847577868*fl[9]+0.1936491673103708*(fr[1]+fl[1]))+vtSq[0]*((-0.3046875*fr[9])+0.3046875*fl[9]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*vtSq[2]*(fr[6]+fl[6]); 
  incr2[10] = vtSq[2]*(0.2118037767457181*(fr[19]+fl[19])+0.1936491673103708*(fr[17]+fl[17])-0.3046875*fr[4]+0.3046875*fl[4]+0.2165063509461096*(fr[0]+fl[0]))+vtSq[5]*((-0.2725207847577868*fr[10])+0.2725207847577868*fl[10]+0.1936491673103708*(fr[2]+fl[2]))+vtSq[0]*((-0.3046875*fr[10])+0.3046875*fl[10]+0.2165063509461096*(fr[2]+fl[2]))+vtSq[3]*((-0.3046875*fr[9])+0.3046875*fl[9]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*vtSq[1]*(fr[6]+fl[6]); 
  incr2[11] = vtSq[0]*((-0.3046875*fr[11])+0.3046875*fl[11]+0.2165063509461096*(fr[3]+fl[3]))+0.2165063509461096*(vtSq[2]*(fr[8]+fl[8])+vtSq[1]*(fr[7]+fl[7])); 
  incr2[15] = vtSq[0]*((-0.3046875*fr[15])+0.3046875*fl[15]+0.2165063509461096*(fr[5]+fl[5]))+0.2165063509461096*(vtSq[2]*(fr[13]+fl[13])+vtSq[1]*(fr[12]+fl[12])); 
  incr2[19] = vtSq[0]*((-0.8203125*(fr[19]+fl[19]))+1.180049613297572*fr[4]-1.180049613297572*fl[4]-0.8385254915624212*(fr[0]+fl[0]))-0.8385254915624212*(vtSq[5]*(fr[17]+fl[17])+vtSq[4]*(fr[16]+fl[16]))+vtSq[2]*(1.180049613297572*fr[10]-1.180049613297572*fl[10]-0.8385254915624212*(fr[2]+fl[2]))+vtSq[1]*(1.180049613297572*fr[9]-1.180049613297572*fl[9]-0.8385254915624212*(fr[1]+fl[1]))-0.8385254915624212*vtSq[3]*(fr[6]+fl[6]); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr1[2]*rdv2nuR; 
  outr[3] += incr1[3]*rdv2nuR; 
  outr[4] += incr2[4]*rdvSq4nuR+incr1[4]*rdv2nuR; 
  outr[5] += incr1[5]*rdv2nuR; 
  outr[6] += incr1[6]*rdv2nuR; 
  outr[7] += incr1[7]*rdv2nuR; 
  outr[8] += incr1[8]*rdv2nuR; 
  outr[9] += incr2[9]*rdvSq4nuR+incr1[9]*rdv2nuR; 
  outr[10] += incr2[10]*rdvSq4nuR+incr1[10]*rdv2nuR; 
  outr[11] += incr2[11]*rdvSq4nuR+incr1[11]*rdv2nuR; 
  outr[12] += incr1[12]*rdv2nuR; 
  outr[13] += incr1[13]*rdv2nuR; 
  outr[14] += incr1[14]*rdv2nuR; 
  outr[15] += incr2[15]*rdvSq4nuR+incr1[15]*rdv2nuR; 
  outr[16] += incr1[16]*rdv2nuR; 
  outr[17] += incr1[17]*rdv2nuR; 
  outr[18] += incr1[18]*rdv2nuR; 
  outr[19] += incr2[19]*rdvSq4nuR+incr1[19]*rdv2nuR; 
  outr[20] += incr1[20]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += -1.0*incr1[2]*rdv2nuL; 
  outl[3] += -1.0*incr1[3]*rdv2nuL; 
  outl[4] += incr1[4]*rdv2nuL-1.0*incr2[4]*rdvSq4nuL; 
  outl[5] += -1.0*incr1[5]*rdv2nuL; 
  outl[6] += -1.0*incr1[6]*rdv2nuL; 
  outl[7] += -1.0*incr1[7]*rdv2nuL; 
  outl[8] += -1.0*incr1[8]*rdv2nuL; 
  outl[9] += incr1[9]*rdv2nuL-1.0*incr2[9]*rdvSq4nuL; 
  outl[10] += incr1[10]*rdv2nuL-1.0*incr2[10]*rdvSq4nuL; 
  outl[11] += incr1[11]*rdv2nuL-1.0*incr2[11]*rdvSq4nuL; 
  outl[12] += -1.0*incr1[12]*rdv2nuL; 
  outl[13] += -1.0*incr1[13]*rdv2nuL; 
  outl[14] += -1.0*incr1[14]*rdv2nuL; 
  outl[15] += incr1[15]*rdv2nuL-1.0*incr2[15]*rdvSq4nuL; 
  outl[16] += -1.0*incr1[16]*rdv2nuL; 
  outl[17] += -1.0*incr1[17]*rdv2nuL; 
  outl[18] += -1.0*incr1[18]*rdv2nuL; 
  outl[19] += incr2[19]*rdvSq4nuL-1.0*incr1[19]*rdv2nuL; 
  outl[20] += -1.0*incr1[20]*rdv2nuL; 

  return std::abs(0.5590169943749475*Uy[5]+0.5590169943749475*Uy[4]+wl[3]-0.5*Uy[0]); 
} 
double VmLBOconstNuSurf2x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates. 
  // dxv[5]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[9]:       bulk velocity (in 3 directions). 
  // vtSq[3]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[4]; 
  double rdv2nuL = 2.0*nu/dxvl[4]; 
  double rdv2nuR = 2.0*nu/dxvr[4]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[4]*dxvl[4]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[4]*dxvr[4]); 

  const double *Uz = &u[6]; 

  double favg[6]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 

  double fjump[6]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(-1*fr[5])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 2.0*wl[4]+dxvl[4]-1.0*Uz[0]; 
  alphaDrag[1] = -1.0*Uz[1]; 
  alphaDrag[2] = -1.0*Uz[2]; 

  double Ghat[6]; 
  for(unsigned short int i=0; i<6; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.082531754730548*vtSq[0]*fr[5])-1.082531754730548*vtSq[0]*fl[5]+1.125*fr[2]*vtSq[2]-1.125*fl[2]*vtSq[2]+1.125*fr[1]*vtSq[1]-1.125*fl[1]*vtSq[1]+1.125*fr[0]*vtSq[0]-1.125*fl[0]*vtSq[0])*rdv-0.8660254037844386*fjump[5]+alphaDrag[0]*(0.4330127018922193*favg[5]+0.25*favg[0])+0.25*alphaDrag[2]*favg[2]+0.25*alphaDrag[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = ((-1.082531754730548*vtSq[1]*fr[5])-1.082531754730548*vtSq[1]*fl[5]+1.125*fr[0]*vtSq[1]-1.125*fl[0]*vtSq[1]+1.125*vtSq[0]*fr[1]-1.125*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.4330127018922193*favg[5]+0.25*favg[0])-0.5*fjump[1]+0.25*alphaDrag[0]*favg[1]; 
  Ghat[2] = ((-1.082531754730548*vtSq[2]*fr[5])-1.082531754730548*vtSq[2]*fl[5]+1.125*fr[0]*vtSq[2]-1.125*fl[0]*vtSq[2]+1.125*vtSq[0]*fr[2]-1.125*vtSq[0]*fl[2])*rdv+alphaDrag[2]*(0.4330127018922193*favg[5]+0.25*favg[0])-0.5*fjump[2]+0.25*alphaDrag[0]*favg[2]; 
  Ghat[3] = (1.125*vtSq[0]*fr[3]-1.125*vtSq[0]*fl[3])*rdv-0.5*fjump[3]+0.25*alphaDrag[0]*favg[3]; 
  Ghat[4] = (1.125*vtSq[0]*fr[4]-1.125*vtSq[0]*fl[4])*rdv-0.5*fjump[4]+0.25*alphaDrag[0]*favg[4]; 

  double incr1[6]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[0]; 

  double incr2[6]; 
  incr2[5] = vtSq[0]*((-0.25*fr[5])+0.25*fl[5]+0.2165063509461096*(fr[0]+fl[0]))+0.2165063509461096*((fr[2]+fl[2])*vtSq[2]+(fr[1]+fl[1])*vtSq[1]); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr1[2]*rdv2nuR; 
  outr[3] += incr1[3]*rdv2nuR; 
  outr[4] += incr1[4]*rdv2nuR; 
  outr[5] += incr2[5]*rdvSq4nuR+incr1[5]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += -1.0*incr1[2]*rdv2nuL; 
  outl[3] += -1.0*incr1[3]*rdv2nuL; 
  outl[4] += -1.0*incr1[4]*rdv2nuL; 
  outl[5] += incr1[5]*rdv2nuL-1.0*incr2[5]*rdvSq4nuL; 

  return std::abs(wl[4]-0.5*Uz[0]); 
} 
double VmLBOconstNuSurf2x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates. 
  // dxv[5]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[18]:       bulk velocity (in 3 directions). 
  // vtSq[6]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[4]; 
  double rdv2nuL = 2.0*nu/dxvl[4]; 
  double rdv2nuR = 2.0*nu/dxvr[4]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[4]*dxvl[4]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[4]*dxvr[4]); 

  const double *Uz = &u[12]; 

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
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = vMuMidMax*(fl[10]-(1*fr[10])); 
  fjump[11] = vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = vMuMidMax*(fl[12]-(-1*fr[12])); 
  fjump[13] = vMuMidMax*(fl[13]-(-1*fr[13])); 
  fjump[14] = vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = vMuMidMax*(fl[15]-(-1*fr[15])); 
  fjump[16] = vMuMidMax*(fl[16]-(1*fr[16])); 
  fjump[17] = vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = vMuMidMax*(fl[19]-(1*fr[19])); 
  fjump[20] = vMuMidMax*(fl[20]-(1*fr[20])); 

  double alphaDrag[6]; 
  alphaDrag[0] = 2.0*wl[4]+dxvl[4]-1.0*Uz[0]; 
  alphaDrag[1] = -1.0*Uz[1]; 
  alphaDrag[2] = -1.0*Uz[2]; 
  alphaDrag[3] = -1.0*Uz[3]; 
  alphaDrag[4] = -1.0*Uz[4]; 
  alphaDrag[5] = -1.0*Uz[5]; 

  double Ghat[21]; 
  for(unsigned short int i=0; i<21; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.341640786499874*vtSq[0]*fr[20]-1.341640786499874*vtSq[0]*fl[20]+1.875*vtSq[5]*fr[17]-1.875*vtSq[5]*fl[17]+1.875*vtSq[4]*fr[16]-1.875*vtSq[4]*fl[16]-2.381569860407206*vtSq[2]*fr[13]-2.381569860407206*vtSq[2]*fl[13]-2.381569860407206*vtSq[1]*fr[12]-2.381569860407206*vtSq[1]*fl[12]+1.875*vtSq[3]*fr[6]-1.875*vtSq[3]*fl[6]-2.381569860407206*vtSq[0]*fr[5]-2.381569860407206*vtSq[0]*fl[5]+1.875*fr[2]*vtSq[2]-1.875*fl[2]*vtSq[2]+1.875*fr[1]*vtSq[1]-1.875*fl[1]*vtSq[1]+1.875*fr[0]*vtSq[0]-1.875*fl[0]*vtSq[0])*rdv-1.118033988749895*fjump[20]+alphaDrag[0]*(0.5590169943749475*favg[20]+0.4330127018922193*favg[5]+0.25*favg[0])+0.25*alphaDrag[5]*favg[17]+0.25*alphaDrag[4]*favg[16]+alphaDrag[2]*(0.4330127018922193*favg[13]+0.25*favg[2])+alphaDrag[1]*(0.4330127018922193*favg[12]+0.25*favg[1])+0.25*alphaDrag[3]*favg[6]-0.8660254037844386*fjump[5]-0.5*fjump[0]; 
  Ghat[1] = (1.341640786499874*vtSq[1]*fr[20]-1.341640786499874*vtSq[1]*fl[20]+1.677050983124842*vtSq[1]*fr[16]-1.677050983124842*vtSq[1]*fl[16]-2.381569860407206*vtSq[3]*fr[13]-2.381569860407206*vtSq[3]*fl[13]-2.130140840414079*vtSq[4]*fr[12]-2.381569860407206*vtSq[0]*fr[12]-2.130140840414079*vtSq[4]*fl[12]-2.381569860407206*vtSq[0]*fl[12]+1.875*vtSq[2]*fr[6]-1.875*vtSq[2]*fl[6]-2.381569860407206*vtSq[1]*fr[5]-2.381569860407206*vtSq[1]*fl[5]+1.677050983124842*fr[1]*vtSq[4]-1.677050983124842*fl[1]*vtSq[4]+1.875*fr[2]*vtSq[3]-1.875*fl[2]*vtSq[3]+1.875*fr[0]*vtSq[1]-1.875*fl[0]*vtSq[1]+1.875*vtSq[0]*fr[1]-1.875*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.5590169943749475*favg[20]+0.223606797749979*favg[16]+0.4330127018922193*favg[5]+0.25*favg[0])+alphaDrag[3]*(0.4330127018922193*favg[13]+0.25*favg[2])-0.8660254037844386*fjump[12]+alphaDrag[0]*(0.4330127018922193*favg[12]+0.25*favg[1])+alphaDrag[4]*(0.3872983346207416*favg[12]+0.223606797749979*favg[1])+0.25*alphaDrag[2]*favg[6]-0.5*fjump[1]; 
  Ghat[2] = (1.341640786499874*vtSq[2]*fr[20]-1.341640786499874*vtSq[2]*fl[20]+1.677050983124842*vtSq[2]*fr[17]-1.677050983124842*vtSq[2]*fl[17]-2.130140840414079*vtSq[5]*fr[13]-2.381569860407206*vtSq[0]*fr[13]-2.130140840414079*vtSq[5]*fl[13]-2.381569860407206*vtSq[0]*fl[13]-2.381569860407206*vtSq[3]*fr[12]-2.381569860407206*vtSq[3]*fl[12]+1.875*vtSq[1]*fr[6]-1.875*vtSq[1]*fl[6]+1.677050983124842*fr[2]*vtSq[5]-1.677050983124842*fl[2]*vtSq[5]-2.381569860407206*vtSq[2]*fr[5]-2.381569860407206*vtSq[2]*fl[5]+1.875*fr[1]*vtSq[3]-1.875*fl[1]*vtSq[3]+1.875*fr[0]*vtSq[2]-1.875*fl[0]*vtSq[2]+1.875*vtSq[0]*fr[2]-1.875*vtSq[0]*fl[2])*rdv+alphaDrag[2]*(0.5590169943749475*favg[20]+0.223606797749979*favg[17]+0.4330127018922193*favg[5]+0.25*favg[0])-0.8660254037844386*fjump[13]+alphaDrag[0]*(0.4330127018922193*favg[13]+0.25*favg[2])+alphaDrag[5]*(0.3872983346207416*favg[13]+0.223606797749979*favg[2])+alphaDrag[3]*(0.4330127018922193*favg[12]+0.25*favg[1])+0.25*alphaDrag[1]*favg[6]-0.5*fjump[2]; 
  Ghat[3] = ((-2.381569860407206*vtSq[0]*fr[14])-2.381569860407206*vtSq[0]*fl[14]+1.875*vtSq[2]*fr[8]-1.875*vtSq[2]*fl[8]+1.875*vtSq[1]*fr[7]-1.875*vtSq[1]*fl[7]+1.875*vtSq[0]*fr[3]-1.875*vtSq[0]*fl[3])*rdv-0.8660254037844386*fjump[14]+alphaDrag[0]*(0.4330127018922193*favg[14]+0.25*favg[3])+0.25*alphaDrag[2]*favg[8]+0.25*alphaDrag[1]*favg[7]-0.5*fjump[3]; 
  Ghat[4] = ((-2.381569860407206*vtSq[0]*fr[15])-2.381569860407206*vtSq[0]*fl[15]+1.875*vtSq[2]*fr[10]-1.875*vtSq[2]*fl[10]+1.875*vtSq[1]*fr[9]-1.875*vtSq[1]*fl[9]+1.875*vtSq[0]*fr[4]-1.875*vtSq[0]*fl[4])*rdv-0.8660254037844386*fjump[15]+alphaDrag[0]*(0.4330127018922193*favg[15]+0.25*favg[4])+0.25*alphaDrag[2]*favg[10]+0.25*alphaDrag[1]*favg[9]-0.5*fjump[4]; 
  Ghat[6] = (1.341640786499874*vtSq[3]*fr[20]-1.341640786499874*vtSq[3]*fl[20]+1.677050983124842*vtSq[3]*fr[17]-1.677050983124842*vtSq[3]*fl[17]+1.677050983124842*vtSq[3]*fr[16]-1.677050983124842*vtSq[3]*fl[16]-2.381569860407206*vtSq[1]*fr[13]-2.381569860407206*vtSq[1]*fl[13]-2.381569860407206*vtSq[2]*fr[12]-2.381569860407206*vtSq[2]*fl[12]+1.677050983124842*vtSq[5]*fr[6]+1.677050983124842*vtSq[4]*fr[6]+1.875*vtSq[0]*fr[6]-1.677050983124842*vtSq[5]*fl[6]-1.677050983124842*vtSq[4]*fl[6]-1.875*vtSq[0]*fl[6]-2.381569860407206*vtSq[3]*fr[5]-2.381569860407206*vtSq[3]*fl[5]+1.875*fr[0]*vtSq[3]-1.875*fl[0]*vtSq[3]+1.875*fr[1]*vtSq[2]-1.875*fl[1]*vtSq[2]+1.875*vtSq[1]*fr[2]-1.875*vtSq[1]*fl[2])*rdv+alphaDrag[3]*(0.5590169943749475*favg[20]+0.223606797749979*favg[17]+0.223606797749979*favg[16]+0.4330127018922193*favg[5]+0.25*favg[0])+alphaDrag[1]*(0.4330127018922193*favg[13]+0.25*favg[2])+alphaDrag[2]*(0.4330127018922193*favg[12]+0.25*favg[1])-0.5*fjump[6]+0.223606797749979*alphaDrag[5]*favg[6]+0.223606797749979*alphaDrag[4]*favg[6]+0.25*alphaDrag[0]*favg[6]; 
  Ghat[7] = ((-2.381569860407206*vtSq[1]*fr[14])-2.381569860407206*vtSq[1]*fl[14]+1.875*vtSq[3]*fr[8]-1.875*vtSq[3]*fl[8]+1.677050983124842*vtSq[4]*fr[7]+1.875*vtSq[0]*fr[7]-1.677050983124842*vtSq[4]*fl[7]-1.875*vtSq[0]*fl[7]+1.875*vtSq[1]*fr[3]-1.875*vtSq[1]*fl[3])*rdv+alphaDrag[1]*(0.4330127018922193*favg[14]+0.25*favg[3])+0.25*alphaDrag[3]*favg[8]-0.5*fjump[7]+0.223606797749979*alphaDrag[4]*favg[7]+0.25*alphaDrag[0]*favg[7]; 
  Ghat[8] = ((-2.381569860407206*vtSq[2]*fr[14])-2.381569860407206*vtSq[2]*fl[14]+1.677050983124842*vtSq[5]*fr[8]+1.875*vtSq[0]*fr[8]-1.677050983124842*vtSq[5]*fl[8]-1.875*vtSq[0]*fl[8]+1.875*vtSq[3]*fr[7]-1.875*vtSq[3]*fl[7]+1.875*vtSq[2]*fr[3]-1.875*vtSq[2]*fl[3])*rdv+alphaDrag[2]*(0.4330127018922193*favg[14]+0.25*favg[3])-0.5*fjump[8]+0.223606797749979*alphaDrag[5]*favg[8]+0.25*alphaDrag[0]*favg[8]+0.25*alphaDrag[3]*favg[7]; 
  Ghat[9] = ((-2.381569860407206*vtSq[1]*fr[15])-2.381569860407206*vtSq[1]*fl[15]+1.875*vtSq[3]*fr[10]-1.875*vtSq[3]*fl[10]+1.677050983124842*vtSq[4]*fr[9]+1.875*vtSq[0]*fr[9]-1.677050983124842*vtSq[4]*fl[9]-1.875*vtSq[0]*fl[9]+1.875*vtSq[1]*fr[4]-1.875*vtSq[1]*fl[4])*rdv+alphaDrag[1]*(0.4330127018922193*favg[15]+0.25*favg[4])+0.25*alphaDrag[3]*favg[10]-0.5*fjump[9]+0.223606797749979*alphaDrag[4]*favg[9]+0.25*alphaDrag[0]*favg[9]; 
  Ghat[10] = ((-2.381569860407206*vtSq[2]*fr[15])-2.381569860407206*vtSq[2]*fl[15]+1.677050983124842*vtSq[5]*fr[10]+1.875*vtSq[0]*fr[10]-1.677050983124842*vtSq[5]*fl[10]-1.875*vtSq[0]*fl[10]+1.875*vtSq[3]*fr[9]-1.875*vtSq[3]*fl[9]+1.875*vtSq[2]*fr[4]-1.875*vtSq[2]*fl[4])*rdv+alphaDrag[2]*(0.4330127018922193*favg[15]+0.25*favg[4])-0.5*fjump[10]+0.223606797749979*alphaDrag[5]*favg[10]+0.25*alphaDrag[0]*favg[10]+0.25*alphaDrag[3]*favg[9]; 
  Ghat[11] = (1.875*vtSq[0]*fr[11]-1.875*vtSq[0]*fl[11])*rdv-0.5*fjump[11]+0.25*alphaDrag[0]*favg[11]; 
  Ghat[16] = (1.341640786499874*vtSq[4]*fr[20]-1.341640786499874*vtSq[4]*fl[20]+1.197893559374888*vtSq[4]*fr[16]+1.875*vtSq[0]*fr[16]-1.197893559374888*vtSq[4]*fl[16]-1.875*vtSq[0]*fl[16]-2.130140840414079*vtSq[1]*fr[12]-2.130140840414079*vtSq[1]*fl[12]+1.677050983124842*vtSq[3]*fr[6]-1.677050983124842*vtSq[3]*fl[6]-2.381569860407206*vtSq[4]*fr[5]-2.381569860407206*vtSq[4]*fl[5]+1.875*fr[0]*vtSq[4]-1.875*fl[0]*vtSq[4]+1.677050983124842*fr[1]*vtSq[1]-1.677050983124842*fl[1]*vtSq[1])*rdv+alphaDrag[4]*(0.5590169943749475*favg[20]+0.159719141249985*favg[16]+0.4330127018922193*favg[5]+0.25*favg[0])-0.5*fjump[16]+0.25*alphaDrag[0]*favg[16]+alphaDrag[1]*(0.3872983346207416*favg[12]+0.223606797749979*favg[1])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[17] = (1.341640786499874*vtSq[5]*fr[20]-1.341640786499874*vtSq[5]*fl[20]+1.197893559374888*vtSq[5]*fr[17]+1.875*vtSq[0]*fr[17]-1.197893559374888*vtSq[5]*fl[17]-1.875*vtSq[0]*fl[17]-2.130140840414079*vtSq[2]*fr[13]-2.130140840414079*vtSq[2]*fl[13]+1.677050983124842*vtSq[3]*fr[6]-1.677050983124842*vtSq[3]*fl[6]-2.381569860407206*fr[5]*vtSq[5]-2.381569860407206*fl[5]*vtSq[5]+1.875*fr[0]*vtSq[5]-1.875*fl[0]*vtSq[5]+1.677050983124842*fr[2]*vtSq[2]-1.677050983124842*fl[2]*vtSq[2])*rdv+alphaDrag[5]*(0.5590169943749475*favg[20]+0.159719141249985*favg[17]+0.4330127018922193*favg[5]+0.25*favg[0])-0.5*fjump[17]+0.25*alphaDrag[0]*favg[17]+alphaDrag[2]*(0.3872983346207416*favg[13]+0.223606797749979*favg[2])+0.223606797749979*alphaDrag[3]*favg[6]; 
  Ghat[18] = (1.875*vtSq[0]*fr[18]-1.875*vtSq[0]*fl[18])*rdv-0.5*fjump[18]+0.25*alphaDrag[0]*favg[18]; 
  Ghat[19] = (1.875*vtSq[0]*fr[19]-1.875*vtSq[0]*fl[19])*rdv-0.5*fjump[19]+0.25*alphaDrag[0]*favg[19]; 

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
  incr2[5] = vtSq[0]*(0.2118037767457181*(fr[20]+fl[20])-0.3046875*fr[5]+0.3046875*fl[5]+0.2165063509461096*(fr[0]+fl[0]))+0.2165063509461096*(vtSq[5]*(fr[17]+fl[17])+vtSq[4]*(fr[16]+fl[16]))+vtSq[2]*((-0.3046875*fr[13])+0.3046875*fl[13]+0.2165063509461096*(fr[2]+fl[2]))+vtSq[1]*((-0.3046875*fr[12])+0.3046875*fl[12]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*vtSq[3]*(fr[6]+fl[6]); 
  incr2[12] = vtSq[1]*(0.2118037767457181*(fr[20]+fl[20])+0.1936491673103708*(fr[16]+fl[16])-0.3046875*fr[5]+0.3046875*fl[5]+0.2165063509461096*(fr[0]+fl[0]))+vtSq[3]*((-0.3046875*fr[13])+0.3046875*fl[13]+0.2165063509461096*(fr[2]+fl[2]))+vtSq[4]*((-0.2725207847577868*fr[12])+0.2725207847577868*fl[12]+0.1936491673103708*(fr[1]+fl[1]))+vtSq[0]*((-0.3046875*fr[12])+0.3046875*fl[12]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*vtSq[2]*(fr[6]+fl[6]); 
  incr2[13] = vtSq[2]*(0.2118037767457181*(fr[20]+fl[20])+0.1936491673103708*(fr[17]+fl[17])-0.3046875*fr[5]+0.3046875*fl[5]+0.2165063509461096*(fr[0]+fl[0]))+vtSq[5]*((-0.2725207847577868*fr[13])+0.2725207847577868*fl[13]+0.1936491673103708*(fr[2]+fl[2]))+vtSq[0]*((-0.3046875*fr[13])+0.3046875*fl[13]+0.2165063509461096*(fr[2]+fl[2]))+vtSq[3]*((-0.3046875*fr[12])+0.3046875*fl[12]+0.2165063509461096*(fr[1]+fl[1]))+0.2165063509461096*vtSq[1]*(fr[6]+fl[6]); 
  incr2[14] = vtSq[0]*((-0.3046875*fr[14])+0.3046875*fl[14]+0.2165063509461096*(fr[3]+fl[3]))+0.2165063509461096*(vtSq[2]*(fr[8]+fl[8])+vtSq[1]*(fr[7]+fl[7])); 
  incr2[15] = vtSq[0]*((-0.3046875*fr[15])+0.3046875*fl[15]+0.2165063509461096*(fr[4]+fl[4]))+0.2165063509461096*(vtSq[2]*(fr[10]+fl[10])+vtSq[1]*(fr[9]+fl[9])); 
  incr2[20] = vtSq[0]*((-0.8203125*(fr[20]+fl[20]))+1.180049613297572*fr[5]-1.180049613297572*fl[5]-0.8385254915624212*(fr[0]+fl[0]))-0.8385254915624212*(vtSq[5]*(fr[17]+fl[17])+vtSq[4]*(fr[16]+fl[16]))+vtSq[2]*(1.180049613297572*fr[13]-1.180049613297572*fl[13]-0.8385254915624212*(fr[2]+fl[2]))+vtSq[1]*(1.180049613297572*fr[12]-1.180049613297572*fl[12]-0.8385254915624212*(fr[1]+fl[1]))-0.8385254915624212*vtSq[3]*(fr[6]+fl[6]); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr1[2]*rdv2nuR; 
  outr[3] += incr1[3]*rdv2nuR; 
  outr[4] += incr1[4]*rdv2nuR; 
  outr[5] += incr2[5]*rdvSq4nuR+incr1[5]*rdv2nuR; 
  outr[6] += incr1[6]*rdv2nuR; 
  outr[7] += incr1[7]*rdv2nuR; 
  outr[8] += incr1[8]*rdv2nuR; 
  outr[9] += incr1[9]*rdv2nuR; 
  outr[10] += incr1[10]*rdv2nuR; 
  outr[11] += incr1[11]*rdv2nuR; 
  outr[12] += incr2[12]*rdvSq4nuR+incr1[12]*rdv2nuR; 
  outr[13] += incr2[13]*rdvSq4nuR+incr1[13]*rdv2nuR; 
  outr[14] += incr2[14]*rdvSq4nuR+incr1[14]*rdv2nuR; 
  outr[15] += incr2[15]*rdvSq4nuR+incr1[15]*rdv2nuR; 
  outr[16] += incr1[16]*rdv2nuR; 
  outr[17] += incr1[17]*rdv2nuR; 
  outr[18] += incr1[18]*rdv2nuR; 
  outr[19] += incr1[19]*rdv2nuR; 
  outr[20] += incr2[20]*rdvSq4nuR+incr1[20]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += -1.0*incr1[2]*rdv2nuL; 
  outl[3] += -1.0*incr1[3]*rdv2nuL; 
  outl[4] += -1.0*incr1[4]*rdv2nuL; 
  outl[5] += incr1[5]*rdv2nuL-1.0*incr2[5]*rdvSq4nuL; 
  outl[6] += -1.0*incr1[6]*rdv2nuL; 
  outl[7] += -1.0*incr1[7]*rdv2nuL; 
  outl[8] += -1.0*incr1[8]*rdv2nuL; 
  outl[9] += -1.0*incr1[9]*rdv2nuL; 
  outl[10] += -1.0*incr1[10]*rdv2nuL; 
  outl[11] += -1.0*incr1[11]*rdv2nuL; 
  outl[12] += incr1[12]*rdv2nuL-1.0*incr2[12]*rdvSq4nuL; 
  outl[13] += incr1[13]*rdv2nuL-1.0*incr2[13]*rdvSq4nuL; 
  outl[14] += incr1[14]*rdv2nuL-1.0*incr2[14]*rdvSq4nuL; 
  outl[15] += incr1[15]*rdv2nuL-1.0*incr2[15]*rdvSq4nuL; 
  outl[16] += -1.0*incr1[16]*rdv2nuL; 
  outl[17] += -1.0*incr1[17]*rdv2nuL; 
  outl[18] += -1.0*incr1[18]*rdv2nuL; 
  outl[19] += -1.0*incr1[19]*rdv2nuL; 
  outl[20] += incr2[20]*rdvSq4nuL-1.0*incr1[20]*rdv2nuL; 

  return std::abs(0.5590169943749475*Uz[5]+wl[4]+0.5590169943749475*Uz[4]-0.5*Uz[0]); 
} 
