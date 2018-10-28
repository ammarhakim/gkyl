#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf1x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates. 
  // dxv[2]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[1*2]:   bulk velocity (in 1 directions). 
  // vtSq[2]:   thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nuL = 2.0*nu/dxvl[1]; 
  double rdv2nuR = 2.0*nu/dxvr[1]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double favg[4]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 

  double fjump[4]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(-1*fr[3])); 

  double drBar[2]; 
  drBar[0] = 1.414213562373095*wl[1]-1.0*u[0]; 
  drBar[1] = -1.0*u[1]; 

  double alpha[4]; 
  alpha[0] = 1.414213562373095*vtSq[0]; 
  alpha[1] = 1.414213562373095*vtSq[1]; 

  double Gdiff[4]; 
  Gdiff[0] = (-1.082531754730548*(alpha[1]*(fr[3]+fl[3])+alpha[0]*(fr[2]+fl[2])))+alpha[1]*(1.125*fr[1]-1.125*fl[1])+alpha[0]*(1.125*fr[0]-1.125*fl[0]); 
  Gdiff[1] = (-1.082531754730548*(alpha[0]*(fr[3]+fl[3])+alpha[1]*(fr[2]+fl[2])))+alpha[0]*(1.125*fr[1]-1.125*fl[1])+(1.125*fr[0]-1.125*fl[0])*alpha[1]; 

  double Ghat[4]; 
  for(unsigned short int i=0; i<4; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+drBar[1]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]+0.4330127018922193*dxvl[1]*favg[2]+drBar[0]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.25*favg[0]*dxvl[1]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv-0.8660254037844386*fjump[3]+0.4330127018922193*dxvl[1]*favg[3]+drBar[0]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[1]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[1]+0.25*dxvl[1]*favg[1]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 

  double incr2[4]; 
  incr2[2] = alpha[1]*(0.25*fl[3]-0.25*fr[3])+alpha[0]*(0.25*fl[2]-0.25*fr[2])+0.2165063509461096*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[3] = alpha[0]*(0.25*fl[3]-0.25*fr[3])+alpha[1]*(0.25*fl[2]-0.25*fr[2])+0.2165063509461096*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr2[2]*rdvSq4nuR+incr1[2]*rdv2nuR; 
  outr[3] += incr2[3]*rdvSq4nuR+incr1[3]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += incr1[2]*rdv2nuL-1.0*incr2[2]*rdvSq4nuL; 
  outl[3] += incr1[3]*rdv2nuL-1.0*incr2[3]*rdvSq4nuL; 

  const double vMuMid = wl[1]-0.7071067811865475*u[0]; 
  return std::abs(vMuMid); 
} 
double GkLBOconstNuSurf1x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates. 
  // dxv[2]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[1*3]:   bulk velocity (in 1 directions). 
  // vtSq[3]:   thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nuL = 2.0*nu/dxvl[1]; 
  double rdv2nuR = 2.0*nu/dxvr[1]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double favg[8]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 

  double fjump[8]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = vMuMidMax*(fl[7]-(1*fr[7])); 

  double drBar[3]; 
  drBar[0] = 1.414213562373095*wl[1]-1.0*u[0]; 
  drBar[1] = -1.0*u[1]; 
  drBar[2] = -1.0*u[2]; 

  double alpha[8]; 
  alpha[0] = 1.414213562373095*vtSq[0]; 
  alpha[1] = 1.414213562373095*vtSq[1]; 
  alpha[4] = 1.414213562373095*vtSq[2]; 

  double Gdiff[8]; 
  Gdiff[0] = alpha[1]*(1.341640786499874*fr[7]-1.341640786499874*fl[7])-2.381569860407206*alpha[4]*(fr[6]+fl[6])+alpha[0]*(1.341640786499874*fr[5]-1.341640786499874*fl[5])+alpha[4]*(1.875*fr[4]-1.875*fl[4])-2.381569860407206*(alpha[1]*(fr[3]+fl[3])+alpha[0]*(fr[2]+fl[2]))+alpha[1]*(1.875*fr[1]-1.875*fl[1])+alpha[0]*(1.875*fr[0]-1.875*fl[0]); 
  Gdiff[1] = (1.2*alpha[4]+1.341640786499874*alpha[0])*fr[7]+((-1.2*alpha[4])-1.341640786499874*alpha[0])*fl[7]+alpha[1]*((-2.130140840414079*(fr[6]+fl[6]))+1.341640786499874*fr[5]-1.341640786499874*fl[5]+1.677050983124842*fr[4]-1.677050983124842*fl[4])+((-2.130140840414079*(fr[3]+fl[3]))+1.677050983124842*fr[1]-1.677050983124842*fl[1])*alpha[4]-2.381569860407206*(alpha[0]*(fr[3]+fl[3])+alpha[1]*(fr[2]+fl[2]))+alpha[0]*(1.875*fr[1]-1.875*fl[1])+(1.875*fr[0]-1.875*fl[0])*alpha[1]; 
  Gdiff[4] = alpha[1]*(1.2*fr[7]-1.2*fl[7])+((-1.521529171724342*alpha[4])-2.381569860407206*alpha[0])*fr[6]+((-1.521529171724342*alpha[4])-2.381569860407206*alpha[0])*fl[6]+alpha[4]*(1.341640786499874*fr[5]-1.341640786499874*fl[5])+(1.197893559374888*alpha[4]+1.875*alpha[0])*fr[4]+((-1.197893559374888*alpha[4])-1.875*alpha[0])*fl[4]+((-2.381569860407206*(fr[2]+fl[2]))+1.875*fr[0]-1.875*fl[0])*alpha[4]+alpha[1]*((-2.130140840414079*(fr[3]+fl[3]))+1.677050983124842*fr[1]-1.677050983124842*fl[1]); 

  double Ghat[8]; 
  for(unsigned short int i=0; i<8; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+drBar[1]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[2]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])-1.118033988749895*fjump[5]+0.5590169943749475*dxvl[1]*favg[5]+drBar[0]*(0.7905694150420947*favg[5]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[2]+0.4330127018922193*dxvl[1]*favg[2]+0.25*favg[0]*dxvl[1]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv-1.118033988749895*fjump[7]+0.5590169943749476*dxvl[1]*favg[7]+drBar[0]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[2]*(0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])+drBar[1]*(0.5477225575051661*favg[6]+0.7905694150420947*favg[5]+0.3162277660168379*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[3]+0.4330127018922193*dxvl[1]*favg[3]-0.5*fjump[1]+0.25*dxvl[1]*favg[1]; 
  Ghat[4] = Gdiff[4]*rdv+drBar[1]*(0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])-0.8660254037844387*fjump[6]+0.4330127018922194*dxvl[1]*favg[6]+drBar[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])+drBar[2]*(0.3912303982179757*favg[6]+0.7905694150420947*favg[5]+0.2258769757263128*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[4]+0.25*dxvl[1]*favg[4]; 

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -1.118033988749895*Ghat[0]; 
  incr1[6] = 0.8660254037844387*Ghat[4]; 
  incr1[7] = -1.118033988749895*Ghat[1]; 

  double incr2[8]; 
  incr2[2] = 0.2118037767457181*alpha[1]*(fr[7]+fl[7])+alpha[4]*(0.3046875*fl[6]-0.3046875*fr[6])+0.2118037767457181*alpha[0]*(fr[5]+fl[5])+0.2165063509461096*alpha[4]*(fr[4]+fl[4])+alpha[1]*(0.3046875*fl[3]-0.3046875*fr[3])+alpha[0]*(0.3046875*fl[2]-0.3046875*fr[2])+0.2165063509461096*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[3] = (0.1894430570778459*alpha[4]+0.2118037767457181*alpha[0])*fr[7]+(0.1894430570778459*alpha[4]+0.2118037767457181*alpha[0])*fl[7]+alpha[1]*((-0.2725207847577868*fr[6])+0.2725207847577868*fl[6]+0.2118037767457181*(fr[5]+fl[5])+0.1936491673103708*(fr[4]+fl[4]))+((-0.2725207847577868*fr[3])+0.2725207847577868*fl[3]+0.1936491673103708*(fr[1]+fl[1]))*alpha[4]+alpha[0]*(0.3046875*fl[3]-0.3046875*fr[3])+alpha[1]*(0.3046875*fl[2]-0.3046875*fr[2])+0.2165063509461096*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 
  incr2[5] = (-0.8203125*alpha[1]*(fr[7]+fl[7]))+alpha[4]*(1.180049613297572*fr[6]-1.180049613297572*fl[6])-0.8203125*alpha[0]*(fr[5]+fl[5])-0.8385254915624212*alpha[4]*(fr[4]+fl[4])+alpha[1]*(1.180049613297572*fr[3]-1.180049613297572*fl[3])+alpha[0]*(1.180049613297572*fr[2]-1.180049613297572*fl[2])-0.8385254915624212*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[6] = 0.1894430570778459*alpha[1]*(fr[7]+fl[7])+((-0.1946577033984192*alpha[4])-0.3046875*alpha[0])*fr[6]+(0.1946577033984192*alpha[4]+0.3046875*alpha[0])*fl[6]+0.2118037767457181*alpha[4]*(fr[5]+fl[5])+(0.138320833793122*alpha[4]+0.2165063509461096*alpha[0])*fr[4]+(0.138320833793122*alpha[4]+0.2165063509461096*alpha[0])*fl[4]+((-0.3046875*fr[2])+0.3046875*fl[2]+0.2165063509461096*(fr[0]+fl[0]))*alpha[4]+alpha[1]*((-0.2725207847577868*fr[3])+0.2725207847577868*fl[3]+0.1936491673103708*(fr[1]+fl[1])); 
  incr2[7] = ((-0.7337098051171185*alpha[4])-0.8203125*alpha[0])*fr[7]+((-0.7337098051171185*alpha[4])-0.8203125*alpha[0])*fl[7]+alpha[1]*(1.055468460862284*fr[6]-1.055468460862284*fl[6]-0.8203125*(fr[5]+fl[5])-0.75*(fr[4]+fl[4]))+(1.055468460862284*fr[3]-1.055468460862284*fl[3]-0.75*(fr[1]+fl[1]))*alpha[4]+alpha[0]*(1.180049613297572*fr[3]-1.180049613297572*fl[3])+alpha[1]*(1.180049613297572*fr[2]-1.180049613297572*fl[2])-0.8385254915624212*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr2[2]*rdvSq4nuR+incr1[2]*rdv2nuR; 
  outr[3] += incr2[3]*rdvSq4nuR+incr1[3]*rdv2nuR; 
  outr[4] += incr1[4]*rdv2nuR; 
  outr[5] += incr2[5]*rdvSq4nuR+incr1[5]*rdv2nuR; 
  outr[6] += incr2[6]*rdvSq4nuR+incr1[6]*rdv2nuR; 
  outr[7] += incr2[7]*rdvSq4nuR+incr1[7]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += incr1[2]*rdv2nuL-1.0*incr2[2]*rdvSq4nuL; 
  outl[3] += incr1[3]*rdv2nuL-1.0*incr2[3]*rdvSq4nuL; 
  outl[4] += -1.0*incr1[4]*rdv2nuL; 
  outl[5] += incr2[5]*rdvSq4nuL-1.0*incr1[5]*rdv2nuL; 
  outl[6] += incr1[6]*rdv2nuL-1.0*incr2[6]*rdvSq4nuL; 
  outl[7] += incr2[7]*rdvSq4nuL-1.0*incr1[7]*rdv2nuL; 

  const double vMuMid = 0.7905694150420947*u[2]+wl[1]-0.7071067811865475*u[0]; 
  return std::abs(vMuMid); 
} 
double GkLBOconstNuSurf1x1vSer_Vpar_P3(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates. 
  // dxv[2]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[1*4]:   bulk velocity (in 1 directions). 
  // vtSq[4]:   thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nuL = 2.0*nu/dxvl[1]; 
  double rdv2nuR = 2.0*nu/dxvr[1]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double favg[12]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 

  double fjump[12]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = vMuMidMax*(fl[9]-(-1*fr[9])); 
  fjump[10] = vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = vMuMidMax*(fl[11]-(-1*fr[11])); 

  double drBar[4]; 
  drBar[0] = 1.414213562373095*wl[1]-1.0*u[0]; 
  drBar[1] = -1.0*u[1]; 
  drBar[2] = -1.0*u[2]; 
  drBar[3] = -1.0*u[3]; 

  double alpha[12]; 
  alpha[0] = 1.414213562373095*vtSq[0]; 
  alpha[1] = 1.414213562373095*vtSq[1]; 
  alpha[4] = 1.414213562373095*vtSq[2]; 
  alpha[8] = 1.414213562373095*vtSq[3]; 

  double Gdiff[12]; 
  Gdiff[0] = (-1.36421551976768*alpha[1]*(fr[11]+fl[11]))-3.87005102316171*alpha[8]*(fr[10]+fl[10])-1.36421551976768*alpha[0]*(fr[9]+fl[9])+alpha[8]*(2.734375*fr[8]-2.734375*fl[8])+alpha[1]*(3.109532031210645*fr[7]-3.109532031210645*fl[7])-3.87005102316171*alpha[4]*(fr[6]+fl[6])+alpha[0]*(3.109532031210645*fr[5]-3.109532031210645*fl[5])+alpha[4]*(2.734375*fr[4]-2.734375*fl[4])-3.87005102316171*(alpha[1]*(fr[3]+fl[3])+alpha[0]*(fr[2]+fl[2]))+alpha[1]*(2.734375*fr[1]-2.734375*fl[1])+alpha[0]*(2.734375*fr[0]-2.734375*fl[0]); 
  Gdiff[1] = ((-1.220191455264296*alpha[4])-1.36421551976768*alpha[0])*fr[11]+((-1.220191455264296*alpha[4])-1.36421551976768*alpha[0])*fl[11]-3.399104768236252*alpha[4]*(fr[10]+fl[10])-1.36421551976768*alpha[1]*(fr[9]+fl[9])+alpha[4]*(2.401629085771781*fr[8]-2.401629085771781*fl[8])+((-3.399104768236252*(fr[6]+fl[6]))+2.401629085771781*fr[4]-2.401629085771781*fl[4])*alpha[8]+(2.78125*alpha[4]+3.109532031210645*alpha[0])*fr[7]+((-2.78125*alpha[4])-3.109532031210645*alpha[0])*fl[7]+alpha[1]*((-3.461478865672879*(fr[6]+fl[6]))+3.109532031210645*fr[5]-3.109532031210645*fl[5]+2.445699350390396*fr[4]-2.445699350390396*fl[4])+((-3.461478865672879*(fr[3]+fl[3]))+2.445699350390396*fr[1]-2.445699350390396*fl[1])*alpha[4]-3.87005102316171*(alpha[0]*(fr[3]+fl[3])+alpha[1]*(fr[2]+fl[2]))+alpha[0]*(2.734375*fr[1]-2.734375*fl[1])+(2.734375*fr[0]-2.734375*fl[0])*alpha[1]; 
  Gdiff[4] = ((-1.198204222732919*alpha[8])-1.220191455264296*alpha[1])*fr[11]+((-1.198204222732919*alpha[8])-1.220191455264296*alpha[1])*fl[11]+((-2.307652577115253*alpha[8])-3.399104768236252*alpha[1])*fr[10]+((-2.307652577115253*alpha[8])-3.399104768236252*alpha[1])*fl[10]-1.36421551976768*alpha[4]*(fr[9]+fl[9])+(1.630466233593597*alpha[8]+2.401629085771781*alpha[1])*fr[8]+((-1.630466233593597*alpha[8])-2.401629085771781*alpha[1])*fl[8]+(2.731133282484842*fr[7]-2.731133282484842*fl[7]-3.399104768236252*(fr[3]+fl[3])+2.401629085771781*fr[1]-2.401629085771781*fl[1])*alpha[8]+alpha[1]*(2.78125*fr[7]-2.78125*fl[7])+((-2.472484904052056*alpha[4])-3.87005102316171*alpha[0])*fr[6]+((-2.472484904052056*alpha[4])-3.87005102316171*alpha[0])*fl[6]+alpha[4]*(3.109532031210645*fr[5]-3.109532031210645*fl[5])+(1.746928107421711*alpha[4]+2.734375*alpha[0])*fr[4]+((-1.746928107421711*alpha[4])-2.734375*alpha[0])*fl[4]+((-3.87005102316171*(fr[2]+fl[2]))+2.734375*fr[0]-2.734375*fl[0])*alpha[4]+alpha[1]*((-3.461478865672879*(fr[3]+fl[3]))+2.445699350390396*fr[1]-2.445699350390396*fl[1]); 
  Gdiff[8] = (-1.198204222732919*alpha[4]*(fr[11]+fl[11]))+((-2.307652577115253*alpha[4])-3.87005102316171*alpha[0])*fr[10]+((-2.307652577115253*alpha[4])-3.87005102316171*alpha[0])*fl[10]-1.36421551976768*alpha[8]*(fr[9]+fl[9])+(1.630466233593597*alpha[4]+2.734375*alpha[0])*fr[8]+((-1.630466233593597*alpha[4])-2.734375*alpha[0])*fl[8]+((-2.307652577115253*(fr[6]+fl[6]))+3.109532031210645*fr[5]-3.109532031210645*fl[5]+1.630466233593597*fr[4]-1.630466233593597*fl[4]-3.87005102316171*(fr[2]+fl[2])+2.734375*fr[0]-2.734375*fl[0])*alpha[8]+alpha[4]*(2.731133282484842*fr[7]-2.731133282484842*fl[7])+alpha[1]*((-3.399104768236252*(fr[6]+fl[6]))+2.401629085771781*fr[4]-2.401629085771781*fl[4])+((-3.399104768236252*(fr[3]+fl[3]))+2.401629085771781*fr[1]-2.401629085771781*fl[1])*alpha[4]; 

  double Ghat[12]; 
  for(unsigned short int i=0; i<12; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+drBar[1]*(0.9354143466934852*favg[11]+0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[3]*(0.6123724356957942*favg[10]+0.3535533905932737*favg[8])-1.322875655532295*fjump[9]+0.6614378277661477*dxvl[1]*favg[9]+drBar[0]*(0.9354143466934851*favg[9]+0.7905694150420947*favg[5]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+drBar[2]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])-1.118033988749895*fjump[5]+0.5590169943749475*dxvl[1]*favg[5]-0.8660254037844386*fjump[2]+0.4330127018922193*dxvl[1]*favg[2]+0.25*favg[0]*dxvl[1]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv-1.322875655532295*fjump[11]+0.6614378277661477*dxvl[1]*favg[11]+drBar[0]*(0.9354143466934852*favg[11]+0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[2]*(0.8366600265340755*favg[11]+0.5378528742004768*favg[10]+0.3105295017040592*favg[8]+0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])+drBar[1]*(0.9354143466934851*favg[9]+0.5477225575051661*favg[6]+0.7905694150420947*favg[5]+0.3162277660168379*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-1.118033988749895*fjump[7]+0.5590169943749476*dxvl[1]*favg[7]+drBar[3]*(0.5378528742004769*favg[6]+0.3105295017040592*favg[4])-0.8660254037844386*fjump[3]+0.4330127018922193*dxvl[1]*favg[3]-0.5*fjump[1]+0.25*dxvl[1]*favg[1]; 
  Ghat[4] = Gdiff[4]*rdv+drBar[1]*(0.8366600265340755*favg[11]+0.5378528742004768*favg[10]+0.3105295017040592*favg[8]+0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])+drBar[3]*(0.8215838362577488*favg[11]+0.3651483716701107*favg[10]+0.210818510677892*favg[8]+0.6943650748294133*favg[7]+0.537852874200477*favg[3]+0.3105295017040592*favg[1])+drBar[2]*(0.9354143466934851*favg[9]+0.3912303982179757*favg[6]+0.7905694150420947*favg[5]+0.2258769757263128*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844387*fjump[6]+0.4330127018922194*dxvl[1]*favg[6]+drBar[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])-0.5*fjump[4]+0.25*dxvl[1]*favg[4]; 
  Ghat[8] = Gdiff[8]*rdv+drBar[2]*(0.8215838362577488*favg[11]+0.3651483716701107*favg[10]+0.210818510677892*favg[8]+0.6943650748294133*favg[7]+0.537852874200477*favg[3]+0.3105295017040592*favg[1])-0.8660254037844386*fjump[10]+0.4330127018922193*dxvl[1]*favg[10]+drBar[0]*(0.6123724356957942*favg[10]+0.3535533905932737*favg[8])+drBar[3]*(0.9354143466934851*favg[9]+0.3651483716701107*favg[6]+0.7905694150420947*favg[5]+0.210818510677892*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[8]+0.25*dxvl[1]*favg[8]+drBar[1]*(0.5378528742004769*favg[6]+0.3105295017040592*favg[4]); 

  double incr1[12]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -1.118033988749895*Ghat[0]; 
  incr1[6] = 0.8660254037844387*Ghat[4]; 
  incr1[7] = -1.118033988749895*Ghat[1]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = 1.322875655532295*Ghat[0]; 
  incr1[10] = 0.8660254037844386*Ghat[8]; 
  incr1[11] = 1.322875655532295*Ghat[1]; 

  double incr2[12]; 
  incr2[2] = alpha[1]*(0.1636634176769943*fl[11]-0.1636634176769943*fr[11])+alpha[8]*(0.328125*fl[10]-0.328125*fr[10])+alpha[0]*(0.1636634176769943*fl[9]-0.1636634176769943*fr[9])+0.2165063509461096*alpha[8]*(fr[8]+fl[8])+0.3025768239224545*alpha[1]*(fr[7]+fl[7])+alpha[4]*(0.328125*fl[6]-0.328125*fr[6])+0.3025768239224545*alpha[0]*(fr[5]+fl[5])+0.2165063509461096*alpha[4]*(fr[4]+fl[4])+alpha[1]*(0.328125*fl[3]-0.328125*fr[3])+alpha[0]*(0.328125*fl[2]-0.328125*fr[2])+0.2165063509461096*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[3] = ((-0.14638501094228*alpha[4])-0.1636634176769943*alpha[0])*fr[11]+(0.14638501094228*alpha[4]+0.1636634176769943*alpha[0])*fl[11]+alpha[4]*(0.2881954902926136*fl[10]-0.2881954902926136*fr[10])+alpha[1]*(0.1636634176769943*fl[9]-0.1636634176769943*fr[9])+0.1901597073139162*alpha[4]*(fr[8]+fl[8])+((-0.2881954902926136*fr[6])+0.2881954902926136*fl[6]+0.1901597073139162*(fr[4]+fl[4]))*alpha[8]+(0.270632938682637*alpha[4]+0.3025768239224545*alpha[0])*fr[7]+(0.270632938682637*alpha[4]+0.3025768239224545*alpha[0])*fl[7]+alpha[1]*((-0.2934839220468474*fr[6])+0.2934839220468474*fl[6]+0.3025768239224545*(fr[5]+fl[5])+0.1936491673103708*(fr[4]+fl[4]))+((-0.2934839220468474*fr[3])+0.2934839220468474*fl[3]+0.1936491673103708*(fr[1]+fl[1]))*alpha[4]+alpha[0]*(0.328125*fl[3]-0.328125*fr[3])+alpha[1]*(0.328125*fl[2]-0.328125*fr[2])+0.2165063509461096*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 
  incr2[5] = alpha[1]*(0.6338656910463875*fr[11]-0.6338656910463875*fl[11])+alpha[8]*(1.270822660474308*fr[10]-1.270822660474308*fl[10])+alpha[0]*(0.6338656910463875*fr[9]-0.6338656910463875*fl[9])-0.8385254915624212*alpha[8]*(fr[8]+fl[8])-1.171875*alpha[1]*(fr[7]+fl[7])+alpha[4]*(1.270822660474308*fr[6]-1.270822660474308*fl[6])-1.171875*alpha[0]*(fr[5]+fl[5])-0.8385254915624212*alpha[4]*(fr[4]+fl[4])+alpha[1]*(1.270822660474308*fr[3]-1.270822660474308*fl[3])+alpha[0]*(1.270822660474308*fr[2]-1.270822660474308*fl[2])-0.8385254915624212*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[6] = ((-0.1437472271249865*alpha[8])-0.14638501094228*alpha[1])*fr[11]+(0.1437472271249865*alpha[8]+0.14638501094228*alpha[1])*fl[11]+((-0.1956559480312316*alpha[8])-0.2881954902926136*alpha[1])*fr[10]+(0.1956559480312316*alpha[8]+0.2881954902926136*alpha[1])*fl[10]+alpha[4]*(0.1636634176769943*fl[9]-0.1636634176769943*fr[9])+(0.1290994448735806*alpha[8]+0.1901597073139162*alpha[1])*fr[8]+(0.1290994448735806*alpha[8]+0.1901597073139162*alpha[1])*fl[8]+(0.2657562700846129*(fr[7]+fl[7])-0.2881954902926136*fr[3]+0.2881954902926136*fl[3]+0.1901597073139162*(fr[1]+fl[1]))*alpha[8]+0.270632938682637*alpha[1]*(fr[7]+fl[7])+((-0.2096313728906053*alpha[4])-0.328125*alpha[0])*fr[6]+(0.2096313728906053*alpha[4]+0.328125*alpha[0])*fl[6]+0.3025768239224545*alpha[4]*(fr[5]+fl[5])+(0.138320833793122*alpha[4]+0.2165063509461096*alpha[0])*fr[4]+(0.138320833793122*alpha[4]+0.2165063509461096*alpha[0])*fl[4]+((-0.328125*fr[2])+0.328125*fl[2]+0.2165063509461096*(fr[0]+fl[0]))*alpha[4]+alpha[1]*((-0.2934839220468474*fr[3])+0.2934839220468474*fl[3]+0.1936491673103708*(fr[1]+fl[1])); 
  incr2[7] = (0.5669467095138407*alpha[4]+0.6338656910463875*alpha[0])*fr[11]+((-0.5669467095138407*alpha[4])-0.6338656910463875*alpha[0])*fl[11]+alpha[4]*(1.116176334355374*fr[10]-1.116176334355374*fl[10])+alpha[1]*(0.6338656910463875*fr[9]-0.6338656910463875*fl[9])-0.736485379546474*alpha[4]*(fr[8]+fl[8])+(1.116176334355374*fr[6]-1.116176334355374*fl[6]-0.736485379546474*(fr[4]+fl[4]))*alpha[8]+((-1.048156864453027*alpha[4])-1.171875*alpha[0])*fr[7]+((-1.048156864453027*alpha[4])-1.171875*alpha[0])*fl[7]+alpha[1]*(1.136658342467076*fr[6]-1.136658342467076*fl[6]-1.171875*(fr[5]+fl[5])-0.75*(fr[4]+fl[4]))+(1.136658342467076*fr[3]-1.136658342467076*fl[3]-0.75*(fr[1]+fl[1]))*alpha[4]+alpha[0]*(1.270822660474308*fr[3]-1.270822660474308*fl[3])+alpha[1]*(1.270822660474308*fr[2]-1.270822660474308*fl[2])-0.8385254915624212*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 
  incr2[9] = alpha[1]*(1.5*fl[11]-1.5*fr[11])+alpha[8]*(3.00731529981477*fl[10]-3.00731529981477*fr[10])+alpha[0]*(1.5*fl[9]-1.5*fr[9])+1.984313483298443*alpha[8]*(fr[8]+fl[8])+2.773162398327945*alpha[1]*(fr[7]+fl[7])+alpha[4]*(3.00731529981477*fl[6]-3.00731529981477*fr[6])+2.773162398327945*alpha[0]*(fr[5]+fl[5])+1.984313483298443*alpha[4]*(fr[4]+fl[4])+alpha[1]*(3.00731529981477*fl[3]-3.00731529981477*fr[3])+alpha[0]*(3.00731529981477*fl[2]-3.00731529981477*fr[2])+1.984313483298443*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[10] = alpha[4]*(0.1437472271249865*fl[11]-0.1437472271249865*fr[11])+((-0.1956559480312316*alpha[4])-0.328125*alpha[0])*fr[10]+(0.1956559480312316*alpha[4]+0.328125*alpha[0])*fl[10]+alpha[8]*(0.1636634176769943*fl[9]-0.1636634176769943*fr[9])+(0.1290994448735806*alpha[4]+0.2165063509461096*alpha[0])*fr[8]+(0.1290994448735806*alpha[4]+0.2165063509461096*alpha[0])*fl[8]+((-0.1956559480312316*fr[6])+0.1956559480312316*fl[6]+0.3025768239224545*(fr[5]+fl[5])+0.1290994448735806*(fr[4]+fl[4])-0.328125*fr[2]+0.328125*fl[2]+0.2165063509461096*(fr[0]+fl[0]))*alpha[8]+0.2657562700846129*alpha[4]*(fr[7]+fl[7])+alpha[1]*((-0.2881954902926136*fr[6])+0.2881954902926136*fl[6]+0.1901597073139162*(fr[4]+fl[4]))+((-0.2881954902926136*fr[3])+0.2881954902926136*fl[3]+0.1901597073139162*(fr[1]+fl[1]))*alpha[4]; 
  incr2[11] = ((-1.341640786499874*alpha[4])-1.5*alpha[0])*fr[11]+(1.341640786499874*alpha[4]+1.5*alpha[0])*fl[11]+alpha[4]*(2.641355298421626*fl[10]-2.641355298421626*fr[10])+alpha[1]*(1.5*fl[9]-1.5*fr[9])+1.742842505793337*alpha[4]*(fr[8]+fl[8])+((-2.641355298421626*fr[6])+2.641355298421626*fl[6]+1.742842505793337*(fr[4]+fl[4]))*alpha[8]+(2.480391854123054*alpha[4]+2.773162398327945*alpha[0])*fr[7]+(2.480391854123054*alpha[4]+2.773162398327945*alpha[0])*fl[7]+alpha[1]*((-2.689824576064395*fr[6])+2.689824576064395*fl[6]+2.773162398327945*(fr[5]+fl[5])+1.774823934929885*(fr[4]+fl[4]))+((-2.689824576064395*fr[3])+2.689824576064395*fl[3]+1.774823934929885*(fr[1]+fl[1]))*alpha[4]+alpha[0]*(3.00731529981477*fl[3]-3.00731529981477*fr[3])+alpha[1]*(3.00731529981477*fl[2]-3.00731529981477*fr[2])+1.984313483298443*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr2[2]*rdvSq4nuR+incr1[2]*rdv2nuR; 
  outr[3] += incr2[3]*rdvSq4nuR+incr1[3]*rdv2nuR; 
  outr[4] += incr1[4]*rdv2nuR; 
  outr[5] += incr2[5]*rdvSq4nuR+incr1[5]*rdv2nuR; 
  outr[6] += incr2[6]*rdvSq4nuR+incr1[6]*rdv2nuR; 
  outr[7] += incr2[7]*rdvSq4nuR+incr1[7]*rdv2nuR; 
  outr[8] += incr1[8]*rdv2nuR; 
  outr[9] += incr2[9]*rdvSq4nuR+incr1[9]*rdv2nuR; 
  outr[10] += incr2[10]*rdvSq4nuR+incr1[10]*rdv2nuR; 
  outr[11] += incr2[11]*rdvSq4nuR+incr1[11]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += incr1[2]*rdv2nuL-1.0*incr2[2]*rdvSq4nuL; 
  outl[3] += incr1[3]*rdv2nuL-1.0*incr2[3]*rdvSq4nuL; 
  outl[4] += -1.0*incr1[4]*rdv2nuL; 
  outl[5] += incr2[5]*rdvSq4nuL-1.0*incr1[5]*rdv2nuL; 
  outl[6] += incr1[6]*rdv2nuL-1.0*incr2[6]*rdvSq4nuL; 
  outl[7] += incr2[7]*rdvSq4nuL-1.0*incr1[7]*rdv2nuL; 
  outl[8] += -1.0*incr1[8]*rdv2nuL; 
  outl[9] += incr1[9]*rdv2nuL-1.0*incr2[9]*rdvSq4nuL; 
  outl[10] += incr1[10]*rdv2nuL-1.0*incr2[10]*rdvSq4nuL; 
  outl[11] += incr1[11]*rdv2nuL-1.0*incr2[11]*rdvSq4nuL; 

  const double vMuMid = 0.7905694150420947*u[2]+wl[1]-0.7071067811865475*u[0]; 
  return std::abs(vMuMid); 
} 
