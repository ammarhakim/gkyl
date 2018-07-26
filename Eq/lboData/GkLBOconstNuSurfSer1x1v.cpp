#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf1x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[1*2]: bulk velocity (in 1 directions) times nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2l = 2.0/dxvl[1]; 
  double rdv2r = 2.0/dxvr[1]; 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double favg[4]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 

  double fjump[4]; 
  fjump[0] = nu*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nu*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nu*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nu*vMuMidMax*(fl[3]-(-1*fr[3])); 

  double drBar[2]; 
  drBar[0] = 1.414213562373095*wl[1]*nu-1.0*nuU[0]; 
  drBar[1] = -1.0*nuU[1]; 

  double alpha[4]; 
  alpha[0] = 1.414213562373095*nuVtSq[0]; 
  alpha[1] = 1.414213562373095*nuVtSq[1]; 

  double Gdiff[4]; 
  Gdiff[0] = (-1.082531754730548*(alpha[1]*(fr[3]+fl[3])+alpha[0]*(fr[2]+fl[2])))+alpha[1]*(1.125*fr[1]-1.125*fl[1])+alpha[0]*(1.125*fr[0]-1.125*fl[0]); 
  Gdiff[1] = (-1.082531754730548*(alpha[0]*(fr[3]+fl[3])+alpha[1]*(fr[2]+fl[2])))+alpha[0]*(1.125*fr[1]-1.125*fl[1])+(1.125*fr[0]-1.125*fl[0])*alpha[1]; 

  double Ghat[4]; 
  for(unsigned short int i=0; i<4; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(0.4330127018922193*dxvl[1]*favg[2]+0.25*favg[0]*dxvl[1])*nu+drBar[1]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]+drBar[0]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+(0.4330127018922193*dxvl[1]*favg[3]+0.25*dxvl[1]*favg[1])*nu-0.8660254037844386*fjump[3]+drBar[0]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[1]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[1]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 

  double incr2[4]; 
  incr2[2] = alpha[1]*(0.25*fl[3]-0.25*fr[3])+alpha[0]*(0.25*fl[2]-0.25*fr[2])+0.2165063509461096*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[3] = alpha[0]*(0.25*fl[3]-0.25*fr[3])+alpha[1]*(0.25*fl[2]-0.25*fr[2])+0.2165063509461096*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 

  outr[0] += incr1[0]*rdv2r; 
  outr[1] += incr1[1]*rdv2r; 
  outr[2] += incr2[2]*rdvSq4r+incr1[2]*rdv2r; 
  outr[3] += incr2[3]*rdvSq4r+incr1[3]*rdv2r; 

  outl[0] += -1.0*incr1[0]*rdv2l; 
  outl[1] += -1.0*incr1[1]*rdv2l; 
  outl[2] += incr1[2]*rdv2l-1.0*incr2[2]*rdvSq4l; 
  outl[3] += incr1[3]*rdv2l-1.0*incr2[3]*rdvSq4l; 

  const double vMuMid = wl[1]-(0.7071067811865475*nuU[0])/nu; 
  return std::abs(vMuMid); 
} 
double GkLBOconstNuSurf1x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[1*3]: bulk velocity (in 1 directions) times nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2l = 2.0/dxvl[1]; 
  double rdv2r = 2.0/dxvr[1]; 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

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
  fjump[0] = nu*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nu*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nu*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nu*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nu*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nu*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nu*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nu*vMuMidMax*(fl[7]-(1*fr[7])); 

  double drBar[3]; 
  drBar[0] = 1.414213562373095*wl[1]*nu-1.0*nuU[0]; 
  drBar[1] = -1.0*nuU[1]; 
  drBar[2] = -1.0*nuU[2]; 

  double alpha[8]; 
  alpha[0] = 1.414213562373095*nuVtSq[0]; 
  alpha[1] = 1.414213562373095*nuVtSq[1]; 
  alpha[4] = 1.414213562373095*nuVtSq[2]; 

  double Gdiff[8]; 
  Gdiff[0] = alpha[1]*(1.341640786499874*fr[7]-1.341640786499874*fl[7])-2.381569860407206*alpha[4]*(fr[6]+fl[6])+alpha[0]*(1.341640786499874*fr[5]-1.341640786499874*fl[5])+alpha[4]*(1.875*fr[4]-1.875*fl[4])-2.381569860407206*(alpha[1]*(fr[3]+fl[3])+alpha[0]*(fr[2]+fl[2]))+alpha[1]*(1.875*fr[1]-1.875*fl[1])+alpha[0]*(1.875*fr[0]-1.875*fl[0]); 
  Gdiff[1] = (1.2*alpha[4]+1.341640786499874*alpha[0])*fr[7]+((-1.2*alpha[4])-1.341640786499874*alpha[0])*fl[7]+alpha[1]*((-2.130140840414079*(fr[6]+fl[6]))+1.341640786499874*fr[5]-1.341640786499874*fl[5]+1.677050983124842*fr[4]-1.677050983124842*fl[4])+((-2.130140840414079*(fr[3]+fl[3]))+1.677050983124842*fr[1]-1.677050983124842*fl[1])*alpha[4]-2.381569860407206*(alpha[0]*(fr[3]+fl[3])+alpha[1]*(fr[2]+fl[2]))+alpha[0]*(1.875*fr[1]-1.875*fl[1])+(1.875*fr[0]-1.875*fl[0])*alpha[1]; 
  Gdiff[4] = alpha[1]*(1.2*fr[7]-1.2*fl[7])+((-1.521529171724342*alpha[4])-2.381569860407206*alpha[0])*fr[6]+((-1.521529171724342*alpha[4])-2.381569860407206*alpha[0])*fl[6]+alpha[4]*(1.341640786499874*fr[5]-1.341640786499874*fl[5])+(1.197893559374888*alpha[4]+1.875*alpha[0])*fr[4]+((-1.197893559374888*alpha[4])-1.875*alpha[0])*fl[4]+((-2.381569860407206*(fr[2]+fl[2]))+1.875*fr[0]-1.875*fl[0])*alpha[4]+alpha[1]*((-2.130140840414079*(fr[3]+fl[3]))+1.677050983124842*fr[1]-1.677050983124842*fl[1]); 

  double Ghat[8]; 
  for(unsigned short int i=0; i<8; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(0.5590169943749475*dxvl[1]*favg[5]+0.4330127018922193*dxvl[1]*favg[2]+0.25*favg[0]*dxvl[1])*nu+drBar[1]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[2]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])-1.118033988749895*fjump[5]+drBar[0]*(0.7905694150420947*favg[5]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+(0.5590169943749476*dxvl[1]*favg[7]+0.4330127018922193*dxvl[1]*favg[3]+0.25*dxvl[1]*favg[1])*nu-1.118033988749895*fjump[7]+drBar[0]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[2]*(0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])+drBar[1]*(0.5477225575051661*favg[6]+0.7905694150420947*favg[5]+0.3162277660168379*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[4] = Gdiff[4]*rdv+(0.4330127018922194*dxvl[1]*favg[6]+0.25*dxvl[1]*favg[4])*nu+drBar[1]*(0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])-0.8660254037844387*fjump[6]+drBar[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])+drBar[2]*(0.3912303982179757*favg[6]+0.7905694150420947*favg[5]+0.2258769757263128*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[4]; 

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

  outr[0] += incr1[0]*rdv2r; 
  outr[1] += incr1[1]*rdv2r; 
  outr[2] += incr2[2]*rdvSq4r+incr1[2]*rdv2r; 
  outr[3] += incr2[3]*rdvSq4r+incr1[3]*rdv2r; 
  outr[4] += incr1[4]*rdv2r; 
  outr[5] += incr2[5]*rdvSq4r+incr1[5]*rdv2r; 
  outr[6] += incr2[6]*rdvSq4r+incr1[6]*rdv2r; 
  outr[7] += incr2[7]*rdvSq4r+incr1[7]*rdv2r; 

  outl[0] += -1.0*incr1[0]*rdv2l; 
  outl[1] += -1.0*incr1[1]*rdv2l; 
  outl[2] += incr1[2]*rdv2l-1.0*incr2[2]*rdvSq4l; 
  outl[3] += incr1[3]*rdv2l-1.0*incr2[3]*rdvSq4l; 
  outl[4] += -1.0*incr1[4]*rdv2l; 
  outl[5] += incr2[5]*rdvSq4l-1.0*incr1[5]*rdv2l; 
  outl[6] += incr1[6]*rdv2l-1.0*incr2[6]*rdvSq4l; 
  outl[7] += incr2[7]*rdvSq4l-1.0*incr1[7]*rdv2l; 

  const double vMuMid = (0.7905694150420947*nuU[2])/nu-(0.7071067811865475*nuU[0])/nu+wl[1]; 
  return std::abs(vMuMid); 
} 
