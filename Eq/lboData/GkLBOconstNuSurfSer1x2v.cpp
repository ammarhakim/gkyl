#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf1x2vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*2]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
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
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 

  double fjump[8]; 
  fjump[0] = nu*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nu*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nu*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nu*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nu*vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = nu*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nu*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nu*vMuMidMax*(fl[7]-(-1*fr[7])); 

  double drBar[2]; 
  drBar[0] = 1.414213562373095*wl[1]*nu-1.0*nuU[0]; 
  drBar[1] = -1.0*nuU[1]; 

  double alpha[8]; 
  alpha[0] = 2.0*nuVtSq[0]; 
  alpha[1] = 2.0*nuVtSq[1]; 

  double Gdiff[8]; 
  Gdiff[0] = (-0.7654655446197428*(alpha[1]*(fr[4]+fl[4])+alpha[0]*(fr[2]+fl[2])))+alpha[1]*(0.7954951288348656*fr[1]-0.7954951288348656*fl[1])+alpha[0]*(0.7954951288348656*fr[0]-0.7954951288348656*fl[0]); 
  Gdiff[1] = (-0.7654655446197428*(alpha[0]*(fr[4]+fl[4])+alpha[1]*(fr[2]+fl[2])))+alpha[0]*(0.7954951288348656*fr[1]-0.7954951288348656*fl[1])+(0.7954951288348656*fr[0]-0.7954951288348656*fl[0])*alpha[1]; 
  Gdiff[3] = (-0.7654655446197428*(alpha[1]*(fr[7]+fl[7])+alpha[0]*(fr[6]+fl[6])))+alpha[1]*(0.7954951288348656*fr[5]-0.7954951288348656*fl[5])+alpha[0]*(0.7954951288348656*fr[3]-0.7954951288348656*fl[3]); 
  Gdiff[5] = (-0.7654655446197428*(alpha[0]*(fr[7]+fl[7])+alpha[1]*(fr[6]+fl[6])))+alpha[0]*(0.7954951288348656*fr[5]-0.7954951288348656*fl[5])+alpha[1]*(0.7954951288348656*fr[3]-0.7954951288348656*fl[3]); 

  double Ghat[8]; 
  for(unsigned short int i=0; i<8; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(0.4330127018922193*dxvl[1]*favg[2]+0.25*favg[0]*dxvl[1])*nu+drBar[1]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]+drBar[0]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+(0.4330127018922193*dxvl[1]*favg[4]+0.25*dxvl[1]*favg[1])*nu-0.8660254037844386*fjump[4]+drBar[0]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+drBar[1]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[1]; 
  Ghat[3] = Gdiff[3]*rdv+(0.4330127018922193*dxvl[1]*favg[6]+0.25*dxvl[1]*favg[3])*nu+drBar[1]*(0.6123724356957944*favg[7]+0.3535533905932737*favg[5])-0.8660254037844386*fjump[6]+drBar[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[3]; 
  Ghat[5] = Gdiff[5]*rdv+(0.4330127018922193*dxvl[1]*favg[7]+0.25*dxvl[1]*favg[5])*nu-0.8660254037844386*fjump[7]+drBar[0]*(0.6123724356957944*favg[7]+0.3535533905932737*favg[5])+drBar[1]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[5]; 

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[1]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[3]; 
  incr1[7] = 0.8660254037844386*Ghat[5]; 

  double incr2[8]; 
  incr2[2] = alpha[1]*(0.1767766952966368*fl[4]-0.1767766952966368*fr[4])+alpha[0]*(0.1767766952966368*fl[2]-0.1767766952966368*fr[2])+0.1530931089239486*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[4] = alpha[0]*(0.1767766952966368*fl[4]-0.1767766952966368*fr[4])+alpha[1]*(0.1767766952966368*fl[2]-0.1767766952966368*fr[2])+0.1530931089239486*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 
  incr2[6] = alpha[1]*(0.1767766952966368*fl[7]-0.1767766952966368*fr[7])+alpha[0]*(0.1767766952966368*fl[6]-0.1767766952966368*fr[6])+0.1530931089239486*(alpha[1]*(fr[5]+fl[5])+alpha[0]*(fr[3]+fl[3])); 
  incr2[7] = alpha[0]*(0.1767766952966368*fl[7]-0.1767766952966368*fr[7])+alpha[1]*(0.1767766952966368*fl[6]-0.1767766952966368*fr[6])+0.1530931089239486*(alpha[0]*(fr[5]+fl[5])+alpha[1]*(fr[3]+fl[3])); 

  outr[0] += incr1[0]*rdv2r; 
  outr[1] += incr1[1]*rdv2r; 
  outr[2] += incr2[2]*rdvSq4r+incr1[2]*rdv2r; 
  outr[3] += incr1[3]*rdv2r; 
  outr[4] += incr2[4]*rdvSq4r+incr1[4]*rdv2r; 
  outr[5] += incr1[5]*rdv2r; 
  outr[6] += incr2[6]*rdvSq4r+incr1[6]*rdv2r; 
  outr[7] += incr2[7]*rdvSq4r+incr1[7]*rdv2r; 

  outl[0] += -1.0*incr1[0]*rdv2l; 
  outl[1] += -1.0*incr1[1]*rdv2l; 
  outl[2] += incr1[2]*rdv2l-1.0*incr2[2]*rdvSq4l; 
  outl[3] += -1.0*incr1[3]*rdv2l; 
  outl[4] += incr1[4]*rdv2l-1.0*incr2[4]*rdvSq4l; 
  outl[5] += -1.0*incr1[5]*rdv2l; 
  outl[6] += incr1[6]*rdv2l-1.0*incr2[6]*rdvSq4l; 
  outl[7] += incr1[7]*rdv2l-1.0*incr2[7]*rdvSq4l; 

  const double vMuMid = wl[1]-(0.7071067811865475*nuU[0])/nu; 
  return std::abs(vMuMid); 
} 
double GkLBOconstNuSurf1x2vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*3]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2l = 2.0/dxvl[1]; 
  double rdv2r = 2.0/dxvr[1]; 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

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
  fjump[0] = nu*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nu*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nu*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nu*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nu*vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = nu*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nu*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nu*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nu*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nu*vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = nu*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nu*vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = nu*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nu*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nu*vMuMidMax*(fl[14]-(1*fr[14])); 
  fjump[15] = nu*vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = nu*vMuMidMax*(fl[16]-(-1*fr[16])); 
  fjump[17] = nu*vMuMidMax*(fl[17]-(-1*fr[17])); 
  fjump[18] = nu*vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = nu*vMuMidMax*(fl[19]-(-1*fr[19])); 

  double drBar[3]; 
  drBar[0] = 1.414213562373095*wl[1]*nu-1.0*nuU[0]; 
  drBar[1] = -1.0*nuU[1]; 
  drBar[2] = -1.0*nuU[2]; 

  double alpha[20]; 
  alpha[0] = 2.0*nuVtSq[0]; 
  alpha[1] = 2.0*nuVtSq[1]; 
  alpha[7] = 2.0*nuVtSq[2]; 

  double Gdiff[20]; 
  Gdiff[0] = alpha[1]*(0.9486832980505137*fr[12]-0.9486832980505137*fl[12])-1.684024198163434*alpha[7]*(fr[11]+fl[11])+alpha[0]*(0.9486832980505137*fr[8]-0.9486832980505137*fl[8])+alpha[7]*(1.325825214724776*fr[7]-1.325825214724776*fl[7])-1.684024198163434*(alpha[1]*(fr[4]+fl[4])+alpha[0]*(fr[2]+fl[2]))+alpha[1]*(1.325825214724776*fr[1]-1.325825214724776*fl[1])+alpha[0]*(1.325825214724776*fr[0]-1.325825214724776*fl[0]); 
  Gdiff[1] = (0.848528137423857*alpha[7]+0.9486832980505137*alpha[0])*fr[12]+((-0.848528137423857*alpha[7])-0.9486832980505137*alpha[0])*fl[12]+alpha[1]*((-1.506237033139206*(fr[11]+fl[11]))+0.9486832980505137*fr[8]-0.9486832980505137*fl[8]+1.185854122563142*fr[7]-1.185854122563142*fl[7])+((-1.506237033139206*(fr[4]+fl[4]))+1.185854122563142*fr[1]-1.185854122563142*fl[1])*alpha[7]-1.684024198163434*(alpha[0]*(fr[4]+fl[4])+alpha[1]*(fr[2]+fl[2]))+alpha[0]*(1.325825214724776*fr[1]-1.325825214724776*fl[1])+(1.325825214724776*fr[0]-1.325825214724776*fl[0])*alpha[1]; 
  Gdiff[3] = alpha[1]*(0.9486832980505137*fr[18]-0.9486832980505137*fl[18])-1.684024198163434*alpha[7]*(fr[17]+fl[17])+alpha[0]*(0.9486832980505137*fr[14]-0.9486832980505137*fl[14])+alpha[7]*(1.325825214724776*fr[13]-1.325825214724776*fl[13])-1.684024198163434*(alpha[1]*(fr[10]+fl[10])+alpha[0]*(fr[6]+fl[6]))+alpha[1]*(1.325825214724776*fr[5]-1.325825214724776*fl[5])+alpha[0]*(1.325825214724776*fr[3]-1.325825214724776*fl[3]); 
  Gdiff[5] = (0.848528137423857*alpha[7]+0.9486832980505137*alpha[0])*fr[18]+((-0.848528137423857*alpha[7])-0.9486832980505137*alpha[0])*fl[18]+alpha[1]*((-1.506237033139206*(fr[17]+fl[17]))+0.9486832980505137*fr[14]-0.9486832980505137*fl[14]+1.185854122563142*fr[13]-1.185854122563142*fl[13])+((-1.506237033139206*alpha[7])-1.684024198163434*alpha[0])*fr[10]+((-1.506237033139206*alpha[7])-1.684024198163434*alpha[0])*fl[10]+(1.185854122563142*fr[5]-1.185854122563142*fl[5])*alpha[7]-1.684024198163434*alpha[1]*(fr[6]+fl[6])+alpha[0]*(1.325825214724776*fr[5]-1.325825214724776*fl[5])+alpha[1]*(1.325825214724776*fr[3]-1.325825214724776*fl[3]); 
  Gdiff[7] = alpha[1]*(0.848528137423857*fr[12]-0.848528137423857*fl[12])+((-1.075883595099433*alpha[7])-1.684024198163434*alpha[0])*fr[11]+((-1.075883595099433*alpha[7])-1.684024198163434*alpha[0])*fl[11]+alpha[7]*(0.9486832980505137*fr[8]-0.9486832980505137*fl[8])+(0.8470386589736728*alpha[7]+1.325825214724776*alpha[0])*fr[7]+((-0.8470386589736728*alpha[7])-1.325825214724776*alpha[0])*fl[7]+((-1.684024198163434*(fr[2]+fl[2]))+1.325825214724776*fr[0]-1.325825214724776*fl[0])*alpha[7]+alpha[1]*((-1.506237033139206*(fr[4]+fl[4]))+1.185854122563142*fr[1]-1.185854122563142*fl[1]); 
  Gdiff[9] = (-1.684024198163434*(alpha[1]*(fr[19]+fl[19])+alpha[0]*(fr[16]+fl[16])))+alpha[1]*(1.325825214724776*fr[15]-1.325825214724776*fl[15])+alpha[0]*(1.325825214724776*fr[9]-1.325825214724776*fl[9]); 
  Gdiff[13] = alpha[1]*(0.848528137423857*fr[18]-0.848528137423857*fl[18])+((-1.075883595099433*alpha[7])-1.684024198163434*alpha[0])*fr[17]+((-1.075883595099433*alpha[7])-1.684024198163434*alpha[0])*fl[17]+alpha[7]*(0.9486832980505137*fr[14]-0.9486832980505137*fl[14])+(0.8470386589736728*alpha[7]+1.325825214724776*alpha[0])*fr[13]+((-0.8470386589736728*alpha[7])-1.325825214724776*alpha[0])*fl[13]-1.506237033139206*alpha[1]*(fr[10]+fl[10])+((-1.684024198163434*(fr[6]+fl[6]))+1.325825214724776*fr[3]-1.325825214724776*fl[3])*alpha[7]+alpha[1]*(1.185854122563142*fr[5]-1.185854122563142*fl[5]); 
  Gdiff[15] = ((-1.506237033139206*alpha[7])-1.684024198163434*alpha[0])*fr[19]-1.506237033139206*alpha[7]*fl[19]-1.684024198163434*(alpha[0]*fl[19]+alpha[1]*(fr[16]+fl[16]))+(1.185854122563142*alpha[7]+1.325825214724776*alpha[0])*fr[15]+((-1.185854122563142*alpha[7])-1.325825214724776*alpha[0])*fl[15]+alpha[1]*(1.325825214724776*fr[9]-1.325825214724776*fl[9]); 

  double Ghat[20]; 
  for(unsigned short int i=0; i<20; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(0.5590169943749475*dxvl[1]*favg[8]+0.4330127018922193*dxvl[1]*favg[2]+0.25*favg[0]*dxvl[1])*nu+drBar[1]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+drBar[2]*(0.6123724356957944*favg[11]+0.3535533905932737*favg[7])-1.118033988749895*fjump[8]+drBar[0]*(0.7905694150420947*favg[8]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+(0.5590169943749476*dxvl[1]*favg[12]+0.4330127018922193*dxvl[1]*favg[4]+0.25*dxvl[1]*favg[1])*nu-1.118033988749895*fjump[12]+drBar[0]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+drBar[2]*(0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])+drBar[1]*(0.5477225575051661*favg[11]+0.7905694150420947*favg[8]+0.3162277660168379*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[4]-0.5*fjump[1]; 
  Ghat[3] = Gdiff[3]*rdv+(0.5590169943749476*dxvl[1]*favg[14]+0.4330127018922193*dxvl[1]*favg[6]+0.25*dxvl[1]*favg[3])*nu+drBar[1]*(0.7905694150420947*favg[18]+0.6123724356957944*favg[10]+0.3535533905932737*favg[5])+drBar[2]*(0.6123724356957944*favg[17]+0.3535533905932737*favg[13])-1.118033988749895*fjump[14]+drBar[0]*(0.7905694150420948*favg[14]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.8660254037844386*fjump[6]-0.5*fjump[3]; 
  Ghat[5] = Gdiff[5]*rdv+(0.5590169943749475*dxvl[1]*favg[18]+0.4330127018922193*dxvl[1]*favg[10]+0.25*dxvl[1]*favg[5])*nu-1.118033988749895*fjump[18]+drBar[0]*(0.7905694150420947*favg[18]+0.6123724356957944*favg[10]+0.3535533905932737*favg[5])+drBar[2]*(0.7071067811865475*favg[18]+0.5477225575051661*favg[10]+0.3162277660168379*favg[5])+drBar[1]*(0.5477225575051661*favg[17]+0.7905694150420948*favg[14]+0.3162277660168379*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.8660254037844386*fjump[10]-0.5*fjump[5]; 
  Ghat[7] = Gdiff[7]*rdv+(0.4330127018922194*dxvl[1]*favg[11]+0.25*dxvl[1]*favg[7])*nu+drBar[1]*(0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])-0.8660254037844387*fjump[11]+drBar[0]*(0.6123724356957944*favg[11]+0.3535533905932737*favg[7])+drBar[2]*(0.3912303982179757*favg[11]+0.7905694150420947*favg[8]+0.2258769757263128*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[7]; 
  Ghat[9] = Gdiff[9]*rdv+(0.4330127018922194*dxvl[1]*favg[16]+0.25*dxvl[1]*favg[9])*nu+drBar[1]*(0.6123724356957944*favg[19]+0.3535533905932737*favg[15])-0.8660254037844387*fjump[16]+drBar[0]*(0.6123724356957944*favg[16]+0.3535533905932737*favg[9])-0.5*fjump[9]; 
  Ghat[13] = Gdiff[13]*rdv+(0.4330127018922194*dxvl[1]*favg[17]+0.25*dxvl[1]*favg[13])*nu+drBar[1]*(0.7071067811865475*favg[18]+0.5477225575051661*favg[10]+0.3162277660168379*favg[5])-0.8660254037844387*fjump[17]+drBar[0]*(0.6123724356957944*favg[17]+0.3535533905932737*favg[13])+drBar[2]*(0.3912303982179757*favg[17]+0.7905694150420947*favg[14]+0.2258769757263128*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[13]; 
  Ghat[15] = Gdiff[15]*rdv+(0.4330127018922194*dxvl[1]*favg[19]+0.25*dxvl[1]*favg[15])*nu-0.8660254037844387*fjump[19]+drBar[0]*(0.6123724356957944*favg[19]+0.3535533905932737*favg[15])+drBar[2]*(0.5477225575051661*favg[19]+0.3162277660168379*favg[15])+drBar[1]*(0.6123724356957944*favg[16]+0.3535533905932737*favg[9])-0.5*fjump[15]; 

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
  incr2[2] = 0.1497678868178187*alpha[1]*(fr[12]+fl[12])+alpha[7]*(0.215446597392776*fl[11]-0.215446597392776*fr[11])+0.1497678868178187*alpha[0]*(fr[8]+fl[8])+0.1530931089239486*alpha[7]*(fr[7]+fl[7])+alpha[1]*(0.215446597392776*fl[4]-0.215446597392776*fr[4])+alpha[0]*(0.215446597392776*fl[2]-0.215446597392776*fr[2])+0.1530931089239486*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[4] = (0.1339564703084549*alpha[7]+0.1497678868178187*alpha[0])*fr[12]+(0.1339564703084549*alpha[7]+0.1497678868178187*alpha[0])*fl[12]+alpha[1]*((-0.1927012949165104*fr[11])+0.1927012949165104*fl[11]+0.1497678868178187*(fr[8]+fl[8])+0.1369306393762915*(fr[7]+fl[7]))+((-0.1927012949165104*fr[4])+0.1927012949165104*fl[4]+0.1369306393762915*(fr[1]+fl[1]))*alpha[7]+alpha[0]*(0.215446597392776*fl[4]-0.215446597392776*fr[4])+alpha[1]*(0.215446597392776*fl[2]-0.215446597392776*fr[2])+0.1530931089239486*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 
  incr2[6] = 0.1497678868178187*alpha[1]*(fr[18]+fl[18])+alpha[7]*(0.215446597392776*fl[17]-0.215446597392776*fr[17])+0.1497678868178187*alpha[0]*(fr[14]+fl[14])+0.1530931089239486*alpha[7]*(fr[13]+fl[13])+alpha[1]*(0.215446597392776*fl[10]-0.215446597392776*fr[10])+alpha[0]*(0.215446597392776*fl[6]-0.215446597392776*fr[6])+0.1530931089239486*(alpha[1]*(fr[5]+fl[5])+alpha[0]*(fr[3]+fl[3])); 
  incr2[8] = (-0.5800485314420891*alpha[1]*(fr[12]+fl[12]))+alpha[7]*(0.8344210836992756*fr[11]-0.8344210836992756*fl[11])-0.5800485314420891*alpha[0]*(fr[8]+fl[8])-0.592927061281571*alpha[7]*(fr[7]+fl[7])+alpha[1]*(0.8344210836992756*fr[4]-0.8344210836992756*fl[4])+alpha[0]*(0.8344210836992756*fr[2]-0.8344210836992756*fl[2])-0.592927061281571*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[10] = (0.1339564703084549*alpha[7]+0.1497678868178187*alpha[0])*fr[18]+(0.1339564703084549*alpha[7]+0.1497678868178187*alpha[0])*fl[18]+alpha[1]*((-0.1927012949165104*fr[17])+0.1927012949165104*fl[17]+0.1497678868178187*(fr[14]+fl[14])+0.1369306393762915*(fr[13]+fl[13]))+((-0.1927012949165104*alpha[7])-0.215446597392776*alpha[0])*fr[10]+(0.1927012949165104*alpha[7]+0.215446597392776*alpha[0])*fl[10]+0.1369306393762915*(fr[5]+fl[5])*alpha[7]+alpha[1]*(0.215446597392776*fl[6]-0.215446597392776*fr[6])+0.1530931089239486*(alpha[0]*(fr[5]+fl[5])+alpha[1]*(fr[3]+fl[3])); 
  incr2[11] = 0.1339564703084549*alpha[1]*(fr[12]+fl[12])+((-0.1376437820832217*alpha[7])-0.215446597392776*alpha[0])*fr[11]+(0.1376437820832217*alpha[7]+0.215446597392776*alpha[0])*fl[11]+0.1497678868178187*alpha[7]*(fr[8]+fl[8])+(0.09780759955449389*alpha[7]+0.1530931089239486*alpha[0])*fr[7]+(0.09780759955449389*alpha[7]+0.1530931089239486*alpha[0])*fl[7]+((-0.215446597392776*fr[2])+0.215446597392776*fl[2]+0.1530931089239486*(fr[0]+fl[0]))*alpha[7]+alpha[1]*((-0.1927012949165104*fr[4])+0.1927012949165104*fl[4]+0.1369306393762915*(fr[1]+fl[1])); 
  incr2[12] = ((-0.5188111786213743*alpha[7])-0.5800485314420891*alpha[0])*fr[12]+((-0.5188111786213743*alpha[7])-0.5800485314420891*alpha[0])*fl[12]+alpha[1]*(0.7463289060042488*fr[11]-0.7463289060042488*fl[11]-0.5800485314420891*(fr[8]+fl[8])-0.5303300858899105*(fr[7]+fl[7]))+(0.7463289060042488*fr[4]-0.7463289060042488*fl[4]-0.5303300858899105*(fr[1]+fl[1]))*alpha[7]+alpha[0]*(0.8344210836992756*fr[4]-0.8344210836992756*fl[4])+alpha[1]*(0.8344210836992756*fr[2]-0.8344210836992756*fl[2])-0.592927061281571*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 
  incr2[14] = (-0.5800485314420891*alpha[1]*(fr[18]+fl[18]))+alpha[7]*(0.8344210836992756*fr[17]-0.8344210836992756*fl[17])-0.5800485314420891*alpha[0]*(fr[14]+fl[14])-0.592927061281571*alpha[7]*(fr[13]+fl[13])+alpha[1]*(0.8344210836992756*fr[10]-0.8344210836992756*fl[10])+alpha[0]*(0.8344210836992756*fr[6]-0.8344210836992756*fl[6])-0.592927061281571*(alpha[1]*(fr[5]+fl[5])+alpha[0]*(fr[3]+fl[3])); 
  incr2[16] = alpha[1]*(0.215446597392776*fl[19]-0.215446597392776*fr[19])+alpha[0]*(0.215446597392776*fl[16]-0.215446597392776*fr[16])+0.1530931089239486*(alpha[1]*(fr[15]+fl[15])+alpha[0]*(fr[9]+fl[9])); 
  incr2[17] = 0.1339564703084549*alpha[1]*(fr[18]+fl[18])+((-0.1376437820832217*alpha[7])-0.215446597392776*alpha[0])*fr[17]+(0.1376437820832217*alpha[7]+0.215446597392776*alpha[0])*fl[17]+0.1497678868178187*alpha[7]*(fr[14]+fl[14])+(0.09780759955449389*alpha[7]+0.1530931089239486*alpha[0])*fr[13]+(0.09780759955449389*alpha[7]+0.1530931089239486*alpha[0])*fl[13]+alpha[1]*(0.1927012949165104*fl[10]-0.1927012949165104*fr[10])+((-0.215446597392776*fr[6])+0.215446597392776*fl[6]+0.1530931089239486*(fr[3]+fl[3]))*alpha[7]+0.1369306393762915*alpha[1]*(fr[5]+fl[5]); 
  incr2[18] = ((-0.5188111786213743*alpha[7])-0.5800485314420891*alpha[0])*fr[18]+((-0.5188111786213743*alpha[7])-0.5800485314420891*alpha[0])*fl[18]+alpha[1]*(0.7463289060042488*fr[17]-0.7463289060042488*fl[17]-0.5800485314420891*(fr[14]+fl[14])-0.5303300858899105*(fr[13]+fl[13]))+(0.7463289060042488*alpha[7]+0.8344210836992756*alpha[0])*fr[10]+((-0.7463289060042488*alpha[7])-0.8344210836992756*alpha[0])*fl[10]-0.5303300858899105*(fr[5]+fl[5])*alpha[7]+alpha[1]*(0.8344210836992756*fr[6]-0.8344210836992756*fl[6])-0.592927061281571*(alpha[0]*(fr[5]+fl[5])+alpha[1]*(fr[3]+fl[3])); 
  incr2[19] = ((-0.1927012949165104*alpha[7])-0.215446597392776*alpha[0])*fr[19]+(0.1927012949165104*alpha[7]+0.215446597392776*alpha[0])*fl[19]+alpha[1]*(0.215446597392776*fl[16]-0.215446597392776*fr[16])+(0.1369306393762915*alpha[7]+0.1530931089239486*alpha[0])*fr[15]+0.1369306393762915*alpha[7]*fl[15]+0.1530931089239486*(alpha[0]*fl[15]+alpha[1]*(fr[9]+fl[9])); 

  outr[0] += incr1[0]*rdv2r; 
  outr[1] += incr1[1]*rdv2r; 
  outr[2] += incr2[2]*rdvSq4r+incr1[2]*rdv2r; 
  outr[3] += incr1[3]*rdv2r; 
  outr[4] += incr2[4]*rdvSq4r+incr1[4]*rdv2r; 
  outr[5] += incr1[5]*rdv2r; 
  outr[6] += incr2[6]*rdvSq4r+incr1[6]*rdv2r; 
  outr[7] += incr1[7]*rdv2r; 
  outr[8] += incr2[8]*rdvSq4r+incr1[8]*rdv2r; 
  outr[9] += incr1[9]*rdv2r; 
  outr[10] += incr2[10]*rdvSq4r+incr1[10]*rdv2r; 
  outr[11] += incr2[11]*rdvSq4r+incr1[11]*rdv2r; 
  outr[12] += incr2[12]*rdvSq4r+incr1[12]*rdv2r; 
  outr[13] += incr1[13]*rdv2r; 
  outr[14] += incr2[14]*rdvSq4r+incr1[14]*rdv2r; 
  outr[15] += incr1[15]*rdv2r; 
  outr[16] += incr2[16]*rdvSq4r+incr1[16]*rdv2r; 
  outr[17] += incr2[17]*rdvSq4r+incr1[17]*rdv2r; 
  outr[18] += incr2[18]*rdvSq4r+incr1[18]*rdv2r; 
  outr[19] += incr2[19]*rdvSq4r+incr1[19]*rdv2r; 

  outl[0] += -1.0*incr1[0]*rdv2l; 
  outl[1] += -1.0*incr1[1]*rdv2l; 
  outl[2] += incr1[2]*rdv2l-1.0*incr2[2]*rdvSq4l; 
  outl[3] += -1.0*incr1[3]*rdv2l; 
  outl[4] += incr1[4]*rdv2l-1.0*incr2[4]*rdvSq4l; 
  outl[5] += -1.0*incr1[5]*rdv2l; 
  outl[6] += incr1[6]*rdv2l-1.0*incr2[6]*rdvSq4l; 
  outl[7] += -1.0*incr1[7]*rdv2l; 
  outl[8] += incr2[8]*rdvSq4l-1.0*incr1[8]*rdv2l; 
  outl[9] += -1.0*incr1[9]*rdv2l; 
  outl[10] += incr1[10]*rdv2l-1.0*incr2[10]*rdvSq4l; 
  outl[11] += incr1[11]*rdv2l-1.0*incr2[11]*rdvSq4l; 
  outl[12] += incr2[12]*rdvSq4l-1.0*incr1[12]*rdv2l; 
  outl[13] += -1.0*incr1[13]*rdv2l; 
  outl[14] += incr2[14]*rdvSq4l-1.0*incr1[14]*rdv2l; 
  outl[15] += -1.0*incr1[15]*rdv2l; 
  outl[16] += incr1[16]*rdv2l-1.0*incr2[16]*rdvSq4l; 
  outl[17] += incr1[17]*rdv2l-1.0*incr2[17]*rdvSq4l; 
  outl[18] += incr2[18]*rdvSq4l-1.0*incr1[18]*rdv2l; 
  outl[19] += incr1[19]*rdv2l-1.0*incr2[19]*rdvSq4l; 

  const double vMuMid = (0.7905694150420947*nuU[2])/nu-(0.7071067811865475*nuU[0])/nu+wl[1]; 
  return std::abs(vMuMid); 
} 
double GkLBOconstNuSurf1x2vSer_Mu_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*2]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2l = 2.0/dxvl[2]; 
  double rdv2r = 2.0/dxvr[2]; 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  double favg[8]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 

  double fjump[8]; 
  fjump[0] = nu*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nu*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nu*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nu*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nu*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nu*vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = nu*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nu*vMuMidMax*(fl[7]-(-1*fr[7])); 

  double alpha[8]; 
  alpha[0] = (BmagInv[1]*nuVtSq[1]+BmagInv[0]*nuVtSq[0])*(2.828427124746191*wl[2]+1.414213562373095*dxvl[2])*m_; 
  alpha[1] = (BmagInv[0]*nuVtSq[1]+nuVtSq[0]*BmagInv[1])*(2.828427124746191*wl[2]+1.414213562373095*dxvl[2])*m_; 

  double Gdiff[8]; 
  Gdiff[0] = (-0.7654655446197428*(alpha[1]*(fr[5]+fl[5])+alpha[0]*(fr[3]+fl[3])))+alpha[1]*(0.7954951288348656*fr[1]-0.7954951288348656*fl[1])+alpha[0]*(0.7954951288348656*fr[0]-0.7954951288348656*fl[0]); 
  Gdiff[1] = (-0.7654655446197428*(alpha[0]*(fr[5]+fl[5])+alpha[1]*(fr[3]+fl[3])))+alpha[0]*(0.7954951288348656*fr[1]-0.7954951288348656*fl[1])+(0.7954951288348656*fr[0]-0.7954951288348656*fl[0])*alpha[1]; 
  Gdiff[2] = (-0.7654655446197428*(alpha[1]*(fr[7]+fl[7])+alpha[0]*(fr[6]+fl[6])))+alpha[1]*(0.7954951288348656*fr[4]-0.7954951288348656*fl[4])+alpha[0]*(0.7954951288348656*fr[2]-0.7954951288348656*fl[2]); 
  Gdiff[4] = (-0.7654655446197428*(alpha[0]*(fr[7]+fl[7])+alpha[1]*(fr[6]+fl[6])))+alpha[0]*(0.7954951288348656*fr[4]-0.7954951288348656*fl[4])+alpha[1]*(0.7954951288348656*fr[2]-0.7954951288348656*fl[2]); 

  double Ghat[8]; 
  for(unsigned short int i=0; i<8; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(1.732050807568877*wl[2]*favg[3]+0.8660254037844386*dxvl[2]*favg[3]+favg[0]*wl[2]+0.5*favg[0]*dxvl[2])*nu-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+(1.732050807568877*wl[2]*favg[5]+0.8660254037844386*dxvl[2]*favg[5]+favg[1]*wl[2]+0.5*favg[1]*dxvl[2])*nu-0.8660254037844386*fjump[5]-0.5*fjump[1]; 
  Ghat[2] = Gdiff[2]*rdv+(1.732050807568877*wl[2]*favg[6]+0.8660254037844386*dxvl[2]*favg[6]+favg[2]*wl[2]+0.5*dxvl[2]*favg[2])*nu-0.8660254037844386*fjump[6]-0.5*fjump[2]; 
  Ghat[4] = Gdiff[4]*rdv+(1.732050807568877*wl[2]*favg[7]+0.8660254037844386*dxvl[2]*favg[7]+wl[2]*favg[4]+0.5*dxvl[2]*favg[4])*nu-0.8660254037844386*fjump[7]-0.5*fjump[4]; 

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = 0.8660254037844386*Ghat[2]; 
  incr1[7] = 0.8660254037844386*Ghat[4]; 

  double incr2[8]; 
  incr2[3] = alpha[1]*(0.1767766952966368*fl[5]-0.1767766952966368*fr[5])+alpha[0]*(0.1767766952966368*fl[3]-0.1767766952966368*fr[3])+0.1530931089239486*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[5] = alpha[0]*(0.1767766952966368*fl[5]-0.1767766952966368*fr[5])+alpha[1]*(0.1767766952966368*fl[3]-0.1767766952966368*fr[3])+0.1530931089239486*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 
  incr2[6] = alpha[1]*(0.1767766952966368*fl[7]-0.1767766952966368*fr[7])+alpha[0]*(0.1767766952966368*fl[6]-0.1767766952966368*fr[6])+0.1530931089239486*(alpha[1]*(fr[4]+fl[4])+alpha[0]*(fr[2]+fl[2])); 
  incr2[7] = alpha[0]*(0.1767766952966368*fl[7]-0.1767766952966368*fr[7])+alpha[1]*(0.1767766952966368*fl[6]-0.1767766952966368*fr[6])+0.1530931089239486*(alpha[0]*(fr[4]+fl[4])+alpha[1]*(fr[2]+fl[2])); 

  outr[0] += incr1[0]*rdv2r; 
  outr[1] += incr1[1]*rdv2r; 
  outr[2] += incr1[2]*rdv2r; 
  outr[3] += incr2[3]*rdvSq4r+incr1[3]*rdv2r; 
  outr[4] += incr1[4]*rdv2r; 
  outr[5] += incr2[5]*rdvSq4r+incr1[5]*rdv2r; 
  outr[6] += incr2[6]*rdvSq4r+incr1[6]*rdv2r; 
  outr[7] += incr2[7]*rdvSq4r+incr1[7]*rdv2r; 

  outl[0] += -1.0*incr1[0]*rdv2l; 
  outl[1] += -1.0*incr1[1]*rdv2l; 
  outl[2] += -1.0*incr1[2]*rdv2l; 
  outl[3] += incr1[3]*rdv2l-1.0*incr2[3]*rdvSq4l; 
  outl[4] += -1.0*incr1[4]*rdv2l; 
  outl[5] += incr1[5]*rdv2l-1.0*incr2[5]*rdvSq4l; 
  outl[6] += incr1[6]*rdv2l-1.0*incr2[6]*rdvSq4l; 
  outl[7] += incr1[7]*rdv2l-1.0*incr2[7]*rdvSq4l; 

  const double vMuMid = 2.0*wl[2]; 
  return std::abs(vMuMid); 
} 
double GkLBOconstNuSurf1x2vSer_Mu_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*3]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2l = 2.0/dxvl[2]; 
  double rdv2r = 2.0/dxvr[2]; 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

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
  fjump[0] = nu*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nu*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nu*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nu*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nu*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nu*vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = nu*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nu*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nu*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nu*vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = nu*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nu*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nu*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nu*vMuMidMax*(fl[13]-(-1*fr[13])); 
  fjump[14] = nu*vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = nu*vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = nu*vMuMidMax*(fl[16]-(1*fr[16])); 
  fjump[17] = nu*vMuMidMax*(fl[17]-(-1*fr[17])); 
  fjump[18] = nu*vMuMidMax*(fl[18]-(-1*fr[18])); 
  fjump[19] = nu*vMuMidMax*(fl[19]-(1*fr[19])); 

  double alpha[20]; 
  alpha[0] = (2.828427124746191*wl[2]+1.414213562373095*dxvl[2])*(BmagInv[2]*nuVtSq[2]+BmagInv[1]*nuVtSq[1]+BmagInv[0]*nuVtSq[0])*m_; 
  alpha[1] = (BmagInv[1]*(2.529822128134704*wl[2]+1.264911064067352*dxvl[2])*nuVtSq[2]+nuVtSq[1]*((2.529822128134704*BmagInv[2]+2.828427124746191*BmagInv[0])*wl[2]+(1.264911064067352*BmagInv[2]+1.414213562373095*BmagInv[0])*dxvl[2])+nuVtSq[0]*BmagInv[1]*(2.828427124746191*wl[2]+1.414213562373095*dxvl[2]))*m_; 
  alpha[7] = (((1.807015805810503*BmagInv[2]+2.828427124746191*BmagInv[0])*wl[2]+(0.9035079029052515*BmagInv[2]+1.414213562373095*BmagInv[0])*dxvl[2])*nuVtSq[2]+nuVtSq[0]*BmagInv[2]*(2.828427124746191*wl[2]+1.414213562373095*dxvl[2])+BmagInv[1]*nuVtSq[1]*(2.529822128134704*wl[2]+1.264911064067352*dxvl[2]))*m_; 

  double Gdiff[20]; 
  Gdiff[0] = alpha[1]*(0.9486832980505137*fr[15]-0.9486832980505137*fl[15])-1.684024198163434*alpha[7]*(fr[13]+fl[13])+alpha[0]*(0.9486832980505137*fr[9]-0.9486832980505137*fl[9])+alpha[7]*(1.325825214724776*fr[7]-1.325825214724776*fl[7])-1.684024198163434*(alpha[1]*(fr[5]+fl[5])+alpha[0]*(fr[3]+fl[3]))+alpha[1]*(1.325825214724776*fr[1]-1.325825214724776*fl[1])+alpha[0]*(1.325825214724776*fr[0]-1.325825214724776*fl[0]); 
  Gdiff[1] = (0.848528137423857*alpha[7]+0.9486832980505137*alpha[0])*fr[15]+((-0.848528137423857*alpha[7])-0.9486832980505137*alpha[0])*fl[15]+alpha[1]*((-1.506237033139206*(fr[13]+fl[13]))+0.9486832980505137*fr[9]-0.9486832980505137*fl[9]+1.185854122563142*fr[7]-1.185854122563142*fl[7])+((-1.506237033139206*(fr[5]+fl[5]))+1.185854122563142*fr[1]-1.185854122563142*fl[1])*alpha[7]-1.684024198163434*(alpha[0]*(fr[5]+fl[5])+alpha[1]*(fr[3]+fl[3]))+alpha[0]*(1.325825214724776*fr[1]-1.325825214724776*fl[1])+(1.325825214724776*fr[0]-1.325825214724776*fl[0])*alpha[1]; 
  Gdiff[2] = alpha[1]*(0.9486832980505137*fr[19]-0.9486832980505137*fl[19])-1.684024198163434*alpha[7]*(fr[17]+fl[17])+alpha[0]*(0.9486832980505137*fr[16]-0.9486832980505137*fl[16])+alpha[7]*(1.325825214724776*fr[11]-1.325825214724776*fl[11])-1.684024198163434*(alpha[1]*(fr[10]+fl[10])+alpha[0]*(fr[6]+fl[6]))+alpha[1]*(1.325825214724776*fr[4]-1.325825214724776*fl[4])+alpha[0]*(1.325825214724776*fr[2]-1.325825214724776*fl[2]); 
  Gdiff[4] = (0.848528137423857*alpha[7]+0.9486832980505137*alpha[0])*fr[19]+((-0.848528137423857*alpha[7])-0.9486832980505137*alpha[0])*fl[19]+alpha[1]*((-1.506237033139206*(fr[17]+fl[17]))+0.9486832980505137*fr[16]-0.9486832980505137*fl[16]+1.185854122563142*fr[11]-1.185854122563142*fl[11])+((-1.506237033139206*alpha[7])-1.684024198163434*alpha[0])*fr[10]+((-1.506237033139206*alpha[7])-1.684024198163434*alpha[0])*fl[10]+(1.185854122563142*fr[4]-1.185854122563142*fl[4])*alpha[7]-1.684024198163434*alpha[1]*(fr[6]+fl[6])+alpha[0]*(1.325825214724776*fr[4]-1.325825214724776*fl[4])+alpha[1]*(1.325825214724776*fr[2]-1.325825214724776*fl[2]); 
  Gdiff[7] = alpha[1]*(0.848528137423857*fr[15]-0.848528137423857*fl[15])+((-1.075883595099433*alpha[7])-1.684024198163434*alpha[0])*fr[13]+((-1.075883595099433*alpha[7])-1.684024198163434*alpha[0])*fl[13]+alpha[7]*(0.9486832980505137*fr[9]-0.9486832980505137*fl[9])+(0.8470386589736728*alpha[7]+1.325825214724776*alpha[0])*fr[7]+((-0.8470386589736728*alpha[7])-1.325825214724776*alpha[0])*fl[7]+((-1.684024198163434*(fr[3]+fl[3]))+1.325825214724776*fr[0]-1.325825214724776*fl[0])*alpha[7]+alpha[1]*((-1.506237033139206*(fr[5]+fl[5]))+1.185854122563142*fr[1]-1.185854122563142*fl[1]); 
  Gdiff[8] = (-1.684024198163434*(alpha[1]*(fr[18]+fl[18])+alpha[0]*(fr[14]+fl[14])))+alpha[1]*(1.325825214724776*fr[12]-1.325825214724776*fl[12])+alpha[0]*(1.325825214724776*fr[8]-1.325825214724776*fl[8]); 
  Gdiff[11] = alpha[1]*(0.848528137423857*fr[19]-0.848528137423857*fl[19])+((-1.075883595099433*alpha[7])-1.684024198163434*alpha[0])*fr[17]+((-1.075883595099433*alpha[7])-1.684024198163434*alpha[0])*fl[17]+alpha[7]*(0.9486832980505137*fr[16]-0.9486832980505137*fl[16])+(0.8470386589736728*alpha[7]+1.325825214724776*alpha[0])*fr[11]+((-0.8470386589736728*alpha[7])-1.325825214724776*alpha[0])*fl[11]-1.506237033139206*alpha[1]*(fr[10]+fl[10])+((-1.684024198163434*(fr[6]+fl[6]))+1.325825214724776*fr[2]-1.325825214724776*fl[2])*alpha[7]+alpha[1]*(1.185854122563142*fr[4]-1.185854122563142*fl[4]); 
  Gdiff[12] = ((-1.506237033139206*alpha[7])-1.684024198163434*alpha[0])*fr[18]-1.506237033139206*alpha[7]*fl[18]-1.684024198163434*(alpha[0]*fl[18]+alpha[1]*(fr[14]+fl[14]))+(1.185854122563142*alpha[7]+1.325825214724776*alpha[0])*fr[12]+((-1.185854122563142*alpha[7])-1.325825214724776*alpha[0])*fl[12]+alpha[1]*(1.325825214724776*fr[8]-1.325825214724776*fl[8]); 

  double Ghat[20]; 
  for(unsigned short int i=0; i<20; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(2.23606797749979*wl[2]*favg[9]+1.118033988749895*dxvl[2]*favg[9]+1.732050807568877*wl[2]*favg[3]+0.8660254037844386*dxvl[2]*favg[3]+favg[0]*wl[2]+0.5*favg[0]*dxvl[2])*nu-1.118033988749895*fjump[9]-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+(2.23606797749979*wl[2]*favg[15]+1.118033988749895*dxvl[2]*favg[15]+1.732050807568877*wl[2]*favg[5]+0.8660254037844386*dxvl[2]*favg[5]+favg[1]*wl[2]+0.5*favg[1]*dxvl[2])*nu-1.118033988749895*fjump[15]-0.8660254037844386*fjump[5]-0.5*fjump[1]; 
  Ghat[2] = Gdiff[2]*rdv+(2.23606797749979*wl[2]*favg[16]+1.118033988749895*dxvl[2]*favg[16]+1.732050807568877*wl[2]*favg[6]+0.8660254037844386*dxvl[2]*favg[6]+favg[2]*wl[2]+0.5*dxvl[2]*favg[2])*nu-1.118033988749895*fjump[16]-0.8660254037844386*fjump[6]-0.5*fjump[2]; 
  Ghat[4] = Gdiff[4]*rdv+(2.23606797749979*wl[2]*favg[19]+1.118033988749895*dxvl[2]*favg[19]+1.732050807568877*wl[2]*favg[10]+0.8660254037844386*dxvl[2]*favg[10]+wl[2]*favg[4]+0.5*dxvl[2]*favg[4])*nu-1.118033988749895*fjump[19]-0.8660254037844386*fjump[10]-0.5*fjump[4]; 
  Ghat[7] = Gdiff[7]*rdv+(1.732050807568877*wl[2]*favg[13]+0.8660254037844387*dxvl[2]*favg[13]+wl[2]*favg[7]+0.5*dxvl[2]*favg[7])*nu-0.8660254037844387*fjump[13]-0.5*fjump[7]; 
  Ghat[8] = Gdiff[8]*rdv+(1.732050807568877*wl[2]*favg[14]+0.8660254037844387*dxvl[2]*favg[14]+wl[2]*favg[8]+0.5*dxvl[2]*favg[8])*nu-0.8660254037844387*fjump[14]-0.5*fjump[8]; 
  Ghat[11] = Gdiff[11]*rdv+(1.732050807568877*wl[2]*favg[17]+0.8660254037844387*dxvl[2]*favg[17]+wl[2]*favg[11]+0.5*dxvl[2]*favg[11])*nu-0.8660254037844387*fjump[17]-0.5*fjump[11]; 
  Ghat[12] = Gdiff[12]*rdv+(1.732050807568877*wl[2]*favg[18]+0.8660254037844387*dxvl[2]*favg[18]+wl[2]*favg[12]+0.5*dxvl[2]*favg[12])*nu-0.8660254037844387*fjump[18]-0.5*fjump[12]; 

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
  incr2[3] = 0.1497678868178187*alpha[1]*(fr[15]+fl[15])+alpha[7]*(0.215446597392776*fl[13]-0.215446597392776*fr[13])+0.1497678868178187*alpha[0]*(fr[9]+fl[9])+0.1530931089239486*alpha[7]*(fr[7]+fl[7])+alpha[1]*(0.215446597392776*fl[5]-0.215446597392776*fr[5])+alpha[0]*(0.215446597392776*fl[3]-0.215446597392776*fr[3])+0.1530931089239486*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[5] = (0.1339564703084549*alpha[7]+0.1497678868178187*alpha[0])*fr[15]+(0.1339564703084549*alpha[7]+0.1497678868178187*alpha[0])*fl[15]+alpha[1]*((-0.1927012949165104*fr[13])+0.1927012949165104*fl[13]+0.1497678868178187*(fr[9]+fl[9])+0.1369306393762915*(fr[7]+fl[7]))+((-0.1927012949165104*fr[5])+0.1927012949165104*fl[5]+0.1369306393762915*(fr[1]+fl[1]))*alpha[7]+alpha[0]*(0.215446597392776*fl[5]-0.215446597392776*fr[5])+alpha[1]*(0.215446597392776*fl[3]-0.215446597392776*fr[3])+0.1530931089239486*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 
  incr2[6] = 0.1497678868178187*alpha[1]*(fr[19]+fl[19])+alpha[7]*(0.215446597392776*fl[17]-0.215446597392776*fr[17])+0.1497678868178187*alpha[0]*(fr[16]+fl[16])+0.1530931089239486*alpha[7]*(fr[11]+fl[11])+alpha[1]*(0.215446597392776*fl[10]-0.215446597392776*fr[10])+alpha[0]*(0.215446597392776*fl[6]-0.215446597392776*fr[6])+0.1530931089239486*(alpha[1]*(fr[4]+fl[4])+alpha[0]*(fr[2]+fl[2])); 
  incr2[9] = (-0.5800485314420891*alpha[1]*(fr[15]+fl[15]))+alpha[7]*(0.8344210836992756*fr[13]-0.8344210836992756*fl[13])-0.5800485314420891*alpha[0]*(fr[9]+fl[9])-0.592927061281571*alpha[7]*(fr[7]+fl[7])+alpha[1]*(0.8344210836992756*fr[5]-0.8344210836992756*fl[5])+alpha[0]*(0.8344210836992756*fr[3]-0.8344210836992756*fl[3])-0.592927061281571*(alpha[1]*(fr[1]+fl[1])+alpha[0]*(fr[0]+fl[0])); 
  incr2[10] = (0.1339564703084549*alpha[7]+0.1497678868178187*alpha[0])*fr[19]+(0.1339564703084549*alpha[7]+0.1497678868178187*alpha[0])*fl[19]+alpha[1]*((-0.1927012949165104*fr[17])+0.1927012949165104*fl[17]+0.1497678868178187*(fr[16]+fl[16])+0.1369306393762915*(fr[11]+fl[11]))+((-0.1927012949165104*alpha[7])-0.215446597392776*alpha[0])*fr[10]+(0.1927012949165104*alpha[7]+0.215446597392776*alpha[0])*fl[10]+0.1369306393762915*(fr[4]+fl[4])*alpha[7]+alpha[1]*(0.215446597392776*fl[6]-0.215446597392776*fr[6])+0.1530931089239486*(alpha[0]*(fr[4]+fl[4])+alpha[1]*(fr[2]+fl[2])); 
  incr2[13] = 0.1339564703084549*alpha[1]*(fr[15]+fl[15])+((-0.1376437820832217*alpha[7])-0.215446597392776*alpha[0])*fr[13]+(0.1376437820832217*alpha[7]+0.215446597392776*alpha[0])*fl[13]+0.1497678868178187*alpha[7]*(fr[9]+fl[9])+(0.09780759955449389*alpha[7]+0.1530931089239486*alpha[0])*fr[7]+(0.09780759955449389*alpha[7]+0.1530931089239486*alpha[0])*fl[7]+((-0.215446597392776*fr[3])+0.215446597392776*fl[3]+0.1530931089239486*(fr[0]+fl[0]))*alpha[7]+alpha[1]*((-0.1927012949165104*fr[5])+0.1927012949165104*fl[5]+0.1369306393762915*(fr[1]+fl[1])); 
  incr2[14] = alpha[1]*(0.215446597392776*fl[18]-0.215446597392776*fr[18])+alpha[0]*(0.215446597392776*fl[14]-0.215446597392776*fr[14])+0.1530931089239486*(alpha[1]*(fr[12]+fl[12])+alpha[0]*(fr[8]+fl[8])); 
  incr2[15] = ((-0.5188111786213743*alpha[7])-0.5800485314420891*alpha[0])*fr[15]+((-0.5188111786213743*alpha[7])-0.5800485314420891*alpha[0])*fl[15]+alpha[1]*(0.7463289060042488*fr[13]-0.7463289060042488*fl[13]-0.5800485314420891*(fr[9]+fl[9])-0.5303300858899105*(fr[7]+fl[7]))+(0.7463289060042488*fr[5]-0.7463289060042488*fl[5]-0.5303300858899105*(fr[1]+fl[1]))*alpha[7]+alpha[0]*(0.8344210836992756*fr[5]-0.8344210836992756*fl[5])+alpha[1]*(0.8344210836992756*fr[3]-0.8344210836992756*fl[3])-0.592927061281571*(alpha[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*alpha[1]); 
  incr2[16] = (-0.5800485314420891*alpha[1]*(fr[19]+fl[19]))+alpha[7]*(0.8344210836992756*fr[17]-0.8344210836992756*fl[17])-0.5800485314420891*alpha[0]*(fr[16]+fl[16])-0.592927061281571*alpha[7]*(fr[11]+fl[11])+alpha[1]*(0.8344210836992756*fr[10]-0.8344210836992756*fl[10])+alpha[0]*(0.8344210836992756*fr[6]-0.8344210836992756*fl[6])-0.592927061281571*(alpha[1]*(fr[4]+fl[4])+alpha[0]*(fr[2]+fl[2])); 
  incr2[17] = 0.1339564703084549*alpha[1]*(fr[19]+fl[19])+((-0.1376437820832217*alpha[7])-0.215446597392776*alpha[0])*fr[17]+(0.1376437820832217*alpha[7]+0.215446597392776*alpha[0])*fl[17]+0.1497678868178187*alpha[7]*(fr[16]+fl[16])+(0.09780759955449389*alpha[7]+0.1530931089239486*alpha[0])*fr[11]+(0.09780759955449389*alpha[7]+0.1530931089239486*alpha[0])*fl[11]+alpha[1]*(0.1927012949165104*fl[10]-0.1927012949165104*fr[10])+((-0.215446597392776*fr[6])+0.215446597392776*fl[6]+0.1530931089239486*(fr[2]+fl[2]))*alpha[7]+0.1369306393762915*alpha[1]*(fr[4]+fl[4]); 
  incr2[18] = ((-0.1927012949165104*alpha[7])-0.215446597392776*alpha[0])*fr[18]+(0.1927012949165104*alpha[7]+0.215446597392776*alpha[0])*fl[18]+alpha[1]*(0.215446597392776*fl[14]-0.215446597392776*fr[14])+(0.1369306393762915*alpha[7]+0.1530931089239486*alpha[0])*fr[12]+0.1369306393762915*alpha[7]*fl[12]+0.1530931089239486*(alpha[0]*fl[12]+alpha[1]*(fr[8]+fl[8])); 
  incr2[19] = ((-0.5188111786213743*alpha[7])-0.5800485314420891*alpha[0])*fr[19]+((-0.5188111786213743*alpha[7])-0.5800485314420891*alpha[0])*fl[19]+alpha[1]*(0.7463289060042488*fr[17]-0.7463289060042488*fl[17]-0.5800485314420891*(fr[16]+fl[16])-0.5303300858899105*(fr[11]+fl[11]))+(0.7463289060042488*alpha[7]+0.8344210836992756*alpha[0])*fr[10]+((-0.7463289060042488*alpha[7])-0.8344210836992756*alpha[0])*fl[10]-0.5303300858899105*(fr[4]+fl[4])*alpha[7]+alpha[1]*(0.8344210836992756*fr[6]-0.8344210836992756*fl[6])-0.592927061281571*(alpha[0]*(fr[4]+fl[4])+alpha[1]*(fr[2]+fl[2])); 

  outr[0] += incr1[0]*rdv2r; 
  outr[1] += incr1[1]*rdv2r; 
  outr[2] += incr1[2]*rdv2r; 
  outr[3] += incr2[3]*rdvSq4r+incr1[3]*rdv2r; 
  outr[4] += incr1[4]*rdv2r; 
  outr[5] += incr2[5]*rdvSq4r+incr1[5]*rdv2r; 
  outr[6] += incr2[6]*rdvSq4r+incr1[6]*rdv2r; 
  outr[7] += incr1[7]*rdv2r; 
  outr[8] += incr1[8]*rdv2r; 
  outr[9] += incr2[9]*rdvSq4r+incr1[9]*rdv2r; 
  outr[10] += incr2[10]*rdvSq4r+incr1[10]*rdv2r; 
  outr[11] += incr1[11]*rdv2r; 
  outr[12] += incr1[12]*rdv2r; 
  outr[13] += incr2[13]*rdvSq4r+incr1[13]*rdv2r; 
  outr[14] += incr2[14]*rdvSq4r+incr1[14]*rdv2r; 
  outr[15] += incr2[15]*rdvSq4r+incr1[15]*rdv2r; 
  outr[16] += incr2[16]*rdvSq4r+incr1[16]*rdv2r; 
  outr[17] += incr2[17]*rdvSq4r+incr1[17]*rdv2r; 
  outr[18] += incr2[18]*rdvSq4r+incr1[18]*rdv2r; 
  outr[19] += incr2[19]*rdvSq4r+incr1[19]*rdv2r; 

  outl[0] += -1.0*incr1[0]*rdv2l; 
  outl[1] += -1.0*incr1[1]*rdv2l; 
  outl[2] += -1.0*incr1[2]*rdv2l; 
  outl[3] += incr1[3]*rdv2l-1.0*incr2[3]*rdvSq4l; 
  outl[4] += -1.0*incr1[4]*rdv2l; 
  outl[5] += incr1[5]*rdv2l-1.0*incr2[5]*rdvSq4l; 
  outl[6] += incr1[6]*rdv2l-1.0*incr2[6]*rdvSq4l; 
  outl[7] += -1.0*incr1[7]*rdv2l; 
  outl[8] += -1.0*incr1[8]*rdv2l; 
  outl[9] += incr2[9]*rdvSq4l-1.0*incr1[9]*rdv2l; 
  outl[10] += incr1[10]*rdv2l-1.0*incr2[10]*rdvSq4l; 
  outl[11] += -1.0*incr1[11]*rdv2l; 
  outl[12] += -1.0*incr1[12]*rdv2l; 
  outl[13] += incr1[13]*rdv2l-1.0*incr2[13]*rdvSq4l; 
  outl[14] += incr1[14]*rdv2l-1.0*incr2[14]*rdvSq4l; 
  outl[15] += incr2[15]*rdvSq4l-1.0*incr1[15]*rdv2l; 
  outl[16] += incr2[16]*rdvSq4l-1.0*incr1[16]*rdv2l; 
  outl[17] += incr1[17]*rdv2l-1.0*incr2[17]*rdvSq4l; 
  outl[18] += incr1[18]*rdv2l-1.0*incr2[18]*rdvSq4l; 
  outl[19] += incr2[19]*rdvSq4l-1.0*incr1[19]*rdv2l; 

  const double vMuMid = 2.0*wl[2]; 
  return std::abs(vMuMid); 
} 
