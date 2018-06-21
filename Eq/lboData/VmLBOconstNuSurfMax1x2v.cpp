#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*2]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2l = 2.0/dxvl[1]; 
  double rdv2r = 2.0/dxvr[1]; 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  const double *nuUx = &nuU[0]; 

  double favg[4]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 

  double fjump[4]; 
  fjump[0] = nu*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nu*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nu*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nu*vMuMidMax*(fl[3]-(1*fr[3])); 

  double drBar[2]; 
  drBar[0] = 1.414213562373095*wl[1]*nu-1.0*nuUx[0]; 
  drBar[1] = -1.0*nuUx[1]; 

  double Gdiff[4]; 
  Gdiff[0] = nuVtSq[0]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0])+(1.590990257669731*fr[1]-1.590990257669731*fl[1])*nuVtSq[1]; 
  Gdiff[1] = nuVtSq[1]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0])+nuVtSq[0]*(1.590990257669731*fr[1]-1.590990257669731*fl[1]); 
  Gdiff[3] = nuVtSq[0]*(1.590990257669731*fr[3]-1.590990257669731*fl[3]); 

  double Ghat[4]; 
  for(unsigned short int i=0; i<4; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(0.4330127018922193*dxvl[1]*favg[2]+0.25*favg[0]*dxvl[1])*nu-0.8660254037844386*fjump[2]+drBar[0]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*drBar[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+0.25*dxvl[1]*favg[1]*nu+drBar[1]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[1]+0.3535533905932737*drBar[0]*favg[1]; 
  Ghat[3] = Gdiff[3]*rdv+0.25*dxvl[1]*favg[3]*nu-0.5*fjump[3]+0.3535533905932737*drBar[0]*favg[3]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 

  double incr2[4]; 
  incr2[2] = nuVtSq[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*(fr[1]+fl[1])*nuVtSq[1]; 

  outr[0] += incr1[0]*rdv2r; 
  outr[1] += incr1[1]*rdv2r; 
  outr[2] += incr2[2]*rdvSq4r+incr1[2]*rdv2r; 
  outr[3] += incr1[3]*rdv2r; 

  outl[0] += -1.0*incr1[0]*rdv2l; 
  outl[1] += -1.0*incr1[1]*rdv2l; 
  outl[2] += incr1[2]*rdv2l-1.0*incr2[2]*rdvSq4l; 
  outl[3] += -1.0*incr1[3]*rdv2l; 

  const double vMuMid = wl[1]-(0.7071067811865475*nuUx[0])/nu; 
  return std::abs(vMuMid); 
} 
double VmLBOconstNuSurf1x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*3]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2l = 2.0/dxvl[1]; 
  double rdv2r = 2.0/dxvr[1]; 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  const double *nuUx = &nuU[0]; 

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

  double drBar[3]; 
  drBar[0] = 1.414213562373095*wl[1]*nu-1.0*nuUx[0]; 
  drBar[1] = -1.0*nuUx[1]; 
  drBar[2] = -1.0*nuUx[2]; 

  double Gdiff[10]; 
  Gdiff[0] = nuVtSq[0]*(1.897366596101028*fr[8]-1.897366596101028*fl[8]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+nuVtSq[2]*(2.651650429449552*fr[7]-2.651650429449552*fl[7])+nuVtSq[1]*((-3.368048396326869*(fr[4]+fl[4]))+2.651650429449552*fr[1]-2.651650429449552*fl[1]); 
  Gdiff[1] = nuVtSq[1]*(1.897366596101028*fr[8]-1.897366596101028*fl[8]+2.371708245126284*fr[7]-2.371708245126284*fl[7]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+nuVtSq[2]*((-3.012474066278413*(fr[4]+fl[4]))+2.371708245126284*fr[1]-2.371708245126284*fl[1])+nuVtSq[0]*((-3.368048396326869*(fr[4]+fl[4]))+2.651650429449552*fr[1]-2.651650429449552*fl[1]); 
  Gdiff[3] = nuVtSq[0]*((-3.368048396326869*(fr[6]+fl[6]))+2.651650429449552*fr[3]-2.651650429449552*fl[3])+nuVtSq[1]*(2.651650429449552*fr[5]-2.651650429449552*fl[5]); 
  Gdiff[5] = nuVtSq[1]*((-3.368048396326869*(fr[6]+fl[6]))+2.651650429449552*fr[3]-2.651650429449552*fl[3])+nuVtSq[0]*(2.651650429449552*fr[5]-2.651650429449552*fl[5])+nuVtSq[2]*(2.371708245126284*fr[5]-2.371708245126284*fl[5]); 
  Gdiff[7] = nuVtSq[2]*(1.897366596101028*fr[8]-1.897366596101028*fl[8]+1.694077317947346*fr[7]-1.694077317947346*fl[7]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+nuVtSq[0]*(2.651650429449552*fr[7]-2.651650429449552*fl[7])+nuVtSq[1]*((-3.012474066278413*(fr[4]+fl[4]))+2.371708245126284*fr[1]-2.371708245126284*fl[1]); 
  Gdiff[9] = nuVtSq[0]*(2.651650429449552*fr[9]-2.651650429449552*fl[9]); 

  double Ghat[10]; 
  for(unsigned short int i=0; i<10; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(0.5590169943749475*dxvl[1]*favg[8]+0.4330127018922193*dxvl[1]*favg[2]+0.25*favg[0]*dxvl[1])*nu-1.118033988749895*fjump[8]+drBar[0]*(0.7905694150420947*favg[8]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*drBar[2]*favg[7]+drBar[1]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+(0.4330127018922193*dxvl[1]*favg[4]+0.25*dxvl[1]*favg[1])*nu+drBar[1]*(0.7905694150420947*favg[8]+0.3162277660168379*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[4]+drBar[0]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+drBar[2]*(0.5477225575051661*favg[4]+0.3162277660168379*favg[1])-0.5*fjump[1]; 
  Ghat[3] = Gdiff[3]*rdv+(0.4330127018922193*dxvl[1]*favg[6]+0.25*dxvl[1]*favg[3])*nu-0.8660254037844386*fjump[6]+drBar[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[3])+0.3535533905932737*drBar[1]*favg[5]-0.5*fjump[3]; 
  Ghat[5] = Gdiff[5]*rdv+0.25*dxvl[1]*favg[5]*nu+drBar[1]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[5]+0.3162277660168379*drBar[2]*favg[5]+0.3535533905932737*drBar[0]*favg[5]; 
  Ghat[7] = Gdiff[7]*rdv+0.25*dxvl[1]*favg[7]*nu+drBar[2]*(0.7905694150420947*favg[8]+0.2258769757263128*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[7]+0.3535533905932737*drBar[0]*favg[7]+drBar[1]*(0.5477225575051661*favg[4]+0.3162277660168379*favg[1]); 
  Ghat[9] = Gdiff[9]*rdv+0.25*dxvl[1]*favg[9]*nu-0.5*fjump[9]+0.3535533905932737*drBar[0]*favg[9]; 

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
  incr2[2] = nuVtSq[0]*(0.2995357736356374*(fr[8]+fl[8])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*nuVtSq[2]*(fr[7]+fl[7])+nuVtSq[1]*((-0.430893194785552*fr[4])+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[4] = nuVtSq[1]*(0.2995357736356374*(fr[8]+fl[8])+0.273861278752583*(fr[7]+fl[7])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSq[2]*((-0.3854025898330209*fr[4])+0.3854025898330209*fl[4]+0.273861278752583*(fr[1]+fl[1]))+nuVtSq[0]*((-0.430893194785552*fr[4])+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[6] = nuVtSq[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3]))+0.3061862178478971*nuVtSq[1]*(fr[5]+fl[5]); 
  incr2[8] = nuVtSq[0]*((-1.160097062884178*(fr[8]+fl[8]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*nuVtSq[2]*(fr[7]+fl[7])+nuVtSq[1]*(1.668842167398551*fr[4]-1.668842167398551*fl[4]-1.185854122563142*(fr[1]+fl[1])); 

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

  const double vMuMid = (0.7905694150420947*nuUx[2])/nu-(0.7071067811865475*nuUx[0])/nu+wl[1]; 
  return std::abs(vMuMid); 
} 
double VmLBOconstNuSurf1x2vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*4]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[4]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2l = 2.0/dxvl[1]; 
  double rdv2r = 2.0/dxvr[1]; 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  const double *nuUx = &nuU[0]; 

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
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = -1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 

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
  fjump[17] = nu*vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = nu*vMuMidMax*(fl[18]-(-1*fr[18])); 
  fjump[19] = nu*vMuMidMax*(fl[19]-(1*fr[19])); 

  double drBar[4]; 
  drBar[0] = 1.414213562373095*wl[1]*nu-1.0*nuUx[0]; 
  drBar[1] = -1.0*nuUx[1]; 
  drBar[2] = -1.0*nuUx[2]; 
  drBar[3] = -1.0*nuUx[3]; 

  double Gdiff[20]; 
  Gdiff[0] = nuVtSq[0]*((-1.929292090055312*(fr[18]+fl[18]))+4.397542371171649*fr[8]-4.397542371171649*fl[8]-5.473078644031159*(fr[2]+fl[2])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+nuVtSq[3]*(3.866990209613929*fr[17]-3.866990209613929*fl[17])+nuVtSq[1]*(4.39754237117165*fr[12]-4.39754237117165*fl[12]-5.473078644031159*(fr[4]+fl[4])+3.866990209613929*fr[1]-3.866990209613929*fl[1])+nuVtSq[2]*((-5.473078644031159*(fr[11]+fl[11]))+3.866990209613929*fr[7]-3.866990209613929*fl[7]); 
  Gdiff[1] = nuVtSq[1]*((-1.929292090055312*(fr[18]+fl[18]))-4.89527035770242*(fr[11]+fl[11])+4.397542371171649*fr[8]-4.397542371171649*fl[8]+3.458741190809164*fr[7]-3.458741190809164*fl[7]-5.473078644031159*(fr[2]+fl[2])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+nuVtSq[2]*(3.396416424888148*fr[17]-3.396416424888148*fl[17]+3.933281470350169*fr[12]-3.933281470350169*fl[12]-4.89527035770242*(fr[4]+fl[4])+3.458741190809164*fr[1]-3.458741190809164*fl[1])+nuVtSq[0]*(4.39754237117165*fr[12]-4.39754237117165*fl[12]-5.473078644031159*(fr[4]+fl[4])+3.866990209613929*fr[1]-3.866990209613929*fl[1])+nuVtSq[3]*((-4.807060063166761*(fr[11]+fl[11]))+3.396416424888148*fr[7]-3.396416424888148*fl[7]); 
  Gdiff[3] = nuVtSq[0]*(4.39754237117165*fr[14]-4.39754237117165*fl[14]-5.473078644031159*(fr[6]+fl[6])+3.866990209613929*fr[3]-3.866990209613929*fl[3])+nuVtSq[2]*(3.86699020961393*fr[13]-3.86699020961393*fl[13])+nuVtSq[1]*((-5.473078644031159*(fr[10]+fl[10]))+3.866990209613929*fr[5]-3.866990209613929*fl[5]); 
  Gdiff[5] = nuVtSq[1]*(4.39754237117165*fr[14]-4.39754237117165*fl[14]+3.458741190809164*fr[13]-3.458741190809164*fl[13]-5.473078644031159*(fr[6]+fl[6])+3.866990209613929*fr[3]-3.866990209613929*fl[3])+nuVtSq[3]*(3.396416424888148*fr[13]-3.396416424888148*fl[13])+nuVtSq[2]*((-4.89527035770242*(fr[10]+fl[10]))+3.458741190809164*fr[5]-3.458741190809164*fl[5])+nuVtSq[0]*((-5.473078644031159*(fr[10]+fl[10]))+3.866990209613929*fr[5]-3.866990209613929*fl[5]); 
  Gdiff[7] = nuVtSq[2]*((-1.929292090055312*(fr[18]+fl[18]))-3.496621684073157*(fr[11]+fl[11])+4.397542371171649*fr[8]-4.397542371171649*fl[8]+2.470529422006546*fr[7]-2.470529422006546*fl[7]-5.473078644031159*(fr[2]+fl[2])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+nuVtSq[1]*(3.396416424888148*fr[17]-3.396416424888148*fl[17]+3.933281470350169*fr[12]-3.933281470350169*fl[12]-4.89527035770242*(fr[4]+fl[4])+3.458741190809164*fr[1]-3.458741190809164*fl[1])+nuVtSq[3]*(2.305827460539443*fr[17]-2.305827460539443*fl[17]+3.86240572873861*fr[12]-3.86240572873861*fl[12]-4.807060063166761*(fr[4]+fl[4])+3.396416424888148*fr[1]-3.396416424888148*fl[1])+nuVtSq[0]*((-5.473078644031159*(fr[11]+fl[11]))+3.866990209613929*fr[7]-3.866990209613929*fl[7]); 
  Gdiff[9] = nuVtSq[0]*((-5.473078644031159*(fr[16]+fl[16]))+3.866990209613929*fr[9]-3.866990209613929*fl[9])+nuVtSq[1]*(3.86699020961393*fr[15]-3.86699020961393*fl[15]); 
  Gdiff[13] = nuVtSq[2]*(4.397542371171649*fr[14]-4.397542371171649*fl[14]+2.470529422006546*fr[13]-2.470529422006546*fl[13]-5.473078644031159*(fr[6]+fl[6])+3.86699020961393*fr[3]-3.86699020961393*fl[3])+nuVtSq[0]*(3.866990209613929*fr[13]-3.866990209613929*fl[13])+nuVtSq[3]*((-4.807060063166761*(fr[10]+fl[10]))+3.396416424888148*fr[5]-3.396416424888148*fl[5])+nuVtSq[1]*((-4.89527035770242*(fr[10]+fl[10]))+3.458741190809164*fr[5]-3.458741190809164*fl[5]); 
  Gdiff[15] = nuVtSq[1]*((-5.473078644031159*(fr[16]+fl[16]))+3.86699020961393*fr[9]-3.86699020961393*fl[9])+nuVtSq[0]*(3.866990209613929*fr[15]-3.866990209613929*fl[15])+nuVtSq[2]*(3.458741190809164*fr[15]-3.458741190809164*fl[15]); 
  Gdiff[17] = nuVtSq[3]*((-1.929292090055312*(fr[18]+fl[18]))-3.263513571801613*(fr[11]+fl[11])+4.397542371171649*fr[8]-4.397542371171649*fl[8]+2.305827460539443*fr[7]-2.305827460539443*fl[7]-5.473078644031159*(fr[2]+fl[2])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+nuVtSq[0]*(3.866990209613929*fr[17]-3.866990209613929*fl[17])+nuVtSq[2]*(2.305827460539443*fr[17]-2.305827460539443*fl[17]+3.86240572873861*fr[12]-3.86240572873861*fl[12]-4.807060063166761*(fr[4]+fl[4])+3.396416424888148*fr[1]-3.396416424888148*fl[1])+nuVtSq[1]*((-4.807060063166761*(fr[11]+fl[11]))+3.396416424888148*fr[7]-3.396416424888148*fl[7]); 
  Gdiff[19] = nuVtSq[0]*(3.866990209613929*fr[19]-3.866990209613929*fl[19]); 

  double Ghat[20]; 
  for(unsigned short int i=0; i<20; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(0.6614378277661477*dxvl[1]*favg[18]+0.5590169943749475*dxvl[1]*favg[8]+0.4330127018922193*dxvl[1]*favg[2]+0.25*favg[0]*dxvl[1])*nu-1.322875655532295*fjump[18]+drBar[0]*(0.9354143466934851*favg[18]+0.7905694150420947*favg[8]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*drBar[3]*favg[17]+drBar[1]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+drBar[2]*(0.6123724356957944*favg[11]+0.3535533905932737*favg[7])-1.118033988749895*fjump[8]-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+(0.5590169943749476*dxvl[1]*favg[12]+0.4330127018922193*dxvl[1]*favg[4]+0.25*dxvl[1]*favg[1])*nu+drBar[1]*(0.9354143466934851*favg[18]+0.5477225575051661*favg[11]+0.7905694150420947*favg[8]+0.3162277660168379*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+drBar[2]*(0.3105295017040592*favg[17]+0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])-1.118033988749895*fjump[12]+drBar[0]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+drBar[3]*(0.5378528742004769*favg[11]+0.3105295017040592*favg[7])-0.8660254037844386*fjump[4]-0.5*fjump[1]; 
  Ghat[3] = Gdiff[3]*rdv+(0.5590169943749476*dxvl[1]*favg[14]+0.4330127018922193*dxvl[1]*favg[6]+0.25*dxvl[1]*favg[3])*nu-1.118033988749895*fjump[14]+drBar[0]*(0.7905694150420948*favg[14]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])+0.3535533905932737*drBar[2]*favg[13]+drBar[1]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[5])-0.8660254037844386*fjump[6]-0.5*fjump[3]; 
  Ghat[5] = Gdiff[5]*rdv+(0.4330127018922193*dxvl[1]*favg[10]+0.25*dxvl[1]*favg[5])*nu+drBar[1]*(0.7905694150420948*favg[14]+0.3162277660168379*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])+0.3105295017040593*drBar[3]*favg[13]-0.8660254037844386*fjump[10]+drBar[0]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[5])+drBar[2]*(0.5477225575051661*favg[10]+0.3162277660168379*favg[5])-0.5*fjump[5]; 
  Ghat[7] = Gdiff[7]*rdv+(0.4330127018922194*dxvl[1]*favg[11]+0.25*dxvl[1]*favg[7])*nu+drBar[2]*(0.9354143466934851*favg[18]+0.3912303982179757*favg[11]+0.7905694150420947*favg[8]+0.2258769757263128*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+drBar[1]*(0.3105295017040592*favg[17]+0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])+drBar[3]*(0.210818510677892*favg[17]+0.6943650748294133*favg[12]+0.537852874200477*favg[4]+0.3105295017040592*favg[1])-0.8660254037844387*fjump[11]+drBar[0]*(0.6123724356957944*favg[11]+0.3535533905932737*favg[7])-0.5*fjump[7]; 
  Ghat[9] = Gdiff[9]*rdv+(0.4330127018922194*dxvl[1]*favg[16]+0.25*dxvl[1]*favg[9])*nu-0.8660254037844387*fjump[16]+drBar[0]*(0.6123724356957944*favg[16]+0.3535533905932737*favg[9])+0.3535533905932737*drBar[1]*favg[15]-0.5*fjump[9]; 
  Ghat[13] = Gdiff[13]*rdv+0.25*dxvl[1]*favg[13]*nu+drBar[2]*(0.7905694150420947*favg[14]+0.2258769757263128*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[13]+0.3535533905932737*drBar[0]*favg[13]+drBar[1]*(0.5477225575051661*favg[10]+0.3162277660168379*favg[5])+drBar[3]*(0.5378528742004769*favg[10]+0.3105295017040593*favg[5]); 
  Ghat[15] = Gdiff[15]*rdv+0.25*dxvl[1]*favg[15]*nu+drBar[1]*(0.6123724356957944*favg[16]+0.3535533905932737*favg[9])-0.5*fjump[15]+0.3162277660168379*drBar[2]*favg[15]+0.3535533905932737*drBar[0]*favg[15]; 
  Ghat[17] = Gdiff[17]*rdv+0.25*dxvl[1]*favg[17]*nu+drBar[3]*(0.9354143466934851*favg[18]+0.3651483716701107*favg[11]+0.7905694150420947*favg[8]+0.210818510677892*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[17]+0.3535533905932737*drBar[0]*favg[17]+drBar[2]*(0.210818510677892*favg[17]+0.6943650748294133*favg[12]+0.537852874200477*favg[4]+0.3105295017040592*favg[1])+drBar[1]*(0.5378528742004769*favg[11]+0.3105295017040592*favg[7]); 
  Ghat[19] = Gdiff[19]*rdv+0.25*dxvl[1]*favg[19]*nu-0.5*fjump[19]+0.3535533905932737*drBar[0]*favg[19]; 

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
  incr1[17] = -0.5*Ghat[17]; 
  incr1[18] = 1.322875655532295*Ghat[0]; 
  incr1[19] = -0.5*Ghat[19]; 

  double incr2[20]; 
  incr2[2] = nuVtSq[0]*((-0.2314550249431378*fr[18])+0.2314550249431378*fl[18]+0.4279082480509108*(fr[8]+fl[8])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*nuVtSq[3]*(fr[17]+fl[17])+nuVtSq[1]*(0.4279082480509107*(fr[12]+fl[12])-0.4640388251536715*fr[4]+0.4640388251536715*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSq[2]*((-0.4640388251536715*fr[11])+0.4640388251536715*fl[11]+0.3061862178478971*(fr[7]+fl[7])); 
  incr2[4] = nuVtSq[1]*((-0.2314550249431378*fr[18])+0.2314550249431378*fl[18]-0.4150489428970995*fr[11]+0.4150489428970995*fl[11]+0.4279082480509108*(fr[8]+fl[8])+0.273861278752583*(fr[7]+fl[7])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSq[2]*(0.2689264371002384*(fr[17]+fl[17])+0.3827327723098714*(fr[12]+fl[12])-0.4150489428970995*fr[4]+0.4150489428970995*fl[4]+0.273861278752583*(fr[1]+fl[1]))+nuVtSq[0]*(0.4279082480509107*(fr[12]+fl[12])-0.4640388251536715*fr[4]+0.4640388251536715*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSq[3]*((-0.4075699709865777*fr[11])+0.4075699709865777*fl[11]+0.2689264371002384*(fr[7]+fl[7])); 
  incr2[6] = nuVtSq[0]*(0.4279082480509107*(fr[14]+fl[14])-0.4640388251536715*fr[6]+0.4640388251536715*fl[6]+0.3061862178478971*(fr[3]+fl[3]))+0.3061862178478971*nuVtSq[2]*(fr[13]+fl[13])+nuVtSq[1]*((-0.4640388251536715*fr[10])+0.4640388251536715*fl[10]+0.3061862178478971*(fr[5]+fl[5])); 
  incr2[8] = nuVtSq[0]*(0.896421457000795*fr[18]-0.896421457000795*fl[18]-1.657281518405969*(fr[8]+fl[8])+1.797214641813825*fr[2]-1.797214641813825*fl[2]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*nuVtSq[3]*(fr[17]+fl[17])+nuVtSq[1]*((-1.65728151840597*(fr[12]+fl[12]))+1.797214641813825*fr[4]-1.797214641813825*fl[4]-1.185854122563142*(fr[1]+fl[1]))+nuVtSq[2]*(1.797214641813825*fr[11]-1.797214641813825*fl[11]-1.185854122563142*(fr[7]+fl[7])); 
  incr2[10] = nuVtSq[1]*(0.4279082480509107*(fr[14]+fl[14])+0.273861278752583*(fr[13]+fl[13])-0.4640388251536715*fr[6]+0.4640388251536715*fl[6]+0.3061862178478971*(fr[3]+fl[3]))+0.2689264371002384*nuVtSq[3]*(fr[13]+fl[13])+nuVtSq[2]*((-0.4150489428970995*fr[10])+0.4150489428970995*fl[10]+0.273861278752583*(fr[5]+fl[5]))+nuVtSq[0]*((-0.4640388251536715*fr[10])+0.4640388251536715*fl[10]+0.3061862178478971*(fr[5]+fl[5])); 
  incr2[11] = nuVtSq[2]*((-0.2314550249431378*fr[18])+0.2314550249431378*fl[18]-0.2964635306407854*fr[11]+0.2964635306407854*fl[11]+0.4279082480509108*(fr[8]+fl[8])+0.1956151991089878*(fr[7]+fl[7])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSq[1]*(0.2689264371002384*(fr[17]+fl[17])+0.3827327723098713*(fr[12]+fl[12])-0.4150489428970995*fr[4]+0.4150489428970995*fl[4]+0.273861278752583*(fr[1]+fl[1]))+nuVtSq[3]*(0.1825741858350553*(fr[17]+fl[17])+0.3758361214393465*(fr[12]+fl[12])-0.4075699709865777*fr[4]+0.4075699709865777*fl[4]+0.2689264371002384*(fr[1]+fl[1]))+nuVtSq[0]*((-0.4640388251536715*fr[11])+0.4640388251536715*fl[11]+0.3061862178478971*(fr[7]+fl[7])); 
  incr2[12] = nuVtSq[1]*(0.8964214570007949*fr[18]-0.8964214570007949*fl[18]+1.60747764370146*fr[11]-1.60747764370146*fl[11]-1.65728151840597*(fr[8]+fl[8])-1.060660171779821*(fr[7]+fl[7])+1.797214641813825*fr[2]-1.797214641813825*fl[2]-1.185854122563142*(fr[0]+fl[0]))+nuVtSq[2]*((-1.04154761224412*(fr[17]+fl[17]))-1.482317653203927*(fr[12]+fl[12])+1.60747764370146*fr[4]-1.60747764370146*fl[4]-1.060660171779821*(fr[1]+fl[1]))+nuVtSq[0]*((-1.657281518405969*(fr[12]+fl[12]))+1.797214641813825*fr[4]-1.797214641813825*fl[4]-1.185854122563142*(fr[1]+fl[1]))+nuVtSq[3]*(1.578511710045255*fr[11]-1.578511710045255*fl[11]-1.04154761224412*(fr[7]+fl[7])); 
  incr2[14] = nuVtSq[0]*((-1.657281518405969*(fr[14]+fl[14]))+1.797214641813825*fr[6]-1.797214641813825*fl[6]-1.185854122563142*(fr[3]+fl[3]))-1.185854122563142*nuVtSq[2]*(fr[13]+fl[13])+nuVtSq[1]*(1.797214641813825*fr[10]-1.797214641813825*fl[10]-1.185854122563142*(fr[5]+fl[5])); 
  incr2[16] = nuVtSq[0]*((-0.4640388251536715*fr[16])+0.4640388251536715*fl[16]+0.3061862178478971*(fr[9]+fl[9]))+0.3061862178478971*nuVtSq[1]*(fr[15]+fl[15]); 
  incr2[18] = nuVtSq[0]*((-2.121320343559642*fr[18])+2.121320343559642*fl[18]+3.921843874378478*(fr[8]+fl[8])-4.252986083330156*fr[2]+4.252986083330156*fl[2]+2.806243040080455*(fr[0]+fl[0]))+2.806243040080455*nuVtSq[3]*(fr[17]+fl[17])+nuVtSq[1]*(3.921843874378477*(fr[12]+fl[12])-4.252986083330156*fr[4]+4.252986083330156*fl[4]+2.806243040080455*(fr[1]+fl[1]))+nuVtSq[2]*((-4.252986083330157*fr[11])+4.252986083330157*fl[11]+2.806243040080455*(fr[7]+fl[7])); 

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
  outr[17] += incr1[17]*rdv2r; 
  outr[18] += incr2[18]*rdvSq4r+incr1[18]*rdv2r; 
  outr[19] += incr1[19]*rdv2r; 

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
  outl[17] += -1.0*incr1[17]*rdv2l; 
  outl[18] += incr1[18]*rdv2l-1.0*incr2[18]*rdvSq4l; 
  outl[19] += -1.0*incr1[19]*rdv2l; 

  const double vMuMid = (0.7905694150420947*nuUx[2])/nu-(0.7071067811865475*nuUx[0])/nu+wl[1]; 
  return std::abs(vMuMid); 
} 
double VmLBOconstNuSurf1x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*2]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2l = 2.0/dxvl[2]; 
  double rdv2r = 2.0/dxvr[2]; 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  const double *nuUy = &nuU[2]; 

  double favg[4]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 

  double fjump[4]; 
  fjump[0] = nu*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nu*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nu*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nu*vMuMidMax*(fl[3]-(-1*fr[3])); 

  double drBar[2]; 
  drBar[0] = 1.414213562373095*wl[2]*nu-1.0*nuUy[0]; 
  drBar[1] = -1.0*nuUy[1]; 

  double Gdiff[4]; 
  Gdiff[0] = nuVtSq[0]*((-1.530931089239486*(fr[3]+fl[3]))+1.590990257669731*fr[0]-1.590990257669731*fl[0])+(1.590990257669731*fr[1]-1.590990257669731*fl[1])*nuVtSq[1]; 
  Gdiff[1] = nuVtSq[1]*((-1.530931089239486*(fr[3]+fl[3]))+1.590990257669731*fr[0]-1.590990257669731*fl[0])+nuVtSq[0]*(1.590990257669731*fr[1]-1.590990257669731*fl[1]); 
  Gdiff[2] = nuVtSq[0]*(1.590990257669731*fr[2]-1.590990257669731*fl[2]); 

  double Ghat[4]; 
  for(unsigned short int i=0; i<4; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(0.4330127018922193*dxvl[2]*favg[3]+0.25*favg[0]*dxvl[2])*nu-0.8660254037844386*fjump[3]+drBar[0]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+0.3535533905932737*drBar[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+0.25*favg[1]*dxvl[2]*nu+drBar[1]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.5*fjump[1]+0.3535533905932737*drBar[0]*favg[1]; 
  Ghat[2] = Gdiff[2]*rdv+0.25*dxvl[2]*favg[2]*nu-0.5*fjump[2]+0.3535533905932737*drBar[0]*favg[2]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 

  double incr2[4]; 
  incr2[3] = nuVtSq[0]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*(fr[1]+fl[1])*nuVtSq[1]; 

  outr[0] += incr1[0]*rdv2r; 
  outr[1] += incr1[1]*rdv2r; 
  outr[2] += incr1[2]*rdv2r; 
  outr[3] += incr2[3]*rdvSq4r+incr1[3]*rdv2r; 

  outl[0] += -1.0*incr1[0]*rdv2l; 
  outl[1] += -1.0*incr1[1]*rdv2l; 
  outl[2] += -1.0*incr1[2]*rdv2l; 
  outl[3] += incr1[3]*rdv2l-1.0*incr2[3]*rdvSq4l; 

  const double vMuMid = wl[2]-(0.7071067811865475*nuUy[0])/nu; 
  return std::abs(vMuMid); 
} 
double VmLBOconstNuSurf1x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*3]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2l = 2.0/dxvl[2]; 
  double rdv2r = 2.0/dxvr[2]; 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  const double *nuUy = &nuU[3]; 

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

  double drBar[3]; 
  drBar[0] = 1.414213562373095*wl[2]*nu-1.0*nuUy[0]; 
  drBar[1] = -1.0*nuUy[1]; 
  drBar[2] = -1.0*nuUy[2]; 

  double Gdiff[10]; 
  Gdiff[0] = nuVtSq[0]*(1.897366596101028*fr[9]-1.897366596101028*fl[9]-3.368048396326869*(fr[3]+fl[3])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+nuVtSq[2]*(2.651650429449552*fr[7]-2.651650429449552*fl[7])+nuVtSq[1]*((-3.368048396326869*(fr[5]+fl[5]))+2.651650429449552*fr[1]-2.651650429449552*fl[1]); 
  Gdiff[1] = nuVtSq[1]*(1.897366596101028*fr[9]-1.897366596101028*fl[9]+2.371708245126284*fr[7]-2.371708245126284*fl[7]-3.368048396326869*(fr[3]+fl[3])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+nuVtSq[2]*((-3.012474066278413*(fr[5]+fl[5]))+2.371708245126284*fr[1]-2.371708245126284*fl[1])+nuVtSq[0]*((-3.368048396326869*(fr[5]+fl[5]))+2.651650429449552*fr[1]-2.651650429449552*fl[1]); 
  Gdiff[2] = nuVtSq[0]*((-3.368048396326869*(fr[6]+fl[6]))+2.651650429449552*fr[2]-2.651650429449552*fl[2])+nuVtSq[1]*(2.651650429449552*fr[4]-2.651650429449552*fl[4]); 
  Gdiff[4] = nuVtSq[1]*((-3.368048396326869*(fr[6]+fl[6]))+2.651650429449552*fr[2]-2.651650429449552*fl[2])+nuVtSq[0]*(2.651650429449552*fr[4]-2.651650429449552*fl[4])+nuVtSq[2]*(2.371708245126284*fr[4]-2.371708245126284*fl[4]); 
  Gdiff[7] = nuVtSq[2]*(1.897366596101028*fr[9]-1.897366596101028*fl[9]+1.694077317947346*fr[7]-1.694077317947346*fl[7]-3.368048396326869*(fr[3]+fl[3])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+nuVtSq[0]*(2.651650429449552*fr[7]-2.651650429449552*fl[7])+nuVtSq[1]*((-3.012474066278413*(fr[5]+fl[5]))+2.371708245126284*fr[1]-2.371708245126284*fl[1]); 
  Gdiff[8] = nuVtSq[0]*(2.651650429449552*fr[8]-2.651650429449552*fl[8]); 

  double Ghat[10]; 
  for(unsigned short int i=0; i<10; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(0.5590169943749475*dxvl[2]*favg[9]+0.4330127018922193*dxvl[2]*favg[3]+0.25*favg[0]*dxvl[2])*nu-1.118033988749895*fjump[9]+drBar[0]*(0.7905694150420947*favg[9]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+0.3535533905932737*drBar[2]*favg[7]+drBar[1]*(0.6123724356957944*favg[5]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+(0.4330127018922193*dxvl[2]*favg[5]+0.25*favg[1]*dxvl[2])*nu+drBar[1]*(0.7905694150420947*favg[9]+0.3162277660168379*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[5]+drBar[0]*(0.6123724356957944*favg[5]+0.3535533905932737*favg[1])+drBar[2]*(0.5477225575051661*favg[5]+0.3162277660168379*favg[1])-0.5*fjump[1]; 
  Ghat[2] = Gdiff[2]*rdv+(0.4330127018922193*dxvl[2]*favg[6]+0.25*dxvl[2]*favg[2])*nu-0.8660254037844386*fjump[6]+drBar[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[2])+0.3535533905932737*drBar[1]*favg[4]-0.5*fjump[2]; 
  Ghat[4] = Gdiff[4]*rdv+0.25*dxvl[2]*favg[4]*nu+drBar[1]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-0.5*fjump[4]+0.3162277660168379*drBar[2]*favg[4]+0.3535533905932737*drBar[0]*favg[4]; 
  Ghat[7] = Gdiff[7]*rdv+0.25*dxvl[2]*favg[7]*nu+drBar[2]*(0.7905694150420947*favg[9]+0.2258769757263128*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.5*fjump[7]+0.3535533905932737*drBar[0]*favg[7]+drBar[1]*(0.5477225575051661*favg[5]+0.3162277660168379*favg[1]); 
  Ghat[8] = Gdiff[8]*rdv+0.25*dxvl[2]*favg[8]*nu-0.5*fjump[8]+0.3535533905932737*drBar[0]*favg[8]; 

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
  incr2[3] = nuVtSq[0]*(0.2995357736356374*(fr[9]+fl[9])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*nuVtSq[2]*(fr[7]+fl[7])+nuVtSq[1]*((-0.430893194785552*fr[5])+0.430893194785552*fl[5]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[5] = nuVtSq[1]*(0.2995357736356374*(fr[9]+fl[9])+0.273861278752583*(fr[7]+fl[7])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSq[2]*((-0.3854025898330209*fr[5])+0.3854025898330209*fl[5]+0.273861278752583*(fr[1]+fl[1]))+nuVtSq[0]*((-0.430893194785552*fr[5])+0.430893194785552*fl[5]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[6] = nuVtSq[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[2]+fl[2]))+0.3061862178478971*nuVtSq[1]*(fr[4]+fl[4]); 
  incr2[9] = nuVtSq[0]*((-1.160097062884178*(fr[9]+fl[9]))+1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*nuVtSq[2]*(fr[7]+fl[7])+nuVtSq[1]*(1.668842167398551*fr[5]-1.668842167398551*fl[5]-1.185854122563142*(fr[1]+fl[1])); 

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

  const double vMuMid = (0.7905694150420947*nuUy[2])/nu-(0.7071067811865475*nuUy[0])/nu+wl[2]; 
  return std::abs(vMuMid); 
} 
double VmLBOconstNuSurf1x2vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*4]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[4]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2l = 2.0/dxvl[2]; 
  double rdv2r = 2.0/dxvr[2]; 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  const double *nuUy = &nuU[4]; 

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
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = -1*fr[19]+fl[19]; 

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
  fjump[17] = nu*vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = nu*vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = nu*vMuMidMax*(fl[19]-(-1*fr[19])); 

  double drBar[4]; 
  drBar[0] = 1.414213562373095*wl[2]*nu-1.0*nuUy[0]; 
  drBar[1] = -1.0*nuUy[1]; 
  drBar[2] = -1.0*nuUy[2]; 
  drBar[3] = -1.0*nuUy[3]; 

  double Gdiff[20]; 
  Gdiff[0] = nuVtSq[0]*((-1.929292090055312*(fr[19]+fl[19]))+4.397542371171649*fr[9]-4.397542371171649*fl[9]-5.473078644031159*(fr[3]+fl[3])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+nuVtSq[3]*(3.866990209613929*fr[17]-3.866990209613929*fl[17])+nuVtSq[1]*(4.39754237117165*fr[15]-4.39754237117165*fl[15]-5.473078644031159*(fr[5]+fl[5])+3.866990209613929*fr[1]-3.866990209613929*fl[1])+nuVtSq[2]*((-5.473078644031159*(fr[13]+fl[13]))+3.866990209613929*fr[7]-3.866990209613929*fl[7]); 
  Gdiff[1] = nuVtSq[1]*((-1.929292090055312*(fr[19]+fl[19]))-4.89527035770242*(fr[13]+fl[13])+4.397542371171649*fr[9]-4.397542371171649*fl[9]+3.458741190809164*fr[7]-3.458741190809164*fl[7]-5.473078644031159*(fr[3]+fl[3])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+nuVtSq[2]*(3.396416424888148*fr[17]-3.396416424888148*fl[17]+3.933281470350169*fr[15]-3.933281470350169*fl[15]-4.89527035770242*(fr[5]+fl[5])+3.458741190809164*fr[1]-3.458741190809164*fl[1])+nuVtSq[0]*(4.39754237117165*fr[15]-4.39754237117165*fl[15]-5.473078644031159*(fr[5]+fl[5])+3.866990209613929*fr[1]-3.866990209613929*fl[1])+nuVtSq[3]*((-4.807060063166761*(fr[13]+fl[13]))+3.396416424888148*fr[7]-3.396416424888148*fl[7]); 
  Gdiff[2] = nuVtSq[0]*(4.39754237117165*fr[16]-4.39754237117165*fl[16]-5.473078644031159*(fr[6]+fl[6])+3.866990209613929*fr[2]-3.866990209613929*fl[2])+nuVtSq[2]*(3.86699020961393*fr[11]-3.86699020961393*fl[11])+nuVtSq[1]*((-5.473078644031159*(fr[10]+fl[10]))+3.866990209613929*fr[4]-3.866990209613929*fl[4]); 
  Gdiff[4] = nuVtSq[1]*(4.39754237117165*fr[16]-4.39754237117165*fl[16]+3.458741190809164*fr[11]-3.458741190809164*fl[11]-5.473078644031159*(fr[6]+fl[6])+3.866990209613929*fr[2]-3.866990209613929*fl[2])+nuVtSq[3]*(3.396416424888148*fr[11]-3.396416424888148*fl[11])+nuVtSq[2]*((-4.89527035770242*(fr[10]+fl[10]))+3.458741190809164*fr[4]-3.458741190809164*fl[4])+nuVtSq[0]*((-5.473078644031159*(fr[10]+fl[10]))+3.866990209613929*fr[4]-3.866990209613929*fl[4]); 
  Gdiff[7] = nuVtSq[2]*((-1.929292090055312*(fr[19]+fl[19]))-3.496621684073157*(fr[13]+fl[13])+4.397542371171649*fr[9]-4.397542371171649*fl[9]+2.470529422006546*fr[7]-2.470529422006546*fl[7]-5.473078644031159*(fr[3]+fl[3])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+nuVtSq[1]*(3.396416424888148*fr[17]-3.396416424888148*fl[17]+3.933281470350169*fr[15]-3.933281470350169*fl[15]-4.89527035770242*(fr[5]+fl[5])+3.458741190809164*fr[1]-3.458741190809164*fl[1])+nuVtSq[3]*(2.305827460539443*fr[17]-2.305827460539443*fl[17]+3.86240572873861*fr[15]-3.86240572873861*fl[15]-4.807060063166761*(fr[5]+fl[5])+3.396416424888148*fr[1]-3.396416424888148*fl[1])+nuVtSq[0]*((-5.473078644031159*(fr[13]+fl[13]))+3.866990209613929*fr[7]-3.866990209613929*fl[7]); 
  Gdiff[8] = nuVtSq[0]*((-5.473078644031159*(fr[14]+fl[14]))+3.866990209613929*fr[8]-3.866990209613929*fl[8])+nuVtSq[1]*(3.86699020961393*fr[12]-3.86699020961393*fl[12]); 
  Gdiff[11] = nuVtSq[2]*(4.397542371171649*fr[16]-4.397542371171649*fl[16]+2.470529422006546*fr[11]-2.470529422006546*fl[11]-5.473078644031159*(fr[6]+fl[6])+3.86699020961393*fr[2]-3.86699020961393*fl[2])+nuVtSq[0]*(3.866990209613929*fr[11]-3.866990209613929*fl[11])+nuVtSq[3]*((-4.807060063166761*(fr[10]+fl[10]))+3.396416424888148*fr[4]-3.396416424888148*fl[4])+nuVtSq[1]*((-4.89527035770242*(fr[10]+fl[10]))+3.458741190809164*fr[4]-3.458741190809164*fl[4]); 
  Gdiff[12] = nuVtSq[1]*((-5.473078644031159*(fr[14]+fl[14]))+3.86699020961393*fr[8]-3.86699020961393*fl[8])+nuVtSq[0]*(3.866990209613929*fr[12]-3.866990209613929*fl[12])+nuVtSq[2]*(3.458741190809164*fr[12]-3.458741190809164*fl[12]); 
  Gdiff[17] = nuVtSq[3]*((-1.929292090055312*(fr[19]+fl[19]))-3.263513571801613*(fr[13]+fl[13])+4.397542371171649*fr[9]-4.397542371171649*fl[9]+2.305827460539443*fr[7]-2.305827460539443*fl[7]-5.473078644031159*(fr[3]+fl[3])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+nuVtSq[0]*(3.866990209613929*fr[17]-3.866990209613929*fl[17])+nuVtSq[2]*(2.305827460539443*fr[17]-2.305827460539443*fl[17]+3.86240572873861*fr[15]-3.86240572873861*fl[15]-4.807060063166761*(fr[5]+fl[5])+3.396416424888148*fr[1]-3.396416424888148*fl[1])+nuVtSq[1]*((-4.807060063166761*(fr[13]+fl[13]))+3.396416424888148*fr[7]-3.396416424888148*fl[7]); 
  Gdiff[18] = nuVtSq[0]*(3.866990209613929*fr[18]-3.866990209613929*fl[18]); 

  double Ghat[20]; 
  for(unsigned short int i=0; i<20; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+(0.6614378277661477*dxvl[2]*favg[19]+0.5590169943749475*dxvl[2]*favg[9]+0.4330127018922193*dxvl[2]*favg[3]+0.25*favg[0]*dxvl[2])*nu-1.322875655532295*fjump[19]+drBar[0]*(0.9354143466934851*favg[19]+0.7905694150420947*favg[9]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+0.3535533905932737*drBar[3]*favg[17]+drBar[1]*(0.7905694150420948*favg[15]+0.6123724356957944*favg[5]+0.3535533905932737*favg[1])+drBar[2]*(0.6123724356957944*favg[13]+0.3535533905932737*favg[7])-1.118033988749895*fjump[9]-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+(0.5590169943749476*dxvl[2]*favg[15]+0.4330127018922193*dxvl[2]*favg[5]+0.25*favg[1]*dxvl[2])*nu+drBar[1]*(0.9354143466934851*favg[19]+0.5477225575051661*favg[13]+0.7905694150420947*favg[9]+0.3162277660168379*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+drBar[2]*(0.3105295017040592*favg[17]+0.7071067811865475*favg[15]+0.5477225575051661*favg[5]+0.3162277660168379*favg[1])-1.118033988749895*fjump[15]+drBar[0]*(0.7905694150420948*favg[15]+0.6123724356957944*favg[5]+0.3535533905932737*favg[1])+drBar[3]*(0.5378528742004769*favg[13]+0.3105295017040592*favg[7])-0.8660254037844386*fjump[5]-0.5*fjump[1]; 
  Ghat[2] = Gdiff[2]*rdv+(0.5590169943749476*dxvl[2]*favg[16]+0.4330127018922193*dxvl[2]*favg[6]+0.25*dxvl[2]*favg[2])*nu-1.118033988749895*fjump[16]+drBar[0]*(0.7905694150420948*favg[16]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])+0.3535533905932737*drBar[2]*favg[11]+drBar[1]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[4])-0.8660254037844386*fjump[6]-0.5*fjump[2]; 
  Ghat[4] = Gdiff[4]*rdv+(0.4330127018922193*dxvl[2]*favg[10]+0.25*dxvl[2]*favg[4])*nu+drBar[1]*(0.7905694150420948*favg[16]+0.3162277660168379*favg[11]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])+0.3105295017040593*drBar[3]*favg[11]-0.8660254037844386*fjump[10]+drBar[0]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[4])+drBar[2]*(0.5477225575051661*favg[10]+0.3162277660168379*favg[4])-0.5*fjump[4]; 
  Ghat[7] = Gdiff[7]*rdv+(0.4330127018922194*dxvl[2]*favg[13]+0.25*dxvl[2]*favg[7])*nu+drBar[2]*(0.9354143466934851*favg[19]+0.3912303982179757*favg[13]+0.7905694150420947*favg[9]+0.2258769757263128*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+drBar[1]*(0.3105295017040592*favg[17]+0.7071067811865475*favg[15]+0.5477225575051661*favg[5]+0.3162277660168379*favg[1])+drBar[3]*(0.210818510677892*favg[17]+0.6943650748294133*favg[15]+0.537852874200477*favg[5]+0.3105295017040592*favg[1])-0.8660254037844387*fjump[13]+drBar[0]*(0.6123724356957944*favg[13]+0.3535533905932737*favg[7])-0.5*fjump[7]; 
  Ghat[8] = Gdiff[8]*rdv+(0.4330127018922194*dxvl[2]*favg[14]+0.25*dxvl[2]*favg[8])*nu-0.8660254037844387*fjump[14]+drBar[0]*(0.6123724356957944*favg[14]+0.3535533905932737*favg[8])+0.3535533905932737*drBar[1]*favg[12]-0.5*fjump[8]; 
  Ghat[11] = Gdiff[11]*rdv+0.25*dxvl[2]*favg[11]*nu+drBar[2]*(0.7905694150420947*favg[16]+0.2258769757263128*favg[11]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-0.5*fjump[11]+0.3535533905932737*drBar[0]*favg[11]+drBar[1]*(0.5477225575051661*favg[10]+0.3162277660168379*favg[4])+drBar[3]*(0.5378528742004769*favg[10]+0.3105295017040593*favg[4]); 
  Ghat[12] = Gdiff[12]*rdv+0.25*dxvl[2]*favg[12]*nu+drBar[1]*(0.6123724356957944*favg[14]+0.3535533905932737*favg[8])-0.5*fjump[12]+0.3162277660168379*drBar[2]*favg[12]+0.3535533905932737*drBar[0]*favg[12]; 
  Ghat[17] = Gdiff[17]*rdv+0.25*dxvl[2]*favg[17]*nu+drBar[3]*(0.9354143466934851*favg[19]+0.3651483716701107*favg[13]+0.7905694150420947*favg[9]+0.210818510677892*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.5*fjump[17]+0.3535533905932737*drBar[0]*favg[17]+drBar[2]*(0.210818510677892*favg[17]+0.6943650748294133*favg[15]+0.537852874200477*favg[5]+0.3105295017040592*favg[1])+drBar[1]*(0.5378528742004769*favg[13]+0.3105295017040592*favg[7]); 
  Ghat[18] = Gdiff[18]*rdv+0.25*dxvl[2]*favg[18]*nu-0.5*fjump[18]+0.3535533905932737*drBar[0]*favg[18]; 

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
  incr1[17] = -0.5*Ghat[17]; 
  incr1[18] = -0.5*Ghat[18]; 
  incr1[19] = 1.322875655532295*Ghat[0]; 

  double incr2[20]; 
  incr2[3] = nuVtSq[0]*((-0.2314550249431378*fr[19])+0.2314550249431378*fl[19]+0.4279082480509108*(fr[9]+fl[9])-0.4640388251536715*fr[3]+0.4640388251536715*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*nuVtSq[3]*(fr[17]+fl[17])+nuVtSq[1]*(0.4279082480509107*(fr[15]+fl[15])-0.4640388251536715*fr[5]+0.4640388251536715*fl[5]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSq[2]*((-0.4640388251536715*fr[13])+0.4640388251536715*fl[13]+0.3061862178478971*(fr[7]+fl[7])); 
  incr2[5] = nuVtSq[1]*((-0.2314550249431378*fr[19])+0.2314550249431378*fl[19]-0.4150489428970995*fr[13]+0.4150489428970995*fl[13]+0.4279082480509108*(fr[9]+fl[9])+0.273861278752583*(fr[7]+fl[7])-0.4640388251536715*fr[3]+0.4640388251536715*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSq[2]*(0.2689264371002384*(fr[17]+fl[17])+0.3827327723098714*(fr[15]+fl[15])-0.4150489428970995*fr[5]+0.4150489428970995*fl[5]+0.273861278752583*(fr[1]+fl[1]))+nuVtSq[0]*(0.4279082480509107*(fr[15]+fl[15])-0.4640388251536715*fr[5]+0.4640388251536715*fl[5]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSq[3]*((-0.4075699709865777*fr[13])+0.4075699709865777*fl[13]+0.2689264371002384*(fr[7]+fl[7])); 
  incr2[6] = nuVtSq[0]*(0.4279082480509107*(fr[16]+fl[16])-0.4640388251536715*fr[6]+0.4640388251536715*fl[6]+0.3061862178478971*(fr[2]+fl[2]))+0.3061862178478971*nuVtSq[2]*(fr[11]+fl[11])+nuVtSq[1]*((-0.4640388251536715*fr[10])+0.4640388251536715*fl[10]+0.3061862178478971*(fr[4]+fl[4])); 
  incr2[9] = nuVtSq[0]*(0.896421457000795*fr[19]-0.896421457000795*fl[19]-1.657281518405969*(fr[9]+fl[9])+1.797214641813825*fr[3]-1.797214641813825*fl[3]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*nuVtSq[3]*(fr[17]+fl[17])+nuVtSq[1]*((-1.65728151840597*(fr[15]+fl[15]))+1.797214641813825*fr[5]-1.797214641813825*fl[5]-1.185854122563142*(fr[1]+fl[1]))+nuVtSq[2]*(1.797214641813825*fr[13]-1.797214641813825*fl[13]-1.185854122563142*(fr[7]+fl[7])); 
  incr2[10] = nuVtSq[1]*(0.4279082480509107*(fr[16]+fl[16])+0.273861278752583*(fr[11]+fl[11])-0.4640388251536715*fr[6]+0.4640388251536715*fl[6]+0.3061862178478971*(fr[2]+fl[2]))+0.2689264371002384*nuVtSq[3]*(fr[11]+fl[11])+nuVtSq[2]*((-0.4150489428970995*fr[10])+0.4150489428970995*fl[10]+0.273861278752583*(fr[4]+fl[4]))+nuVtSq[0]*((-0.4640388251536715*fr[10])+0.4640388251536715*fl[10]+0.3061862178478971*(fr[4]+fl[4])); 
  incr2[13] = nuVtSq[2]*((-0.2314550249431378*fr[19])+0.2314550249431378*fl[19]-0.2964635306407854*fr[13]+0.2964635306407854*fl[13]+0.4279082480509108*(fr[9]+fl[9])+0.1956151991089878*(fr[7]+fl[7])-0.4640388251536715*fr[3]+0.4640388251536715*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSq[1]*(0.2689264371002384*(fr[17]+fl[17])+0.3827327723098713*(fr[15]+fl[15])-0.4150489428970995*fr[5]+0.4150489428970995*fl[5]+0.273861278752583*(fr[1]+fl[1]))+nuVtSq[3]*(0.1825741858350553*(fr[17]+fl[17])+0.3758361214393465*(fr[15]+fl[15])-0.4075699709865777*fr[5]+0.4075699709865777*fl[5]+0.2689264371002384*(fr[1]+fl[1]))+nuVtSq[0]*((-0.4640388251536715*fr[13])+0.4640388251536715*fl[13]+0.3061862178478971*(fr[7]+fl[7])); 
  incr2[14] = nuVtSq[0]*((-0.4640388251536715*fr[14])+0.4640388251536715*fl[14]+0.3061862178478971*(fr[8]+fl[8]))+0.3061862178478971*nuVtSq[1]*(fr[12]+fl[12]); 
  incr2[15] = nuVtSq[1]*(0.8964214570007949*fr[19]-0.8964214570007949*fl[19]+1.60747764370146*fr[13]-1.60747764370146*fl[13]-1.65728151840597*(fr[9]+fl[9])-1.060660171779821*(fr[7]+fl[7])+1.797214641813825*fr[3]-1.797214641813825*fl[3]-1.185854122563142*(fr[0]+fl[0]))+nuVtSq[2]*((-1.04154761224412*(fr[17]+fl[17]))-1.482317653203927*(fr[15]+fl[15])+1.60747764370146*fr[5]-1.60747764370146*fl[5]-1.060660171779821*(fr[1]+fl[1]))+nuVtSq[0]*((-1.657281518405969*(fr[15]+fl[15]))+1.797214641813825*fr[5]-1.797214641813825*fl[5]-1.185854122563142*(fr[1]+fl[1]))+nuVtSq[3]*(1.578511710045255*fr[13]-1.578511710045255*fl[13]-1.04154761224412*(fr[7]+fl[7])); 
  incr2[16] = nuVtSq[0]*((-1.657281518405969*(fr[16]+fl[16]))+1.797214641813825*fr[6]-1.797214641813825*fl[6]-1.185854122563142*(fr[2]+fl[2]))-1.185854122563142*nuVtSq[2]*(fr[11]+fl[11])+nuVtSq[1]*(1.797214641813825*fr[10]-1.797214641813825*fl[10]-1.185854122563142*(fr[4]+fl[4])); 
  incr2[19] = nuVtSq[0]*((-2.121320343559642*fr[19])+2.121320343559642*fl[19]+3.921843874378478*(fr[9]+fl[9])-4.252986083330156*fr[3]+4.252986083330156*fl[3]+2.806243040080455*(fr[0]+fl[0]))+2.806243040080455*nuVtSq[3]*(fr[17]+fl[17])+nuVtSq[1]*(3.921843874378477*(fr[15]+fl[15])-4.252986083330156*fr[5]+4.252986083330156*fl[5]+2.806243040080455*(fr[1]+fl[1]))+nuVtSq[2]*((-4.252986083330157*fr[13])+4.252986083330157*fl[13]+2.806243040080455*(fr[7]+fl[7])); 

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
  outr[17] += incr1[17]*rdv2r; 
  outr[18] += incr1[18]*rdv2r; 
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
  outl[17] += -1.0*incr1[17]*rdv2l; 
  outl[18] += -1.0*incr1[18]*rdv2l; 
  outl[19] += incr1[19]*rdv2l-1.0*incr2[19]*rdvSq4l; 

  const double vMuMid = (0.7905694150420947*nuUy[2])/nu-(0.7071067811865475*nuUy[0])/nu+wl[2]; 
  return std::abs(vMuMid); 
} 
