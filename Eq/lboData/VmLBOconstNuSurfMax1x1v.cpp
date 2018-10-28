#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x1vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates. 
  // dxv[2]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[2]:       bulk velocity (in 1 directions). 
  // vtSq[2]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nuL = 2.0*nu/dxvl[1]; 
  double rdv2nuR = 2.0*nu/dxvr[1]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  const double *Ux = &u[0]; 

  double favg[3]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 

  double fjump[3]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(-1*fr[2])); 

  double drBar[2]; 
  drBar[0] = 1.414213562373095*wl[1]-1.0*Ux[0]; 
  drBar[1] = -1.0*Ux[1]; 

  double Gdiff[3]; 
  Gdiff[0] = vtSq[0]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0])+(1.590990257669731*fr[1]-1.590990257669731*fl[1])*vtSq[1]; 
  Gdiff[1] = vtSq[1]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0])+vtSq[0]*(1.590990257669731*fr[1]-1.590990257669731*fl[1]); 

  double Ghat[3]; 
  for(unsigned short int i=0; i<3; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv-0.8660254037844386*fjump[2]+0.4330127018922193*dxvl[1]*favg[2]+drBar[0]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*drBar[1]*favg[1]+0.25*favg[0]*dxvl[1]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+drBar[1]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[1]+0.25*dxvl[1]*favg[1]+0.3535533905932737*drBar[0]*favg[1]; 

  double incr1[3]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 

  double incr2[3]; 
  incr2[2] = vtSq[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*(fr[1]+fl[1])*vtSq[1]; 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr2[2]*rdvSq4nuR+incr1[2]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += incr1[2]*rdv2nuL-1.0*incr2[2]*rdvSq4nuL; 

  const double vMuMid = wl[1]-0.7071067811865475*Ux[0]; 
  return std::abs(vMuMid); 
} 
double VmLBOconstNuSurf1x1vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates. 
  // dxv[2]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[3]:       bulk velocity (in 1 directions). 
  // vtSq[3]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nuL = 2.0*nu/dxvl[1]; 
  double rdv2nuR = 2.0*nu/dxvr[1]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  const double *Ux = &u[0]; 

  double favg[6]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 

  double fjump[6]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(1*fr[5])); 

  double drBar[3]; 
  drBar[0] = 1.414213562373095*wl[1]-1.0*Ux[0]; 
  drBar[1] = -1.0*Ux[1]; 
  drBar[2] = -1.0*Ux[2]; 

  double Gdiff[6]; 
  Gdiff[0] = vtSq[0]*(1.897366596101028*fr[5]-1.897366596101028*fl[5]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+vtSq[2]*(2.651650429449552*fr[4]-2.651650429449552*fl[4])+vtSq[1]*((-3.368048396326869*(fr[3]+fl[3]))+2.651650429449552*fr[1]-2.651650429449552*fl[1]); 
  Gdiff[1] = vtSq[1]*(1.897366596101028*fr[5]-1.897366596101028*fl[5]+2.371708245126284*fr[4]-2.371708245126284*fl[4]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+vtSq[2]*((-3.012474066278413*(fr[3]+fl[3]))+2.371708245126284*fr[1]-2.371708245126284*fl[1])+vtSq[0]*((-3.368048396326869*(fr[3]+fl[3]))+2.651650429449552*fr[1]-2.651650429449552*fl[1]); 
  Gdiff[4] = vtSq[2]*(1.897366596101028*fr[5]-1.897366596101028*fl[5]+1.694077317947346*fr[4]-1.694077317947346*fl[4]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+vtSq[0]*(2.651650429449552*fr[4]-2.651650429449552*fl[4])+vtSq[1]*((-3.012474066278413*(fr[3]+fl[3]))+2.371708245126284*fr[1]-2.371708245126284*fl[1]); 

  double Ghat[6]; 
  for(unsigned short int i=0; i<6; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv-1.118033988749895*fjump[5]+0.5590169943749475*dxvl[1]*favg[5]+drBar[0]*(0.7905694150420947*favg[5]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*drBar[2]*favg[4]+drBar[1]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]+0.4330127018922193*dxvl[1]*favg[2]+0.25*favg[0]*dxvl[1]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+drBar[1]*(0.7905694150420947*favg[5]+0.3162277660168379*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[3]+0.4330127018922193*dxvl[1]*favg[3]+drBar[0]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[2]*(0.5477225575051661*favg[3]+0.3162277660168379*favg[1])-0.5*fjump[1]+0.25*dxvl[1]*favg[1]; 
  Ghat[4] = Gdiff[4]*rdv+drBar[2]*(0.7905694150420947*favg[5]+0.2258769757263128*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[4]+0.25*dxvl[1]*favg[4]+0.3535533905932737*drBar[0]*favg[4]+drBar[1]*(0.5477225575051661*favg[3]+0.3162277660168379*favg[1]); 

  double incr1[6]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -1.118033988749895*Ghat[0]; 

  double incr2[6]; 
  incr2[2] = vtSq[0]*(0.2995357736356374*(fr[5]+fl[5])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*vtSq[2]*(fr[4]+fl[4])+vtSq[1]*((-0.430893194785552*fr[3])+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[3] = vtSq[1]*(0.2995357736356374*(fr[5]+fl[5])+0.273861278752583*(fr[4]+fl[4])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[2]*((-0.3854025898330209*fr[3])+0.3854025898330209*fl[3]+0.273861278752583*(fr[1]+fl[1]))+vtSq[0]*((-0.430893194785552*fr[3])+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[5] = vtSq[0]*((-1.160097062884178*(fr[5]+fl[5]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*vtSq[2]*(fr[4]+fl[4])+vtSq[1]*(1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[1]+fl[1])); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr2[2]*rdvSq4nuR+incr1[2]*rdv2nuR; 
  outr[3] += incr2[3]*rdvSq4nuR+incr1[3]*rdv2nuR; 
  outr[4] += incr1[4]*rdv2nuR; 
  outr[5] += incr2[5]*rdvSq4nuR+incr1[5]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += incr1[2]*rdv2nuL-1.0*incr2[2]*rdvSq4nuL; 
  outl[3] += incr1[3]*rdv2nuL-1.0*incr2[3]*rdvSq4nuL; 
  outl[4] += -1.0*incr1[4]*rdv2nuL; 
  outl[5] += incr2[5]*rdvSq4nuL-1.0*incr1[5]*rdv2nuL; 

  const double vMuMid = 0.7905694150420947*Ux[2]+wl[1]-0.7071067811865475*Ux[0]; 
  return std::abs(vMuMid); 
} 
double VmLBOconstNuSurf1x1vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates. 
  // dxv[2]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[4]:       bulk velocity (in 1 directions). 
  // vtSq[4]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nuL = 2.0*nu/dxvl[1]; 
  double rdv2nuR = 2.0*nu/dxvr[1]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  const double *Ux = &u[0]; 

  double favg[10]; 
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

  double fjump[10]; 
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

  double drBar[4]; 
  drBar[0] = 1.414213562373095*wl[1]-1.0*Ux[0]; 
  drBar[1] = -1.0*Ux[1]; 
  drBar[2] = -1.0*Ux[2]; 
  drBar[3] = -1.0*Ux[3]; 

  double Gdiff[10]; 
  Gdiff[0] = vtSq[0]*((-1.929292090055312*(fr[9]+fl[9]))+4.397542371171649*fr[5]-4.397542371171649*fl[5]-5.473078644031159*(fr[2]+fl[2])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+vtSq[3]*(3.866990209613929*fr[8]-3.866990209613929*fl[8])+vtSq[1]*(4.39754237117165*fr[7]-4.39754237117165*fl[7]-5.473078644031159*(fr[3]+fl[3])+3.866990209613929*fr[1]-3.866990209613929*fl[1])+vtSq[2]*((-5.473078644031159*(fr[6]+fl[6]))+3.866990209613929*fr[4]-3.866990209613929*fl[4]); 
  Gdiff[1] = vtSq[1]*((-1.929292090055312*(fr[9]+fl[9]))-4.89527035770242*(fr[6]+fl[6])+4.397542371171649*fr[5]-4.397542371171649*fl[5]+3.458741190809164*fr[4]-3.458741190809164*fl[4]-5.473078644031159*(fr[2]+fl[2])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+vtSq[2]*(3.396416424888148*fr[8]-3.396416424888148*fl[8]+3.933281470350169*fr[7]-3.933281470350169*fl[7]-4.89527035770242*(fr[3]+fl[3])+3.458741190809164*fr[1]-3.458741190809164*fl[1])+vtSq[0]*(4.39754237117165*fr[7]-4.39754237117165*fl[7]-5.473078644031159*(fr[3]+fl[3])+3.866990209613929*fr[1]-3.866990209613929*fl[1])+vtSq[3]*((-4.807060063166761*(fr[6]+fl[6]))+3.396416424888148*fr[4]-3.396416424888148*fl[4]); 
  Gdiff[4] = vtSq[2]*((-1.929292090055312*(fr[9]+fl[9]))-3.496621684073157*(fr[6]+fl[6])+4.397542371171649*fr[5]-4.397542371171649*fl[5]+2.470529422006546*fr[4]-2.470529422006546*fl[4]-5.473078644031159*(fr[2]+fl[2])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+vtSq[1]*(3.396416424888148*fr[8]-3.396416424888148*fl[8]+3.933281470350169*fr[7]-3.933281470350169*fl[7]-4.89527035770242*(fr[3]+fl[3])+3.458741190809164*fr[1]-3.458741190809164*fl[1])+vtSq[3]*(2.305827460539443*fr[8]-2.305827460539443*fl[8]+3.86240572873861*fr[7]-3.86240572873861*fl[7]-4.807060063166761*(fr[3]+fl[3])+3.396416424888148*fr[1]-3.396416424888148*fl[1])+vtSq[0]*((-5.473078644031159*(fr[6]+fl[6]))+3.866990209613929*fr[4]-3.866990209613929*fl[4]); 
  Gdiff[8] = vtSq[3]*((-1.929292090055312*(fr[9]+fl[9]))-3.263513571801613*(fr[6]+fl[6])+4.397542371171649*fr[5]-4.397542371171649*fl[5]+2.305827460539443*fr[4]-2.305827460539443*fl[4]-5.473078644031159*(fr[2]+fl[2])+3.866990209613929*fr[0]-3.866990209613929*fl[0])+vtSq[0]*(3.866990209613929*fr[8]-3.866990209613929*fl[8])+vtSq[2]*(2.305827460539443*fr[8]-2.305827460539443*fl[8]+3.86240572873861*fr[7]-3.86240572873861*fl[7]-4.807060063166761*(fr[3]+fl[3])+3.396416424888148*fr[1]-3.396416424888148*fl[1])+vtSq[1]*((-4.807060063166761*(fr[6]+fl[6]))+3.396416424888148*fr[4]-3.396416424888148*fl[4]); 

  double Ghat[10]; 
  for(unsigned short int i=0; i<10; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv-1.322875655532295*fjump[9]+0.6614378277661477*dxvl[1]*favg[9]+drBar[0]*(0.9354143466934851*favg[9]+0.7905694150420947*favg[5]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*drBar[3]*favg[8]+drBar[1]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[2]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])-1.118033988749895*fjump[5]+0.5590169943749475*dxvl[1]*favg[5]-0.8660254037844386*fjump[2]+0.4330127018922193*dxvl[1]*favg[2]+0.25*favg[0]*dxvl[1]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv+drBar[1]*(0.9354143466934851*favg[9]+0.5477225575051661*favg[6]+0.7905694150420947*favg[5]+0.3162277660168379*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+drBar[2]*(0.3105295017040592*favg[8]+0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])-1.118033988749895*fjump[7]+0.5590169943749476*dxvl[1]*favg[7]+drBar[0]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[3]*(0.5378528742004769*favg[6]+0.3105295017040592*favg[4])-0.8660254037844386*fjump[3]+0.4330127018922193*dxvl[1]*favg[3]-0.5*fjump[1]+0.25*dxvl[1]*favg[1]; 
  Ghat[4] = Gdiff[4]*rdv+drBar[2]*(0.9354143466934851*favg[9]+0.3912303982179757*favg[6]+0.7905694150420947*favg[5]+0.2258769757263128*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+drBar[1]*(0.3105295017040592*favg[8]+0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])+drBar[3]*(0.210818510677892*favg[8]+0.6943650748294133*favg[7]+0.537852874200477*favg[3]+0.3105295017040592*favg[1])-0.8660254037844387*fjump[6]+0.4330127018922194*dxvl[1]*favg[6]+drBar[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])-0.5*fjump[4]+0.25*dxvl[1]*favg[4]; 
  Ghat[8] = Gdiff[8]*rdv+drBar[3]*(0.9354143466934851*favg[9]+0.3651483716701107*favg[6]+0.7905694150420947*favg[5]+0.210818510677892*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[8]+0.25*dxvl[1]*favg[8]+0.3535533905932737*drBar[0]*favg[8]+drBar[2]*(0.210818510677892*favg[8]+0.6943650748294133*favg[7]+0.537852874200477*favg[3]+0.3105295017040592*favg[1])+drBar[1]*(0.5378528742004769*favg[6]+0.3105295017040592*favg[4]); 

  double incr1[10]; 
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

  double incr2[10]; 
  incr2[2] = vtSq[0]*((-0.2314550249431378*fr[9])+0.2314550249431378*fl[9]+0.4279082480509108*(fr[5]+fl[5])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*vtSq[3]*(fr[8]+fl[8])+vtSq[1]*(0.4279082480509107*(fr[7]+fl[7])-0.4640388251536715*fr[3]+0.4640388251536715*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[2]*((-0.4640388251536715*fr[6])+0.4640388251536715*fl[6]+0.3061862178478971*(fr[4]+fl[4])); 
  incr2[3] = vtSq[1]*((-0.2314550249431378*fr[9])+0.2314550249431378*fl[9]-0.4150489428970995*fr[6]+0.4150489428970995*fl[6]+0.4279082480509108*(fr[5]+fl[5])+0.273861278752583*(fr[4]+fl[4])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[2]*(0.2689264371002384*(fr[8]+fl[8])+0.3827327723098714*(fr[7]+fl[7])-0.4150489428970995*fr[3]+0.4150489428970995*fl[3]+0.273861278752583*(fr[1]+fl[1]))+vtSq[0]*(0.4279082480509107*(fr[7]+fl[7])-0.4640388251536715*fr[3]+0.4640388251536715*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[3]*((-0.4075699709865777*fr[6])+0.4075699709865777*fl[6]+0.2689264371002384*(fr[4]+fl[4])); 
  incr2[5] = vtSq[0]*(0.896421457000795*fr[9]-0.896421457000795*fl[9]-1.657281518405969*(fr[5]+fl[5])+1.797214641813825*fr[2]-1.797214641813825*fl[2]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*vtSq[3]*(fr[8]+fl[8])+vtSq[1]*((-1.65728151840597*(fr[7]+fl[7]))+1.797214641813825*fr[3]-1.797214641813825*fl[3]-1.185854122563142*(fr[1]+fl[1]))+vtSq[2]*(1.797214641813825*fr[6]-1.797214641813825*fl[6]-1.185854122563142*(fr[4]+fl[4])); 
  incr2[6] = vtSq[2]*((-0.2314550249431378*fr[9])+0.2314550249431378*fl[9]-0.2964635306407854*fr[6]+0.2964635306407854*fl[6]+0.4279082480509108*(fr[5]+fl[5])+0.1956151991089878*(fr[4]+fl[4])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[1]*(0.2689264371002384*(fr[8]+fl[8])+0.3827327723098713*(fr[7]+fl[7])-0.4150489428970995*fr[3]+0.4150489428970995*fl[3]+0.273861278752583*(fr[1]+fl[1]))+vtSq[3]*(0.1825741858350553*(fr[8]+fl[8])+0.3758361214393465*(fr[7]+fl[7])-0.4075699709865777*fr[3]+0.4075699709865777*fl[3]+0.2689264371002384*(fr[1]+fl[1]))+vtSq[0]*((-0.4640388251536715*fr[6])+0.4640388251536715*fl[6]+0.3061862178478971*(fr[4]+fl[4])); 
  incr2[7] = vtSq[1]*(0.8964214570007949*fr[9]-0.8964214570007949*fl[9]+1.60747764370146*fr[6]-1.60747764370146*fl[6]-1.65728151840597*(fr[5]+fl[5])-1.060660171779821*(fr[4]+fl[4])+1.797214641813825*fr[2]-1.797214641813825*fl[2]-1.185854122563142*(fr[0]+fl[0]))+vtSq[2]*((-1.04154761224412*(fr[8]+fl[8]))-1.482317653203927*(fr[7]+fl[7])+1.60747764370146*fr[3]-1.60747764370146*fl[3]-1.060660171779821*(fr[1]+fl[1]))+vtSq[0]*((-1.657281518405969*(fr[7]+fl[7]))+1.797214641813825*fr[3]-1.797214641813825*fl[3]-1.185854122563142*(fr[1]+fl[1]))+vtSq[3]*(1.578511710045255*fr[6]-1.578511710045255*fl[6]-1.04154761224412*(fr[4]+fl[4])); 
  incr2[9] = vtSq[0]*((-2.121320343559642*fr[9])+2.121320343559642*fl[9]+3.921843874378478*(fr[5]+fl[5])-4.252986083330156*fr[2]+4.252986083330156*fl[2]+2.806243040080455*(fr[0]+fl[0]))+2.806243040080455*vtSq[3]*(fr[8]+fl[8])+vtSq[1]*(3.921843874378477*(fr[7]+fl[7])-4.252986083330156*fr[3]+4.252986083330156*fl[3]+2.806243040080455*(fr[1]+fl[1]))+vtSq[2]*((-4.252986083330157*fr[6])+4.252986083330157*fl[6]+2.806243040080455*(fr[4]+fl[4])); 

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

  const double vMuMid = 0.7905694150420947*Ux[2]+wl[1]-0.7071067811865475*Ux[0]; 
  return std::abs(vMuMid); 
} 
