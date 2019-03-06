#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:       Cell-center coordinates. 
  // dxv[3]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[4]:       bulk velocity (in 2 directions). 
  // vtSq[2]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nuL = 2.0*nu/dxvl[1]; 
  double rdv2nuR = 2.0*nu/dxvr[1]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  const double *Ux = &u[0]; 

  double favg[4]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 

  double fjump[4]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(1*fr[3])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[1]+0.7071067811865475*dxvl[1]-1.0*Ux[0]; 
  alphaDrag[1] = -1.0*Ux[1]; 

  double Ghat[4]; 
  for(unsigned short int i=0; i<4; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.530931089239486*vtSq[0]*fr[2])-1.530931089239486*vtSq[0]*fl[2]+1.590990257669731*fr[1]*vtSq[1]-1.590990257669731*fl[1]*vtSq[1]+1.590990257669731*fr[0]*vtSq[0]-1.590990257669731*fl[0]*vtSq[0])*rdv-0.8660254037844386*fjump[2]+alphaDrag[0]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = ((-1.530931089239486*vtSq[1]*fr[2])-1.530931089239486*vtSq[1]*fl[2]+1.590990257669731*fr[0]*vtSq[1]-1.590990257669731*fl[0]*vtSq[1]+1.590990257669731*vtSq[0]*fr[1]-1.590990257669731*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[1]+0.3535533905932737*alphaDrag[0]*favg[1]; 
  Ghat[3] = (1.590990257669731*vtSq[0]*fr[3]-1.590990257669731*vtSq[0]*fl[3])*rdv-0.5*fjump[3]+0.3535533905932737*alphaDrag[0]*favg[3]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 

  double incr2[4]; 
  incr2[2] = vtSq[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*(fr[1]+fl[1])*vtSq[1]; 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr2[2]*rdvSq4nuR+incr1[2]*rdv2nuR; 
  outr[3] += incr1[3]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += incr1[2]*rdv2nuL-1.0*incr2[2]*rdvSq4nuL; 
  outl[3] += -1.0*incr1[3]*rdv2nuL; 

  return std::abs(wl[1]-0.7071067811865475*Ux[0]); 
} 
double VmLBOconstNuSurf1x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:       Cell-center coordinates. 
  // dxv[3]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[6]:       bulk velocity (in 2 directions). 
  // vtSq[3]:    thermal speed squared. 
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
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 

  double fjump[10]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = vMuMidMax*(fl[9]-(1*fr[9])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[1]+0.7071067811865475*dxvl[1]-1.0*Ux[0]; 
  alphaDrag[1] = -1.0*Ux[1]; 
  alphaDrag[2] = -1.0*Ux[2]; 

  double Ghat[10]; 
  for(unsigned short int i=0; i<10; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*vtSq[0]*fr[8]-1.897366596101028*vtSq[0]*fl[8]+2.651650429449552*vtSq[2]*fr[7]-2.651650429449552*vtSq[2]*fl[7]-3.368048396326869*vtSq[1]*fr[4]-3.368048396326869*vtSq[1]*fl[4]-3.368048396326869*vtSq[0]*fr[2]-3.368048396326869*vtSq[0]*fl[2]+2.651650429449552*fr[1]*vtSq[1]-2.651650429449552*fl[1]*vtSq[1]+2.651650429449552*fr[0]*vtSq[0]-2.651650429449552*fl[0]*vtSq[0])*rdv-1.118033988749895*fjump[8]+alphaDrag[0]*(0.7905694150420947*favg[8]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[2]*favg[7]+alphaDrag[1]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (1.897366596101028*vtSq[1]*fr[8]-1.897366596101028*vtSq[1]*fl[8]+2.371708245126284*vtSq[1]*fr[7]-2.371708245126284*vtSq[1]*fl[7]-3.012474066278413*vtSq[2]*fr[4]-3.368048396326869*vtSq[0]*fr[4]-3.012474066278413*vtSq[2]*fl[4]-3.368048396326869*vtSq[0]*fl[4]+2.371708245126284*fr[1]*vtSq[2]-2.371708245126284*fl[1]*vtSq[2]-3.368048396326869*vtSq[1]*fr[2]-3.368048396326869*vtSq[1]*fl[2]+2.651650429449552*fr[0]*vtSq[1]-2.651650429449552*fl[0]*vtSq[1]+2.651650429449552*vtSq[0]*fr[1]-2.651650429449552*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.7905694150420947*favg[8]+0.3162277660168379*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[4]+alphaDrag[0]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.5477225575051661*favg[4]+0.3162277660168379*favg[1])-0.5*fjump[1]; 
  Ghat[3] = ((-3.368048396326869*vtSq[0]*fr[6])-3.368048396326869*vtSq[0]*fl[6]+2.651650429449552*vtSq[1]*fr[5]-2.651650429449552*vtSq[1]*fl[5]+2.651650429449552*vtSq[0]*fr[3]-2.651650429449552*vtSq[0]*fl[3])*rdv-0.8660254037844386*fjump[6]+alphaDrag[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[3])+0.3535533905932737*alphaDrag[1]*favg[5]-0.5*fjump[3]; 
  Ghat[5] = ((-3.368048396326869*vtSq[1]*fr[6])-3.368048396326869*vtSq[1]*fl[6]+2.371708245126284*vtSq[2]*fr[5]+2.651650429449552*vtSq[0]*fr[5]-2.371708245126284*vtSq[2]*fl[5]-2.651650429449552*vtSq[0]*fl[5]+2.651650429449552*vtSq[1]*fr[3]-2.651650429449552*vtSq[1]*fl[3])*rdv+alphaDrag[1]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[5]+0.3162277660168379*alphaDrag[2]*favg[5]+0.3535533905932737*alphaDrag[0]*favg[5]; 
  Ghat[7] = (1.897366596101028*vtSq[2]*fr[8]-1.897366596101028*vtSq[2]*fl[8]+1.694077317947346*vtSq[2]*fr[7]+2.651650429449552*vtSq[0]*fr[7]-1.694077317947346*vtSq[2]*fl[7]-2.651650429449552*vtSq[0]*fl[7]-3.012474066278413*vtSq[1]*fr[4]-3.012474066278413*vtSq[1]*fl[4]-3.368048396326869*fr[2]*vtSq[2]-3.368048396326869*fl[2]*vtSq[2]+2.651650429449552*fr[0]*vtSq[2]-2.651650429449552*fl[0]*vtSq[2]+2.371708245126284*fr[1]*vtSq[1]-2.371708245126284*fl[1]*vtSq[1])*rdv+alphaDrag[2]*(0.7905694150420947*favg[8]+0.2258769757263128*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[7]+0.3535533905932737*alphaDrag[0]*favg[7]+alphaDrag[1]*(0.5477225575051661*favg[4]+0.3162277660168379*favg[1]); 
  Ghat[9] = (2.651650429449552*vtSq[0]*fr[9]-2.651650429449552*vtSq[0]*fl[9])*rdv-0.5*fjump[9]+0.3535533905932737*alphaDrag[0]*favg[9]; 

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
  incr2[2] = vtSq[0]*(0.2995357736356374*(fr[8]+fl[8])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*vtSq[2]*(fr[7]+fl[7])+vtSq[1]*((-0.430893194785552*fr[4])+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[4] = vtSq[1]*(0.2995357736356374*(fr[8]+fl[8])+0.273861278752583*(fr[7]+fl[7])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[2]*((-0.3854025898330209*fr[4])+0.3854025898330209*fl[4]+0.273861278752583*(fr[1]+fl[1]))+vtSq[0]*((-0.430893194785552*fr[4])+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[6] = vtSq[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3]))+0.3061862178478971*vtSq[1]*(fr[5]+fl[5]); 
  incr2[8] = vtSq[0]*((-1.160097062884178*(fr[8]+fl[8]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*vtSq[2]*(fr[7]+fl[7])+vtSq[1]*(1.668842167398551*fr[4]-1.668842167398551*fl[4]-1.185854122563142*(fr[1]+fl[1])); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr2[2]*rdvSq4nuR+incr1[2]*rdv2nuR; 
  outr[3] += incr1[3]*rdv2nuR; 
  outr[4] += incr2[4]*rdvSq4nuR+incr1[4]*rdv2nuR; 
  outr[5] += incr1[5]*rdv2nuR; 
  outr[6] += incr2[6]*rdvSq4nuR+incr1[6]*rdv2nuR; 
  outr[7] += incr1[7]*rdv2nuR; 
  outr[8] += incr2[8]*rdvSq4nuR+incr1[8]*rdv2nuR; 
  outr[9] += incr1[9]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += incr1[2]*rdv2nuL-1.0*incr2[2]*rdvSq4nuL; 
  outl[3] += -1.0*incr1[3]*rdv2nuL; 
  outl[4] += incr1[4]*rdv2nuL-1.0*incr2[4]*rdvSq4nuL; 
  outl[5] += -1.0*incr1[5]*rdv2nuL; 
  outl[6] += incr1[6]*rdv2nuL-1.0*incr2[6]*rdvSq4nuL; 
  outl[7] += -1.0*incr1[7]*rdv2nuL; 
  outl[8] += incr2[8]*rdvSq4nuL-1.0*incr1[8]*rdv2nuL; 
  outl[9] += -1.0*incr1[9]*rdv2nuL; 

  return std::abs(0.7905694150420947*Ux[2]+wl[1]-0.7071067811865475*Ux[0]); 
} 
double VmLBOconstNuSurf1x2vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:       Cell-center coordinates. 
  // dxv[3]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[8]:       bulk velocity (in 2 directions). 
  // vtSq[4]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nuL = 2.0*nu/dxvl[1]; 
  double rdv2nuR = 2.0*nu/dxvr[1]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  const double *Ux = &u[0]; 

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
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = vMuMidMax*(fl[14]-(1*fr[14])); 
  fjump[15] = vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = vMuMidMax*(fl[16]-(-1*fr[16])); 
  fjump[17] = vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = vMuMidMax*(fl[18]-(-1*fr[18])); 
  fjump[19] = vMuMidMax*(fl[19]-(1*fr[19])); 

  double alphaDrag[4]; 
  alphaDrag[0] = 1.414213562373095*wl[1]+0.7071067811865475*dxvl[1]-1.0*Ux[0]; 
  alphaDrag[1] = -1.0*Ux[1]; 
  alphaDrag[2] = -1.0*Ux[2]; 
  alphaDrag[3] = -1.0*Ux[3]; 

  double Ghat[20]; 
  for(unsigned short int i=0; i<20; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.929292090055312*vtSq[0]*fr[18])-1.929292090055312*vtSq[0]*fl[18]+3.866990209613929*vtSq[3]*fr[17]-3.866990209613929*vtSq[3]*fl[17]+4.39754237117165*vtSq[1]*fr[12]-4.39754237117165*vtSq[1]*fl[12]-5.473078644031159*vtSq[2]*fr[11]-5.473078644031159*vtSq[2]*fl[11]+4.397542371171649*vtSq[0]*fr[8]-4.397542371171649*vtSq[0]*fl[8]+3.866990209613929*vtSq[2]*fr[7]-3.866990209613929*vtSq[2]*fl[7]-5.473078644031159*vtSq[1]*fr[4]-5.473078644031159*vtSq[1]*fl[4]-5.473078644031159*vtSq[0]*fr[2]-5.473078644031159*vtSq[0]*fl[2]+3.866990209613929*fr[1]*vtSq[1]-3.866990209613929*fl[1]*vtSq[1]+3.866990209613929*fr[0]*vtSq[0]-3.866990209613929*fl[0]*vtSq[0])*rdv-1.322875655532295*fjump[18]+alphaDrag[0]*(0.9354143466934851*favg[18]+0.7905694150420947*favg[8]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[3]*favg[17]+alphaDrag[1]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.6123724356957944*favg[11]+0.3535533905932737*favg[7])-1.118033988749895*fjump[8]-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = ((-1.929292090055312*vtSq[1]*fr[18])-1.929292090055312*vtSq[1]*fl[18]+3.396416424888148*vtSq[2]*fr[17]-3.396416424888148*vtSq[2]*fl[17]+3.933281470350169*vtSq[2]*fr[12]+4.39754237117165*vtSq[0]*fr[12]-3.933281470350169*vtSq[2]*fl[12]-4.39754237117165*vtSq[0]*fl[12]-4.807060063166761*vtSq[3]*fr[11]-4.89527035770242*vtSq[1]*fr[11]-4.807060063166761*vtSq[3]*fl[11]-4.89527035770242*vtSq[1]*fl[11]+4.397542371171649*vtSq[1]*fr[8]-4.397542371171649*vtSq[1]*fl[8]+3.396416424888148*vtSq[3]*fr[7]+3.458741190809164*vtSq[1]*fr[7]-3.396416424888148*vtSq[3]*fl[7]-3.458741190809164*vtSq[1]*fl[7]-4.89527035770242*vtSq[2]*fr[4]-5.473078644031159*vtSq[0]*fr[4]-4.89527035770242*vtSq[2]*fl[4]-5.473078644031159*vtSq[0]*fl[4]+3.458741190809164*fr[1]*vtSq[2]-3.458741190809164*fl[1]*vtSq[2]-5.473078644031159*vtSq[1]*fr[2]-5.473078644031159*vtSq[1]*fl[2]+3.866990209613929*fr[0]*vtSq[1]-3.866990209613929*fl[0]*vtSq[1]+3.866990209613929*vtSq[0]*fr[1]-3.866990209613929*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.9354143466934851*favg[18]+0.5477225575051661*favg[11]+0.7905694150420947*favg[8]+0.3162277660168379*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+alphaDrag[2]*(0.3105295017040592*favg[17]+0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])-1.118033988749895*fjump[12]+alphaDrag[0]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+alphaDrag[3]*(0.5378528742004769*favg[11]+0.3105295017040592*favg[7])-0.8660254037844386*fjump[4]-0.5*fjump[1]; 
  Ghat[3] = (4.39754237117165*vtSq[0]*fr[14]-4.39754237117165*vtSq[0]*fl[14]+3.86699020961393*vtSq[2]*fr[13]-3.86699020961393*vtSq[2]*fl[13]-5.473078644031159*vtSq[1]*fr[10]-5.473078644031159*vtSq[1]*fl[10]-5.473078644031159*vtSq[0]*fr[6]-5.473078644031159*vtSq[0]*fl[6]+3.866990209613929*vtSq[1]*fr[5]-3.866990209613929*vtSq[1]*fl[5]+3.866990209613929*vtSq[0]*fr[3]-3.866990209613929*vtSq[0]*fl[3])*rdv-1.118033988749895*fjump[14]+alphaDrag[0]*(0.7905694150420948*favg[14]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])+0.3535533905932737*alphaDrag[2]*favg[13]+alphaDrag[1]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[5])-0.8660254037844386*fjump[6]-0.5*fjump[3]; 
  Ghat[5] = (4.39754237117165*vtSq[1]*fr[14]-4.39754237117165*vtSq[1]*fl[14]+3.396416424888148*vtSq[3]*fr[13]+3.458741190809164*vtSq[1]*fr[13]-3.396416424888148*vtSq[3]*fl[13]-3.458741190809164*vtSq[1]*fl[13]-4.89527035770242*vtSq[2]*fr[10]-5.473078644031159*vtSq[0]*fr[10]-4.89527035770242*vtSq[2]*fl[10]-5.473078644031159*vtSq[0]*fl[10]-5.473078644031159*vtSq[1]*fr[6]-5.473078644031159*vtSq[1]*fl[6]+3.458741190809164*vtSq[2]*fr[5]+3.866990209613929*vtSq[0]*fr[5]-3.458741190809164*vtSq[2]*fl[5]-3.866990209613929*vtSq[0]*fl[5]+3.866990209613929*vtSq[1]*fr[3]-3.866990209613929*vtSq[1]*fl[3])*rdv+alphaDrag[1]*(0.7905694150420948*favg[14]+0.3162277660168379*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])+0.3105295017040593*alphaDrag[3]*favg[13]-0.8660254037844386*fjump[10]+alphaDrag[0]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[5])+alphaDrag[2]*(0.5477225575051661*favg[10]+0.3162277660168379*favg[5])-0.5*fjump[5]; 
  Ghat[7] = ((-1.929292090055312*vtSq[2]*fr[18])-1.929292090055312*vtSq[2]*fl[18]+2.305827460539443*vtSq[3]*fr[17]+3.396416424888148*vtSq[1]*fr[17]-2.305827460539443*vtSq[3]*fl[17]-3.396416424888148*vtSq[1]*fl[17]+3.86240572873861*vtSq[3]*fr[12]+3.933281470350169*vtSq[1]*fr[12]-3.86240572873861*vtSq[3]*fl[12]-3.933281470350169*vtSq[1]*fl[12]-3.496621684073157*vtSq[2]*fr[11]-5.473078644031159*vtSq[0]*fr[11]-3.496621684073157*vtSq[2]*fl[11]-5.473078644031159*vtSq[0]*fl[11]+4.397542371171649*vtSq[2]*fr[8]-4.397542371171649*vtSq[2]*fl[8]+2.470529422006546*vtSq[2]*fr[7]+3.866990209613929*vtSq[0]*fr[7]-2.470529422006546*vtSq[2]*fl[7]-3.866990209613929*vtSq[0]*fl[7]-4.807060063166761*vtSq[3]*fr[4]-4.89527035770242*vtSq[1]*fr[4]-4.807060063166761*vtSq[3]*fl[4]-4.89527035770242*vtSq[1]*fl[4]+3.396416424888148*fr[1]*vtSq[3]-3.396416424888148*fl[1]*vtSq[3]-5.473078644031159*fr[2]*vtSq[2]-5.473078644031159*fl[2]*vtSq[2]+3.866990209613929*fr[0]*vtSq[2]-3.866990209613929*fl[0]*vtSq[2]+3.458741190809164*fr[1]*vtSq[1]-3.458741190809164*fl[1]*vtSq[1])*rdv+alphaDrag[2]*(0.9354143466934851*favg[18]+0.3912303982179757*favg[11]+0.7905694150420947*favg[8]+0.2258769757263128*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+alphaDrag[1]*(0.3105295017040592*favg[17]+0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])+alphaDrag[3]*(0.210818510677892*favg[17]+0.6943650748294133*favg[12]+0.537852874200477*favg[4]+0.3105295017040592*favg[1])-0.8660254037844387*fjump[11]+alphaDrag[0]*(0.6123724356957944*favg[11]+0.3535533905932737*favg[7])-0.5*fjump[7]; 
  Ghat[9] = ((-5.473078644031159*vtSq[0]*fr[16])-5.473078644031159*vtSq[0]*fl[16]+3.86699020961393*vtSq[1]*fr[15]-3.86699020961393*vtSq[1]*fl[15]+3.866990209613929*vtSq[0]*fr[9]-3.866990209613929*vtSq[0]*fl[9])*rdv-0.8660254037844387*fjump[16]+alphaDrag[0]*(0.6123724356957944*favg[16]+0.3535533905932737*favg[9])+0.3535533905932737*alphaDrag[1]*favg[15]-0.5*fjump[9]; 
  Ghat[13] = (4.397542371171649*vtSq[2]*fr[14]-4.397542371171649*vtSq[2]*fl[14]+2.470529422006546*vtSq[2]*fr[13]+3.866990209613929*vtSq[0]*fr[13]-2.470529422006546*vtSq[2]*fl[13]-3.866990209613929*vtSq[0]*fl[13]-4.807060063166761*vtSq[3]*fr[10]-4.89527035770242*vtSq[1]*fr[10]-4.807060063166761*vtSq[3]*fl[10]-4.89527035770242*vtSq[1]*fl[10]-5.473078644031159*vtSq[2]*fr[6]-5.473078644031159*vtSq[2]*fl[6]+3.396416424888148*vtSq[3]*fr[5]+3.458741190809164*vtSq[1]*fr[5]-3.396416424888148*vtSq[3]*fl[5]-3.458741190809164*vtSq[1]*fl[5]+3.86699020961393*vtSq[2]*fr[3]-3.86699020961393*vtSq[2]*fl[3])*rdv+alphaDrag[2]*(0.7905694150420947*favg[14]+0.2258769757263128*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[13]+0.3535533905932737*alphaDrag[0]*favg[13]+alphaDrag[1]*(0.5477225575051661*favg[10]+0.3162277660168379*favg[5])+alphaDrag[3]*(0.5378528742004769*favg[10]+0.3105295017040593*favg[5]); 
  Ghat[15] = ((-5.473078644031159*vtSq[1]*fr[16])-5.473078644031159*vtSq[1]*fl[16]+3.458741190809164*vtSq[2]*fr[15]+3.866990209613929*vtSq[0]*fr[15]-3.458741190809164*vtSq[2]*fl[15]-3.866990209613929*vtSq[0]*fl[15]+3.86699020961393*vtSq[1]*fr[9]-3.86699020961393*vtSq[1]*fl[9])*rdv+alphaDrag[1]*(0.6123724356957944*favg[16]+0.3535533905932737*favg[9])-0.5*fjump[15]+0.3162277660168379*alphaDrag[2]*favg[15]+0.3535533905932737*alphaDrag[0]*favg[15]; 
  Ghat[17] = ((-1.929292090055312*vtSq[3]*fr[18])-1.929292090055312*vtSq[3]*fl[18]+2.305827460539443*vtSq[2]*fr[17]+3.866990209613929*vtSq[0]*fr[17]-2.305827460539443*vtSq[2]*fl[17]-3.866990209613929*vtSq[0]*fl[17]+3.86240572873861*vtSq[2]*fr[12]-3.86240572873861*vtSq[2]*fl[12]-3.263513571801613*vtSq[3]*fr[11]-4.807060063166761*vtSq[1]*fr[11]-3.263513571801613*vtSq[3]*fl[11]-4.807060063166761*vtSq[1]*fl[11]+4.397542371171649*vtSq[3]*fr[8]-4.397542371171649*vtSq[3]*fl[8]+2.305827460539443*vtSq[3]*fr[7]+3.396416424888148*vtSq[1]*fr[7]-2.305827460539443*vtSq[3]*fl[7]-3.396416424888148*vtSq[1]*fl[7]-4.807060063166761*vtSq[2]*fr[4]-4.807060063166761*vtSq[2]*fl[4]-5.473078644031159*fr[2]*vtSq[3]-5.473078644031159*fl[2]*vtSq[3]+3.866990209613929*fr[0]*vtSq[3]-3.866990209613929*fl[0]*vtSq[3]+3.396416424888148*fr[1]*vtSq[2]-3.396416424888148*fl[1]*vtSq[2])*rdv+alphaDrag[3]*(0.9354143466934851*favg[18]+0.3651483716701107*favg[11]+0.7905694150420947*favg[8]+0.210818510677892*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[17]+0.3535533905932737*alphaDrag[0]*favg[17]+alphaDrag[2]*(0.210818510677892*favg[17]+0.6943650748294133*favg[12]+0.537852874200477*favg[4]+0.3105295017040592*favg[1])+alphaDrag[1]*(0.5378528742004769*favg[11]+0.3105295017040592*favg[7]); 
  Ghat[19] = (3.866990209613929*vtSq[0]*fr[19]-3.866990209613929*vtSq[0]*fl[19])*rdv-0.5*fjump[19]+0.3535533905932737*alphaDrag[0]*favg[19]; 

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
  incr2[2] = vtSq[0]*((-0.2314550249431378*fr[18])+0.2314550249431378*fl[18]+0.4279082480509108*(fr[8]+fl[8])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*vtSq[3]*(fr[17]+fl[17])+vtSq[1]*(0.4279082480509107*(fr[12]+fl[12])-0.4640388251536715*fr[4]+0.4640388251536715*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[2]*((-0.4640388251536715*fr[11])+0.4640388251536715*fl[11]+0.3061862178478971*(fr[7]+fl[7])); 
  incr2[4] = vtSq[1]*((-0.2314550249431378*fr[18])+0.2314550249431378*fl[18]-0.4150489428970995*fr[11]+0.4150489428970995*fl[11]+0.4279082480509108*(fr[8]+fl[8])+0.273861278752583*(fr[7]+fl[7])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[2]*(0.2689264371002384*(fr[17]+fl[17])+0.3827327723098714*(fr[12]+fl[12])-0.4150489428970995*fr[4]+0.4150489428970995*fl[4]+0.273861278752583*(fr[1]+fl[1]))+vtSq[0]*(0.4279082480509107*(fr[12]+fl[12])-0.4640388251536715*fr[4]+0.4640388251536715*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[3]*((-0.4075699709865777*fr[11])+0.4075699709865777*fl[11]+0.2689264371002384*(fr[7]+fl[7])); 
  incr2[6] = vtSq[0]*(0.4279082480509107*(fr[14]+fl[14])-0.4640388251536715*fr[6]+0.4640388251536715*fl[6]+0.3061862178478971*(fr[3]+fl[3]))+0.3061862178478971*vtSq[2]*(fr[13]+fl[13])+vtSq[1]*((-0.4640388251536715*fr[10])+0.4640388251536715*fl[10]+0.3061862178478971*(fr[5]+fl[5])); 
  incr2[8] = vtSq[0]*(0.896421457000795*fr[18]-0.896421457000795*fl[18]-1.657281518405969*(fr[8]+fl[8])+1.797214641813825*fr[2]-1.797214641813825*fl[2]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*vtSq[3]*(fr[17]+fl[17])+vtSq[1]*((-1.65728151840597*(fr[12]+fl[12]))+1.797214641813825*fr[4]-1.797214641813825*fl[4]-1.185854122563142*(fr[1]+fl[1]))+vtSq[2]*(1.797214641813825*fr[11]-1.797214641813825*fl[11]-1.185854122563142*(fr[7]+fl[7])); 
  incr2[10] = vtSq[1]*(0.4279082480509107*(fr[14]+fl[14])+0.273861278752583*(fr[13]+fl[13])-0.4640388251536715*fr[6]+0.4640388251536715*fl[6]+0.3061862178478971*(fr[3]+fl[3]))+0.2689264371002384*vtSq[3]*(fr[13]+fl[13])+vtSq[2]*((-0.4150489428970995*fr[10])+0.4150489428970995*fl[10]+0.273861278752583*(fr[5]+fl[5]))+vtSq[0]*((-0.4640388251536715*fr[10])+0.4640388251536715*fl[10]+0.3061862178478971*(fr[5]+fl[5])); 
  incr2[11] = vtSq[2]*((-0.2314550249431378*fr[18])+0.2314550249431378*fl[18]-0.2964635306407854*fr[11]+0.2964635306407854*fl[11]+0.4279082480509108*(fr[8]+fl[8])+0.1956151991089878*(fr[7]+fl[7])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[1]*(0.2689264371002384*(fr[17]+fl[17])+0.3827327723098713*(fr[12]+fl[12])-0.4150489428970995*fr[4]+0.4150489428970995*fl[4]+0.273861278752583*(fr[1]+fl[1]))+vtSq[3]*(0.1825741858350553*(fr[17]+fl[17])+0.3758361214393465*(fr[12]+fl[12])-0.4075699709865777*fr[4]+0.4075699709865777*fl[4]+0.2689264371002384*(fr[1]+fl[1]))+vtSq[0]*((-0.4640388251536715*fr[11])+0.4640388251536715*fl[11]+0.3061862178478971*(fr[7]+fl[7])); 
  incr2[12] = vtSq[1]*(0.8964214570007949*fr[18]-0.8964214570007949*fl[18]+1.60747764370146*fr[11]-1.60747764370146*fl[11]-1.65728151840597*(fr[8]+fl[8])-1.060660171779821*(fr[7]+fl[7])+1.797214641813825*fr[2]-1.797214641813825*fl[2]-1.185854122563142*(fr[0]+fl[0]))+vtSq[2]*((-1.04154761224412*(fr[17]+fl[17]))-1.482317653203927*(fr[12]+fl[12])+1.60747764370146*fr[4]-1.60747764370146*fl[4]-1.060660171779821*(fr[1]+fl[1]))+vtSq[0]*((-1.657281518405969*(fr[12]+fl[12]))+1.797214641813825*fr[4]-1.797214641813825*fl[4]-1.185854122563142*(fr[1]+fl[1]))+vtSq[3]*(1.578511710045255*fr[11]-1.578511710045255*fl[11]-1.04154761224412*(fr[7]+fl[7])); 
  incr2[14] = vtSq[0]*((-1.657281518405969*(fr[14]+fl[14]))+1.797214641813825*fr[6]-1.797214641813825*fl[6]-1.185854122563142*(fr[3]+fl[3]))-1.185854122563142*vtSq[2]*(fr[13]+fl[13])+vtSq[1]*(1.797214641813825*fr[10]-1.797214641813825*fl[10]-1.185854122563142*(fr[5]+fl[5])); 
  incr2[16] = vtSq[0]*((-0.4640388251536715*fr[16])+0.4640388251536715*fl[16]+0.3061862178478971*(fr[9]+fl[9]))+0.3061862178478971*vtSq[1]*(fr[15]+fl[15]); 
  incr2[18] = vtSq[0]*((-2.121320343559642*fr[18])+2.121320343559642*fl[18]+3.921843874378478*(fr[8]+fl[8])-4.252986083330156*fr[2]+4.252986083330156*fl[2]+2.806243040080455*(fr[0]+fl[0]))+2.806243040080455*vtSq[3]*(fr[17]+fl[17])+vtSq[1]*(3.921843874378477*(fr[12]+fl[12])-4.252986083330156*fr[4]+4.252986083330156*fl[4]+2.806243040080455*(fr[1]+fl[1]))+vtSq[2]*((-4.252986083330157*fr[11])+4.252986083330157*fl[11]+2.806243040080455*(fr[7]+fl[7])); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr2[2]*rdvSq4nuR+incr1[2]*rdv2nuR; 
  outr[3] += incr1[3]*rdv2nuR; 
  outr[4] += incr2[4]*rdvSq4nuR+incr1[4]*rdv2nuR; 
  outr[5] += incr1[5]*rdv2nuR; 
  outr[6] += incr2[6]*rdvSq4nuR+incr1[6]*rdv2nuR; 
  outr[7] += incr1[7]*rdv2nuR; 
  outr[8] += incr2[8]*rdvSq4nuR+incr1[8]*rdv2nuR; 
  outr[9] += incr1[9]*rdv2nuR; 
  outr[10] += incr2[10]*rdvSq4nuR+incr1[10]*rdv2nuR; 
  outr[11] += incr2[11]*rdvSq4nuR+incr1[11]*rdv2nuR; 
  outr[12] += incr2[12]*rdvSq4nuR+incr1[12]*rdv2nuR; 
  outr[13] += incr1[13]*rdv2nuR; 
  outr[14] += incr2[14]*rdvSq4nuR+incr1[14]*rdv2nuR; 
  outr[15] += incr1[15]*rdv2nuR; 
  outr[16] += incr2[16]*rdvSq4nuR+incr1[16]*rdv2nuR; 
  outr[17] += incr1[17]*rdv2nuR; 
  outr[18] += incr2[18]*rdvSq4nuR+incr1[18]*rdv2nuR; 
  outr[19] += incr1[19]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += incr1[2]*rdv2nuL-1.0*incr2[2]*rdvSq4nuL; 
  outl[3] += -1.0*incr1[3]*rdv2nuL; 
  outl[4] += incr1[4]*rdv2nuL-1.0*incr2[4]*rdvSq4nuL; 
  outl[5] += -1.0*incr1[5]*rdv2nuL; 
  outl[6] += incr1[6]*rdv2nuL-1.0*incr2[6]*rdvSq4nuL; 
  outl[7] += -1.0*incr1[7]*rdv2nuL; 
  outl[8] += incr2[8]*rdvSq4nuL-1.0*incr1[8]*rdv2nuL; 
  outl[9] += -1.0*incr1[9]*rdv2nuL; 
  outl[10] += incr1[10]*rdv2nuL-1.0*incr2[10]*rdvSq4nuL; 
  outl[11] += incr1[11]*rdv2nuL-1.0*incr2[11]*rdvSq4nuL; 
  outl[12] += incr2[12]*rdvSq4nuL-1.0*incr1[12]*rdv2nuL; 
  outl[13] += -1.0*incr1[13]*rdv2nuL; 
  outl[14] += incr2[14]*rdvSq4nuL-1.0*incr1[14]*rdv2nuL; 
  outl[15] += -1.0*incr1[15]*rdv2nuL; 
  outl[16] += incr1[16]*rdv2nuL-1.0*incr2[16]*rdvSq4nuL; 
  outl[17] += -1.0*incr1[17]*rdv2nuL; 
  outl[18] += incr1[18]*rdv2nuL-1.0*incr2[18]*rdvSq4nuL; 
  outl[19] += -1.0*incr1[19]*rdv2nuL; 

  return std::abs(0.7905694150420947*Ux[2]+wl[1]-0.7071067811865475*Ux[0]); 
} 
double VmLBOconstNuSurf1x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:       Cell-center coordinates. 
  // dxv[3]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[4]:       bulk velocity (in 2 directions). 
  // vtSq[2]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2nuL = 2.0*nu/dxvl[2]; 
  double rdv2nuR = 2.0*nu/dxvr[2]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[2]*dxvl[2]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[2]*dxvr[2]); 

  const double *Uy = &u[2]; 

  double favg[4]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 

  double fjump[4]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(-1*fr[3])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[2]+0.7071067811865475*dxvl[2]-1.0*Uy[0]; 
  alphaDrag[1] = -1.0*Uy[1]; 

  double Ghat[4]; 
  for(unsigned short int i=0; i<4; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.530931089239486*vtSq[0]*fr[3])-1.530931089239486*vtSq[0]*fl[3]+1.590990257669731*fr[1]*vtSq[1]-1.590990257669731*fl[1]*vtSq[1]+1.590990257669731*fr[0]*vtSq[0]-1.590990257669731*fl[0]*vtSq[0])*rdv-0.8660254037844386*fjump[3]+alphaDrag[0]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = ((-1.530931089239486*vtSq[1]*fr[3])-1.530931089239486*vtSq[1]*fl[3]+1.590990257669731*fr[0]*vtSq[1]-1.590990257669731*fl[0]*vtSq[1]+1.590990257669731*vtSq[0]*fr[1]-1.590990257669731*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.5*fjump[1]+0.3535533905932737*alphaDrag[0]*favg[1]; 
  Ghat[2] = (1.590990257669731*vtSq[0]*fr[2]-1.590990257669731*vtSq[0]*fl[2])*rdv-0.5*fjump[2]+0.3535533905932737*alphaDrag[0]*favg[2]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 

  double incr2[4]; 
  incr2[3] = vtSq[0]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*(fr[1]+fl[1])*vtSq[1]; 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr1[2]*rdv2nuR; 
  outr[3] += incr2[3]*rdvSq4nuR+incr1[3]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += -1.0*incr1[2]*rdv2nuL; 
  outl[3] += incr1[3]*rdv2nuL-1.0*incr2[3]*rdvSq4nuL; 

  return std::abs(wl[2]-0.7071067811865475*Uy[0]); 
} 
double VmLBOconstNuSurf1x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:       Cell-center coordinates. 
  // dxv[3]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[6]:       bulk velocity (in 2 directions). 
  // vtSq[3]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2nuL = 2.0*nu/dxvl[2]; 
  double rdv2nuR = 2.0*nu/dxvr[2]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[2]*dxvl[2]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[2]*dxvr[2]); 

  const double *Uy = &u[3]; 

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
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = vMuMidMax*(fl[9]-(1*fr[9])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[2]+0.7071067811865475*dxvl[2]-1.0*Uy[0]; 
  alphaDrag[1] = -1.0*Uy[1]; 
  alphaDrag[2] = -1.0*Uy[2]; 

  double Ghat[10]; 
  for(unsigned short int i=0; i<10; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*vtSq[0]*fr[9]-1.897366596101028*vtSq[0]*fl[9]+2.651650429449552*vtSq[2]*fr[7]-2.651650429449552*vtSq[2]*fl[7]-3.368048396326869*vtSq[1]*fr[5]-3.368048396326869*vtSq[1]*fl[5]-3.368048396326869*vtSq[0]*fr[3]-3.368048396326869*vtSq[0]*fl[3]+2.651650429449552*fr[1]*vtSq[1]-2.651650429449552*fl[1]*vtSq[1]+2.651650429449552*fr[0]*vtSq[0]-2.651650429449552*fl[0]*vtSq[0])*rdv-1.118033988749895*fjump[9]+alphaDrag[0]*(0.7905694150420947*favg[9]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[2]*favg[7]+alphaDrag[1]*(0.6123724356957944*favg[5]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = (1.897366596101028*vtSq[1]*fr[9]-1.897366596101028*vtSq[1]*fl[9]+2.371708245126284*vtSq[1]*fr[7]-2.371708245126284*vtSq[1]*fl[7]-3.012474066278413*vtSq[2]*fr[5]-3.368048396326869*vtSq[0]*fr[5]-3.012474066278413*vtSq[2]*fl[5]-3.368048396326869*vtSq[0]*fl[5]-3.368048396326869*vtSq[1]*fr[3]-3.368048396326869*vtSq[1]*fl[3]+2.371708245126284*fr[1]*vtSq[2]-2.371708245126284*fl[1]*vtSq[2]+2.651650429449552*fr[0]*vtSq[1]-2.651650429449552*fl[0]*vtSq[1]+2.651650429449552*vtSq[0]*fr[1]-2.651650429449552*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.7905694150420947*favg[9]+0.3162277660168379*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[5]+alphaDrag[0]*(0.6123724356957944*favg[5]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.5477225575051661*favg[5]+0.3162277660168379*favg[1])-0.5*fjump[1]; 
  Ghat[2] = ((-3.368048396326869*vtSq[0]*fr[6])-3.368048396326869*vtSq[0]*fl[6]+2.651650429449552*vtSq[1]*fr[4]-2.651650429449552*vtSq[1]*fl[4]+2.651650429449552*vtSq[0]*fr[2]-2.651650429449552*vtSq[0]*fl[2])*rdv-0.8660254037844386*fjump[6]+alphaDrag[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[2])+0.3535533905932737*alphaDrag[1]*favg[4]-0.5*fjump[2]; 
  Ghat[4] = ((-3.368048396326869*vtSq[1]*fr[6])-3.368048396326869*vtSq[1]*fl[6]+2.371708245126284*vtSq[2]*fr[4]+2.651650429449552*vtSq[0]*fr[4]-2.371708245126284*vtSq[2]*fl[4]-2.651650429449552*vtSq[0]*fl[4]+2.651650429449552*vtSq[1]*fr[2]-2.651650429449552*vtSq[1]*fl[2])*rdv+alphaDrag[1]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-0.5*fjump[4]+0.3162277660168379*alphaDrag[2]*favg[4]+0.3535533905932737*alphaDrag[0]*favg[4]; 
  Ghat[7] = (1.897366596101028*vtSq[2]*fr[9]-1.897366596101028*vtSq[2]*fl[9]+1.694077317947346*vtSq[2]*fr[7]+2.651650429449552*vtSq[0]*fr[7]-1.694077317947346*vtSq[2]*fl[7]-2.651650429449552*vtSq[0]*fl[7]-3.012474066278413*vtSq[1]*fr[5]-3.012474066278413*vtSq[1]*fl[5]-3.368048396326869*vtSq[2]*fr[3]-3.368048396326869*vtSq[2]*fl[3]+2.651650429449552*fr[0]*vtSq[2]-2.651650429449552*fl[0]*vtSq[2]+2.371708245126284*fr[1]*vtSq[1]-2.371708245126284*fl[1]*vtSq[1])*rdv+alphaDrag[2]*(0.7905694150420947*favg[9]+0.2258769757263128*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.5*fjump[7]+0.3535533905932737*alphaDrag[0]*favg[7]+alphaDrag[1]*(0.5477225575051661*favg[5]+0.3162277660168379*favg[1]); 
  Ghat[8] = (2.651650429449552*vtSq[0]*fr[8]-2.651650429449552*vtSq[0]*fl[8])*rdv-0.5*fjump[8]+0.3535533905932737*alphaDrag[0]*favg[8]; 

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
  incr2[3] = vtSq[0]*(0.2995357736356374*(fr[9]+fl[9])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*vtSq[2]*(fr[7]+fl[7])+vtSq[1]*((-0.430893194785552*fr[5])+0.430893194785552*fl[5]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[5] = vtSq[1]*(0.2995357736356374*(fr[9]+fl[9])+0.273861278752583*(fr[7]+fl[7])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[2]*((-0.3854025898330209*fr[5])+0.3854025898330209*fl[5]+0.273861278752583*(fr[1]+fl[1]))+vtSq[0]*((-0.430893194785552*fr[5])+0.430893194785552*fl[5]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[6] = vtSq[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[2]+fl[2]))+0.3061862178478971*vtSq[1]*(fr[4]+fl[4]); 
  incr2[9] = vtSq[0]*((-1.160097062884178*(fr[9]+fl[9]))+1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*vtSq[2]*(fr[7]+fl[7])+vtSq[1]*(1.668842167398551*fr[5]-1.668842167398551*fl[5]-1.185854122563142*(fr[1]+fl[1])); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr1[2]*rdv2nuR; 
  outr[3] += incr2[3]*rdvSq4nuR+incr1[3]*rdv2nuR; 
  outr[4] += incr1[4]*rdv2nuR; 
  outr[5] += incr2[5]*rdvSq4nuR+incr1[5]*rdv2nuR; 
  outr[6] += incr2[6]*rdvSq4nuR+incr1[6]*rdv2nuR; 
  outr[7] += incr1[7]*rdv2nuR; 
  outr[8] += incr1[8]*rdv2nuR; 
  outr[9] += incr2[9]*rdvSq4nuR+incr1[9]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += -1.0*incr1[2]*rdv2nuL; 
  outl[3] += incr1[3]*rdv2nuL-1.0*incr2[3]*rdvSq4nuL; 
  outl[4] += -1.0*incr1[4]*rdv2nuL; 
  outl[5] += incr1[5]*rdv2nuL-1.0*incr2[5]*rdvSq4nuL; 
  outl[6] += incr1[6]*rdv2nuL-1.0*incr2[6]*rdvSq4nuL; 
  outl[7] += -1.0*incr1[7]*rdv2nuL; 
  outl[8] += -1.0*incr1[8]*rdv2nuL; 
  outl[9] += incr2[9]*rdvSq4nuL-1.0*incr1[9]*rdv2nuL; 

  return std::abs(wl[2]+0.7905694150420947*Uy[2]-0.7071067811865475*Uy[0]); 
} 
double VmLBOconstNuSurf1x2vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:       Cell-center coordinates. 
  // dxv[3]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[8]:       bulk velocity (in 2 directions). 
  // vtSq[4]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2nuL = 2.0*nu/dxvl[2]; 
  double rdv2nuR = 2.0*nu/dxvr[2]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[2]*dxvl[2]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[2]*dxvr[2]); 

  const double *Uy = &u[4]; 

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
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = vMuMidMax*(fl[13]-(-1*fr[13])); 
  fjump[14] = vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = vMuMidMax*(fl[16]-(1*fr[16])); 
  fjump[17] = vMuMidMax*(fl[17]-(1*fr[17])); 
  fjump[18] = vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = vMuMidMax*(fl[19]-(-1*fr[19])); 

  double alphaDrag[4]; 
  alphaDrag[0] = 1.414213562373095*wl[2]+0.7071067811865475*dxvl[2]-1.0*Uy[0]; 
  alphaDrag[1] = -1.0*Uy[1]; 
  alphaDrag[2] = -1.0*Uy[2]; 
  alphaDrag[3] = -1.0*Uy[3]; 

  double Ghat[20]; 
  for(unsigned short int i=0; i<20; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.929292090055312*vtSq[0]*fr[19])-1.929292090055312*vtSq[0]*fl[19]+3.866990209613929*vtSq[3]*fr[17]-3.866990209613929*vtSq[3]*fl[17]+4.39754237117165*vtSq[1]*fr[15]-4.39754237117165*vtSq[1]*fl[15]-5.473078644031159*vtSq[2]*fr[13]-5.473078644031159*vtSq[2]*fl[13]+4.397542371171649*vtSq[0]*fr[9]-4.397542371171649*vtSq[0]*fl[9]+3.866990209613929*vtSq[2]*fr[7]-3.866990209613929*vtSq[2]*fl[7]-5.473078644031159*vtSq[1]*fr[5]-5.473078644031159*vtSq[1]*fl[5]-5.473078644031159*vtSq[0]*fr[3]-5.473078644031159*vtSq[0]*fl[3]+3.866990209613929*fr[1]*vtSq[1]-3.866990209613929*fl[1]*vtSq[1]+3.866990209613929*fr[0]*vtSq[0]-3.866990209613929*fl[0]*vtSq[0])*rdv-1.322875655532295*fjump[19]+alphaDrag[0]*(0.9354143466934851*favg[19]+0.7905694150420947*favg[9]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[3]*favg[17]+alphaDrag[1]*(0.7905694150420948*favg[15]+0.6123724356957944*favg[5]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.6123724356957944*favg[13]+0.3535533905932737*favg[7])-1.118033988749895*fjump[9]-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = ((-1.929292090055312*vtSq[1]*fr[19])-1.929292090055312*vtSq[1]*fl[19]+3.396416424888148*vtSq[2]*fr[17]-3.396416424888148*vtSq[2]*fl[17]+3.933281470350169*vtSq[2]*fr[15]+4.39754237117165*vtSq[0]*fr[15]-3.933281470350169*vtSq[2]*fl[15]-4.39754237117165*vtSq[0]*fl[15]-4.807060063166761*vtSq[3]*fr[13]-4.89527035770242*vtSq[1]*fr[13]-4.807060063166761*vtSq[3]*fl[13]-4.89527035770242*vtSq[1]*fl[13]+4.397542371171649*vtSq[1]*fr[9]-4.397542371171649*vtSq[1]*fl[9]+3.396416424888148*vtSq[3]*fr[7]+3.458741190809164*vtSq[1]*fr[7]-3.396416424888148*vtSq[3]*fl[7]-3.458741190809164*vtSq[1]*fl[7]-4.89527035770242*vtSq[2]*fr[5]-5.473078644031159*vtSq[0]*fr[5]-4.89527035770242*vtSq[2]*fl[5]-5.473078644031159*vtSq[0]*fl[5]-5.473078644031159*vtSq[1]*fr[3]-5.473078644031159*vtSq[1]*fl[3]+3.458741190809164*fr[1]*vtSq[2]-3.458741190809164*fl[1]*vtSq[2]+3.866990209613929*fr[0]*vtSq[1]-3.866990209613929*fl[0]*vtSq[1]+3.866990209613929*vtSq[0]*fr[1]-3.866990209613929*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.9354143466934851*favg[19]+0.5477225575051661*favg[13]+0.7905694150420947*favg[9]+0.3162277660168379*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+alphaDrag[2]*(0.3105295017040592*favg[17]+0.7071067811865475*favg[15]+0.5477225575051661*favg[5]+0.3162277660168379*favg[1])-1.118033988749895*fjump[15]+alphaDrag[0]*(0.7905694150420948*favg[15]+0.6123724356957944*favg[5]+0.3535533905932737*favg[1])+alphaDrag[3]*(0.5378528742004769*favg[13]+0.3105295017040592*favg[7])-0.8660254037844386*fjump[5]-0.5*fjump[1]; 
  Ghat[2] = (4.39754237117165*vtSq[0]*fr[16]-4.39754237117165*vtSq[0]*fl[16]+3.86699020961393*vtSq[2]*fr[11]-3.86699020961393*vtSq[2]*fl[11]-5.473078644031159*vtSq[1]*fr[10]-5.473078644031159*vtSq[1]*fl[10]-5.473078644031159*vtSq[0]*fr[6]-5.473078644031159*vtSq[0]*fl[6]+3.866990209613929*vtSq[1]*fr[4]-3.866990209613929*vtSq[1]*fl[4]+3.866990209613929*vtSq[0]*fr[2]-3.866990209613929*vtSq[0]*fl[2])*rdv-1.118033988749895*fjump[16]+alphaDrag[0]*(0.7905694150420948*favg[16]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])+0.3535533905932737*alphaDrag[2]*favg[11]+alphaDrag[1]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[4])-0.8660254037844386*fjump[6]-0.5*fjump[2]; 
  Ghat[4] = (4.39754237117165*vtSq[1]*fr[16]-4.39754237117165*vtSq[1]*fl[16]+3.396416424888148*vtSq[3]*fr[11]+3.458741190809164*vtSq[1]*fr[11]-3.396416424888148*vtSq[3]*fl[11]-3.458741190809164*vtSq[1]*fl[11]-4.89527035770242*vtSq[2]*fr[10]-5.473078644031159*vtSq[0]*fr[10]-4.89527035770242*vtSq[2]*fl[10]-5.473078644031159*vtSq[0]*fl[10]-5.473078644031159*vtSq[1]*fr[6]-5.473078644031159*vtSq[1]*fl[6]+3.458741190809164*vtSq[2]*fr[4]+3.866990209613929*vtSq[0]*fr[4]-3.458741190809164*vtSq[2]*fl[4]-3.866990209613929*vtSq[0]*fl[4]+3.866990209613929*vtSq[1]*fr[2]-3.866990209613929*vtSq[1]*fl[2])*rdv+alphaDrag[1]*(0.7905694150420948*favg[16]+0.3162277660168379*favg[11]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])+0.3105295017040593*alphaDrag[3]*favg[11]-0.8660254037844386*fjump[10]+alphaDrag[0]*(0.6123724356957944*favg[10]+0.3535533905932737*favg[4])+alphaDrag[2]*(0.5477225575051661*favg[10]+0.3162277660168379*favg[4])-0.5*fjump[4]; 
  Ghat[7] = ((-1.929292090055312*vtSq[2]*fr[19])-1.929292090055312*vtSq[2]*fl[19]+2.305827460539443*vtSq[3]*fr[17]+3.396416424888148*vtSq[1]*fr[17]-2.305827460539443*vtSq[3]*fl[17]-3.396416424888148*vtSq[1]*fl[17]+3.86240572873861*vtSq[3]*fr[15]+3.933281470350169*vtSq[1]*fr[15]-3.86240572873861*vtSq[3]*fl[15]-3.933281470350169*vtSq[1]*fl[15]-3.496621684073157*vtSq[2]*fr[13]-5.473078644031159*vtSq[0]*fr[13]-3.496621684073157*vtSq[2]*fl[13]-5.473078644031159*vtSq[0]*fl[13]+4.397542371171649*vtSq[2]*fr[9]-4.397542371171649*vtSq[2]*fl[9]+2.470529422006546*vtSq[2]*fr[7]+3.866990209613929*vtSq[0]*fr[7]-2.470529422006546*vtSq[2]*fl[7]-3.866990209613929*vtSq[0]*fl[7]-4.807060063166761*vtSq[3]*fr[5]-4.89527035770242*vtSq[1]*fr[5]-4.807060063166761*vtSq[3]*fl[5]-4.89527035770242*vtSq[1]*fl[5]+3.396416424888148*fr[1]*vtSq[3]-3.396416424888148*fl[1]*vtSq[3]-5.473078644031159*vtSq[2]*fr[3]-5.473078644031159*vtSq[2]*fl[3]+3.866990209613929*fr[0]*vtSq[2]-3.866990209613929*fl[0]*vtSq[2]+3.458741190809164*fr[1]*vtSq[1]-3.458741190809164*fl[1]*vtSq[1])*rdv+alphaDrag[2]*(0.9354143466934851*favg[19]+0.3912303982179757*favg[13]+0.7905694150420947*favg[9]+0.2258769757263128*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+alphaDrag[1]*(0.3105295017040592*favg[17]+0.7071067811865475*favg[15]+0.5477225575051661*favg[5]+0.3162277660168379*favg[1])+alphaDrag[3]*(0.210818510677892*favg[17]+0.6943650748294133*favg[15]+0.537852874200477*favg[5]+0.3105295017040592*favg[1])-0.8660254037844387*fjump[13]+alphaDrag[0]*(0.6123724356957944*favg[13]+0.3535533905932737*favg[7])-0.5*fjump[7]; 
  Ghat[8] = ((-5.473078644031159*vtSq[0]*fr[14])-5.473078644031159*vtSq[0]*fl[14]+3.86699020961393*vtSq[1]*fr[12]-3.86699020961393*vtSq[1]*fl[12]+3.866990209613929*vtSq[0]*fr[8]-3.866990209613929*vtSq[0]*fl[8])*rdv-0.8660254037844387*fjump[14]+alphaDrag[0]*(0.6123724356957944*favg[14]+0.3535533905932737*favg[8])+0.3535533905932737*alphaDrag[1]*favg[12]-0.5*fjump[8]; 
  Ghat[11] = (4.397542371171649*vtSq[2]*fr[16]-4.397542371171649*vtSq[2]*fl[16]+2.470529422006546*vtSq[2]*fr[11]+3.866990209613929*vtSq[0]*fr[11]-2.470529422006546*vtSq[2]*fl[11]-3.866990209613929*vtSq[0]*fl[11]-4.807060063166761*vtSq[3]*fr[10]-4.89527035770242*vtSq[1]*fr[10]-4.807060063166761*vtSq[3]*fl[10]-4.89527035770242*vtSq[1]*fl[10]-5.473078644031159*vtSq[2]*fr[6]-5.473078644031159*vtSq[2]*fl[6]+3.396416424888148*vtSq[3]*fr[4]+3.458741190809164*vtSq[1]*fr[4]-3.396416424888148*vtSq[3]*fl[4]-3.458741190809164*vtSq[1]*fl[4]+3.86699020961393*fr[2]*vtSq[2]-3.86699020961393*fl[2]*vtSq[2])*rdv+alphaDrag[2]*(0.7905694150420947*favg[16]+0.2258769757263128*favg[11]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-0.5*fjump[11]+0.3535533905932737*alphaDrag[0]*favg[11]+alphaDrag[1]*(0.5477225575051661*favg[10]+0.3162277660168379*favg[4])+alphaDrag[3]*(0.5378528742004769*favg[10]+0.3105295017040593*favg[4]); 
  Ghat[12] = ((-5.473078644031159*vtSq[1]*fr[14])-5.473078644031159*vtSq[1]*fl[14]+3.458741190809164*vtSq[2]*fr[12]+3.866990209613929*vtSq[0]*fr[12]-3.458741190809164*vtSq[2]*fl[12]-3.866990209613929*vtSq[0]*fl[12]+3.86699020961393*vtSq[1]*fr[8]-3.86699020961393*vtSq[1]*fl[8])*rdv+alphaDrag[1]*(0.6123724356957944*favg[14]+0.3535533905932737*favg[8])-0.5*fjump[12]+0.3162277660168379*alphaDrag[2]*favg[12]+0.3535533905932737*alphaDrag[0]*favg[12]; 
  Ghat[17] = ((-1.929292090055312*vtSq[3]*fr[19])-1.929292090055312*vtSq[3]*fl[19]+2.305827460539443*vtSq[2]*fr[17]+3.866990209613929*vtSq[0]*fr[17]-2.305827460539443*vtSq[2]*fl[17]-3.866990209613929*vtSq[0]*fl[17]+3.86240572873861*vtSq[2]*fr[15]-3.86240572873861*vtSq[2]*fl[15]-3.263513571801613*vtSq[3]*fr[13]-4.807060063166761*vtSq[1]*fr[13]-3.263513571801613*vtSq[3]*fl[13]-4.807060063166761*vtSq[1]*fl[13]+4.397542371171649*vtSq[3]*fr[9]-4.397542371171649*vtSq[3]*fl[9]+2.305827460539443*vtSq[3]*fr[7]+3.396416424888148*vtSq[1]*fr[7]-2.305827460539443*vtSq[3]*fl[7]-3.396416424888148*vtSq[1]*fl[7]-4.807060063166761*vtSq[2]*fr[5]-4.807060063166761*vtSq[2]*fl[5]-5.473078644031159*fr[3]*vtSq[3]-5.473078644031159*fl[3]*vtSq[3]+3.866990209613929*fr[0]*vtSq[3]-3.866990209613929*fl[0]*vtSq[3]+3.396416424888148*fr[1]*vtSq[2]-3.396416424888148*fl[1]*vtSq[2])*rdv+alphaDrag[3]*(0.9354143466934851*favg[19]+0.3651483716701107*favg[13]+0.7905694150420947*favg[9]+0.210818510677892*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.5*fjump[17]+0.3535533905932737*alphaDrag[0]*favg[17]+alphaDrag[2]*(0.210818510677892*favg[17]+0.6943650748294133*favg[15]+0.537852874200477*favg[5]+0.3105295017040592*favg[1])+alphaDrag[1]*(0.5378528742004769*favg[13]+0.3105295017040592*favg[7]); 
  Ghat[18] = (3.866990209613929*vtSq[0]*fr[18]-3.866990209613929*vtSq[0]*fl[18])*rdv-0.5*fjump[18]+0.3535533905932737*alphaDrag[0]*favg[18]; 

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
  incr2[3] = vtSq[0]*((-0.2314550249431378*fr[19])+0.2314550249431378*fl[19]+0.4279082480509108*(fr[9]+fl[9])-0.4640388251536715*fr[3]+0.4640388251536715*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*vtSq[3]*(fr[17]+fl[17])+vtSq[1]*(0.4279082480509107*(fr[15]+fl[15])-0.4640388251536715*fr[5]+0.4640388251536715*fl[5]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[2]*((-0.4640388251536715*fr[13])+0.4640388251536715*fl[13]+0.3061862178478971*(fr[7]+fl[7])); 
  incr2[5] = vtSq[1]*((-0.2314550249431378*fr[19])+0.2314550249431378*fl[19]-0.4150489428970995*fr[13]+0.4150489428970995*fl[13]+0.4279082480509108*(fr[9]+fl[9])+0.273861278752583*(fr[7]+fl[7])-0.4640388251536715*fr[3]+0.4640388251536715*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[2]*(0.2689264371002384*(fr[17]+fl[17])+0.3827327723098714*(fr[15]+fl[15])-0.4150489428970995*fr[5]+0.4150489428970995*fl[5]+0.273861278752583*(fr[1]+fl[1]))+vtSq[0]*(0.4279082480509107*(fr[15]+fl[15])-0.4640388251536715*fr[5]+0.4640388251536715*fl[5]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[3]*((-0.4075699709865777*fr[13])+0.4075699709865777*fl[13]+0.2689264371002384*(fr[7]+fl[7])); 
  incr2[6] = vtSq[0]*(0.4279082480509107*(fr[16]+fl[16])-0.4640388251536715*fr[6]+0.4640388251536715*fl[6]+0.3061862178478971*(fr[2]+fl[2]))+0.3061862178478971*vtSq[2]*(fr[11]+fl[11])+vtSq[1]*((-0.4640388251536715*fr[10])+0.4640388251536715*fl[10]+0.3061862178478971*(fr[4]+fl[4])); 
  incr2[9] = vtSq[0]*(0.896421457000795*fr[19]-0.896421457000795*fl[19]-1.657281518405969*(fr[9]+fl[9])+1.797214641813825*fr[3]-1.797214641813825*fl[3]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*vtSq[3]*(fr[17]+fl[17])+vtSq[1]*((-1.65728151840597*(fr[15]+fl[15]))+1.797214641813825*fr[5]-1.797214641813825*fl[5]-1.185854122563142*(fr[1]+fl[1]))+vtSq[2]*(1.797214641813825*fr[13]-1.797214641813825*fl[13]-1.185854122563142*(fr[7]+fl[7])); 
  incr2[10] = vtSq[1]*(0.4279082480509107*(fr[16]+fl[16])+0.273861278752583*(fr[11]+fl[11])-0.4640388251536715*fr[6]+0.4640388251536715*fl[6]+0.3061862178478971*(fr[2]+fl[2]))+0.2689264371002384*vtSq[3]*(fr[11]+fl[11])+vtSq[2]*((-0.4150489428970995*fr[10])+0.4150489428970995*fl[10]+0.273861278752583*(fr[4]+fl[4]))+vtSq[0]*((-0.4640388251536715*fr[10])+0.4640388251536715*fl[10]+0.3061862178478971*(fr[4]+fl[4])); 
  incr2[13] = vtSq[2]*((-0.2314550249431378*fr[19])+0.2314550249431378*fl[19]-0.2964635306407854*fr[13]+0.2964635306407854*fl[13]+0.4279082480509108*(fr[9]+fl[9])+0.1956151991089878*(fr[7]+fl[7])-0.4640388251536715*fr[3]+0.4640388251536715*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[1]*(0.2689264371002384*(fr[17]+fl[17])+0.3827327723098713*(fr[15]+fl[15])-0.4150489428970995*fr[5]+0.4150489428970995*fl[5]+0.273861278752583*(fr[1]+fl[1]))+vtSq[3]*(0.1825741858350553*(fr[17]+fl[17])+0.3758361214393465*(fr[15]+fl[15])-0.4075699709865777*fr[5]+0.4075699709865777*fl[5]+0.2689264371002384*(fr[1]+fl[1]))+vtSq[0]*((-0.4640388251536715*fr[13])+0.4640388251536715*fl[13]+0.3061862178478971*(fr[7]+fl[7])); 
  incr2[14] = vtSq[0]*((-0.4640388251536715*fr[14])+0.4640388251536715*fl[14]+0.3061862178478971*(fr[8]+fl[8]))+0.3061862178478971*vtSq[1]*(fr[12]+fl[12]); 
  incr2[15] = vtSq[1]*(0.8964214570007949*fr[19]-0.8964214570007949*fl[19]+1.60747764370146*fr[13]-1.60747764370146*fl[13]-1.65728151840597*(fr[9]+fl[9])-1.060660171779821*(fr[7]+fl[7])+1.797214641813825*fr[3]-1.797214641813825*fl[3]-1.185854122563142*(fr[0]+fl[0]))+vtSq[2]*((-1.04154761224412*(fr[17]+fl[17]))-1.482317653203927*(fr[15]+fl[15])+1.60747764370146*fr[5]-1.60747764370146*fl[5]-1.060660171779821*(fr[1]+fl[1]))+vtSq[0]*((-1.657281518405969*(fr[15]+fl[15]))+1.797214641813825*fr[5]-1.797214641813825*fl[5]-1.185854122563142*(fr[1]+fl[1]))+vtSq[3]*(1.578511710045255*fr[13]-1.578511710045255*fl[13]-1.04154761224412*(fr[7]+fl[7])); 
  incr2[16] = vtSq[0]*((-1.657281518405969*(fr[16]+fl[16]))+1.797214641813825*fr[6]-1.797214641813825*fl[6]-1.185854122563142*(fr[2]+fl[2]))-1.185854122563142*vtSq[2]*(fr[11]+fl[11])+vtSq[1]*(1.797214641813825*fr[10]-1.797214641813825*fl[10]-1.185854122563142*(fr[4]+fl[4])); 
  incr2[19] = vtSq[0]*((-2.121320343559642*fr[19])+2.121320343559642*fl[19]+3.921843874378478*(fr[9]+fl[9])-4.252986083330156*fr[3]+4.252986083330156*fl[3]+2.806243040080455*(fr[0]+fl[0]))+2.806243040080455*vtSq[3]*(fr[17]+fl[17])+vtSq[1]*(3.921843874378477*(fr[15]+fl[15])-4.252986083330156*fr[5]+4.252986083330156*fl[5]+2.806243040080455*(fr[1]+fl[1]))+vtSq[2]*((-4.252986083330157*fr[13])+4.252986083330157*fl[13]+2.806243040080455*(fr[7]+fl[7])); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr1[2]*rdv2nuR; 
  outr[3] += incr2[3]*rdvSq4nuR+incr1[3]*rdv2nuR; 
  outr[4] += incr1[4]*rdv2nuR; 
  outr[5] += incr2[5]*rdvSq4nuR+incr1[5]*rdv2nuR; 
  outr[6] += incr2[6]*rdvSq4nuR+incr1[6]*rdv2nuR; 
  outr[7] += incr1[7]*rdv2nuR; 
  outr[8] += incr1[8]*rdv2nuR; 
  outr[9] += incr2[9]*rdvSq4nuR+incr1[9]*rdv2nuR; 
  outr[10] += incr2[10]*rdvSq4nuR+incr1[10]*rdv2nuR; 
  outr[11] += incr1[11]*rdv2nuR; 
  outr[12] += incr1[12]*rdv2nuR; 
  outr[13] += incr2[13]*rdvSq4nuR+incr1[13]*rdv2nuR; 
  outr[14] += incr2[14]*rdvSq4nuR+incr1[14]*rdv2nuR; 
  outr[15] += incr2[15]*rdvSq4nuR+incr1[15]*rdv2nuR; 
  outr[16] += incr2[16]*rdvSq4nuR+incr1[16]*rdv2nuR; 
  outr[17] += incr1[17]*rdv2nuR; 
  outr[18] += incr1[18]*rdv2nuR; 
  outr[19] += incr2[19]*rdvSq4nuR+incr1[19]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += -1.0*incr1[2]*rdv2nuL; 
  outl[3] += incr1[3]*rdv2nuL-1.0*incr2[3]*rdvSq4nuL; 
  outl[4] += -1.0*incr1[4]*rdv2nuL; 
  outl[5] += incr1[5]*rdv2nuL-1.0*incr2[5]*rdvSq4nuL; 
  outl[6] += incr1[6]*rdv2nuL-1.0*incr2[6]*rdvSq4nuL; 
  outl[7] += -1.0*incr1[7]*rdv2nuL; 
  outl[8] += -1.0*incr1[8]*rdv2nuL; 
  outl[9] += incr2[9]*rdvSq4nuL-1.0*incr1[9]*rdv2nuL; 
  outl[10] += incr1[10]*rdv2nuL-1.0*incr2[10]*rdvSq4nuL; 
  outl[11] += -1.0*incr1[11]*rdv2nuL; 
  outl[12] += -1.0*incr1[12]*rdv2nuL; 
  outl[13] += incr1[13]*rdv2nuL-1.0*incr2[13]*rdvSq4nuL; 
  outl[14] += incr1[14]*rdv2nuL-1.0*incr2[14]*rdvSq4nuL; 
  outl[15] += incr2[15]*rdvSq4nuL-1.0*incr1[15]*rdv2nuL; 
  outl[16] += incr2[16]*rdvSq4nuL-1.0*incr1[16]*rdv2nuL; 
  outl[17] += -1.0*incr1[17]*rdv2nuL; 
  outl[18] += -1.0*incr1[18]*rdv2nuL; 
  outl[19] += incr1[19]*rdv2nuL-1.0*incr2[19]*rdvSq4nuL; 

  return std::abs(wl[2]+0.7905694150420947*Uy[2]-0.7071067811865475*Uy[0]); 
} 
