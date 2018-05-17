#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x1vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
// nu: constant collisionality. 
// vMuMidMax: maximum midpoint value of v-u. 
// u[1*2]: mean flow speed in 1 directions. 
// vtSq[2]: thermal speed squared. 
// fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nul = 2.0*nu/dxvl[1]; 
  double rdv2nur = 2.0*nu/dxvr[1]; 
  double rdvSq4nul = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nur = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double favg[4]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 

  double fjump[4]; 
  fjump[0] = vMuMidMax*(1*fr[0]-fl[0]); 
  fjump[1] = vMuMidMax*(1*fr[1]-fl[1]); 
  fjump[2] = vMuMidMax*(-1*fr[2]-fl[2]); 
  fjump[3] = vMuMidMax*(-1*fr[3]-fl[3]); 

  const double *uX = &u[0]; 
  double drBar[2]; 
  drBar[0] = 1.414213562373095*wl[1]-1.0*uX[0]; 
  drBar[1] = -1.0*uX[1]; 

  double Gdiff[4]; 
  Gdiff[0] = vtSq[1]*((-1.530931089239486*(fr[3]+fl[3]))+1.590990257669731*fr[1]-1.590990257669731*fl[1])+vtSq[0]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0]); 
  Gdiff[1] = vtSq[0]*((-1.530931089239486*(fr[3]+fl[3]))+1.590990257669731*fr[1]-1.590990257669731*fl[1])+vtSq[1]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0]); 

  double Ghat[4]; 
  for(unsigned int i=0; i<4; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+drBar[1]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]+drBar[0]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv-0.8660254037844386*fjump[3]+drBar[0]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[1]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[1]; 
  Ghat[2] = 0.25*dxvl[1]*favg[2]+0.1443375672974065*favg[0]*dxvl[1]; 
  Ghat[3] = 0.25*dxvl[1]*favg[3]+0.1443375672974065*dxvl[1]*favg[1]; 

  double incr1[4]; 
  incr1[0] = (-0.8660254037844386*Ghat[2])-0.5*Ghat[0]; 
  incr1[1] = (-0.8660254037844386*Ghat[3])-0.5*Ghat[1]; 
  incr1[2] = 1.5*Ghat[2]+0.8660254037844386*Ghat[0]; 
  incr1[3] = 1.5*Ghat[3]+0.8660254037844386*Ghat[1]; 

  double incr2[4]; 
  incr2[2] = vtSq[1]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[3] = vtSq[0]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[1]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 

  outr[0] += incr1[0]*rdv2nur; 
  outr[1] += incr1[1]*rdv2nur; 
  outr[2] += incr2[2]*rdvSq4nur+incr1[2]*rdv2nur; 
  outr[3] += incr2[3]*rdvSq4nur+incr1[3]*rdv2nur; 

  outl[0] += -1.0*incr1[0]*rdv2nul; 
  outl[1] += -1.0*incr1[1]*rdv2nul; 
  outl[2] += incr1[2]*rdv2nul-1.0*incr2[2]*rdvSq4nul; 
  outl[3] += incr1[3]*rdv2nul-1.0*incr2[3]*rdvSq4nul; 

  const double vMuMid = wl[1]-0.7071067811865475*uX[0]; 
  return std::abs(vMuMid); 
} 
double VmLBOconstNuSurf1x1vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
// nu: constant collisionality. 
// vMuMidMax: maximum midpoint value of v-u. 
// u[1*3]: mean flow speed in 1 directions. 
// vtSq[3]: thermal speed squared. 
// fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nul = 2.0*nu/dxvl[1]; 
  double rdv2nur = 2.0*nu/dxvr[1]; 
  double rdvSq4nul = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nur = 4.0*nu/(dxvr[1]*dxvr[1]); 

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
  fjump[0] = vMuMidMax*(1*fr[0]-fl[0]); 
  fjump[1] = vMuMidMax*(1*fr[1]-fl[1]); 
  fjump[2] = vMuMidMax*(-1*fr[2]-fl[2]); 
  fjump[3] = vMuMidMax*(-1*fr[3]-fl[3]); 
  fjump[4] = vMuMidMax*(1*fr[4]-fl[4]); 
  fjump[5] = vMuMidMax*(1*fr[5]-fl[5]); 
  fjump[6] = vMuMidMax*(-1*fr[6]-fl[6]); 
  fjump[7] = vMuMidMax*(1*fr[7]-fl[7]); 

  const double *uX = &u[0]; 
  double drBar[3]; 
  drBar[0] = 1.414213562373095*wl[1]-1.0*uX[0]; 
  drBar[1] = -1.0*uX[1]; 
  drBar[2] = -1.0*uX[2]; 

  double Gdiff[8]; 
  Gdiff[0] = vtSq[1]*(1.897366596101028*fr[7]-1.897366596101028*fl[7]-3.368048396326869*(fr[3]+fl[3])+2.651650429449552*fr[1]-2.651650429449552*fl[1])+vtSq[2]*((-3.368048396326869*(fr[6]+fl[6]))+2.651650429449552*fr[4]-2.651650429449552*fl[4])+vtSq[0]*(1.897366596101028*fr[5]-1.897366596101028*fl[5]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0]); 
  Gdiff[1] = vtSq[0]*(1.897366596101028*fr[7]-1.897366596101028*fl[7]-3.368048396326869*(fr[3]+fl[3])+2.651650429449552*fr[1]-2.651650429449552*fl[1])+vtSq[2]*(1.697056274847714*fr[7]-1.697056274847714*fl[7]-3.012474066278413*(fr[3]+fl[3])+2.371708245126284*fr[1]-2.371708245126284*fl[1])+vtSq[1]*((-3.012474066278414*(fr[6]+fl[6]))+1.897366596101028*fr[5]-1.897366596101028*fl[5]+2.371708245126284*fr[4]-2.371708245126284*fl[4]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0]); 
  Gdiff[4] = vtSq[1]*(1.697056274847714*fr[7]-1.697056274847714*fl[7]-3.012474066278413*(fr[3]+fl[3])+2.371708245126284*fr[1]-2.371708245126284*fl[1])+vtSq[2]*((-2.151767190198866*(fr[6]+fl[6]))+1.897366596101028*fr[5]-1.897366596101028*fl[5]+1.694077317947346*fr[4]-1.694077317947346*fl[4]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+vtSq[0]*((-3.368048396326869*(fr[6]+fl[6]))+2.651650429449552*fr[4]-2.651650429449552*fl[4]); 

  double Ghat[8]; 
  for(unsigned int i=0; i<8; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdv+drBar[1]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[2]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])-1.118033988749895*fjump[5]+drBar[0]*(0.7905694150420947*favg[5]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv-1.118033988749895*fjump[7]+drBar[0]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+drBar[2]*(0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])+drBar[1]*(0.5477225575051661*favg[6]+0.7905694150420947*favg[5]+0.3162277660168379*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[2] = 0.3227486121839515*dxvl[1]*favg[5]+0.25*dxvl[1]*favg[2]+0.1443375672974065*favg[0]*dxvl[1]; 
  Ghat[3] = 0.3227486121839514*dxvl[1]*favg[7]+0.25*dxvl[1]*favg[3]+0.1443375672974065*dxvl[1]*favg[1]; 
  Ghat[4] = Gdiff[4]*rdv+drBar[1]*(0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])-0.8660254037844387*fjump[6]+drBar[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])+drBar[2]*(0.3912303982179757*favg[6]+0.7905694150420947*favg[5]+0.2258769757263128*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[4]; 
  Ghat[6] = 0.25*dxvl[1]*favg[6]+0.1443375672974064*dxvl[1]*favg[4]; 

  double incr1[8]; 
  incr1[0] = (-0.8660254037844386*Ghat[2])-0.5*Ghat[0]; 
  incr1[1] = (-0.8660254037844386*Ghat[3])-0.5*Ghat[1]; 
  incr1[2] = 1.5*Ghat[2]+0.8660254037844386*Ghat[0]; 
  incr1[3] = 1.5*Ghat[3]+0.8660254037844386*Ghat[1]; 
  incr1[4] = (-0.8660254037844387*Ghat[6])-0.5*Ghat[4]; 
  incr1[5] = (-1.936491673103709*Ghat[2])-1.118033988749895*Ghat[0]; 
  incr1[6] = 1.5*Ghat[6]+0.8660254037844387*Ghat[4]; 
  incr1[7] = (-1.936491673103709*Ghat[3])-1.118033988749895*Ghat[1]; 

  double incr2[8]; 
  incr2[2] = vtSq[1]*(0.2995357736356374*(fr[7]+fl[7])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[2]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[4]+fl[4]))+vtSq[0]*(0.2995357736356374*(fr[5]+fl[5])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[3] = vtSq[0]*(0.2995357736356374*(fr[7]+fl[7])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[2]*(0.2679129406169099*(fr[7]+fl[7])-0.3854025898330209*fr[3]+0.3854025898330209*fl[3]+0.273861278752583*(fr[1]+fl[1]))+vtSq[1]*((-0.3854025898330209*fr[6])+0.3854025898330209*fl[6]+0.2995357736356374*(fr[5]+fl[5])+0.273861278752583*(fr[4]+fl[4])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[5] = vtSq[1]*((-1.160097062884178*(fr[7]+fl[7]))+1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[1]+fl[1]))+vtSq[2]*(1.668842167398552*fr[6]-1.668842167398552*fl[6]-1.185854122563142*(fr[4]+fl[4]))+vtSq[0]*((-1.160097062884178*(fr[5]+fl[5]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0])); 
  incr2[6] = vtSq[1]*(0.2679129406169099*(fr[7]+fl[7])-0.3854025898330209*fr[3]+0.3854025898330209*fl[3]+0.273861278752583*(fr[1]+fl[1]))+vtSq[2]*((-0.2752875641664436*fr[6])+0.2752875641664436*fl[6]+0.2995357736356374*(fr[5]+fl[5])+0.1956151991089878*(fr[4]+fl[4])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[4]+fl[4])); 
  incr2[7] = vtSq[2]*((-1.037622357242749*(fr[7]+fl[7]))+1.492657812008498*fr[3]-1.492657812008498*fl[3]-1.060660171779821*(fr[1]+fl[1]))+vtSq[0]*((-1.160097062884178*(fr[7]+fl[7]))+1.668842167398552*fr[3]-1.668842167398552*fl[3]-1.185854122563142*(fr[1]+fl[1]))+vtSq[1]*(1.492657812008498*fr[6]-1.492657812008498*fl[6]-1.160097062884178*(fr[5]+fl[5])-1.060660171779821*(fr[4]+fl[4])+1.668842167398552*fr[2]-1.668842167398552*fl[2]-1.185854122563142*(fr[0]+fl[0])); 

  outr[0] += incr1[0]*rdv2nur; 
  outr[1] += incr1[1]*rdv2nur; 
  outr[2] += incr2[2]*rdvSq4nur+incr1[2]*rdv2nur; 
  outr[3] += incr2[3]*rdvSq4nur+incr1[3]*rdv2nur; 
  outr[4] += incr1[4]*rdv2nur; 
  outr[5] += incr2[5]*rdvSq4nur+incr1[5]*rdv2nur; 
  outr[6] += incr2[6]*rdvSq4nur+incr1[6]*rdv2nur; 
  outr[7] += incr2[7]*rdvSq4nur+incr1[7]*rdv2nur; 

  outl[0] += -1.0*incr1[0]*rdv2nul; 
  outl[1] += -1.0*incr1[1]*rdv2nul; 
  outl[2] += incr1[2]*rdv2nul-1.0*incr2[2]*rdvSq4nul; 
  outl[3] += incr1[3]*rdv2nul-1.0*incr2[3]*rdvSq4nul; 
  outl[4] += -1.0*incr1[4]*rdv2nul; 
  outl[5] += incr2[5]*rdvSq4nul-1.0*incr1[5]*rdv2nul; 
  outl[6] += incr1[6]*rdv2nul-1.0*incr2[6]*rdvSq4nul; 
  outl[7] += incr2[7]*rdvSq4nul-1.0*incr1[7]*rdv2nul; 

  const double vMuMid = 0.7905694150420947*uX[2]+wl[1]-0.7071067811865475*uX[0]; 
  return std::abs(vMuMid); 
} 
