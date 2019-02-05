#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf1x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates. 
  // dxv[2]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[1*2]:     bulk velocity (in 1 directions). 
  // vtSq[2]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nuL = 2.0*nu/dxvl[1]; 
  double rdv2nuR = 2.0*nu/dxvr[1]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double favg[5]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 

  double fjump[5]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(1*fr[4])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[1]+0.7071067811865475*dxvl[1]-1.0*u[0]; 
  alphaDrag[1] = -1.0*u[1]; 

  double Ghat[5]; 
  for(unsigned short int i=0; i<5; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*vtSq[0]*fr[4]-1.897366596101028*vtSq[0]*fl[4]-3.368048396326869*vtSq[1]*fr[3]-3.368048396326869*vtSq[1]*fl[3]-3.368048396326869*vtSq[0]*fr[2]-3.368048396326869*vtSq[0]*fl[2]+2.651650429449552*fr[1]*vtSq[1]-2.651650429449552*fl[1]*vtSq[1]+2.651650429449552*fr[0]*vtSq[0]-2.651650429449552*fl[0]*vtSq[0])*rdv-1.118033988749895*fjump[4]+alphaDrag[0]*(0.7905694150420947*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+alphaDrag[1]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (1.897366596101028*vtSq[1]*fr[4]-1.897366596101028*vtSq[1]*fl[4]-3.368048396326869*vtSq[0]*fr[3]-3.368048396326869*vtSq[0]*fl[3]-3.368048396326869*vtSq[1]*fr[2]-3.368048396326869*vtSq[1]*fl[2]+2.651650429449552*fr[0]*vtSq[1]-2.651650429449552*fl[0]*vtSq[1]+2.651650429449552*vtSq[0]*fr[1]-2.651650429449552*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.7905694150420947*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[3]+alphaDrag[0]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])-0.5*fjump[1]; 

  double incr1[5]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 
  incr1[4] = -1.118033988749895*Ghat[0]; 

  double incr2[5]; 
  incr2[2] = vtSq[0]*(0.2995357736356374*(fr[4]+fl[4])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[1]*((-0.430893194785552*fr[3])+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[3] = vtSq[1]*(0.2995357736356374*(fr[4]+fl[4])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[0]*((-0.430893194785552*fr[3])+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[4] = vtSq[0]*((-1.160097062884178*(fr[4]+fl[4]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0]))+vtSq[1]*(1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[1]+fl[1])); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr2[2]*rdvSq4nuR+incr1[2]*rdv2nuR; 
  outr[3] += incr2[3]*rdvSq4nuR+incr1[3]*rdv2nuR; 
  outr[4] += incr2[4]*rdvSq4nuR+incr1[4]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += incr1[2]*rdv2nuL-1.0*incr2[2]*rdvSq4nuL; 
  outl[3] += incr1[3]*rdv2nuL-1.0*incr2[3]*rdvSq4nuL; 
  outl[4] += incr2[4]*rdvSq4nuL-1.0*incr1[4]*rdv2nuL; 

  return std::abs(wl[1]-0.7071067811865475*u[0]); 
} 
double GkLBOconstNuSurf1x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates. 
  // dxv[2]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[1*3]:     bulk velocity (in 1 directions). 
  // vtSq[3]:    thermal speed squared. 
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

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[1]+0.7071067811865475*dxvl[1]-1.0*u[0]; 
  alphaDrag[1] = -1.0*u[1]; 
  alphaDrag[2] = -1.0*u[2]; 

  double Ghat[8]; 
  for(unsigned short int i=0; i<8; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*vtSq[1]*fr[7]-1.897366596101028*vtSq[1]*fl[7]-3.368048396326869*vtSq[2]*fr[6]-3.368048396326869*vtSq[2]*fl[6]+1.897366596101028*vtSq[0]*fr[5]-1.897366596101028*vtSq[0]*fl[5]+2.651650429449552*vtSq[2]*fr[4]-2.651650429449552*vtSq[2]*fl[4]-3.368048396326869*vtSq[1]*fr[3]-3.368048396326869*vtSq[1]*fl[3]-3.368048396326869*vtSq[0]*fr[2]-3.368048396326869*vtSq[0]*fl[2]+2.651650429449552*fr[1]*vtSq[1]-2.651650429449552*fl[1]*vtSq[1]+2.651650429449552*fr[0]*vtSq[0]-2.651650429449552*fl[0]*vtSq[0])*rdv+alphaDrag[1]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])-1.118033988749895*fjump[5]+alphaDrag[0]*(0.7905694150420947*favg[5]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (1.697056274847714*vtSq[2]*fr[7]+1.897366596101028*vtSq[0]*fr[7]-1.697056274847714*vtSq[2]*fl[7]-1.897366596101028*vtSq[0]*fl[7]-3.012474066278414*vtSq[1]*fr[6]-3.012474066278414*vtSq[1]*fl[6]+1.897366596101028*vtSq[1]*fr[5]-1.897366596101028*vtSq[1]*fl[5]+2.371708245126284*vtSq[1]*fr[4]-2.371708245126284*vtSq[1]*fl[4]-3.012474066278413*vtSq[2]*fr[3]-3.368048396326869*vtSq[0]*fr[3]-3.012474066278413*vtSq[2]*fl[3]-3.368048396326869*vtSq[0]*fl[3]+2.371708245126284*fr[1]*vtSq[2]-2.371708245126284*fl[1]*vtSq[2]-3.368048396326869*vtSq[1]*fr[2]-3.368048396326869*vtSq[1]*fl[2]+2.651650429449552*fr[0]*vtSq[1]-2.651650429449552*fl[0]*vtSq[1]+2.651650429449552*vtSq[0]*fr[1]-2.651650429449552*vtSq[0]*fl[1])*rdv-1.118033988749895*fjump[7]+alphaDrag[0]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])+alphaDrag[1]*(0.5477225575051661*favg[6]+0.7905694150420947*favg[5]+0.3162277660168379*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[4] = (1.697056274847714*vtSq[1]*fr[7]-1.697056274847714*vtSq[1]*fl[7]-2.151767190198866*vtSq[2]*fr[6]-3.368048396326869*vtSq[0]*fr[6]-2.151767190198866*vtSq[2]*fl[6]-3.368048396326869*vtSq[0]*fl[6]+1.897366596101028*vtSq[2]*fr[5]-1.897366596101028*vtSq[2]*fl[5]+1.694077317947346*vtSq[2]*fr[4]+2.651650429449552*vtSq[0]*fr[4]-1.694077317947346*vtSq[2]*fl[4]-2.651650429449552*vtSq[0]*fl[4]-3.012474066278413*vtSq[1]*fr[3]-3.012474066278413*vtSq[1]*fl[3]-3.368048396326869*fr[2]*vtSq[2]-3.368048396326869*fl[2]*vtSq[2]+2.651650429449552*fr[0]*vtSq[2]-2.651650429449552*fl[0]*vtSq[2]+2.371708245126284*fr[1]*vtSq[1]-2.371708245126284*fl[1]*vtSq[1])*rdv+alphaDrag[1]*(0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])-0.8660254037844387*fjump[6]+alphaDrag[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])+alphaDrag[2]*(0.3912303982179757*favg[6]+0.7905694150420947*favg[5]+0.2258769757263128*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[4]; 

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
  incr2[2] = vtSq[1]*(0.2995357736356374*(fr[7]+fl[7])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[2]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[4]+fl[4]))+vtSq[0]*(0.2995357736356374*(fr[5]+fl[5])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[3] = vtSq[0]*(0.2995357736356374*(fr[7]+fl[7])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[2]*(0.2679129406169099*(fr[7]+fl[7])-0.3854025898330209*fr[3]+0.3854025898330209*fl[3]+0.273861278752583*(fr[1]+fl[1]))+vtSq[1]*((-0.3854025898330209*fr[6])+0.3854025898330209*fl[6]+0.2995357736356374*(fr[5]+fl[5])+0.273861278752583*(fr[4]+fl[4])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[5] = vtSq[1]*((-1.160097062884178*(fr[7]+fl[7]))+1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[1]+fl[1]))+vtSq[2]*(1.668842167398551*fr[6]-1.668842167398551*fl[6]-1.185854122563142*(fr[4]+fl[4]))+vtSq[0]*((-1.160097062884178*(fr[5]+fl[5]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0])); 
  incr2[6] = vtSq[1]*(0.2679129406169099*(fr[7]+fl[7])-0.3854025898330209*fr[3]+0.3854025898330209*fl[3]+0.273861278752583*(fr[1]+fl[1]))+vtSq[2]*((-0.2752875641664436*fr[6])+0.2752875641664436*fl[6]+0.2995357736356374*(fr[5]+fl[5])+0.1956151991089878*(fr[4]+fl[4])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[4]+fl[4])); 
  incr2[7] = vtSq[2]*((-1.037622357242749*(fr[7]+fl[7]))+1.492657812008498*fr[3]-1.492657812008498*fl[3]-1.060660171779821*(fr[1]+fl[1]))+vtSq[0]*((-1.160097062884178*(fr[7]+fl[7]))+1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[1]+fl[1]))+vtSq[1]*(1.492657812008498*fr[6]-1.492657812008498*fl[6]-1.160097062884178*(fr[5]+fl[5])-1.060660171779821*(fr[4]+fl[4])+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0])); 

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

  return std::abs(0.7905694150420947*u[2]+wl[1]-0.7071067811865475*u[0]); 
} 
