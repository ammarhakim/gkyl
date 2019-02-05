#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf1x2vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:       Cell-center coordinates. 
  // dxv[3]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[2*2]:     bulk velocity (in 2 directions). 
  // vtSq[2]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nuL = 2.0*nu/dxvl[1]; 
  double rdv2nuR = 2.0*nu/dxvr[1]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double favg[9]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 

  double fjump[9]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = vMuMidMax*(fl[7]-(-1*fr[7])); 
  fjump[8] = vMuMidMax*(fl[8]-(1*fr[8])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[1]+0.7071067811865475*dxvl[1]-1.0*u[0]; 
  alphaDrag[1] = -1.0*u[1]; 

  double Ghat[9]; 
  for(unsigned short int i=0; i<9; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*vtSq[0]*fr[8]-1.897366596101028*vtSq[0]*fl[8]-3.368048396326869*vtSq[1]*fr[4]-3.368048396326869*vtSq[1]*fl[4]-3.368048396326869*vtSq[0]*fr[2]-3.368048396326869*vtSq[0]*fl[2]+2.651650429449552*fr[1]*vtSq[1]-2.651650429449552*fl[1]*vtSq[1]+2.651650429449552*fr[0]*vtSq[0]-2.651650429449552*fl[0]*vtSq[0])*rdv-1.118033988749895*fjump[8]+alphaDrag[0]*(0.7905694150420947*favg[8]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+alphaDrag[1]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (1.897366596101028*vtSq[1]*fr[8]-1.897366596101028*vtSq[1]*fl[8]-3.368048396326869*vtSq[0]*fr[4]-3.368048396326869*vtSq[0]*fl[4]-3.368048396326869*vtSq[1]*fr[2]-3.368048396326869*vtSq[1]*fl[2]+2.651650429449552*fr[0]*vtSq[1]-2.651650429449552*fl[0]*vtSq[1]+2.651650429449552*vtSq[0]*fr[1]-2.651650429449552*vtSq[0]*fl[1])*rdv+alphaDrag[1]*(0.7905694150420947*favg[8]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[4]+alphaDrag[0]*(0.6123724356957944*favg[4]+0.3535533905932737*favg[1])-0.5*fjump[1]; 
  Ghat[3] = ((-3.368048396326869*vtSq[1]*fr[7])-3.368048396326869*vtSq[1]*fl[7]-3.368048396326869*vtSq[0]*fr[6]-3.368048396326869*vtSq[0]*fl[6]+2.651650429449552*vtSq[1]*fr[5]-2.651650429449552*vtSq[1]*fl[5]+2.651650429449552*vtSq[0]*fr[3]-2.651650429449552*vtSq[0]*fl[3])*rdv+alphaDrag[1]*(0.6123724356957944*favg[7]+0.3535533905932737*favg[5])-0.8660254037844386*fjump[6]+alphaDrag[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[3]; 
  Ghat[5] = ((-3.368048396326869*vtSq[0]*fr[7])-3.368048396326869*vtSq[0]*fl[7]-3.368048396326869*vtSq[1]*fr[6]-3.368048396326869*vtSq[1]*fl[6]+2.651650429449552*vtSq[0]*fr[5]-2.651650429449552*vtSq[0]*fl[5]+2.651650429449552*vtSq[1]*fr[3]-2.651650429449552*vtSq[1]*fl[3])*rdv-0.8660254037844386*fjump[7]+alphaDrag[0]*(0.6123724356957944*favg[7]+0.3535533905932737*favg[5])+alphaDrag[1]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[5]; 

  double incr1[9]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[1]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[3]; 
  incr1[7] = 0.8660254037844386*Ghat[5]; 
  incr1[8] = -1.118033988749895*Ghat[0]; 

  double incr2[9]; 
  incr2[2] = vtSq[0]*(0.2995357736356374*(fr[8]+fl[8])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[1]*((-0.430893194785552*fr[4])+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[4] = vtSq[1]*(0.2995357736356374*(fr[8]+fl[8])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[0]*((-0.430893194785552*fr[4])+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[6] = vtSq[1]*((-0.430893194785552*fr[7])+0.430893194785552*fl[7]+0.3061862178478971*(fr[5]+fl[5]))+vtSq[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 
  incr2[7] = vtSq[0]*((-0.430893194785552*fr[7])+0.430893194785552*fl[7]+0.3061862178478971*(fr[5]+fl[5]))+vtSq[1]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 
  incr2[8] = vtSq[0]*((-1.160097062884178*(fr[8]+fl[8]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0]))+vtSq[1]*(1.668842167398551*fr[4]-1.668842167398551*fl[4]-1.185854122563142*(fr[1]+fl[1])); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr2[2]*rdvSq4nuR+incr1[2]*rdv2nuR; 
  outr[3] += incr1[3]*rdv2nuR; 
  outr[4] += incr2[4]*rdvSq4nuR+incr1[4]*rdv2nuR; 
  outr[5] += incr1[5]*rdv2nuR; 
  outr[6] += incr2[6]*rdvSq4nuR+incr1[6]*rdv2nuR; 
  outr[7] += incr2[7]*rdvSq4nuR+incr1[7]*rdv2nuR; 
  outr[8] += incr2[8]*rdvSq4nuR+incr1[8]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += incr1[2]*rdv2nuL-1.0*incr2[2]*rdvSq4nuL; 
  outl[3] += -1.0*incr1[3]*rdv2nuL; 
  outl[4] += incr1[4]*rdv2nuL-1.0*incr2[4]*rdvSq4nuL; 
  outl[5] += -1.0*incr1[5]*rdv2nuL; 
  outl[6] += incr1[6]*rdv2nuL-1.0*incr2[6]*rdvSq4nuL; 
  outl[7] += incr1[7]*rdv2nuL-1.0*incr2[7]*rdvSq4nuL; 
  outl[8] += incr2[8]*rdvSq4nuL-1.0*incr1[8]*rdv2nuL; 

  return std::abs(wl[1]-0.7071067811865475*u[0]); 
} 
double GkLBOconstNuSurf1x2vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:       Cell-center coordinates. 
  // dxv[3]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[2*3]:     bulk velocity (in 2 directions). 
  // vtSq[3]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2nuL = 2.0*nu/dxvl[1]; 
  double rdv2nuR = 2.0*nu/dxvr[1]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

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
  fjump[17] = vMuMidMax*(fl[17]-(-1*fr[17])); 
  fjump[18] = vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = vMuMidMax*(fl[19]-(-1*fr[19])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[1]+0.7071067811865475*dxvl[1]-1.0*u[0]; 
  alphaDrag[1] = -1.0*u[1]; 
  alphaDrag[2] = -1.0*u[2]; 

  double Ghat[20]; 
  for(unsigned short int i=0; i<20; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*vtSq[1]*fr[12]-1.897366596101028*vtSq[1]*fl[12]-3.368048396326869*vtSq[2]*fr[11]-3.368048396326869*vtSq[2]*fl[11]+1.897366596101028*vtSq[0]*fr[8]-1.897366596101028*vtSq[0]*fl[8]+2.651650429449552*vtSq[2]*fr[7]-2.651650429449552*vtSq[2]*fl[7]-3.368048396326869*vtSq[1]*fr[4]-3.368048396326869*vtSq[1]*fl[4]-3.368048396326869*vtSq[0]*fr[2]-3.368048396326869*vtSq[0]*fl[2]+2.651650429449552*fr[1]*vtSq[1]-2.651650429449552*fl[1]*vtSq[1]+2.651650429449552*fr[0]*vtSq[0]-2.651650429449552*fl[0]*vtSq[0])*rdv+alphaDrag[1]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.6123724356957944*favg[11]+0.3535533905932737*favg[7])-1.118033988749895*fjump[8]+alphaDrag[0]*(0.7905694150420947*favg[8]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (1.697056274847714*vtSq[2]*fr[12]+1.897366596101028*vtSq[0]*fr[12]-1.697056274847714*vtSq[2]*fl[12]-1.897366596101028*vtSq[0]*fl[12]-3.012474066278414*vtSq[1]*fr[11]-3.012474066278414*vtSq[1]*fl[11]+1.897366596101028*vtSq[1]*fr[8]-1.897366596101028*vtSq[1]*fl[8]+2.371708245126284*vtSq[1]*fr[7]-2.371708245126284*vtSq[1]*fl[7]-3.012474066278413*vtSq[2]*fr[4]-3.368048396326869*vtSq[0]*fr[4]-3.012474066278413*vtSq[2]*fl[4]-3.368048396326869*vtSq[0]*fl[4]+2.371708245126284*fr[1]*vtSq[2]-2.371708245126284*fl[1]*vtSq[2]-3.368048396326869*vtSq[1]*fr[2]-3.368048396326869*vtSq[1]*fl[2]+2.651650429449552*fr[0]*vtSq[1]-2.651650429449552*fl[0]*vtSq[1]+2.651650429449552*vtSq[0]*fr[1]-2.651650429449552*vtSq[0]*fl[1])*rdv-1.118033988749895*fjump[12]+alphaDrag[0]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])+alphaDrag[1]*(0.5477225575051661*favg[11]+0.7905694150420947*favg[8]+0.3162277660168379*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[4]-0.5*fjump[1]; 
  Ghat[3] = (1.897366596101028*vtSq[1]*fr[18]-1.897366596101028*vtSq[1]*fl[18]-3.368048396326869*vtSq[2]*fr[17]-3.368048396326869*vtSq[2]*fl[17]+1.897366596101028*vtSq[0]*fr[14]-1.897366596101028*vtSq[0]*fl[14]+2.651650429449552*vtSq[2]*fr[13]-2.651650429449552*vtSq[2]*fl[13]-3.368048396326869*vtSq[1]*fr[10]-3.368048396326869*vtSq[1]*fl[10]-3.368048396326869*vtSq[0]*fr[6]-3.368048396326869*vtSq[0]*fl[6]+2.651650429449552*vtSq[1]*fr[5]-2.651650429449552*vtSq[1]*fl[5]+2.651650429449552*vtSq[0]*fr[3]-2.651650429449552*vtSq[0]*fl[3])*rdv+alphaDrag[1]*(0.7905694150420947*favg[18]+0.6123724356957944*favg[10]+0.3535533905932737*favg[5])+alphaDrag[2]*(0.6123724356957944*favg[17]+0.3535533905932737*favg[13])-1.118033988749895*fjump[14]+alphaDrag[0]*(0.7905694150420948*favg[14]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.8660254037844386*fjump[6]-0.5*fjump[3]; 
  Ghat[5] = (1.697056274847714*vtSq[2]*fr[18]+1.897366596101028*vtSq[0]*fr[18]-1.697056274847714*vtSq[2]*fl[18]-1.897366596101028*vtSq[0]*fl[18]-3.012474066278413*vtSq[1]*fr[17]-3.012474066278413*vtSq[1]*fl[17]+1.897366596101028*vtSq[1]*fr[14]-1.897366596101028*vtSq[1]*fl[14]+2.371708245126284*vtSq[1]*fr[13]-2.371708245126284*vtSq[1]*fl[13]-3.012474066278413*vtSq[2]*fr[10]-3.368048396326869*vtSq[0]*fr[10]-3.012474066278413*vtSq[2]*fl[10]-3.368048396326869*vtSq[0]*fl[10]-3.368048396326869*vtSq[1]*fr[6]-3.368048396326869*vtSq[1]*fl[6]+2.371708245126284*vtSq[2]*fr[5]+2.651650429449552*vtSq[0]*fr[5]-2.371708245126284*vtSq[2]*fl[5]-2.651650429449552*vtSq[0]*fl[5]+2.651650429449552*vtSq[1]*fr[3]-2.651650429449552*vtSq[1]*fl[3])*rdv-1.118033988749895*fjump[18]+alphaDrag[0]*(0.7905694150420947*favg[18]+0.6123724356957944*favg[10]+0.3535533905932737*favg[5])+alphaDrag[2]*(0.7071067811865475*favg[18]+0.5477225575051661*favg[10]+0.3162277660168379*favg[5])+alphaDrag[1]*(0.5477225575051661*favg[17]+0.7905694150420948*favg[14]+0.3162277660168379*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.8660254037844386*fjump[10]-0.5*fjump[5]; 
  Ghat[7] = (1.697056274847714*vtSq[1]*fr[12]-1.697056274847714*vtSq[1]*fl[12]-2.151767190198866*vtSq[2]*fr[11]-3.368048396326869*vtSq[0]*fr[11]-2.151767190198866*vtSq[2]*fl[11]-3.368048396326869*vtSq[0]*fl[11]+1.897366596101028*vtSq[2]*fr[8]-1.897366596101028*vtSq[2]*fl[8]+1.694077317947346*vtSq[2]*fr[7]+2.651650429449552*vtSq[0]*fr[7]-1.694077317947346*vtSq[2]*fl[7]-2.651650429449552*vtSq[0]*fl[7]-3.012474066278413*vtSq[1]*fr[4]-3.012474066278413*vtSq[1]*fl[4]-3.368048396326869*fr[2]*vtSq[2]-3.368048396326869*fl[2]*vtSq[2]+2.651650429449552*fr[0]*vtSq[2]-2.651650429449552*fl[0]*vtSq[2]+2.371708245126284*fr[1]*vtSq[1]-2.371708245126284*fl[1]*vtSq[1])*rdv+alphaDrag[1]*(0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])-0.8660254037844387*fjump[11]+alphaDrag[0]*(0.6123724356957944*favg[11]+0.3535533905932737*favg[7])+alphaDrag[2]*(0.3912303982179757*favg[11]+0.7905694150420947*favg[8]+0.2258769757263128*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[7]; 
  Ghat[9] = ((-3.368048396326869*vtSq[1]*fr[19])-3.368048396326869*vtSq[1]*fl[19]-3.368048396326869*vtSq[0]*fr[16]-3.368048396326869*vtSq[0]*fl[16]+2.651650429449552*vtSq[1]*fr[15]-2.651650429449552*vtSq[1]*fl[15]+2.651650429449552*vtSq[0]*fr[9]-2.651650429449552*vtSq[0]*fl[9])*rdv+alphaDrag[1]*(0.6123724356957944*favg[19]+0.3535533905932737*favg[15])-0.8660254037844387*fjump[16]+alphaDrag[0]*(0.6123724356957944*favg[16]+0.3535533905932737*favg[9])-0.5*fjump[9]; 
  Ghat[13] = (1.697056274847714*vtSq[1]*fr[18]-1.697056274847714*vtSq[1]*fl[18]-2.151767190198866*vtSq[2]*fr[17]-3.368048396326869*vtSq[0]*fr[17]-2.151767190198866*vtSq[2]*fl[17]-3.368048396326869*vtSq[0]*fl[17]+1.897366596101028*vtSq[2]*fr[14]-1.897366596101028*vtSq[2]*fl[14]+1.694077317947346*vtSq[2]*fr[13]+2.651650429449552*vtSq[0]*fr[13]-1.694077317947346*vtSq[2]*fl[13]-2.651650429449552*vtSq[0]*fl[13]-3.012474066278414*vtSq[1]*fr[10]-3.012474066278414*vtSq[1]*fl[10]-3.368048396326869*vtSq[2]*fr[6]-3.368048396326869*vtSq[2]*fl[6]+2.371708245126284*vtSq[1]*fr[5]-2.371708245126284*vtSq[1]*fl[5]+2.651650429449552*vtSq[2]*fr[3]-2.651650429449552*vtSq[2]*fl[3])*rdv+alphaDrag[1]*(0.7071067811865475*favg[18]+0.5477225575051661*favg[10]+0.3162277660168379*favg[5])-0.8660254037844387*fjump[17]+alphaDrag[0]*(0.6123724356957944*favg[17]+0.3535533905932737*favg[13])+alphaDrag[2]*(0.3912303982179757*favg[17]+0.7905694150420947*favg[14]+0.2258769757263128*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.5*fjump[13]; 
  Ghat[15] = ((-3.012474066278414*vtSq[2]*fr[19])-3.368048396326869*vtSq[0]*fr[19]-3.012474066278414*vtSq[2]*fl[19]-3.368048396326869*vtSq[0]*fl[19]-3.368048396326869*vtSq[1]*fr[16]-3.368048396326869*vtSq[1]*fl[16]+2.371708245126284*vtSq[2]*fr[15]+2.651650429449552*vtSq[0]*fr[15]-2.371708245126284*vtSq[2]*fl[15]-2.651650429449552*vtSq[0]*fl[15]+2.651650429449552*vtSq[1]*fr[9]-2.651650429449552*vtSq[1]*fl[9])*rdv-0.8660254037844387*fjump[19]+alphaDrag[0]*(0.6123724356957944*favg[19]+0.3535533905932737*favg[15])+alphaDrag[2]*(0.5477225575051661*favg[19]+0.3162277660168379*favg[15])+alphaDrag[1]*(0.6123724356957944*favg[16]+0.3535533905932737*favg[9])-0.5*fjump[15]; 

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
  incr2[2] = vtSq[1]*(0.2995357736356374*(fr[12]+fl[12])-0.430893194785552*fr[4]+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[2]*((-0.430893194785552*fr[11])+0.430893194785552*fl[11]+0.3061862178478971*(fr[7]+fl[7]))+vtSq[0]*(0.2995357736356374*(fr[8]+fl[8])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[4] = vtSq[0]*(0.2995357736356374*(fr[12]+fl[12])-0.430893194785552*fr[4]+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+vtSq[2]*(0.2679129406169099*(fr[12]+fl[12])-0.3854025898330209*fr[4]+0.3854025898330209*fl[4]+0.273861278752583*(fr[1]+fl[1]))+vtSq[1]*((-0.3854025898330209*fr[11])+0.3854025898330209*fl[11]+0.2995357736356374*(fr[8]+fl[8])+0.273861278752583*(fr[7]+fl[7])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[6] = vtSq[1]*(0.2995357736356374*(fr[18]+fl[18])-0.430893194785552*fr[10]+0.430893194785552*fl[10]+0.3061862178478971*(fr[5]+fl[5]))+vtSq[2]*((-0.430893194785552*fr[17])+0.430893194785552*fl[17]+0.3061862178478971*(fr[13]+fl[13]))+vtSq[0]*(0.2995357736356374*(fr[14]+fl[14])-0.430893194785552*fr[6]+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 
  incr2[8] = vtSq[1]*((-1.160097062884178*(fr[12]+fl[12]))+1.668842167398551*fr[4]-1.668842167398551*fl[4]-1.185854122563142*(fr[1]+fl[1]))+vtSq[2]*(1.668842167398551*fr[11]-1.668842167398551*fl[11]-1.185854122563142*(fr[7]+fl[7]))+vtSq[0]*((-1.160097062884178*(fr[8]+fl[8]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0])); 
  incr2[10] = vtSq[0]*(0.2995357736356374*(fr[18]+fl[18])-0.430893194785552*fr[10]+0.430893194785552*fl[10]+0.3061862178478971*(fr[5]+fl[5]))+vtSq[2]*(0.2679129406169099*(fr[18]+fl[18])-0.3854025898330209*fr[10]+0.3854025898330209*fl[10]+0.273861278752583*(fr[5]+fl[5]))+vtSq[1]*((-0.3854025898330209*fr[17])+0.3854025898330209*fl[17]+0.2995357736356374*(fr[14]+fl[14])+0.273861278752583*(fr[13]+fl[13])-0.430893194785552*fr[6]+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 
  incr2[11] = vtSq[1]*(0.2679129406169099*(fr[12]+fl[12])-0.3854025898330209*fr[4]+0.3854025898330209*fl[4]+0.273861278752583*(fr[1]+fl[1]))+vtSq[2]*((-0.2752875641664436*fr[11])+0.2752875641664436*fl[11]+0.2995357736356374*(fr[8]+fl[8])+0.1956151991089878*(fr[7]+fl[7])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+vtSq[0]*((-0.430893194785552*fr[11])+0.430893194785552*fl[11]+0.3061862178478971*(fr[7]+fl[7])); 
  incr2[12] = vtSq[2]*((-1.037622357242749*(fr[12]+fl[12]))+1.492657812008498*fr[4]-1.492657812008498*fl[4]-1.060660171779821*(fr[1]+fl[1]))+vtSq[0]*((-1.160097062884178*(fr[12]+fl[12]))+1.668842167398551*fr[4]-1.668842167398551*fl[4]-1.185854122563142*(fr[1]+fl[1]))+vtSq[1]*(1.492657812008498*fr[11]-1.492657812008498*fl[11]-1.160097062884178*(fr[8]+fl[8])-1.060660171779821*(fr[7]+fl[7])+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0])); 
  incr2[14] = vtSq[1]*((-1.160097062884178*(fr[18]+fl[18]))+1.668842167398551*fr[10]-1.668842167398551*fl[10]-1.185854122563142*(fr[5]+fl[5]))+vtSq[2]*(1.668842167398551*fr[17]-1.668842167398551*fl[17]-1.185854122563142*(fr[13]+fl[13]))+vtSq[0]*((-1.160097062884178*(fr[14]+fl[14]))+1.668842167398551*fr[6]-1.668842167398551*fl[6]-1.185854122563142*(fr[3]+fl[3])); 
  incr2[16] = vtSq[1]*((-0.430893194785552*fr[19])+0.430893194785552*fl[19]+0.3061862178478971*(fr[15]+fl[15]))+vtSq[0]*((-0.430893194785552*fr[16])+0.430893194785552*fl[16]+0.3061862178478971*(fr[9]+fl[9])); 
  incr2[17] = vtSq[1]*(0.2679129406169099*(fr[18]+fl[18])-0.3854025898330209*fr[10]+0.3854025898330209*fl[10]+0.273861278752583*(fr[5]+fl[5]))+vtSq[2]*((-0.2752875641664436*fr[17])+0.2752875641664436*fl[17]+0.2995357736356374*(fr[14]+fl[14])+0.1956151991089878*(fr[13]+fl[13])-0.430893194785552*fr[6]+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3]))+vtSq[0]*((-0.430893194785552*fr[17])+0.430893194785552*fl[17]+0.3061862178478971*(fr[13]+fl[13])); 
  incr2[18] = vtSq[2]*((-1.037622357242749*(fr[18]+fl[18]))+1.492657812008498*fr[10]-1.492657812008498*fl[10]-1.060660171779821*(fr[5]+fl[5]))+vtSq[0]*((-1.160097062884178*(fr[18]+fl[18]))+1.668842167398551*fr[10]-1.668842167398551*fl[10]-1.185854122563142*(fr[5]+fl[5]))+vtSq[1]*(1.492657812008498*fr[17]-1.492657812008498*fl[17]-1.160097062884178*(fr[14]+fl[14])-1.060660171779821*(fr[13]+fl[13])+1.668842167398551*fr[6]-1.668842167398551*fl[6]-1.185854122563142*(fr[3]+fl[3])); 
  incr2[19] = vtSq[2]*((-0.3854025898330209*fr[19])+0.3854025898330209*fl[19]+0.273861278752583*(fr[15]+fl[15]))+vtSq[0]*((-0.430893194785552*fr[19])+0.430893194785552*fl[19]+0.3061862178478971*(fr[15]+fl[15]))+vtSq[1]*((-0.430893194785552*fr[16])+0.430893194785552*fl[16]+0.3061862178478971*(fr[9]+fl[9])); 

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
  outr[17] += incr2[17]*rdvSq4nuR+incr1[17]*rdv2nuR; 
  outr[18] += incr2[18]*rdvSq4nuR+incr1[18]*rdv2nuR; 
  outr[19] += incr2[19]*rdvSq4nuR+incr1[19]*rdv2nuR; 

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
  outl[17] += incr1[17]*rdv2nuL-1.0*incr2[17]*rdvSq4nuL; 
  outl[18] += incr2[18]*rdvSq4nuL-1.0*incr1[18]*rdv2nuL; 
  outl[19] += incr1[19]*rdv2nuL-1.0*incr2[19]*rdvSq4nuL; 

  return std::abs(0.7905694150420947*u[2]+wl[1]-0.7071067811865475*u[0]); 
} 
double GkLBOconstNuSurf1x2vSer_Mu_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:       Cell-center coordinates. 
  // dxv[3]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[2*2]:     bulk velocity (in 2 directions). 
  // vtSq[2]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2nuL = 2.0*nu/dxvl[2]; 
  double rdv2nuR = 2.0*nu/dxvr[2]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[2]*dxvl[2]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[2]*dxvr[2]); 

  double favg[9]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 

  double fjump[9]; 
  fjump[0] = vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = vMuMidMax*(fl[7]-(-1*fr[7])); 
  fjump[8] = vMuMidMax*(fl[8]-(1*fr[8])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 2.828427124746191*wl[2]+1.414213562373095*dxvl[2]; 

  double diffFac[2]; 
  diffFac[0] = (BmagInv[1]*vtSq[1]+BmagInv[0]*vtSq[0])*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])*m_; 
  diffFac[1] = (BmagInv[0]*vtSq[1]+vtSq[0]*BmagInv[1])*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])*m_; 

  double Ghat[9]; 
  for(unsigned short int i=0; i<9; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.530931089239486*diffFac[1]*fr[5])-1.530931089239486*diffFac[1]*fl[5]-1.530931089239486*diffFac[0]*fr[3]-1.530931089239486*diffFac[0]*fl[3]+1.590990257669731*diffFac[1]*fr[1]-1.590990257669731*diffFac[1]*fl[1]+1.590990257669731*diffFac[0]*fr[0]-1.590990257669731*diffFac[0]*fl[0])*rdv-0.8660254037844386*fjump[3]+alphaDrag[0]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.5*fjump[0]; 
  Ghat[1] = ((-1.530931089239486*diffFac[0]*fr[5])-1.530931089239486*diffFac[0]*fl[5]-1.530931089239486*diffFac[1]*fr[3]-1.530931089239486*diffFac[1]*fl[3]+1.590990257669731*diffFac[0]*fr[1]-1.590990257669731*diffFac[0]*fl[1]+1.590990257669731*fr[0]*diffFac[1]-1.590990257669731*fl[0]*diffFac[1])*rdv-0.8660254037844386*fjump[5]+alphaDrag[0]*(0.6123724356957944*favg[5]+0.3535533905932737*favg[1])-0.5*fjump[1]; 
  Ghat[2] = ((-1.530931089239486*diffFac[1]*fr[7])-1.530931089239486*diffFac[1]*fl[7]-1.530931089239486*diffFac[0]*fr[6]-1.530931089239486*diffFac[0]*fl[6]+1.590990257669731*diffFac[1]*fr[4]-1.590990257669731*diffFac[1]*fl[4]+1.590990257669731*diffFac[0]*fr[2]-1.590990257669731*diffFac[0]*fl[2])*rdv-0.8660254037844386*fjump[6]+alphaDrag[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-0.5*fjump[2]; 
  Ghat[4] = ((-1.530931089239486*diffFac[0]*fr[7])-1.530931089239486*diffFac[0]*fl[7]-1.530931089239486*diffFac[1]*fr[6]-1.530931089239486*diffFac[1]*fl[6]+1.590990257669731*diffFac[0]*fr[4]-1.590990257669731*diffFac[0]*fl[4]+1.590990257669731*diffFac[1]*fr[2]-1.590990257669731*diffFac[1]*fl[2])*rdv-0.8660254037844386*fjump[7]+alphaDrag[0]*(0.6123724356957944*favg[7]+0.3535533905932737*favg[4])-0.5*fjump[4]; 
  Ghat[8] = (1.590990257669731*diffFac[0]*fr[8]-1.590990257669731*diffFac[0]*fl[8])*rdv-0.5*fjump[8]+0.3535533905932737*alphaDrag[0]*favg[8]; 

  double incr1[9]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = 0.8660254037844386*Ghat[2]; 
  incr1[7] = 0.8660254037844386*Ghat[4]; 
  incr1[8] = -0.5*Ghat[8]; 

  double incr2[9]; 
  incr2[3] = diffFac[1]*(0.3535533905932737*fl[5]-0.3535533905932737*fr[5])+diffFac[0]*(0.3535533905932737*fl[3]-0.3535533905932737*fr[3])+0.3061862178478971*(diffFac[1]*(fr[1]+fl[1])+diffFac[0]*(fr[0]+fl[0])); 
  incr2[5] = diffFac[0]*(0.3535533905932737*fl[5]-0.3535533905932737*fr[5])+diffFac[1]*(0.3535533905932737*fl[3]-0.3535533905932737*fr[3])+0.3061862178478971*(diffFac[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*diffFac[1]); 
  incr2[6] = diffFac[1]*(0.3535533905932737*fl[7]-0.3535533905932737*fr[7])+diffFac[0]*(0.3535533905932737*fl[6]-0.3535533905932737*fr[6])+0.3061862178478971*(diffFac[1]*(fr[4]+fl[4])+diffFac[0]*(fr[2]+fl[2])); 
  incr2[7] = diffFac[0]*(0.3535533905932737*fl[7]-0.3535533905932737*fr[7])+diffFac[1]*(0.3535533905932737*fl[6]-0.3535533905932737*fr[6])+0.3061862178478971*(diffFac[0]*(fr[4]+fl[4])+diffFac[1]*(fr[2]+fl[2])); 

  outr[0] += incr1[0]*rdv2nuR; 
  outr[1] += incr1[1]*rdv2nuR; 
  outr[2] += incr1[2]*rdv2nuR; 
  outr[3] += incr2[3]*rdvSq4nuR+incr1[3]*rdv2nuR; 
  outr[4] += incr1[4]*rdv2nuR; 
  outr[5] += incr2[5]*rdvSq4nuR+incr1[5]*rdv2nuR; 
  outr[6] += incr2[6]*rdvSq4nuR+incr1[6]*rdv2nuR; 
  outr[7] += incr2[7]*rdvSq4nuR+incr1[7]*rdv2nuR; 
  outr[8] += incr1[8]*rdv2nuR; 

  outl[0] += -1.0*incr1[0]*rdv2nuL; 
  outl[1] += -1.0*incr1[1]*rdv2nuL; 
  outl[2] += -1.0*incr1[2]*rdv2nuL; 
  outl[3] += incr1[3]*rdv2nuL-1.0*incr2[3]*rdvSq4nuL; 
  outl[4] += -1.0*incr1[4]*rdv2nuL; 
  outl[5] += incr1[5]*rdv2nuL-1.0*incr2[5]*rdvSq4nuL; 
  outl[6] += incr1[6]*rdv2nuL-1.0*incr2[6]*rdvSq4nuL; 
  outl[7] += incr1[7]*rdv2nuL-1.0*incr2[7]*rdvSq4nuL; 
  outl[8] += -1.0*incr1[8]*rdv2nuL; 

  return std::abs(2.0*wl[2]); 
} 
double GkLBOconstNuSurf1x2vSer_Mu_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:       Cell-center coordinates. 
  // dxv[3]:     Cell spacing. 
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[2*3]:     bulk velocity (in 2 directions). 
  // vtSq[3]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2nuL = 2.0*nu/dxvl[2]; 
  double rdv2nuR = 2.0*nu/dxvr[2]; 
  double rdvSq4nuL = 4.0*nu/(dxvl[2]*dxvl[2]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[2]*dxvr[2]); 

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
  fjump[17] = vMuMidMax*(fl[17]-(-1*fr[17])); 
  fjump[18] = vMuMidMax*(fl[18]-(-1*fr[18])); 
  fjump[19] = vMuMidMax*(fl[19]-(1*fr[19])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 2.828427124746191*wl[2]+1.414213562373095*dxvl[2]; 

  double diffFac[3]; 
  diffFac[0] = (BmagInv[2]*vtSq[2]+BmagInv[1]*vtSq[1]+BmagInv[0]*vtSq[0])*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])*m_; 
  diffFac[1] = (vtSq[1]*((1.264911064067352*BmagInv[2]+1.414213562373095*BmagInv[0])*wl[2]+(0.6324555320336759*BmagInv[2]+0.7071067811865475*BmagInv[0])*dxvl[2])+BmagInv[1]*(vtSq[0]*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])+vtSq[2]*(1.264911064067352*wl[2]+0.6324555320336759*dxvl[2])))*m_; 
  diffFac[2] = (vtSq[2]*((0.9035079029052515*BmagInv[2]+1.414213562373095*BmagInv[0])*wl[2]+(0.4517539514526256*BmagInv[2]+0.7071067811865475*BmagInv[0])*dxvl[2])+vtSq[0]*BmagInv[2]*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])+BmagInv[1]*vtSq[1]*(1.264911064067352*wl[2]+0.6324555320336759*dxvl[2]))*m_; 

  double Ghat[20]; 
  for(unsigned short int i=0; i<20; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*diffFac[1]*fr[15]-1.897366596101028*diffFac[1]*fl[15]-3.368048396326869*diffFac[2]*fr[13]-3.368048396326869*diffFac[2]*fl[13]+1.897366596101028*diffFac[0]*fr[9]-1.897366596101028*diffFac[0]*fl[9]+2.651650429449552*diffFac[2]*fr[7]-2.651650429449552*diffFac[2]*fl[7]-3.368048396326869*diffFac[1]*fr[5]-3.368048396326869*diffFac[1]*fl[5]-3.368048396326869*diffFac[0]*fr[3]-3.368048396326869*diffFac[0]*fl[3]+2.651650429449552*diffFac[1]*fr[1]-2.651650429449552*diffFac[1]*fl[1]+2.651650429449552*diffFac[0]*fr[0]-2.651650429449552*diffFac[0]*fl[0])*rdv-1.118033988749895*fjump[9]+alphaDrag[0]*(0.7905694150420947*favg[9]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = (1.697056274847714*diffFac[2]*fr[15]+1.897366596101028*diffFac[0]*fr[15]-1.697056274847714*diffFac[2]*fl[15]-1.897366596101028*diffFac[0]*fl[15]-3.012474066278414*diffFac[1]*fr[13]-3.012474066278414*diffFac[1]*fl[13]+1.897366596101028*diffFac[1]*fr[9]-1.897366596101028*diffFac[1]*fl[9]+2.371708245126284*diffFac[1]*fr[7]-2.371708245126284*diffFac[1]*fl[7]-3.012474066278413*diffFac[2]*fr[5]-3.368048396326869*diffFac[0]*fr[5]-3.012474066278413*diffFac[2]*fl[5]-3.368048396326869*diffFac[0]*fl[5]-3.368048396326869*diffFac[1]*fr[3]-3.368048396326869*diffFac[1]*fl[3]+2.371708245126284*fr[1]*diffFac[2]-2.371708245126284*fl[1]*diffFac[2]+2.651650429449552*diffFac[0]*fr[1]-2.651650429449552*diffFac[0]*fl[1]+2.651650429449552*fr[0]*diffFac[1]-2.651650429449552*fl[0]*diffFac[1])*rdv-1.118033988749895*fjump[15]+alphaDrag[0]*(0.7905694150420948*favg[15]+0.6123724356957944*favg[5]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[5]-0.5*fjump[1]; 
  Ghat[2] = (1.897366596101028*diffFac[1]*fr[19]-1.897366596101028*diffFac[1]*fl[19]-3.368048396326869*diffFac[2]*fr[17]-3.368048396326869*diffFac[2]*fl[17]+1.897366596101028*diffFac[0]*fr[16]-1.897366596101028*diffFac[0]*fl[16]+2.651650429449552*diffFac[2]*fr[11]-2.651650429449552*diffFac[2]*fl[11]-3.368048396326869*diffFac[1]*fr[10]-3.368048396326869*diffFac[1]*fl[10]-3.368048396326869*diffFac[0]*fr[6]-3.368048396326869*diffFac[0]*fl[6]+2.651650429449552*diffFac[1]*fr[4]-2.651650429449552*diffFac[1]*fl[4]+2.651650429449552*diffFac[0]*fr[2]-2.651650429449552*diffFac[0]*fl[2])*rdv-1.118033988749895*fjump[16]+alphaDrag[0]*(0.7905694150420948*favg[16]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-0.8660254037844386*fjump[6]-0.5*fjump[2]; 
  Ghat[4] = (1.697056274847714*diffFac[2]*fr[19]+1.897366596101028*diffFac[0]*fr[19]-1.697056274847714*diffFac[2]*fl[19]-1.897366596101028*diffFac[0]*fl[19]-3.012474066278413*diffFac[1]*fr[17]-3.012474066278413*diffFac[1]*fl[17]+1.897366596101028*diffFac[1]*fr[16]-1.897366596101028*diffFac[1]*fl[16]+2.371708245126284*diffFac[1]*fr[11]-2.371708245126284*diffFac[1]*fl[11]-3.012474066278413*diffFac[2]*fr[10]-3.368048396326869*diffFac[0]*fr[10]-3.012474066278413*diffFac[2]*fl[10]-3.368048396326869*diffFac[0]*fl[10]-3.368048396326869*diffFac[1]*fr[6]-3.368048396326869*diffFac[1]*fl[6]+2.371708245126284*diffFac[2]*fr[4]+2.651650429449552*diffFac[0]*fr[4]-2.371708245126284*diffFac[2]*fl[4]-2.651650429449552*diffFac[0]*fl[4]+2.651650429449552*diffFac[1]*fr[2]-2.651650429449552*diffFac[1]*fl[2])*rdv-1.118033988749895*fjump[19]+alphaDrag[0]*(0.7905694150420947*favg[19]+0.6123724356957944*favg[10]+0.3535533905932737*favg[4])-0.8660254037844386*fjump[10]-0.5*fjump[4]; 
  Ghat[7] = (1.697056274847714*diffFac[1]*fr[15]-1.697056274847714*diffFac[1]*fl[15]-2.151767190198866*diffFac[2]*fr[13]-3.368048396326869*diffFac[0]*fr[13]-2.151767190198866*diffFac[2]*fl[13]-3.368048396326869*diffFac[0]*fl[13]+1.897366596101028*diffFac[2]*fr[9]-1.897366596101028*diffFac[2]*fl[9]+1.694077317947346*diffFac[2]*fr[7]+2.651650429449552*diffFac[0]*fr[7]-1.694077317947346*diffFac[2]*fl[7]-2.651650429449552*diffFac[0]*fl[7]-3.012474066278413*diffFac[1]*fr[5]-3.012474066278413*diffFac[1]*fl[5]-3.368048396326869*diffFac[2]*fr[3]-3.368048396326869*diffFac[2]*fl[3]+2.651650429449552*fr[0]*diffFac[2]-2.651650429449552*fl[0]*diffFac[2]+2.371708245126284*diffFac[1]*fr[1]-2.371708245126284*diffFac[1]*fl[1])*rdv-0.8660254037844387*fjump[13]+alphaDrag[0]*(0.6123724356957944*favg[13]+0.3535533905932737*favg[7])-0.5*fjump[7]; 
  Ghat[8] = ((-3.368048396326869*diffFac[1]*fr[18])-3.368048396326869*diffFac[1]*fl[18]-3.368048396326869*diffFac[0]*fr[14]-3.368048396326869*diffFac[0]*fl[14]+2.651650429449552*diffFac[1]*fr[12]-2.651650429449552*diffFac[1]*fl[12]+2.651650429449552*diffFac[0]*fr[8]-2.651650429449552*diffFac[0]*fl[8])*rdv-0.8660254037844387*fjump[14]+alphaDrag[0]*(0.6123724356957944*favg[14]+0.3535533905932737*favg[8])-0.5*fjump[8]; 
  Ghat[11] = (1.697056274847714*diffFac[1]*fr[19]-1.697056274847714*diffFac[1]*fl[19]-2.151767190198866*diffFac[2]*fr[17]-3.368048396326869*diffFac[0]*fr[17]-2.151767190198866*diffFac[2]*fl[17]-3.368048396326869*diffFac[0]*fl[17]+1.897366596101028*diffFac[2]*fr[16]-1.897366596101028*diffFac[2]*fl[16]+1.694077317947346*diffFac[2]*fr[11]+2.651650429449552*diffFac[0]*fr[11]-1.694077317947346*diffFac[2]*fl[11]-2.651650429449552*diffFac[0]*fl[11]-3.012474066278414*diffFac[1]*fr[10]-3.012474066278414*diffFac[1]*fl[10]-3.368048396326869*diffFac[2]*fr[6]-3.368048396326869*diffFac[2]*fl[6]+2.371708245126284*diffFac[1]*fr[4]-2.371708245126284*diffFac[1]*fl[4]+2.651650429449552*diffFac[2]*fr[2]-2.651650429449552*diffFac[2]*fl[2])*rdv-0.8660254037844387*fjump[17]+alphaDrag[0]*(0.6123724356957944*favg[17]+0.3535533905932737*favg[11])-0.5*fjump[11]; 
  Ghat[12] = ((-3.012474066278414*diffFac[2]*fr[18])-3.368048396326869*diffFac[0]*fr[18]-3.012474066278414*diffFac[2]*fl[18]-3.368048396326869*diffFac[0]*fl[18]-3.368048396326869*diffFac[1]*fr[14]-3.368048396326869*diffFac[1]*fl[14]+2.371708245126284*diffFac[2]*fr[12]+2.651650429449552*diffFac[0]*fr[12]-2.371708245126284*diffFac[2]*fl[12]-2.651650429449552*diffFac[0]*fl[12]+2.651650429449552*diffFac[1]*fr[8]-2.651650429449552*diffFac[1]*fl[8])*rdv-0.8660254037844387*fjump[18]+alphaDrag[0]*(0.6123724356957944*favg[18]+0.3535533905932737*favg[12])-0.5*fjump[12]; 

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
  incr2[3] = 0.2995357736356374*diffFac[1]*(fr[15]+fl[15])+diffFac[2]*(0.430893194785552*fl[13]-0.430893194785552*fr[13])+0.2995357736356374*diffFac[0]*(fr[9]+fl[9])+0.3061862178478971*diffFac[2]*(fr[7]+fl[7])+diffFac[1]*(0.430893194785552*fl[5]-0.430893194785552*fr[5])+diffFac[0]*(0.430893194785552*fl[3]-0.430893194785552*fr[3])+0.3061862178478971*(diffFac[1]*(fr[1]+fl[1])+diffFac[0]*(fr[0]+fl[0])); 
  incr2[5] = (0.2679129406169099*diffFac[2]+0.2995357736356374*diffFac[0])*fr[15]+(0.2679129406169099*diffFac[2]+0.2995357736356374*diffFac[0])*fl[15]+diffFac[1]*((-0.3854025898330209*fr[13])+0.3854025898330209*fl[13]+0.2995357736356374*(fr[9]+fl[9])+0.273861278752583*(fr[7]+fl[7]))+((-0.3854025898330209*diffFac[2])-0.430893194785552*diffFac[0])*fr[5]+(0.3854025898330209*diffFac[2]+0.430893194785552*diffFac[0])*fl[5]+diffFac[1]*(0.430893194785552*fl[3]-0.430893194785552*fr[3])+0.273861278752583*(fr[1]+fl[1])*diffFac[2]+0.3061862178478971*(diffFac[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*diffFac[1]); 
  incr2[6] = 0.2995357736356374*diffFac[1]*(fr[19]+fl[19])+diffFac[2]*(0.430893194785552*fl[17]-0.430893194785552*fr[17])+0.2995357736356374*diffFac[0]*(fr[16]+fl[16])+0.3061862178478971*diffFac[2]*(fr[11]+fl[11])+diffFac[1]*(0.430893194785552*fl[10]-0.430893194785552*fr[10])+diffFac[0]*(0.430893194785552*fl[6]-0.430893194785552*fr[6])+0.3061862178478971*(diffFac[1]*(fr[4]+fl[4])+diffFac[0]*(fr[2]+fl[2])); 
  incr2[9] = (-1.160097062884178*diffFac[1]*(fr[15]+fl[15]))+diffFac[2]*(1.668842167398551*fr[13]-1.668842167398551*fl[13])-1.160097062884178*diffFac[0]*(fr[9]+fl[9])-1.185854122563142*diffFac[2]*(fr[7]+fl[7])+diffFac[1]*(1.668842167398551*fr[5]-1.668842167398551*fl[5])+diffFac[0]*(1.668842167398551*fr[3]-1.668842167398551*fl[3])-1.185854122563142*(diffFac[1]*(fr[1]+fl[1])+diffFac[0]*(fr[0]+fl[0])); 
  incr2[10] = (0.2679129406169099*diffFac[2]+0.2995357736356374*diffFac[0])*fr[19]+(0.2679129406169099*diffFac[2]+0.2995357736356374*diffFac[0])*fl[19]+diffFac[1]*((-0.3854025898330209*fr[17])+0.3854025898330209*fl[17]+0.2995357736356374*(fr[16]+fl[16])+0.273861278752583*(fr[11]+fl[11]))+((-0.3854025898330209*diffFac[2])-0.430893194785552*diffFac[0])*fr[10]+(0.3854025898330209*diffFac[2]+0.430893194785552*diffFac[0])*fl[10]+diffFac[1]*(0.430893194785552*fl[6]-0.430893194785552*fr[6])+(0.273861278752583*diffFac[2]+0.3061862178478971*diffFac[0])*fr[4]+0.273861278752583*diffFac[2]*fl[4]+0.3061862178478971*(diffFac[0]*fl[4]+diffFac[1]*(fr[2]+fl[2])); 
  incr2[13] = 0.2679129406169099*diffFac[1]*(fr[15]+fl[15])+((-0.2752875641664436*diffFac[2])-0.430893194785552*diffFac[0])*fr[13]+(0.2752875641664436*diffFac[2]+0.430893194785552*diffFac[0])*fl[13]+0.2995357736356374*diffFac[2]*(fr[9]+fl[9])+(0.1956151991089878*diffFac[2]+0.3061862178478971*diffFac[0])*fr[7]+(0.1956151991089878*diffFac[2]+0.3061862178478971*diffFac[0])*fl[7]+diffFac[1]*(0.3854025898330209*fl[5]-0.3854025898330209*fr[5])+diffFac[2]*((-0.430893194785552*fr[3])+0.430893194785552*fl[3]+0.3061862178478971*(fr[0]+fl[0]))+0.273861278752583*diffFac[1]*(fr[1]+fl[1]); 
  incr2[14] = diffFac[1]*(0.430893194785552*fl[18]-0.430893194785552*fr[18])+diffFac[0]*(0.430893194785552*fl[14]-0.430893194785552*fr[14])+0.3061862178478971*(diffFac[1]*(fr[12]+fl[12])+diffFac[0]*(fr[8]+fl[8])); 
  incr2[15] = ((-1.037622357242749*diffFac[2])-1.160097062884178*diffFac[0])*fr[15]+((-1.037622357242749*diffFac[2])-1.160097062884178*diffFac[0])*fl[15]+diffFac[1]*(1.492657812008498*fr[13]-1.492657812008498*fl[13]-1.160097062884178*(fr[9]+fl[9])-1.060660171779821*(fr[7]+fl[7]))+(1.492657812008498*diffFac[2]+1.668842167398551*diffFac[0])*fr[5]+((-1.492657812008498*diffFac[2])-1.668842167398551*diffFac[0])*fl[5]+diffFac[1]*(1.668842167398551*fr[3]-1.668842167398551*fl[3])-1.060660171779821*(fr[1]+fl[1])*diffFac[2]-1.185854122563142*(diffFac[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*diffFac[1]); 
  incr2[16] = (-1.160097062884178*diffFac[1]*(fr[19]+fl[19]))+diffFac[2]*(1.668842167398551*fr[17]-1.668842167398551*fl[17])-1.160097062884178*diffFac[0]*(fr[16]+fl[16])-1.185854122563142*diffFac[2]*(fr[11]+fl[11])+diffFac[1]*(1.668842167398551*fr[10]-1.668842167398551*fl[10])+diffFac[0]*(1.668842167398551*fr[6]-1.668842167398551*fl[6])-1.185854122563142*(diffFac[1]*(fr[4]+fl[4])+diffFac[0]*(fr[2]+fl[2])); 
  incr2[17] = 0.2679129406169099*diffFac[1]*(fr[19]+fl[19])+((-0.2752875641664436*diffFac[2])-0.430893194785552*diffFac[0])*fr[17]+(0.2752875641664436*diffFac[2]+0.430893194785552*diffFac[0])*fl[17]+0.2995357736356374*diffFac[2]*(fr[16]+fl[16])+(0.1956151991089878*diffFac[2]+0.3061862178478971*diffFac[0])*fr[11]+(0.1956151991089878*diffFac[2]+0.3061862178478971*diffFac[0])*fl[11]+diffFac[1]*(0.3854025898330209*fl[10]-0.3854025898330209*fr[10])+diffFac[2]*(0.430893194785552*fl[6]-0.430893194785552*fr[6])+0.273861278752583*diffFac[1]*(fr[4]+fl[4])+0.3061862178478971*diffFac[2]*(fr[2]+fl[2]); 
  incr2[18] = ((-0.3854025898330209*diffFac[2])-0.430893194785552*diffFac[0])*fr[18]+(0.3854025898330209*diffFac[2]+0.430893194785552*diffFac[0])*fl[18]+diffFac[1]*(0.430893194785552*fl[14]-0.430893194785552*fr[14])+(0.273861278752583*diffFac[2]+0.3061862178478971*diffFac[0])*fr[12]+0.273861278752583*diffFac[2]*fl[12]+0.3061862178478971*(diffFac[0]*fl[12]+diffFac[1]*(fr[8]+fl[8])); 
  incr2[19] = ((-1.037622357242749*diffFac[2])-1.160097062884178*diffFac[0])*fr[19]+((-1.037622357242749*diffFac[2])-1.160097062884178*diffFac[0])*fl[19]+diffFac[1]*(1.492657812008498*fr[17]-1.492657812008498*fl[17]-1.160097062884178*(fr[16]+fl[16])-1.060660171779821*(fr[11]+fl[11]))+(1.492657812008498*diffFac[2]+1.668842167398551*diffFac[0])*fr[10]+((-1.492657812008498*diffFac[2])-1.668842167398551*diffFac[0])*fl[10]+diffFac[1]*(1.668842167398551*fr[6]-1.668842167398551*fl[6])+((-1.060660171779821*diffFac[2])-1.185854122563142*diffFac[0])*fr[4]-1.060660171779821*diffFac[2]*fl[4]-1.185854122563142*(diffFac[0]*fl[4]+diffFac[1]*(fr[2]+fl[2])); 

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
  outr[17] += incr2[17]*rdvSq4nuR+incr1[17]*rdv2nuR; 
  outr[18] += incr2[18]*rdvSq4nuR+incr1[18]*rdv2nuR; 
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
  outl[17] += incr1[17]*rdv2nuL-1.0*incr2[17]*rdvSq4nuL; 
  outl[18] += incr1[18]*rdv2nuL-1.0*incr2[18]*rdvSq4nuL; 
  outl[19] += incr2[19]*rdvSq4nuL-1.0*incr1[19]*rdv2nuL; 

  return std::abs(2.0*wl[2]); 
} 
