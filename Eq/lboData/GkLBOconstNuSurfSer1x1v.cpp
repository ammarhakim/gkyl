#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf1x1vSer_Vpar_P1(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:         Cell-center coordinates. 
  // dxv[2]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[2]; 
  alphaDrSurf[0] = 0.7071067811865475*((2.0*wl[1]+dxvl[1])*nuSum-1.414213562373095*nuUSum[0]); 
  alphaDrSurf[1] = -1.0*nuUSum[1]; 

  double fUpwindQuad[2];
  fUpwindQuad[0] = 0.5*(0.8660254037844386*fr[3]-0.8660254037844386*(fl[3]+fr[2])+0.8660254037844386*fl[2]-0.5*(fr[1]+fl[1])+0.5*(fr[0]+fl[0]))-0.5*((-0.8660254037844386*(fr[3]+fl[3]))+0.8660254037844386*(fr[2]+fl[2])+0.5*fr[1]-0.5*(fl[1]+fr[0])+0.5*fl[0])*sgn(alphaDrSurf[0]-alphaDrSurf[1]); 
  fUpwindQuad[1] = 0.5*((-0.8660254037844386*(fr[3]+fr[2]))+0.8660254037844386*(fl[3]+fl[2])+0.5*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.8660254037844386*(fr[3]+fl[3]+fr[2]+fl[2])-0.5*(fr[1]+fr[0])+0.5*(fl[1]+fl[0]))*sgn(alphaDrSurf[1]+alphaDrSurf[0]); 

  double fUpwind[2];
  fUpwind[0] = 0.7071067811865475*(fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.7071067811865475*(fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  double Ghat[4]; 
  for(unsigned short int i=0; i<4; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.530931089239486*nuVtSqSum[1]*fr[3])-1.530931089239486*nuVtSqSum[1]*fl[3]-1.530931089239486*nuVtSqSum[0]*fr[2]-1.530931089239486*nuVtSqSum[0]*fl[2]+1.590990257669731*fr[1]*nuVtSqSum[1]-1.590990257669731*fl[1]*nuVtSqSum[1]+1.590990257669731*fr[0]*nuVtSqSum[0]-1.590990257669731*fl[0]*nuVtSqSum[0])*rdv+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = ((-1.530931089239486*nuVtSqSum[0]*fr[3])-1.530931089239486*nuVtSqSum[0]*fl[3]-1.530931089239486*nuVtSqSum[1]*fr[2]-1.530931089239486*nuVtSqSum[1]*fl[2]+1.590990257669731*fr[0]*nuVtSqSum[1]-1.590990257669731*fl[0]*nuVtSqSum[1]+1.590990257669731*nuVtSqSum[0]*fr[1]-1.590990257669731*nuVtSqSum[0]*fl[1])*rdv+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 

  double incr2[4]; 
  incr2[2] = nuVtSqSum[1]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[3] = nuVtSqSum[0]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 

  return std::abs(wl[1]-(0.7071067811865475*nuUSum[0])/nuSum); 
} 
double GkLBOconstNuSurf1x1vSer_Vpar_P2(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:         Cell-center coordinates. 
  // dxv[2]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[3]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[3]; 
  alphaDrSurf[0] = 0.7071067811865475*((2.0*wl[1]+dxvl[1])*nuSum-1.414213562373095*nuUSum[0]); 
  alphaDrSurf[1] = -1.0*nuUSum[1]; 
  alphaDrSurf[2] = -1.0*nuUSum[2]; 

  double fUpwindQuad[3];
  fUpwindQuad[0] = 0.5*((-1.5*(fr[7]+fl[7]))-0.7745966692414833*fr[6]+0.7745966692414833*fl[6]+1.118033988749895*(fr[5]+fl[5])+0.4472135954999579*(fr[4]+fl[4])+1.161895003862225*fr[3]-1.161895003862225*fl[3]-0.8660254037844386*fr[2]+0.8660254037844386*fl[2]-0.6708203932499369*(fr[1]+fl[1])+0.5*(fr[0]+fl[0]))-0.5*(1.5*fr[7]-1.5*fl[7]+0.7745966692414833*(fr[6]+fl[6])-1.118033988749895*fr[5]+1.118033988749895*fl[5]-0.4472135954999579*fr[4]+0.4472135954999579*fl[4]-1.161895003862225*(fr[3]+fl[3])+0.8660254037844386*(fr[2]+fl[2])+0.6708203932499369*fr[1]-0.6708203932499369*fl[1]-0.5*fr[0]+0.5*fl[0])*sgn(0.6324555320336759*alphaDrSurf[2]-0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]); 
  fUpwindQuad[1] = 0.5*(0.9682458365518543*fr[6]-0.9682458365518543*fl[6]+1.118033988749895*(fr[5]+fl[5])-0.5590169943749475*(fr[4]+fl[4])-0.8660254037844386*fr[2]+0.8660254037844386*fl[2]+0.5*(fr[0]+fl[0]))-0.5*((-0.9682458365518543*(fr[6]+fl[6]))-1.118033988749895*fr[5]+1.118033988749895*fl[5]+0.5590169943749475*fr[4]-0.5590169943749475*fl[4]+0.8660254037844386*(fr[2]+fl[2])-0.5*fr[0]+0.5*fl[0])*sgn(0.7071067811865475*alphaDrSurf[0]-0.7905694150420947*alphaDrSurf[2]); 
  fUpwindQuad[2] = 0.5*(1.5*(fr[7]+fl[7])-0.7745966692414833*fr[6]+0.7745966692414833*fl[6]+1.118033988749895*(fr[5]+fl[5])+0.4472135954999579*(fr[4]+fl[4])-1.161895003862225*fr[3]+1.161895003862225*fl[3]-0.8660254037844386*fr[2]+0.8660254037844386*fl[2]+0.6708203932499369*(fr[1]+fl[1])+0.5*(fr[0]+fl[0]))-0.5*((-1.5*fr[7])+1.5*fl[7]+0.7745966692414833*(fr[6]+fl[6])-1.118033988749895*fr[5]+1.118033988749895*fl[5]-0.4472135954999579*fr[4]+0.4472135954999579*fl[4]+1.161895003862225*(fr[3]+fl[3])+0.8660254037844386*(fr[2]+fl[2])-0.6708203932499369*fr[1]+0.6708203932499369*fl[1]-0.5*fr[0]+0.5*fl[0])*sgn(0.6324555320336759*alphaDrSurf[2]+0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]); 

  double fUpwind[3];
  fUpwind[0] = 0.07856742013183861*(5.0*fUpwindQuad[2]+8.0*fUpwindQuad[1]+5.0*fUpwindQuad[0]); 
  fUpwind[1] = 0.5270462766947298*(fUpwindQuad[2]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.3513641844631533*(fUpwindQuad[2]-2.0*fUpwindQuad[1]+fUpwindQuad[0]); 

  double Ghat[8]; 
  for(unsigned short int i=0; i<8; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = (1.897366596101028*nuVtSqSum[1]*fr[7]-1.897366596101028*nuVtSqSum[1]*fl[7]-3.368048396326869*nuVtSqSum[2]*fr[6]-3.368048396326869*nuVtSqSum[2]*fl[6]+1.897366596101028*nuVtSqSum[0]*fr[5]-1.897366596101028*nuVtSqSum[0]*fl[5]+2.651650429449552*nuVtSqSum[2]*fr[4]-2.651650429449552*nuVtSqSum[2]*fl[4]-3.368048396326869*nuVtSqSum[1]*fr[3]-3.368048396326869*nuVtSqSum[1]*fl[3]-3.368048396326869*nuVtSqSum[0]*fr[2]-3.368048396326869*nuVtSqSum[0]*fl[2]+2.651650429449552*fr[1]*nuVtSqSum[1]-2.651650429449552*fl[1]*nuVtSqSum[1]+2.651650429449552*fr[0]*nuVtSqSum[0]-2.651650429449552*fl[0]*nuVtSqSum[0])*rdv+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = (1.697056274847714*nuVtSqSum[2]*fr[7]+1.897366596101028*nuVtSqSum[0]*fr[7]-1.697056274847714*nuVtSqSum[2]*fl[7]-1.897366596101028*nuVtSqSum[0]*fl[7]-3.012474066278414*nuVtSqSum[1]*fr[6]-3.012474066278414*nuVtSqSum[1]*fl[6]+1.897366596101028*nuVtSqSum[1]*fr[5]-1.897366596101028*nuVtSqSum[1]*fl[5]+2.371708245126284*nuVtSqSum[1]*fr[4]-2.371708245126284*nuVtSqSum[1]*fl[4]-3.012474066278413*nuVtSqSum[2]*fr[3]-3.368048396326869*nuVtSqSum[0]*fr[3]-3.012474066278413*nuVtSqSum[2]*fl[3]-3.368048396326869*nuVtSqSum[0]*fl[3]+2.371708245126284*fr[1]*nuVtSqSum[2]-2.371708245126284*fl[1]*nuVtSqSum[2]-3.368048396326869*nuVtSqSum[1]*fr[2]-3.368048396326869*nuVtSqSum[1]*fl[2]+2.651650429449552*fr[0]*nuVtSqSum[1]-2.651650429449552*fl[0]*nuVtSqSum[1]+2.651650429449552*nuVtSqSum[0]*fr[1]-2.651650429449552*nuVtSqSum[0]*fl[1])*rdv+0.8944271909999159*alphaDrSurf[1]*fUpwind[2]+0.8944271909999159*fUpwind[1]*alphaDrSurf[2]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]; 
  Ghat[4] = (1.697056274847714*nuVtSqSum[1]*fr[7]-1.697056274847714*nuVtSqSum[1]*fl[7]-2.151767190198866*nuVtSqSum[2]*fr[6]-3.368048396326869*nuVtSqSum[0]*fr[6]-2.151767190198866*nuVtSqSum[2]*fl[6]-3.368048396326869*nuVtSqSum[0]*fl[6]+1.897366596101028*nuVtSqSum[2]*fr[5]-1.897366596101028*nuVtSqSum[2]*fl[5]+1.694077317947346*nuVtSqSum[2]*fr[4]+2.651650429449552*nuVtSqSum[0]*fr[4]-1.694077317947346*nuVtSqSum[2]*fl[4]-2.651650429449552*nuVtSqSum[0]*fl[4]-3.012474066278413*nuVtSqSum[1]*fr[3]-3.012474066278413*nuVtSqSum[1]*fl[3]-3.368048396326869*fr[2]*nuVtSqSum[2]-3.368048396326869*fl[2]*nuVtSqSum[2]+2.651650429449552*fr[0]*nuVtSqSum[2]-2.651650429449552*fl[0]*nuVtSqSum[2]+2.371708245126284*fr[1]*nuVtSqSum[1]-2.371708245126284*fl[1]*nuVtSqSum[1])*rdv+0.6388765649999399*alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]+0.8944271909999159*alphaDrSurf[1]*fUpwind[1]; 

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
  incr2[2] = nuVtSqSum[1]*(0.2995357736356374*(fr[7]+fl[7])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[4]+fl[4]))+nuVtSqSum[0]*(0.2995357736356374*(fr[5]+fl[5])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[3] = nuVtSqSum[0]*(0.2995357736356374*(fr[7]+fl[7])-0.430893194785552*fr[3]+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[2]*(0.2679129406169099*(fr[7]+fl[7])-0.3854025898330209*fr[3]+0.3854025898330209*fl[3]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.3854025898330209*fr[6])+0.3854025898330209*fl[6]+0.2995357736356374*(fr[5]+fl[5])+0.273861278752583*(fr[4]+fl[4])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[5] = nuVtSqSum[1]*((-1.160097062884178*(fr[7]+fl[7]))+1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[1]+fl[1]))+nuVtSqSum[2]*(1.668842167398551*fr[6]-1.668842167398551*fl[6]-1.185854122563142*(fr[4]+fl[4]))+nuVtSqSum[0]*((-1.160097062884178*(fr[5]+fl[5]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0])); 
  incr2[6] = nuVtSqSum[1]*(0.2679129406169099*(fr[7]+fl[7])-0.3854025898330209*fr[3]+0.3854025898330209*fl[3]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.2752875641664436*fr[6])+0.2752875641664436*fl[6]+0.2995357736356374*(fr[5]+fl[5])+0.1956151991089878*(fr[4]+fl[4])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[4]+fl[4])); 
  incr2[7] = nuVtSqSum[2]*((-1.037622357242749*(fr[7]+fl[7]))+1.492657812008498*fr[3]-1.492657812008498*fl[3]-1.060660171779821*(fr[1]+fl[1]))+nuVtSqSum[0]*((-1.160097062884178*(fr[7]+fl[7]))+1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[1]+fl[1]))+nuVtSqSum[1]*(1.492657812008498*fr[6]-1.492657812008498*fl[6]-1.160097062884178*(fr[5]+fl[5])-1.060660171779821*(fr[4]+fl[4])+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr2[5]*rdvSq4L-1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr2[7]*rdvSq4L-1.0*incr1[7]*rdv2L; 

  return std::abs((0.7905694150420947*nuUSum[2])/nuSum-(0.7071067811865475*nuUSum[0])/nuSum+wl[1]); 
} 
