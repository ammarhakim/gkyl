#include <GkLBOModDecl.h> 
double GkLBOconstNuSurfNonUniform1x2vTensor_Vpar_P1(const double m_, const double cflL, const double cflR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (2.0*wl[1]+dxvl[1])*nuSum-1.414213562373095*nuUSum[0]; 
  alphaDrSurf[1] = -1.414213562373095*nuUSum[1]; 

  double fUpwindQuad[4];
  fUpwindQuad[0] = 0.5*((-0.6123724356957944*fr[7])+0.6123724356957944*(fl[7]+fr[6])-0.6123724356957944*fl[6]+0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]-alphaDrSurf[1]); 
  fUpwindQuad[1] = 0.5*(0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6])-0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[2] = 0.5*(0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*fl[6]-0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]-alphaDrSurf[1]); 
  fUpwindQuad[3] = 0.5*((-0.6123724356957944*(fr[7]+fr[6]))+0.6123724356957944*(fl[7]+fl[6])+0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaDrSurf[1]+alphaDrSurf[0]); 

  double fUpwind[4];
  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  double Gdiff[8]; 
  double Ghat[8]; 
  double incr2[8]; 

  const double dxvl1R2 = pow(dxvl[1],2);
  const double dxvl1R3 = pow(dxvl[1],3);
  const double dxvr1R2 = pow(dxvr[1],2);
  const double dxvr1R3 = pow(dxvr[1],3);

  incr2[2] = -(1.0*((5.0*dxvl1R2*dxvr[1]+3.0*dxvl1R3)*nuVtSqSum[1]*fr[4]+((-3.0*dxvr1R3)-5.0*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[4]+(5.0*nuVtSqSum[0]*dxvl1R2*dxvr[1]+3.0*nuVtSqSum[0]*dxvl1R3)*fr[2]+((-3.0*nuVtSqSum[0]*dxvr1R3)-5.0*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[2]+(((-5.196152422706631*dxvl1R2*dxvr[1])-1.732050807568877*dxvl1R3)*fr[1]+((-1.732050807568877*dxvr1R3)-5.196152422706631*dxvl[1]*dxvr1R2)*fl[1])*nuVtSqSum[1]-1.732050807568877*fl[0]*nuVtSqSum[0]*dxvr1R3-5.196152422706631*fl[0]*nuVtSqSum[0]*dxvl[1]*dxvr1R2-5.196152422706631*fr[0]*nuVtSqSum[0]*dxvl1R2*dxvr[1]-1.732050807568877*fr[0]*nuVtSqSum[0]*dxvl1R3))/(2.828427124746191*dxvr1R3+8.485281374238571*dxvl[1]*dxvr1R2+8.485281374238571*dxvl1R2*dxvr[1]+2.828427124746191*dxvl1R3); 
  incr2[4] = -(1.0*((5.0*nuVtSqSum[0]*dxvl1R2*dxvr[1]+3.0*nuVtSqSum[0]*dxvl1R3)*fr[4]+((-3.0*nuVtSqSum[0]*dxvr1R3)-5.0*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[4]+(5.0*dxvl1R2*dxvr[1]+3.0*dxvl1R3)*nuVtSqSum[1]*fr[2]+((-3.0*dxvr1R3)-5.0*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[2]+((-1.732050807568877*fl[0]*dxvr1R3)-5.196152422706631*fl[0]*dxvl[1]*dxvr1R2-5.196152422706631*fr[0]*dxvl1R2*dxvr[1]-1.732050807568877*fr[0]*dxvl1R3)*nuVtSqSum[1]+((-5.196152422706631*nuVtSqSum[0]*dxvl1R2*dxvr[1])-1.732050807568877*nuVtSqSum[0]*dxvl1R3)*fr[1]+((-1.732050807568877*nuVtSqSum[0]*dxvr1R3)-5.196152422706631*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[1]))/(2.828427124746191*dxvr1R3+8.485281374238571*dxvl[1]*dxvr1R2+8.485281374238571*dxvl1R2*dxvr[1]+2.828427124746191*dxvl1R3); 
  incr2[6] = -(1.0*((5.0*dxvl1R2*dxvr[1]+3.0*dxvl1R3)*nuVtSqSum[1]*fr[7]+((-3.0*dxvr1R3)-5.0*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[7]+(5.0*nuVtSqSum[0]*dxvl1R2*dxvr[1]+3.0*nuVtSqSum[0]*dxvl1R3)*fr[6]+((-3.0*nuVtSqSum[0]*dxvr1R3)-5.0*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[6]+((-5.196152422706631*dxvl1R2*dxvr[1])-1.732050807568877*dxvl1R3)*nuVtSqSum[1]*fr[5]+((-1.732050807568877*dxvr1R3)-5.196152422706631*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[5]+((-5.196152422706631*nuVtSqSum[0]*dxvl1R2*dxvr[1])-1.732050807568877*nuVtSqSum[0]*dxvl1R3)*fr[3]+((-1.732050807568877*nuVtSqSum[0]*dxvr1R3)-5.196152422706631*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[3]))/(2.828427124746191*dxvr1R3+8.485281374238571*dxvl[1]*dxvr1R2+8.485281374238571*dxvl1R2*dxvr[1]+2.828427124746191*dxvl1R3); 
  incr2[7] = -(1.0*((5.0*nuVtSqSum[0]*dxvl1R2*dxvr[1]+3.0*nuVtSqSum[0]*dxvl1R3)*fr[7]+((-3.0*nuVtSqSum[0]*dxvr1R3)-5.0*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[7]+(5.0*dxvl1R2*dxvr[1]+3.0*dxvl1R3)*nuVtSqSum[1]*fr[6]+((-3.0*dxvr1R3)-5.0*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[6]+((-5.196152422706631*nuVtSqSum[0]*dxvl1R2*dxvr[1])-1.732050807568877*nuVtSqSum[0]*dxvl1R3)*fr[5]+((-1.732050807568877*nuVtSqSum[0]*dxvr1R3)-5.196152422706631*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[5]+((-5.196152422706631*dxvl1R2*dxvr[1])-1.732050807568877*dxvl1R3)*nuVtSqSum[1]*fr[3]+((-1.732050807568877*dxvr1R3)-5.196152422706631*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[3]))/(2.828427124746191*dxvr1R3+8.485281374238571*dxvl[1]*dxvr1R2+8.485281374238571*dxvl1R2*dxvr[1]+2.828427124746191*dxvl1R3); 

  const double dxvl1R4 = pow(dxvl[1],4);
  const double dxvr1R4 = pow(dxvr[1],4);

  Gdiff[0] = -(1.0*((17.32050807568877*dxvl1R2*dxvr1R2+3.464101615137754*dxvl1R3*dxvr[1]-3.464101615137754*dxvl1R4)*nuVtSqSum[1]*fr[4]+((-3.464101615137754*dxvr1R4)+3.464101615137754*dxvl[1]*dxvr1R3+17.32050807568877*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[4]+(17.32050807568877*nuVtSqSum[0]*dxvl1R2*dxvr1R2+3.464101615137754*nuVtSqSum[0]*dxvl1R3*dxvr[1]-3.464101615137754*nuVtSqSum[0]*dxvl1R4)*fr[2]+((-3.464101615137754*nuVtSqSum[0]*dxvr1R4)+3.464101615137754*nuVtSqSum[0]*dxvl[1]*dxvr1R3+17.32050807568877*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[2]+(18.0*dxvl1R2*dxvr1R2*fl[1]-18.0*dxvl1R2*dxvr1R2*fr[1])*nuVtSqSum[1]+(18.0*fl[0]-18.0*fr[0])*nuVtSqSum[0]*dxvl1R2*dxvr1R2))/(1.414213562373095*dxvl[1]*dxvr1R4+4.242640687119286*dxvl1R2*dxvr1R3+4.242640687119286*dxvl1R3*dxvr1R2+1.414213562373095*dxvl1R4*dxvr[1]); 
  Gdiff[1] = -(1.0*((17.32050807568877*nuVtSqSum[0]*dxvl1R2*dxvr1R2+3.464101615137754*nuVtSqSum[0]*dxvl1R3*dxvr[1]-3.464101615137754*nuVtSqSum[0]*dxvl1R4)*fr[4]+((-3.464101615137754*nuVtSqSum[0]*dxvr1R4)+3.464101615137754*nuVtSqSum[0]*dxvl[1]*dxvr1R3+17.32050807568877*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[4]+(17.32050807568877*dxvl1R2*dxvr1R2+3.464101615137754*dxvl1R3*dxvr[1]-3.464101615137754*dxvl1R4)*nuVtSqSum[1]*fr[2]+((-3.464101615137754*dxvr1R4)+3.464101615137754*dxvl[1]*dxvr1R3+17.32050807568877*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[2]+(18.0*fl[0]-18.0*fr[0])*dxvl1R2*dxvr1R2*nuVtSqSum[1]-18.0*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[1]+18.0*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[1]))/(1.414213562373095*dxvl[1]*dxvr1R4+4.242640687119286*dxvl1R2*dxvr1R3+4.242640687119286*dxvl1R3*dxvr1R2+1.414213562373095*dxvl1R4*dxvr[1]); 
  Gdiff[3] = -(1.0*((17.32050807568877*dxvl1R2*dxvr1R2+3.464101615137754*dxvl1R3*dxvr[1]-3.464101615137754*dxvl1R4)*nuVtSqSum[1]*fr[7]+((-3.464101615137754*dxvr1R4)+3.464101615137754*dxvl[1]*dxvr1R3+17.32050807568877*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[7]+(17.32050807568877*nuVtSqSum[0]*dxvl1R2*dxvr1R2+3.464101615137754*nuVtSqSum[0]*dxvl1R3*dxvr[1]-3.464101615137754*nuVtSqSum[0]*dxvl1R4)*fr[6]+((-3.464101615137754*nuVtSqSum[0]*dxvr1R4)+3.464101615137754*nuVtSqSum[0]*dxvl[1]*dxvr1R3+17.32050807568877*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[6]-18.0*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fr[5]+18.0*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fl[5]-18.0*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[3]+18.0*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[3]))/(1.414213562373095*dxvl[1]*dxvr1R4+4.242640687119286*dxvl1R2*dxvr1R3+4.242640687119286*dxvl1R3*dxvr1R2+1.414213562373095*dxvl1R4*dxvr[1]); 
  Gdiff[5] = -(1.0*((17.32050807568877*nuVtSqSum[0]*dxvl1R2*dxvr1R2+3.464101615137754*nuVtSqSum[0]*dxvl1R3*dxvr[1]-3.464101615137754*nuVtSqSum[0]*dxvl1R4)*fr[7]+((-3.464101615137754*nuVtSqSum[0]*dxvr1R4)+3.464101615137754*nuVtSqSum[0]*dxvl[1]*dxvr1R3+17.32050807568877*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[7]+(17.32050807568877*dxvl1R2*dxvr1R2+3.464101615137754*dxvl1R3*dxvr[1]-3.464101615137754*dxvl1R4)*nuVtSqSum[1]*fr[6]+((-3.464101615137754*dxvr1R4)+3.464101615137754*dxvl[1]*dxvr1R3+17.32050807568877*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[6]-18.0*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[5]+18.0*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[5]-18.0*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fr[3]+18.0*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fl[3]))/(1.414213562373095*dxvl[1]*dxvr1R4+4.242640687119286*dxvl1R2*dxvr1R3+4.242640687119286*dxvl1R3*dxvr1R2+1.414213562373095*dxvl1R4*dxvr[1]); 

  Ghat[0] = 0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]+Gdiff[0]; 
  Ghat[1] = 0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]+Gdiff[1]; 
  Ghat[3] = 0.7071067811865475*alphaDrSurf[1]*fUpwind[3]+Gdiff[3]+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]; 
  Ghat[5] = Gdiff[5]+0.7071067811865475*alphaDrSurf[0]*fUpwind[3]+0.7071067811865475*alphaDrSurf[1]*fUpwind[2]; 

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[1]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[3]; 
  incr1[7] = 0.8660254037844386*Ghat[5]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 

  return std::abs(wl[1]-(0.7071067811865475*nuUSum[0])/nuSum); 
} 
double GkLBOconstNuSurfNonUniform1x2vTensor_Mu_P1(const double m_, const double cflL, const double cflR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (4.0*wl[2]+2.0*dxvl[2])*nuSum; 

  double fUpwindQuad[4];
  fUpwindQuad[0] = 0.5*((-0.6123724356957944*fr[7])+0.6123724356957944*(fl[7]+fr[6]+fr[5])-0.6123724356957944*(fl[6]+fl[5])+0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.3535533905932737*(fr[2]+fl[2]+fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6]+fr[5]+fl[5])-0.3535533905932737*fr[4]+0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*(fr[2]+fr[1])-0.3535533905932737*(fl[2]+fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]); 
  fUpwindQuad[1] = 0.5*(0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6]+fr[5])+0.6123724356957944*fl[5]-0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.3535533905932737*(fr[2]+fl[2])+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.6123724356957944*(fr[5]+fl[5])+0.3535533905932737*fr[4]-0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*fr[2]-0.3535533905932737*(fl[2]+fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaDrSurf[0]); 
  fUpwindQuad[2] = 0.5*(0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*(fl[6]+fr[5])-0.6123724356957944*fl[5]-0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.3535533905932737*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])-0.6123724356957944*(fr[5]+fl[5])+0.3535533905932737*fr[4]-0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.3535533905932737*fr[2]+0.3535533905932737*(fl[2]+fr[1])-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]); 
  fUpwindQuad[3] = 0.5*((-0.6123724356957944*(fr[7]+fr[6]+fr[5]))+0.6123724356957944*(fl[7]+fl[6]+fl[5])+0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.3535533905932737*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5])-0.3535533905932737*fr[4]+0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.3535533905932737*(fr[2]+fr[1]+fr[0])+0.3535533905932737*(fl[2]+fl[1]+fl[0]))*sgn(alphaDrSurf[0]); 

  double fUpwind[4];
  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  double diffFac[2]; 
  diffFac[0] = nuVtSqSum[1]*(1.414213562373095*BmagInv[1]*wl[2]*m_+0.7071067811865475*BmagInv[1]*dxvl[2]*m_)+nuVtSqSum[0]*(1.414213562373095*BmagInv[0]*wl[2]*m_+0.7071067811865475*BmagInv[0]*dxvl[2]*m_); 
  diffFac[1] = nuVtSqSum[0]*(1.414213562373095*BmagInv[1]*wl[2]*m_+0.7071067811865475*BmagInv[1]*dxvl[2]*m_)+nuVtSqSum[1]*(1.414213562373095*BmagInv[0]*wl[2]*m_+0.7071067811865475*BmagInv[0]*dxvl[2]*m_); 

  double Gdiff[8]; 
  double Ghat[8]; 
  double incr2[8]; 

  const double dxvl2R2 = pow(dxvl[2],2);
  const double dxvl2R3 = pow(dxvl[2],3);
  const double dxvr2R2 = pow(dxvr[2],2);
  const double dxvr2R3 = pow(dxvr[2],3);

  incr2[3] = -(1.0*((5.0*diffFac[1]*dxvl2R2*dxvr[2]+3.0*diffFac[1]*dxvl2R3)*fr[5]+((-3.0*diffFac[1]*dxvr2R3)-5.0*diffFac[1]*dxvl[2]*dxvr2R2)*fl[5]+(5.0*diffFac[0]*dxvl2R2*dxvr[2]+3.0*diffFac[0]*dxvl2R3)*fr[3]+((-3.0*diffFac[0]*dxvr2R3)-5.0*diffFac[0]*dxvl[2]*dxvr2R2)*fl[3]+((-1.732050807568877*diffFac[1]*fl[1])-1.732050807568877*diffFac[0]*fl[0])*dxvr2R3+((-5.196152422706631*diffFac[1]*fl[1])-5.196152422706631*diffFac[0]*fl[0])*dxvl[2]*dxvr2R2+((-5.196152422706631*diffFac[1]*fr[1])-5.196152422706631*diffFac[0]*fr[0])*dxvl2R2*dxvr[2]+((-1.732050807568877*diffFac[1]*fr[1])-1.732050807568877*diffFac[0]*fr[0])*dxvl2R3))/(2.828427124746191*dxvr2R3+8.485281374238571*dxvl[2]*dxvr2R2+8.485281374238571*dxvl2R2*dxvr[2]+2.828427124746191*dxvl2R3); 
  incr2[5] = -(1.0*((5.0*diffFac[0]*dxvl2R2*dxvr[2]+3.0*diffFac[0]*dxvl2R3)*fr[5]+((-3.0*diffFac[0]*dxvr2R3)-5.0*diffFac[0]*dxvl[2]*dxvr2R2)*fl[5]+(5.0*diffFac[1]*dxvl2R2*dxvr[2]+3.0*diffFac[1]*dxvl2R3)*fr[3]+((-3.0*diffFac[1]*dxvr2R3)-5.0*diffFac[1]*dxvl[2]*dxvr2R2)*fl[3]+((-1.732050807568877*diffFac[0]*fl[1])-1.732050807568877*fl[0]*diffFac[1])*dxvr2R3+((-5.196152422706631*diffFac[0]*fl[1])-5.196152422706631*fl[0]*diffFac[1])*dxvl[2]*dxvr2R2+((-5.196152422706631*diffFac[0]*fr[1])-5.196152422706631*fr[0]*diffFac[1])*dxvl2R2*dxvr[2]+((-1.732050807568877*diffFac[0]*fr[1])-1.732050807568877*fr[0]*diffFac[1])*dxvl2R3))/(2.828427124746191*dxvr2R3+8.485281374238571*dxvl[2]*dxvr2R2+8.485281374238571*dxvl2R2*dxvr[2]+2.828427124746191*dxvl2R3); 
  incr2[6] = -(1.0*((5.0*diffFac[1]*dxvl2R2*dxvr[2]+3.0*diffFac[1]*dxvl2R3)*fr[7]+((-3.0*diffFac[1]*dxvr2R3)-5.0*diffFac[1]*dxvl[2]*dxvr2R2)*fl[7]+(5.0*diffFac[0]*dxvl2R2*dxvr[2]+3.0*diffFac[0]*dxvl2R3)*fr[6]+((-3.0*diffFac[0]*dxvr2R3)-5.0*diffFac[0]*dxvl[2]*dxvr2R2)*fl[6]+((-5.196152422706631*diffFac[1]*dxvl2R2*dxvr[2])-1.732050807568877*diffFac[1]*dxvl2R3)*fr[4]+((-1.732050807568877*diffFac[1]*dxvr2R3)-5.196152422706631*diffFac[1]*dxvl[2]*dxvr2R2)*fl[4]+((-5.196152422706631*diffFac[0]*dxvl2R2*dxvr[2])-1.732050807568877*diffFac[0]*dxvl2R3)*fr[2]+((-1.732050807568877*diffFac[0]*dxvr2R3)-5.196152422706631*diffFac[0]*dxvl[2]*dxvr2R2)*fl[2]))/(2.828427124746191*dxvr2R3+8.485281374238571*dxvl[2]*dxvr2R2+8.485281374238571*dxvl2R2*dxvr[2]+2.828427124746191*dxvl2R3); 
  incr2[7] = -(1.0*((5.0*diffFac[0]*dxvl2R2*dxvr[2]+3.0*diffFac[0]*dxvl2R3)*fr[7]+((-3.0*diffFac[0]*dxvr2R3)-5.0*diffFac[0]*dxvl[2]*dxvr2R2)*fl[7]+(5.0*diffFac[1]*dxvl2R2*dxvr[2]+3.0*diffFac[1]*dxvl2R3)*fr[6]+((-3.0*diffFac[1]*dxvr2R3)-5.0*diffFac[1]*dxvl[2]*dxvr2R2)*fl[6]+((-5.196152422706631*diffFac[0]*dxvl2R2*dxvr[2])-1.732050807568877*diffFac[0]*dxvl2R3)*fr[4]+((-1.732050807568877*diffFac[0]*dxvr2R3)-5.196152422706631*diffFac[0]*dxvl[2]*dxvr2R2)*fl[4]+((-5.196152422706631*diffFac[1]*dxvl2R2*dxvr[2])-1.732050807568877*diffFac[1]*dxvl2R3)*fr[2]+((-1.732050807568877*diffFac[1]*dxvr2R3)-5.196152422706631*diffFac[1]*dxvl[2]*dxvr2R2)*fl[2]))/(2.828427124746191*dxvr2R3+8.485281374238571*dxvl[2]*dxvr2R2+8.485281374238571*dxvl2R2*dxvr[2]+2.828427124746191*dxvl2R3); 

  const double dxvl2R4 = pow(dxvl[2],4);
  const double dxvr2R4 = pow(dxvr[2],4);

  Gdiff[0] = -(1.0*((17.32050807568877*diffFac[1]*dxvl2R2*dxvr2R2+3.464101615137754*diffFac[1]*dxvl2R3*dxvr[2]-3.464101615137754*diffFac[1]*dxvl2R4)*fr[5]+((-3.464101615137754*diffFac[1]*dxvr2R4)+3.464101615137754*diffFac[1]*dxvl[2]*dxvr2R3+17.32050807568877*diffFac[1]*dxvl2R2*dxvr2R2)*fl[5]+(17.32050807568877*diffFac[0]*dxvl2R2*dxvr2R2+3.464101615137754*diffFac[0]*dxvl2R3*dxvr[2]-3.464101615137754*diffFac[0]*dxvl2R4)*fr[3]+((-3.464101615137754*diffFac[0]*dxvr2R4)+3.464101615137754*diffFac[0]*dxvl[2]*dxvr2R3+17.32050807568877*diffFac[0]*dxvl2R2*dxvr2R2)*fl[3]+((-18.0*diffFac[1]*fr[1])+18.0*diffFac[1]*fl[1]-18.0*diffFac[0]*fr[0]+18.0*diffFac[0]*fl[0])*dxvl2R2*dxvr2R2))/(1.414213562373095*dxvl[2]*dxvr2R4+4.242640687119286*dxvl2R2*dxvr2R3+4.242640687119286*dxvl2R3*dxvr2R2+1.414213562373095*dxvl2R4*dxvr[2]); 
  Gdiff[1] = -(1.0*((17.32050807568877*diffFac[0]*dxvl2R2*dxvr2R2+3.464101615137754*diffFac[0]*dxvl2R3*dxvr[2]-3.464101615137754*diffFac[0]*dxvl2R4)*fr[5]+((-3.464101615137754*diffFac[0]*dxvr2R4)+3.464101615137754*diffFac[0]*dxvl[2]*dxvr2R3+17.32050807568877*diffFac[0]*dxvl2R2*dxvr2R2)*fl[5]+(17.32050807568877*diffFac[1]*dxvl2R2*dxvr2R2+3.464101615137754*diffFac[1]*dxvl2R3*dxvr[2]-3.464101615137754*diffFac[1]*dxvl2R4)*fr[3]+((-3.464101615137754*diffFac[1]*dxvr2R4)+3.464101615137754*diffFac[1]*dxvl[2]*dxvr2R3+17.32050807568877*diffFac[1]*dxvl2R2*dxvr2R2)*fl[3]+((-18.0*diffFac[0]*fr[1])+18.0*diffFac[0]*fl[1]+(18.0*fl[0]-18.0*fr[0])*diffFac[1])*dxvl2R2*dxvr2R2))/(1.414213562373095*dxvl[2]*dxvr2R4+4.242640687119286*dxvl2R2*dxvr2R3+4.242640687119286*dxvl2R3*dxvr2R2+1.414213562373095*dxvl2R4*dxvr[2]); 
  Gdiff[2] = -(1.0*((17.32050807568877*diffFac[1]*dxvl2R2*dxvr2R2+3.464101615137754*diffFac[1]*dxvl2R3*dxvr[2]-3.464101615137754*diffFac[1]*dxvl2R4)*fr[7]+((-3.464101615137754*diffFac[1]*dxvr2R4)+3.464101615137754*diffFac[1]*dxvl[2]*dxvr2R3+17.32050807568877*diffFac[1]*dxvl2R2*dxvr2R2)*fl[7]+(17.32050807568877*diffFac[0]*dxvl2R2*dxvr2R2+3.464101615137754*diffFac[0]*dxvl2R3*dxvr[2]-3.464101615137754*diffFac[0]*dxvl2R4)*fr[6]+((-3.464101615137754*diffFac[0]*dxvr2R4)+3.464101615137754*diffFac[0]*dxvl[2]*dxvr2R3+17.32050807568877*diffFac[0]*dxvl2R2*dxvr2R2)*fl[6]-18.0*diffFac[1]*dxvl2R2*dxvr2R2*fr[4]+18.0*diffFac[1]*dxvl2R2*dxvr2R2*fl[4]-18.0*diffFac[0]*dxvl2R2*dxvr2R2*fr[2]+18.0*diffFac[0]*dxvl2R2*dxvr2R2*fl[2]))/(1.414213562373095*dxvl[2]*dxvr2R4+4.242640687119286*dxvl2R2*dxvr2R3+4.242640687119286*dxvl2R3*dxvr2R2+1.414213562373095*dxvl2R4*dxvr[2]); 
  Gdiff[4] = -(1.0*((17.32050807568877*diffFac[0]*dxvl2R2*dxvr2R2+3.464101615137754*diffFac[0]*dxvl2R3*dxvr[2]-3.464101615137754*diffFac[0]*dxvl2R4)*fr[7]+((-3.464101615137754*diffFac[0]*dxvr2R4)+3.464101615137754*diffFac[0]*dxvl[2]*dxvr2R3+17.32050807568877*diffFac[0]*dxvl2R2*dxvr2R2)*fl[7]+(17.32050807568877*diffFac[1]*dxvl2R2*dxvr2R2+3.464101615137754*diffFac[1]*dxvl2R3*dxvr[2]-3.464101615137754*diffFac[1]*dxvl2R4)*fr[6]+((-3.464101615137754*diffFac[1]*dxvr2R4)+3.464101615137754*diffFac[1]*dxvl[2]*dxvr2R3+17.32050807568877*diffFac[1]*dxvl2R2*dxvr2R2)*fl[6]-18.0*diffFac[0]*dxvl2R2*dxvr2R2*fr[4]+18.0*diffFac[0]*dxvl2R2*dxvr2R2*fl[4]-18.0*diffFac[1]*dxvl2R2*dxvr2R2*fr[2]+18.0*diffFac[1]*dxvl2R2*dxvr2R2*fl[2]))/(1.414213562373095*dxvl[2]*dxvr2R4+4.242640687119286*dxvl2R2*dxvr2R3+4.242640687119286*dxvl2R3*dxvr2R2+1.414213562373095*dxvl2R4*dxvr[2]); 

  Ghat[0] = 0.7071067811865475*alphaDrSurf[0]*fUpwind[0]+Gdiff[0]; 
  Ghat[1] = 0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+Gdiff[1]; 
  Ghat[2] = 0.7071067811865475*alphaDrSurf[0]*fUpwind[2]+Gdiff[2]; 
  Ghat[4] = Gdiff[4]+0.7071067811865475*alphaDrSurf[0]*fUpwind[3]; 

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = 0.8660254037844386*Ghat[2]; 
  incr1[7] = 0.8660254037844386*Ghat[4]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 

  return std::abs(2.0*wl[2]); 
} 
