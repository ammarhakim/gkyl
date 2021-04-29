#include <VmLBOModDecl.h> 
double VmLBOconstNuSurfNonUniform2x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[5]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 

  double fjump[5]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 2.0*wl[2]*nuSum+dxvl[2]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 

  double Ghat[5]; 
  for(unsigned short int i=0; i<5; ++i){ 
    Ghat[i]=0.0; 
  }; 

  const double dxvl2R2 = pow(dxvl[2],2);
  const double dxvl2R3 = pow(dxvl[2],3);
  const double dxvl2R4 = pow(dxvl[2],4);
  const double dxvr2R2 = pow(dxvr[2],2);
  const double dxvr2R3 = pow(dxvr[2],3);
  const double dxvr2R4 = pow(dxvr[2],4);

  Ghat[0] = (-(1.0*((8.660254037844386*nuVtSqSum[0]*dxvl2R2*dxvr2R2+1.732050807568877*nuVtSqSum[0]*dxvl2R3*dxvr[2]-1.732050807568877*nuVtSqSum[0]*dxvl2R4)*fr[3]+((-1.732050807568877*nuVtSqSum[0]*dxvr2R4)+1.732050807568877*nuVtSqSum[0]*dxvl[2]*dxvr2R3+8.660254037844386*nuVtSqSum[0]*dxvl2R2*dxvr2R2)*fl[3]+(9.0*dxvl2R2*dxvr2R2*fl[2]-9.0*dxvl2R2*dxvr2R2*fr[2])*nuVtSqSum[2]+((9.0*fl[1]-9.0*fr[1])*nuVtSqSum[1]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[0])*dxvl2R2*dxvr2R2))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))-0.5*(1.732050807568877*fjump[3]+fjump[0])+0.25*(1.732050807568877*alphaDrag[0]*favg[3]+alphaDrag[2]*favg[2]+alphaDrag[1]*favg[1]+alphaDrag[0]*favg[0]); 
  Ghat[1] = (-(1.0*((8.660254037844386*nuVtSqSum[1]*dxvl2R2*dxvr2R2+1.732050807568877*nuVtSqSum[1]*dxvl2R3*dxvr[2]-1.732050807568877*nuVtSqSum[1]*dxvl2R4)*fr[3]+((-1.732050807568877*nuVtSqSum[1]*dxvr2R4)+1.732050807568877*nuVtSqSum[1]*dxvl[2]*dxvr2R3+8.660254037844386*nuVtSqSum[1]*dxvl2R2*dxvr2R2)*fl[3]+((9.0*fl[0]-9.0*fr[0])*nuVtSqSum[1]-9.0*nuVtSqSum[0]*fr[1]+9.0*nuVtSqSum[0]*fl[1])*dxvl2R2*dxvr2R2))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))+0.25*(1.732050807568877*alphaDrag[1]*favg[3]+alphaDrag[0]*favg[1]+favg[0]*alphaDrag[1])-0.5*fjump[1]; 
  Ghat[2] = (-(1.0*((8.660254037844386*dxvl2R2*dxvr2R2+1.732050807568877*dxvl2R3*dxvr[2]-1.732050807568877*dxvl2R4)*nuVtSqSum[2]*fr[3]+((-1.732050807568877*dxvr2R4)+1.732050807568877*dxvl[2]*dxvr2R3+8.660254037844386*dxvl2R2*dxvr2R2)*nuVtSqSum[2]*fl[3]+(9.0*fl[0]-9.0*fr[0])*dxvl2R2*dxvr2R2*nuVtSqSum[2]-9.0*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fr[2]+9.0*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fl[2]))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))+0.25*(1.732050807568877*alphaDrag[2]*favg[3]+alphaDrag[0]*favg[2]+favg[0]*alphaDrag[2])-0.5*fjump[2]; 
  Ghat[4] = (9.0*nuVtSqSum[0]*dxvl[2]*dxvr[2]*fr[4]-9.0*nuVtSqSum[0]*dxvl[2]*dxvr[2]*fl[4])/(dxvr2R3+3.0*dxvl[2]*dxvr2R2+3.0*dxvl2R2*dxvr[2]+dxvl2R3)-0.5*fjump[4]+0.25*alphaDrag[0]*favg[4]; 

  double incr1[5]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 

  double incr2[5]; 

  incr2[3] = -(1.0*((5.0*nuVtSqSum[0]*dxvl2R2*dxvr[2]+3.0*nuVtSqSum[0]*dxvl2R3)*fr[3]+((-3.0*nuVtSqSum[0]*dxvr2R3)-5.0*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[3]+(((-5.196152422706631*dxvl2R2*dxvr[2])-1.732050807568877*dxvl2R3)*fr[2]+((-1.732050807568877*dxvr2R3)-5.196152422706631*dxvl[2]*dxvr2R2)*fl[2])*nuVtSqSum[2]+((-1.732050807568877*fl[1]*nuVtSqSum[1])-1.732050807568877*fl[0]*nuVtSqSum[0])*dxvr2R3+((-5.196152422706631*fl[1]*nuVtSqSum[1])-5.196152422706631*fl[0]*nuVtSqSum[0])*dxvl[2]*dxvr2R2+((-5.196152422706631*fr[1]*nuVtSqSum[1])-5.196152422706631*fr[0]*nuVtSqSum[0])*dxvl2R2*dxvr[2]+((-1.732050807568877*fr[1]*nuVtSqSum[1])-1.732050807568877*fr[0]*nuVtSqSum[0])*dxvl2R3))/(4.0*dxvr2R3+12.0*dxvl[2]*dxvr2R2+12.0*dxvl2R2*dxvr[2]+4.0*dxvl2R3); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 

  return std::abs(wl[2]-(0.5*sumNuUx[0])/nuSum); 
} 
double VmLBOconstNuSurfNonUniform2x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  const double *sumNuUy = &nuUSum[3]; 

  double favg[5]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 

  double fjump[5]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(-1*fr[4])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 2.0*wl[3]*nuSum+dxvl[3]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 

  double Ghat[5]; 
  for(unsigned short int i=0; i<5; ++i){ 
    Ghat[i]=0.0; 
  }; 

  const double dxvl3R2 = pow(dxvl[3],2);
  const double dxvl3R3 = pow(dxvl[3],3);
  const double dxvr3R2 = pow(dxvr[3],2);
  const double dxvr3R3 = pow(dxvr[3],3);

  Ghat[0] = (-(1.0*((8.660254037844386*nuVtSqSum[0]*dxvl[3]*dxvr3R2+1.732050807568877*nuVtSqSum[0]*dxvl3R2*dxvr[3]-1.732050807568877*nuVtSqSum[0]*dxvl3R3)*fr[4]+((-1.732050807568877*nuVtSqSum[0]*dxvr3R3)+1.732050807568877*nuVtSqSum[0]*dxvl[3]*dxvr3R2+8.660254037844386*nuVtSqSum[0]*dxvl3R2*dxvr[3])*fl[4]+((18.0*fl[2]-18.0*fr[2])*nuVtSqSum[2]+(18.0*fl[1]-18.0*fr[1])*nuVtSqSum[1]+(18.0*fl[0]-18.0*fr[0])*nuVtSqSum[0])*dxvl[3]*dxvr[3]))/(2.0*dxvr3R3+6.0*dxvl[3]*dxvr3R2+6.0*dxvl3R2*dxvr[3]+2.0*dxvl3R3))-0.5*(1.732050807568877*fjump[4]+fjump[0])+0.25*(1.732050807568877*alphaDrag[0]*favg[4]+alphaDrag[2]*favg[2]+alphaDrag[1]*favg[1]+alphaDrag[0]*favg[0]); 
  Ghat[1] = (-(1.0*((8.660254037844386*nuVtSqSum[1]*dxvl[3]*dxvr3R2+1.732050807568877*nuVtSqSum[1]*dxvl3R2*dxvr[3]-1.732050807568877*nuVtSqSum[1]*dxvl3R3)*fr[4]+((-1.732050807568877*nuVtSqSum[1]*dxvr3R3)+1.732050807568877*nuVtSqSum[1]*dxvl[3]*dxvr3R2+8.660254037844386*nuVtSqSum[1]*dxvl3R2*dxvr[3])*fl[4]+((18.0*fl[0]-18.0*fr[0])*nuVtSqSum[1]-18.0*nuVtSqSum[0]*fr[1]+18.0*nuVtSqSum[0]*fl[1])*dxvl[3]*dxvr[3]))/(2.0*dxvr3R3+6.0*dxvl[3]*dxvr3R2+6.0*dxvl3R2*dxvr[3]+2.0*dxvl3R3))+0.25*(1.732050807568877*alphaDrag[1]*favg[4]+alphaDrag[0]*favg[1]+favg[0]*alphaDrag[1])-0.5*fjump[1]; 
  Ghat[2] = (-(1.0*((8.660254037844386*nuVtSqSum[2]*dxvl[3]*dxvr3R2+1.732050807568877*nuVtSqSum[2]*dxvl3R2*dxvr[3]-1.732050807568877*nuVtSqSum[2]*dxvl3R3)*fr[4]+((-1.732050807568877*nuVtSqSum[2]*dxvr3R3)+1.732050807568877*nuVtSqSum[2]*dxvl[3]*dxvr3R2+8.660254037844386*nuVtSqSum[2]*dxvl3R2*dxvr[3])*fl[4]+((18.0*fl[0]-18.0*fr[0])*nuVtSqSum[2]-18.0*nuVtSqSum[0]*fr[2]+18.0*nuVtSqSum[0]*fl[2])*dxvl[3]*dxvr[3]))/(2.0*dxvr3R3+6.0*dxvl[3]*dxvr3R2+6.0*dxvl3R2*dxvr[3]+2.0*dxvl3R3))+0.25*(1.732050807568877*alphaDrag[2]*favg[4]+alphaDrag[0]*favg[2]+favg[0]*alphaDrag[2])-0.5*fjump[2]; 
  Ghat[3] = (9.0*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fr[3]-9.0*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fl[3])/(dxvr3R3+3.0*dxvl[3]*dxvr3R2+3.0*dxvl3R2*dxvr[3]+dxvl3R3)-0.5*fjump[3]+0.25*alphaDrag[0]*favg[3]; 

  double incr1[5]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[0]; 

  double incr2[5]; 

  incr2[4] = -(1.0*((5.0*nuVtSqSum[0]*dxvl3R2*dxvr3R2+3.0*nuVtSqSum[0]*dxvl3R3*dxvr[3])*fr[4]+((-3.0*nuVtSqSum[0]*dxvl[3]*dxvr3R3)-5.0*nuVtSqSum[0]*dxvl3R2*dxvr3R2)*fl[4]+((-3.464101615137754*fl[2]*nuVtSqSum[2])-3.464101615137754*fl[1]*nuVtSqSum[1]-3.464101615137754*fl[0]*nuVtSqSum[0])*dxvr3R3+((-10.39230484541326*fl[2]*nuVtSqSum[2])-10.39230484541326*fl[1]*nuVtSqSum[1]-10.39230484541326*fl[0]*nuVtSqSum[0])*dxvl[3]*dxvr3R2+((-10.39230484541326*fr[2]*nuVtSqSum[2])-10.39230484541326*fr[1]*nuVtSqSum[1]-10.39230484541326*fr[0]*nuVtSqSum[0])*dxvl3R2*dxvr[3]+((-3.464101615137754*fr[2]*nuVtSqSum[2])-3.464101615137754*fr[1]*nuVtSqSum[1]-3.464101615137754*fr[0]*nuVtSqSum[0])*dxvl3R3))/(8.0*dxvr3R3+24.0*dxvl[3]*dxvr3R2+24.0*dxvl3R2*dxvr[3]+8.0*dxvl3R3); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 

  return std::abs(wl[3]-(0.5*sumNuUy[0])/nuSum); 
} 
