#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x1vTensor_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:          Cell-center coordinates. 
  // dxv[2]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[2]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[4]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 

  double fjump[4]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 

  double Gdiff[4]; 
  double Ghat[4]; 
  double incr2[4]; 


  incr2[2] = -0.125*(2.828427124746191*nuVtSqSum[1]*fr[3]-2.828427124746191*nuVtSqSum[1]*fl[3]+2.828427124746191*nuVtSqSum[0]*fr[2]-2.828427124746191*nuVtSqSum[0]*fl[2]+((-2.449489742783178*fr[1])-2.449489742783178*fl[1])*nuVtSqSum[1]+((-2.449489742783178*fr[0])-2.449489742783178*fl[0])*nuVtSqSum[0]); 
  incr2[3] = -0.125*(2.828427124746191*nuVtSqSum[0]*fr[3]-2.828427124746191*nuVtSqSum[0]*fl[3]+2.828427124746191*nuVtSqSum[1]*fr[2]-2.828427124746191*nuVtSqSum[1]*fl[2]+((-2.449489742783178*fr[0])-2.449489742783178*fl[0])*nuVtSqSum[1]-2.449489742783178*nuVtSqSum[0]*fr[1]-2.449489742783178*nuVtSqSum[0]*fl[1]); 


  Gdiff[0] = -0.0625*(12.24744871391589*nuVtSqSum[1]*fr[3]+12.24744871391589*nuVtSqSum[1]*fl[3]+12.24744871391589*nuVtSqSum[0]*fr[2]+12.24744871391589*nuVtSqSum[0]*fl[2]+(12.72792206135786*fl[1]-12.72792206135786*fr[1])*nuVtSqSum[1]+(12.72792206135786*fl[0]-12.72792206135786*fr[0])*nuVtSqSum[0]); 
  Gdiff[1] = -0.0625*(12.24744871391589*nuVtSqSum[0]*fr[3]+12.24744871391589*nuVtSqSum[0]*fl[3]+12.24744871391589*nuVtSqSum[1]*fr[2]+12.24744871391589*nuVtSqSum[1]*fl[2]+(12.72792206135786*fl[0]-12.72792206135786*fr[0])*nuVtSqSum[1]-12.72792206135786*nuVtSqSum[0]*fr[1]+12.72792206135786*nuVtSqSum[0]*fl[1]); 

  Ghat[0] = Gdiff[0]*rdv2L+alphaDrag[1]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])-0.8660254037844386*fjump[2]+alphaDrag[0]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv2L-0.8660254037844386*fjump[3]+alphaDrag[0]*(0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+alphaDrag[1]*(0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[1]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 

  return std::abs(wl[1]-(0.7071067811865475*sumNuUx[0])/nuSum); 
} 
