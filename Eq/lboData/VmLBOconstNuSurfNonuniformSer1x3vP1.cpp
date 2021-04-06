#include <VmLBOModDecl.h> 
double VmLBOconstNuSurfNonuniform1x3vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[16]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = -1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 

  double fjump[16]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(-1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(-1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(-1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(-1*fr[15])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 

  double Ghat[16]; 
  for(unsigned short int i=0; i<16; ++i){ 
    Ghat[i]=0.0; 
  }; 

  const double dxvl1R2 = pow(dxvl[1],2);
  const double dxvl1R3 = pow(dxvl[1],3);
  const double dxvl1R4 = pow(dxvl[1],4);
  const double dxvr1R2 = pow(dxvr[1],2);
  const double dxvr1R3 = pow(dxvr[1],3);
  const double dxvr1R4 = pow(dxvr[1],4);

  Ghat[0] = (-(1.0*((12.24744871391589*dxvl1R2*dxvr1R2+2.449489742783178*dxvl1R3*dxvr[1]-2.449489742783178*dxvl1R4)*nuVtSqSum[1]*fr[5]+((-2.449489742783178*dxvr1R4)+2.449489742783178*dxvl[1]*dxvr1R3+12.24744871391589*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[5]+(12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2+2.449489742783178*nuVtSqSum[0]*dxvl1R3*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R4)*fr[2]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R4)+2.449489742783178*nuVtSqSum[0]*dxvl[1]*dxvr1R3+12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[2]+(12.72792206135786*dxvl1R2*dxvr1R2*fl[1]-12.72792206135786*dxvl1R2*dxvr1R2*fr[1])*nuVtSqSum[1]+(12.72792206135786*fl[0]-12.72792206135786*fr[0])*nuVtSqSum[0]*dxvl1R2*dxvr1R2))/(dxvl[1]*dxvr1R4+3.0*dxvl1R2*dxvr1R3+3.0*dxvl1R3*dxvr1R2+dxvl1R4*dxvr[1]))+0.25*(2.449489742783178*alphaDrag[1]*favg[5]+2.449489742783178*alphaDrag[0]*favg[2]+1.414213562373095*alphaDrag[1]*favg[1]+1.414213562373095*alphaDrag[0]*favg[0])-0.5*(1.732050807568877*fjump[2]+fjump[0]); 
  Ghat[1] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2+2.449489742783178*nuVtSqSum[0]*dxvl1R3*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R4)*fr[5]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R4)+2.449489742783178*nuVtSqSum[0]*dxvl[1]*dxvr1R3+12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[5]+(12.24744871391589*dxvl1R2*dxvr1R2+2.449489742783178*dxvl1R3*dxvr[1]-2.449489742783178*dxvl1R4)*nuVtSqSum[1]*fr[2]+((-2.449489742783178*dxvr1R4)+2.449489742783178*dxvl[1]*dxvr1R3+12.24744871391589*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[2]+(12.72792206135786*fl[0]-12.72792206135786*fr[0])*dxvl1R2*dxvr1R2*nuVtSqSum[1]-12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[1]+12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[1]))/(dxvl[1]*dxvr1R4+3.0*dxvl1R2*dxvr1R3+3.0*dxvl1R3*dxvr1R2+dxvl1R4*dxvr[1]))-0.5*(1.732050807568877*fjump[5]+fjump[1])+0.25*(2.449489742783178*alphaDrag[0]*favg[5]+2.449489742783178*alphaDrag[1]*favg[2]+1.414213562373095*alphaDrag[0]*favg[1]+1.414213562373095*favg[0]*alphaDrag[1]); 
  Ghat[3] = (-(1.0*((12.24744871391589*dxvl1R2*dxvr1R2+2.449489742783178*dxvl1R3*dxvr[1]-2.449489742783178*dxvl1R4)*nuVtSqSum[1]*fr[11]+((-2.449489742783178*dxvr1R4)+2.449489742783178*dxvl[1]*dxvr1R3+12.24744871391589*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[11]+(12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2+2.449489742783178*nuVtSqSum[0]*dxvl1R3*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R4)*fr[7]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R4)+2.449489742783178*nuVtSqSum[0]*dxvl[1]*dxvr1R3+12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[7]-12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fr[6]+12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fl[6]-12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[3]+12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[3]))/(dxvl[1]*dxvr1R4+3.0*dxvl1R2*dxvr1R3+3.0*dxvl1R3*dxvr1R2+dxvl1R4*dxvr[1]))+0.25*(2.449489742783178*alphaDrag[1]*favg[11]+2.449489742783178*alphaDrag[0]*favg[7]+1.414213562373095*alphaDrag[1]*favg[6]+1.414213562373095*alphaDrag[0]*favg[3])-0.5*(1.732050807568877*fjump[7]+fjump[3]); 
  Ghat[4] = (-(1.0*((12.24744871391589*dxvl1R2*dxvr1R2+2.449489742783178*dxvl1R3*dxvr[1]-2.449489742783178*dxvl1R4)*nuVtSqSum[1]*fr[12]+((-2.449489742783178*dxvr1R4)+2.449489742783178*dxvl[1]*dxvr1R3+12.24744871391589*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[12]+(12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2+2.449489742783178*nuVtSqSum[0]*dxvl1R3*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R4)*fr[9]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R4)+2.449489742783178*nuVtSqSum[0]*dxvl[1]*dxvr1R3+12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[9]-12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fr[8]+12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fl[8]-12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[4]+12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[4]))/(dxvl[1]*dxvr1R4+3.0*dxvl1R2*dxvr1R3+3.0*dxvl1R3*dxvr1R2+dxvl1R4*dxvr[1]))+0.25*(2.449489742783178*alphaDrag[1]*favg[12]+2.449489742783178*alphaDrag[0]*favg[9]+1.414213562373095*alphaDrag[1]*favg[8]+1.414213562373095*alphaDrag[0]*favg[4])-0.5*(1.732050807568877*fjump[9]+fjump[4]); 
  Ghat[6] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2+2.449489742783178*nuVtSqSum[0]*dxvl1R3*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R4)*fr[11]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R4)+2.449489742783178*nuVtSqSum[0]*dxvl[1]*dxvr1R3+12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[11]+(12.24744871391589*dxvl1R2*dxvr1R2+2.449489742783178*dxvl1R3*dxvr[1]-2.449489742783178*dxvl1R4)*nuVtSqSum[1]*fr[7]+((-2.449489742783178*dxvr1R4)+2.449489742783178*dxvl[1]*dxvr1R3+12.24744871391589*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[7]-12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[6]+12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[6]-12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fr[3]+12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fl[3]))/(dxvl[1]*dxvr1R4+3.0*dxvl1R2*dxvr1R3+3.0*dxvl1R3*dxvr1R2+dxvl1R4*dxvr[1]))-0.5*(1.732050807568877*fjump[11]+fjump[6])+0.25*(2.449489742783178*alphaDrag[0]*favg[11]+2.449489742783178*alphaDrag[1]*favg[7]+1.414213562373095*alphaDrag[0]*favg[6]+1.414213562373095*alphaDrag[1]*favg[3]); 
  Ghat[8] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2+2.449489742783178*nuVtSqSum[0]*dxvl1R3*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R4)*fr[12]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R4)+2.449489742783178*nuVtSqSum[0]*dxvl[1]*dxvr1R3+12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[12]+(12.24744871391589*dxvl1R2*dxvr1R2+2.449489742783178*dxvl1R3*dxvr[1]-2.449489742783178*dxvl1R4)*nuVtSqSum[1]*fr[9]+((-2.449489742783178*dxvr1R4)+2.449489742783178*dxvl[1]*dxvr1R3+12.24744871391589*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[9]-12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[8]+12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[8]-12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fr[4]+12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fl[4]))/(dxvl[1]*dxvr1R4+3.0*dxvl1R2*dxvr1R3+3.0*dxvl1R3*dxvr1R2+dxvl1R4*dxvr[1]))-0.5*(1.732050807568877*fjump[12]+fjump[8])+0.25*(2.449489742783178*alphaDrag[0]*favg[12]+2.449489742783178*alphaDrag[1]*favg[9]+1.414213562373095*alphaDrag[0]*favg[8]+1.414213562373095*alphaDrag[1]*favg[4]); 
  Ghat[10] = (-(1.0*((12.24744871391589*dxvl1R2*dxvr1R2+2.449489742783178*dxvl1R3*dxvr[1]-2.449489742783178*dxvl1R4)*nuVtSqSum[1]*fr[15]+((-2.449489742783178*dxvr1R4)+2.449489742783178*dxvl[1]*dxvr1R3+12.24744871391589*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[15]+(12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2+2.449489742783178*nuVtSqSum[0]*dxvl1R3*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R4)*fr[14]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R4)+2.449489742783178*nuVtSqSum[0]*dxvl[1]*dxvr1R3+12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[14]-12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fr[13]+12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fl[13]-12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[10]+12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[10]))/(dxvl[1]*dxvr1R4+3.0*dxvl1R2*dxvr1R3+3.0*dxvl1R3*dxvr1R2+dxvl1R4*dxvr[1]))+0.25*(2.449489742783178*alphaDrag[1]*favg[15]+2.449489742783178*alphaDrag[0]*favg[14]+1.414213562373095*alphaDrag[1]*favg[13]+1.414213562373095*alphaDrag[0]*favg[10])-0.5*(1.732050807568877*fjump[14]+fjump[10]); 
  Ghat[13] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2+2.449489742783178*nuVtSqSum[0]*dxvl1R3*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R4)*fr[15]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R4)+2.449489742783178*nuVtSqSum[0]*dxvl[1]*dxvr1R3+12.24744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R2)*fl[15]+(12.24744871391589*dxvl1R2*dxvr1R2+2.449489742783178*dxvl1R3*dxvr[1]-2.449489742783178*dxvl1R4)*nuVtSqSum[1]*fr[14]+((-2.449489742783178*dxvr1R4)+2.449489742783178*dxvl[1]*dxvr1R3+12.24744871391589*dxvl1R2*dxvr1R2)*nuVtSqSum[1]*fl[14]-12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[13]+12.72792206135786*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[13]-12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fr[10]+12.72792206135786*dxvl1R2*dxvr1R2*nuVtSqSum[1]*fl[10]))/(dxvl[1]*dxvr1R4+3.0*dxvl1R2*dxvr1R3+3.0*dxvl1R3*dxvr1R2+dxvl1R4*dxvr[1]))-0.5*(1.732050807568877*fjump[15]+fjump[13])+0.25*(2.449489742783178*alphaDrag[0]*favg[15]+2.449489742783178*alphaDrag[1]*favg[14]+1.414213562373095*alphaDrag[0]*favg[13]+1.414213562373095*alphaDrag[1]*favg[10]); 

  double incr1[16]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = -0.5*Ghat[6]; 
  incr1[7] = 0.8660254037844386*Ghat[3]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = 0.8660254037844386*Ghat[4]; 
  incr1[10] = -0.5*Ghat[10]; 
  incr1[11] = 0.8660254037844386*Ghat[6]; 
  incr1[12] = 0.8660254037844386*Ghat[8]; 
  incr1[13] = -0.5*Ghat[13]; 
  incr1[14] = 0.8660254037844386*Ghat[10]; 
  incr1[15] = 0.8660254037844386*Ghat[13]; 

  double incr2[16]; 

  incr2[2] = -(1.0*((7.071067811865476*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)*nuVtSqSum[1]*fr[5]+((-4.242640687119286*dxvr1R3)-7.071067811865476*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[5]+(7.071067811865476*nuVtSqSum[0]*dxvl1R2*dxvr[1]+4.242640687119286*nuVtSqSum[0]*dxvl1R3)*fr[2]+((-4.242640687119286*nuVtSqSum[0]*dxvr1R3)-7.071067811865476*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[2]+(((-7.348469228349534*dxvl1R2*dxvr[1])-2.449489742783178*dxvl1R3)*fr[1]+((-2.449489742783178*dxvr1R3)-7.348469228349534*dxvl[1]*dxvr1R2)*fl[1])*nuVtSqSum[1]-2.449489742783178*fl[0]*nuVtSqSum[0]*dxvr1R3-7.348469228349534*fl[0]*nuVtSqSum[0]*dxvl[1]*dxvr1R2-7.348469228349534*fr[0]*nuVtSqSum[0]*dxvl1R2*dxvr[1]-2.449489742783178*fr[0]*nuVtSqSum[0]*dxvl1R3))/(4.0*dxvr1R3+12.0*dxvl[1]*dxvr1R2+12.0*dxvl1R2*dxvr[1]+4.0*dxvl1R3); 
  incr2[5] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl1R2*dxvr[1]+4.242640687119286*nuVtSqSum[0]*dxvl1R3)*fr[5]+((-4.242640687119286*nuVtSqSum[0]*dxvr1R3)-7.071067811865476*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[5]+(7.071067811865476*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)*nuVtSqSum[1]*fr[2]+((-4.242640687119286*dxvr1R3)-7.071067811865476*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[2]+((-2.449489742783178*fl[0]*dxvr1R3)-7.348469228349534*fl[0]*dxvl[1]*dxvr1R2-7.348469228349534*fr[0]*dxvl1R2*dxvr[1]-2.449489742783178*fr[0]*dxvl1R3)*nuVtSqSum[1]+((-7.348469228349534*nuVtSqSum[0]*dxvl1R2*dxvr[1])-2.449489742783178*nuVtSqSum[0]*dxvl1R3)*fr[1]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R3)-7.348469228349534*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[1]))/(4.0*dxvr1R3+12.0*dxvl[1]*dxvr1R2+12.0*dxvl1R2*dxvr[1]+4.0*dxvl1R3); 
  incr2[7] = -(1.0*((7.071067811865476*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)*nuVtSqSum[1]*fr[11]+((-4.242640687119286*dxvr1R3)-7.071067811865476*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[11]+(7.071067811865476*nuVtSqSum[0]*dxvl1R2*dxvr[1]+4.242640687119286*nuVtSqSum[0]*dxvl1R3)*fr[7]+((-4.242640687119286*nuVtSqSum[0]*dxvr1R3)-7.071067811865476*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[7]+((-7.348469228349534*dxvl1R2*dxvr[1])-2.449489742783178*dxvl1R3)*nuVtSqSum[1]*fr[6]+((-2.449489742783178*dxvr1R3)-7.348469228349534*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[6]+((-7.348469228349534*nuVtSqSum[0]*dxvl1R2*dxvr[1])-2.449489742783178*nuVtSqSum[0]*dxvl1R3)*fr[3]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R3)-7.348469228349534*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[3]))/(4.0*dxvr1R3+12.0*dxvl[1]*dxvr1R2+12.0*dxvl1R2*dxvr[1]+4.0*dxvl1R3); 
  incr2[9] = -(1.0*((7.071067811865476*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)*nuVtSqSum[1]*fr[12]+((-4.242640687119286*dxvr1R3)-7.071067811865476*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[12]+(7.071067811865476*nuVtSqSum[0]*dxvl1R2*dxvr[1]+4.242640687119286*nuVtSqSum[0]*dxvl1R3)*fr[9]+((-4.242640687119286*nuVtSqSum[0]*dxvr1R3)-7.071067811865476*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[9]+((-7.348469228349534*dxvl1R2*dxvr[1])-2.449489742783178*dxvl1R3)*nuVtSqSum[1]*fr[8]+((-2.449489742783178*dxvr1R3)-7.348469228349534*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[8]+((-7.348469228349534*nuVtSqSum[0]*dxvl1R2*dxvr[1])-2.449489742783178*nuVtSqSum[0]*dxvl1R3)*fr[4]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R3)-7.348469228349534*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[4]))/(4.0*dxvr1R3+12.0*dxvl[1]*dxvr1R2+12.0*dxvl1R2*dxvr[1]+4.0*dxvl1R3); 
  incr2[11] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl1R2*dxvr[1]+4.242640687119286*nuVtSqSum[0]*dxvl1R3)*fr[11]+((-4.242640687119286*nuVtSqSum[0]*dxvr1R3)-7.071067811865476*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[11]+(7.071067811865476*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)*nuVtSqSum[1]*fr[7]+((-4.242640687119286*dxvr1R3)-7.071067811865476*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[7]+((-7.348469228349534*nuVtSqSum[0]*dxvl1R2*dxvr[1])-2.449489742783178*nuVtSqSum[0]*dxvl1R3)*fr[6]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R3)-7.348469228349534*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[6]+((-7.348469228349534*dxvl1R2*dxvr[1])-2.449489742783178*dxvl1R3)*nuVtSqSum[1]*fr[3]+((-2.449489742783178*dxvr1R3)-7.348469228349534*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[3]))/(4.0*dxvr1R3+12.0*dxvl[1]*dxvr1R2+12.0*dxvl1R2*dxvr[1]+4.0*dxvl1R3); 
  incr2[12] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl1R2*dxvr[1]+4.242640687119286*nuVtSqSum[0]*dxvl1R3)*fr[12]+((-4.242640687119286*nuVtSqSum[0]*dxvr1R3)-7.071067811865476*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[12]+(7.071067811865476*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)*nuVtSqSum[1]*fr[9]+((-4.242640687119286*dxvr1R3)-7.071067811865476*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[9]+((-7.348469228349534*nuVtSqSum[0]*dxvl1R2*dxvr[1])-2.449489742783178*nuVtSqSum[0]*dxvl1R3)*fr[8]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R3)-7.348469228349534*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[8]+((-7.348469228349534*dxvl1R2*dxvr[1])-2.449489742783178*dxvl1R3)*nuVtSqSum[1]*fr[4]+((-2.449489742783178*dxvr1R3)-7.348469228349534*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[4]))/(4.0*dxvr1R3+12.0*dxvl[1]*dxvr1R2+12.0*dxvl1R2*dxvr[1]+4.0*dxvl1R3); 
  incr2[14] = -(1.0*((7.071067811865476*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)*nuVtSqSum[1]*fr[15]+((-4.242640687119286*dxvr1R3)-7.071067811865476*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[15]+(7.071067811865476*nuVtSqSum[0]*dxvl1R2*dxvr[1]+4.242640687119286*nuVtSqSum[0]*dxvl1R3)*fr[14]+((-4.242640687119286*nuVtSqSum[0]*dxvr1R3)-7.071067811865476*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[14]+((-7.348469228349534*dxvl1R2*dxvr[1])-2.449489742783178*dxvl1R3)*nuVtSqSum[1]*fr[13]+((-2.449489742783178*dxvr1R3)-7.348469228349534*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[13]+((-7.348469228349534*nuVtSqSum[0]*dxvl1R2*dxvr[1])-2.449489742783178*nuVtSqSum[0]*dxvl1R3)*fr[10]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R3)-7.348469228349534*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[10]))/(4.0*dxvr1R3+12.0*dxvl[1]*dxvr1R2+12.0*dxvl1R2*dxvr[1]+4.0*dxvl1R3); 
  incr2[15] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl1R2*dxvr[1]+4.242640687119286*nuVtSqSum[0]*dxvl1R3)*fr[15]+((-4.242640687119286*nuVtSqSum[0]*dxvr1R3)-7.071067811865476*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[15]+(7.071067811865476*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)*nuVtSqSum[1]*fr[14]+((-4.242640687119286*dxvr1R3)-7.071067811865476*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[14]+((-7.348469228349534*nuVtSqSum[0]*dxvl1R2*dxvr[1])-2.449489742783178*nuVtSqSum[0]*dxvl1R3)*fr[13]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R3)-7.348469228349534*nuVtSqSum[0]*dxvl[1]*dxvr1R2)*fl[13]+((-7.348469228349534*dxvl1R2*dxvr[1])-2.449489742783178*dxvl1R3)*nuVtSqSum[1]*fr[10]+((-2.449489742783178*dxvr1R3)-7.348469228349534*dxvl[1]*dxvr1R2)*nuVtSqSum[1]*fl[10]))/(4.0*dxvr1R3+12.0*dxvl[1]*dxvr1R2+12.0*dxvl1R2*dxvr[1]+4.0*dxvl1R3); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr2[9]*rdvSq4R+incr1[9]*rdv2R; 
  outr[10] += incr1[10]*rdv2R; 
  outr[11] += incr2[11]*rdvSq4R+incr1[11]*rdv2R; 
  outr[12] += incr2[12]*rdvSq4R+incr1[12]*rdv2R; 
  outr[13] += incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += incr1[9]*rdv2L-1.0*incr2[9]*rdvSq4L; 
  outl[10] += -1.0*incr1[10]*rdv2L; 
  outl[11] += incr1[11]*rdv2L-1.0*incr2[11]*rdvSq4L; 
  outl[12] += incr1[12]*rdv2L-1.0*incr2[12]*rdvSq4L; 
  outl[13] += -1.0*incr1[13]*rdv2L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += incr1[15]*rdv2L-1.0*incr2[15]*rdvSq4L; 

  return std::abs(wl[1]-(0.7071067811865475*sumNuUx[0])/nuSum); 
} 
double VmLBOconstNuSurfNonuniform1x3vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUy = &nuUSum[2]; 

  double favg[16]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 

  double fjump[16]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(-1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(-1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(-1*fr[15])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[2]*nuSum+0.7071067811865475*dxvl[2]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 

  double Ghat[16]; 
  for(unsigned short int i=0; i<16; ++i){ 
    Ghat[i]=0.0; 
  }; 

  const double dxvl2R2 = pow(dxvl[2],2);
  const double dxvl2R3 = pow(dxvl[2],3);
  const double dxvl2R4 = pow(dxvl[2],4);
  const double dxvr2R2 = pow(dxvr[2],2);
  const double dxvr2R3 = pow(dxvr[2],3);
  const double dxvr2R4 = pow(dxvr[2],4);

  Ghat[0] = (-(1.0*((12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[1]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R4)*fr[6]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R4)+2.449489742783178*nuVtSqSum[1]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2)*fl[6]+(12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[0]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R4)*fr[3]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R4)+2.449489742783178*nuVtSqSum[0]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2)*fl[3]+((12.72792206135786*fl[1]-12.72792206135786*fr[1])*nuVtSqSum[1]+(12.72792206135786*fl[0]-12.72792206135786*fr[0])*nuVtSqSum[0])*dxvl2R2*dxvr2R2))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))+0.25*(2.449489742783178*alphaDrag[1]*favg[6]+2.449489742783178*alphaDrag[0]*favg[3]+1.414213562373095*alphaDrag[1]*favg[1]+1.414213562373095*alphaDrag[0]*favg[0])-0.5*(1.732050807568877*fjump[3]+fjump[0]); 
  Ghat[1] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[0]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R4)*fr[6]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R4)+2.449489742783178*nuVtSqSum[0]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2)*fl[6]+(12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[1]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R4)*fr[3]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R4)+2.449489742783178*nuVtSqSum[1]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2)*fl[3]+((12.72792206135786*fl[0]-12.72792206135786*fr[0])*nuVtSqSum[1]-12.72792206135786*nuVtSqSum[0]*fr[1]+12.72792206135786*nuVtSqSum[0]*fl[1])*dxvl2R2*dxvr2R2))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))-0.5*(1.732050807568877*fjump[6]+fjump[1])+0.25*(2.449489742783178*alphaDrag[0]*favg[6]+2.449489742783178*alphaDrag[1]*favg[3]+1.414213562373095*alphaDrag[0]*favg[1]+1.414213562373095*favg[0]*alphaDrag[1]); 
  Ghat[2] = (-(1.0*((12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[1]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R4)*fr[11]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R4)+2.449489742783178*nuVtSqSum[1]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2)*fl[11]+(12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[0]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R4)*fr[7]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R4)+2.449489742783178*nuVtSqSum[0]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2)*fl[7]-12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fr[5]+12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fl[5]-12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fr[2]+12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fl[2]))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))+0.25*(2.449489742783178*alphaDrag[1]*favg[11]+2.449489742783178*alphaDrag[0]*favg[7]+1.414213562373095*alphaDrag[1]*favg[5]+1.414213562373095*alphaDrag[0]*favg[2])-0.5*(1.732050807568877*fjump[7]+fjump[2]); 
  Ghat[4] = (-(1.0*((12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[1]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R4)*fr[13]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R4)+2.449489742783178*nuVtSqSum[1]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2)*fl[13]+(12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[0]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R4)*fr[10]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R4)+2.449489742783178*nuVtSqSum[0]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2)*fl[10]-12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fr[8]+12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fl[8]-12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fr[4]+12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fl[4]))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))+0.25*(2.449489742783178*alphaDrag[1]*favg[13]+2.449489742783178*alphaDrag[0]*favg[10]+1.414213562373095*alphaDrag[1]*favg[8]+1.414213562373095*alphaDrag[0]*favg[4])-0.5*(1.732050807568877*fjump[10]+fjump[4]); 
  Ghat[5] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[0]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R4)*fr[11]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R4)+2.449489742783178*nuVtSqSum[0]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2)*fl[11]+(12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[1]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R4)*fr[7]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R4)+2.449489742783178*nuVtSqSum[1]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2)*fl[7]-12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fr[5]+12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fl[5]-12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fr[2]+12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fl[2]))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))-0.5*(1.732050807568877*fjump[11]+fjump[5])+0.25*(2.449489742783178*alphaDrag[0]*favg[11]+2.449489742783178*alphaDrag[1]*favg[7]+1.414213562373095*alphaDrag[0]*favg[5]+1.414213562373095*alphaDrag[1]*favg[2]); 
  Ghat[8] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[0]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R4)*fr[13]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R4)+2.449489742783178*nuVtSqSum[0]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2)*fl[13]+(12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[1]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R4)*fr[10]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R4)+2.449489742783178*nuVtSqSum[1]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2)*fl[10]-12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fr[8]+12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fl[8]-12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fr[4]+12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fl[4]))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))-0.5*(1.732050807568877*fjump[13]+fjump[8])+0.25*(2.449489742783178*alphaDrag[0]*favg[13]+2.449489742783178*alphaDrag[1]*favg[10]+1.414213562373095*alphaDrag[0]*favg[8]+1.414213562373095*alphaDrag[1]*favg[4]); 
  Ghat[9] = (-(1.0*((12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[1]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R4)*fr[15]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R4)+2.449489742783178*nuVtSqSum[1]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2)*fl[15]+(12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[0]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R4)*fr[14]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R4)+2.449489742783178*nuVtSqSum[0]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2)*fl[14]-12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fr[12]+12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fl[12]-12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fr[9]+12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fl[9]))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))+0.25*(2.449489742783178*alphaDrag[1]*favg[15]+2.449489742783178*alphaDrag[0]*favg[14]+1.414213562373095*alphaDrag[1]*favg[12]+1.414213562373095*alphaDrag[0]*favg[9])-0.5*(1.732050807568877*fjump[14]+fjump[9]); 
  Ghat[12] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[0]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R4)*fr[15]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R4)+2.449489742783178*nuVtSqSum[0]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[0]*dxvl2R2*dxvr2R2)*fl[15]+(12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2+2.449489742783178*nuVtSqSum[1]*dxvl2R3*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R4)*fr[14]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R4)+2.449489742783178*nuVtSqSum[1]*dxvl[2]*dxvr2R3+12.24744871391589*nuVtSqSum[1]*dxvl2R2*dxvr2R2)*fl[14]-12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fr[12]+12.72792206135786*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fl[12]-12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fr[9]+12.72792206135786*nuVtSqSum[1]*dxvl2R2*dxvr2R2*fl[9]))/(dxvl[2]*dxvr2R4+3.0*dxvl2R2*dxvr2R3+3.0*dxvl2R3*dxvr2R2+dxvl2R4*dxvr[2]))-0.5*(1.732050807568877*fjump[15]+fjump[12])+0.25*(2.449489742783178*alphaDrag[0]*favg[15]+2.449489742783178*alphaDrag[1]*favg[14]+1.414213562373095*alphaDrag[0]*favg[12]+1.414213562373095*alphaDrag[1]*favg[9]); 

  double incr1[16]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[1]; 
  incr1[7] = 0.8660254037844386*Ghat[2]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = -0.5*Ghat[9]; 
  incr1[10] = 0.8660254037844386*Ghat[4]; 
  incr1[11] = 0.8660254037844386*Ghat[5]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = 0.8660254037844386*Ghat[8]; 
  incr1[14] = 0.8660254037844386*Ghat[9]; 
  incr1[15] = 0.8660254037844386*Ghat[12]; 

  double incr2[16]; 

  incr2[3] = -(1.0*((7.071067811865476*nuVtSqSum[1]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[1]*dxvl2R3)*fr[6]+((-4.242640687119286*nuVtSqSum[1]*dxvr2R3)-7.071067811865476*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[6]+(7.071067811865476*nuVtSqSum[0]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[0]*dxvl2R3)*fr[3]+((-4.242640687119286*nuVtSqSum[0]*dxvr2R3)-7.071067811865476*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[3]+((-2.449489742783178*fl[1]*nuVtSqSum[1])-2.449489742783178*fl[0]*nuVtSqSum[0])*dxvr2R3+((-7.348469228349534*fl[1]*nuVtSqSum[1])-7.348469228349534*fl[0]*nuVtSqSum[0])*dxvl[2]*dxvr2R2+((-7.348469228349534*fr[1]*nuVtSqSum[1])-7.348469228349534*fr[0]*nuVtSqSum[0])*dxvl2R2*dxvr[2]+((-2.449489742783178*fr[1]*nuVtSqSum[1])-2.449489742783178*fr[0]*nuVtSqSum[0])*dxvl2R3))/(4.0*dxvr2R3+12.0*dxvl[2]*dxvr2R2+12.0*dxvl2R2*dxvr[2]+4.0*dxvl2R3); 
  incr2[6] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[0]*dxvl2R3)*fr[6]+((-4.242640687119286*nuVtSqSum[0]*dxvr2R3)-7.071067811865476*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[6]+(7.071067811865476*nuVtSqSum[1]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[1]*dxvl2R3)*fr[3]+((-4.242640687119286*nuVtSqSum[1]*dxvr2R3)-7.071067811865476*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[3]+((-2.449489742783178*fl[0]*nuVtSqSum[1])-2.449489742783178*nuVtSqSum[0]*fl[1])*dxvr2R3+((-7.348469228349534*fl[0]*nuVtSqSum[1])-7.348469228349534*nuVtSqSum[0]*fl[1])*dxvl[2]*dxvr2R2+((-7.348469228349534*fr[0]*nuVtSqSum[1])-7.348469228349534*nuVtSqSum[0]*fr[1])*dxvl2R2*dxvr[2]+((-2.449489742783178*fr[0]*nuVtSqSum[1])-2.449489742783178*nuVtSqSum[0]*fr[1])*dxvl2R3))/(4.0*dxvr2R3+12.0*dxvl[2]*dxvr2R2+12.0*dxvl2R2*dxvr[2]+4.0*dxvl2R3); 
  incr2[7] = -(1.0*((7.071067811865476*nuVtSqSum[1]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[1]*dxvl2R3)*fr[11]+((-4.242640687119286*nuVtSqSum[1]*dxvr2R3)-7.071067811865476*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[11]+(7.071067811865476*nuVtSqSum[0]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[0]*dxvl2R3)*fr[7]+((-4.242640687119286*nuVtSqSum[0]*dxvr2R3)-7.071067811865476*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[7]+((-7.348469228349534*nuVtSqSum[1]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[1]*dxvl2R3)*fr[5]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R3)-7.348469228349534*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[5]+((-7.348469228349534*nuVtSqSum[0]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[0]*dxvl2R3)*fr[2]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R3)-7.348469228349534*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[2]))/(4.0*dxvr2R3+12.0*dxvl[2]*dxvr2R2+12.0*dxvl2R2*dxvr[2]+4.0*dxvl2R3); 
  incr2[10] = -(1.0*((7.071067811865476*nuVtSqSum[1]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[1]*dxvl2R3)*fr[13]+((-4.242640687119286*nuVtSqSum[1]*dxvr2R3)-7.071067811865476*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[13]+(7.071067811865476*nuVtSqSum[0]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[0]*dxvl2R3)*fr[10]+((-4.242640687119286*nuVtSqSum[0]*dxvr2R3)-7.071067811865476*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[10]+((-7.348469228349534*nuVtSqSum[1]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[1]*dxvl2R3)*fr[8]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R3)-7.348469228349534*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[8]+((-7.348469228349534*nuVtSqSum[0]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[0]*dxvl2R3)*fr[4]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R3)-7.348469228349534*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[4]))/(4.0*dxvr2R3+12.0*dxvl[2]*dxvr2R2+12.0*dxvl2R2*dxvr[2]+4.0*dxvl2R3); 
  incr2[11] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[0]*dxvl2R3)*fr[11]+((-4.242640687119286*nuVtSqSum[0]*dxvr2R3)-7.071067811865476*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[11]+(7.071067811865476*nuVtSqSum[1]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[1]*dxvl2R3)*fr[7]+((-4.242640687119286*nuVtSqSum[1]*dxvr2R3)-7.071067811865476*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[7]+((-7.348469228349534*nuVtSqSum[0]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[0]*dxvl2R3)*fr[5]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R3)-7.348469228349534*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[5]+((-7.348469228349534*nuVtSqSum[1]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[1]*dxvl2R3)*fr[2]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R3)-7.348469228349534*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[2]))/(4.0*dxvr2R3+12.0*dxvl[2]*dxvr2R2+12.0*dxvl2R2*dxvr[2]+4.0*dxvl2R3); 
  incr2[13] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[0]*dxvl2R3)*fr[13]+((-4.242640687119286*nuVtSqSum[0]*dxvr2R3)-7.071067811865476*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[13]+(7.071067811865476*nuVtSqSum[1]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[1]*dxvl2R3)*fr[10]+((-4.242640687119286*nuVtSqSum[1]*dxvr2R3)-7.071067811865476*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[10]+((-7.348469228349534*nuVtSqSum[0]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[0]*dxvl2R3)*fr[8]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R3)-7.348469228349534*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[8]+((-7.348469228349534*nuVtSqSum[1]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[1]*dxvl2R3)*fr[4]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R3)-7.348469228349534*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[4]))/(4.0*dxvr2R3+12.0*dxvl[2]*dxvr2R2+12.0*dxvl2R2*dxvr[2]+4.0*dxvl2R3); 
  incr2[14] = -(1.0*((7.071067811865476*nuVtSqSum[1]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[1]*dxvl2R3)*fr[15]+((-4.242640687119286*nuVtSqSum[1]*dxvr2R3)-7.071067811865476*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[15]+(7.071067811865476*nuVtSqSum[0]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[0]*dxvl2R3)*fr[14]+((-4.242640687119286*nuVtSqSum[0]*dxvr2R3)-7.071067811865476*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[14]+((-7.348469228349534*nuVtSqSum[1]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[1]*dxvl2R3)*fr[12]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R3)-7.348469228349534*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[12]+((-7.348469228349534*nuVtSqSum[0]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[0]*dxvl2R3)*fr[9]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R3)-7.348469228349534*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[9]))/(4.0*dxvr2R3+12.0*dxvl[2]*dxvr2R2+12.0*dxvl2R2*dxvr[2]+4.0*dxvl2R3); 
  incr2[15] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[0]*dxvl2R3)*fr[15]+((-4.242640687119286*nuVtSqSum[0]*dxvr2R3)-7.071067811865476*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[15]+(7.071067811865476*nuVtSqSum[1]*dxvl2R2*dxvr[2]+4.242640687119286*nuVtSqSum[1]*dxvl2R3)*fr[14]+((-4.242640687119286*nuVtSqSum[1]*dxvr2R3)-7.071067811865476*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[14]+((-7.348469228349534*nuVtSqSum[0]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[0]*dxvl2R3)*fr[12]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R3)-7.348469228349534*nuVtSqSum[0]*dxvl[2]*dxvr2R2)*fl[12]+((-7.348469228349534*nuVtSqSum[1]*dxvl2R2*dxvr[2])-2.449489742783178*nuVtSqSum[1]*dxvl2R3)*fr[9]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R3)-7.348469228349534*nuVtSqSum[1]*dxvl[2]*dxvr2R2)*fl[9]))/(4.0*dxvr2R3+12.0*dxvl[2]*dxvr2R2+12.0*dxvl2R2*dxvr[2]+4.0*dxvl2R3); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr1[9]*rdv2R; 
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr2[11]*rdvSq4R+incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += -1.0*incr1[9]*rdv2L; 
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += incr1[11]*rdv2L-1.0*incr2[11]*rdvSq4L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += incr1[13]*rdv2L-1.0*incr2[13]*rdvSq4L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += incr1[15]*rdv2L-1.0*incr2[15]*rdvSq4L; 

  return std::abs(wl[2]-(0.7071067811865475*sumNuUy[0])/nuSum); 
} 
double VmLBOconstNuSurfNonuniform1x3vSer_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  const double *sumNuUz = &nuUSum[4]; 

  double favg[16]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = -1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = -1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 

  double fjump[16]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(-1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(-1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(-1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(-1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(-1*fr[15])); 

  double alphaDrag[2]; 
  alphaDrag[0] = 1.414213562373095*wl[3]*nuSum+0.7071067811865475*dxvl[3]*nuSum-1.0*sumNuUz[0]; 
  alphaDrag[1] = -1.0*sumNuUz[1]; 

  double Ghat[16]; 
  for(unsigned short int i=0; i<16; ++i){ 
    Ghat[i]=0.0; 
  }; 

  const double dxvl3R2 = pow(dxvl[3],2);
  const double dxvl3R3 = pow(dxvl[3],3);
  const double dxvr3R2 = pow(dxvr[3],2);
  const double dxvr3R3 = pow(dxvr[3],3);

  Ghat[0] = (-(1.0*((12.24744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[1]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[1]*dxvl3R3)*fr[8]+((-2.449489742783178*nuVtSqSum[1]*dxvr3R3)+2.449489742783178*nuVtSqSum[1]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[1]*dxvl3R2*dxvr[3])*fl[8]+(12.24744871391589*nuVtSqSum[0]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[0]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[0]*dxvl3R3)*fr[4]+((-2.449489742783178*nuVtSqSum[0]*dxvr3R3)+2.449489742783178*nuVtSqSum[0]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[0]*dxvl3R2*dxvr[3])*fl[4]+((25.45584412271572*fl[1]-25.45584412271572*fr[1])*nuVtSqSum[1]+(25.45584412271572*fl[0]-25.45584412271572*fr[0])*nuVtSqSum[0])*dxvl[3]*dxvr[3]))/(2.0*dxvr3R3+6.0*dxvl[3]*dxvr3R2+6.0*dxvl3R2*dxvr[3]+2.0*dxvl3R3))+0.25*(2.449489742783178*alphaDrag[1]*favg[8]+2.449489742783178*alphaDrag[0]*favg[4]+1.414213562373095*alphaDrag[1]*favg[1]+1.414213562373095*alphaDrag[0]*favg[0])-0.5*(1.732050807568877*fjump[4]+fjump[0]); 
  Ghat[1] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[0]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[0]*dxvl3R3)*fr[8]+((-2.449489742783178*nuVtSqSum[0]*dxvr3R3)+2.449489742783178*nuVtSqSum[0]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[0]*dxvl3R2*dxvr[3])*fl[8]+(12.24744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[1]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[1]*dxvl3R3)*fr[4]+((-2.449489742783178*nuVtSqSum[1]*dxvr3R3)+2.449489742783178*nuVtSqSum[1]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[1]*dxvl3R2*dxvr[3])*fl[4]+((25.45584412271572*fl[0]-25.45584412271572*fr[0])*nuVtSqSum[1]-25.45584412271572*nuVtSqSum[0]*fr[1]+25.45584412271572*nuVtSqSum[0]*fl[1])*dxvl[3]*dxvr[3]))/(2.0*dxvr3R3+6.0*dxvl[3]*dxvr3R2+6.0*dxvl3R2*dxvr[3]+2.0*dxvl3R3))-0.5*(1.732050807568877*fjump[8]+fjump[1])+0.25*(2.449489742783178*alphaDrag[0]*favg[8]+2.449489742783178*alphaDrag[1]*favg[4]+1.414213562373095*alphaDrag[0]*favg[1]+1.414213562373095*favg[0]*alphaDrag[1]); 
  Ghat[2] = (-(1.0*((12.24744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[1]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[1]*dxvl3R3)*fr[12]+((-2.449489742783178*nuVtSqSum[1]*dxvr3R3)+2.449489742783178*nuVtSqSum[1]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[1]*dxvl3R2*dxvr[3])*fl[12]+(12.24744871391589*nuVtSqSum[0]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[0]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[0]*dxvl3R3)*fr[9]+((-2.449489742783178*nuVtSqSum[0]*dxvr3R3)+2.449489742783178*nuVtSqSum[0]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[0]*dxvl3R2*dxvr[3])*fl[9]-25.45584412271572*nuVtSqSum[1]*dxvl[3]*dxvr[3]*fr[5]+25.45584412271572*nuVtSqSum[1]*dxvl[3]*dxvr[3]*fl[5]+(25.45584412271572*nuVtSqSum[0]*fl[2]-25.45584412271572*nuVtSqSum[0]*fr[2])*dxvl[3]*dxvr[3]))/(2.0*dxvr3R3+6.0*dxvl[3]*dxvr3R2+6.0*dxvl3R2*dxvr[3]+2.0*dxvl3R3))+0.25*(2.449489742783178*alphaDrag[1]*favg[12]+2.449489742783178*alphaDrag[0]*favg[9]+1.414213562373095*alphaDrag[1]*favg[5]+1.414213562373095*alphaDrag[0]*favg[2])-0.5*(1.732050807568877*fjump[9]+fjump[2]); 
  Ghat[3] = (-(1.0*((12.24744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[1]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[1]*dxvl3R3)*fr[13]+((-2.449489742783178*nuVtSqSum[1]*dxvr3R3)+2.449489742783178*nuVtSqSum[1]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[1]*dxvl3R2*dxvr[3])*fl[13]+(12.24744871391589*nuVtSqSum[0]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[0]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[0]*dxvl3R3)*fr[10]+((-2.449489742783178*nuVtSqSum[0]*dxvr3R3)+2.449489742783178*nuVtSqSum[0]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[0]*dxvl3R2*dxvr[3])*fl[10]-25.45584412271572*nuVtSqSum[1]*dxvl[3]*dxvr[3]*fr[6]+25.45584412271572*nuVtSqSum[1]*dxvl[3]*dxvr[3]*fl[6]-25.45584412271572*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fr[3]+25.45584412271572*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fl[3]))/(2.0*dxvr3R3+6.0*dxvl[3]*dxvr3R2+6.0*dxvl3R2*dxvr[3]+2.0*dxvl3R3))+0.25*(2.449489742783178*alphaDrag[1]*favg[13]+2.449489742783178*alphaDrag[0]*favg[10]+1.414213562373095*alphaDrag[1]*favg[6]+1.414213562373095*alphaDrag[0]*favg[3])-0.5*(1.732050807568877*fjump[10]+fjump[3]); 
  Ghat[5] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[0]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[0]*dxvl3R3)*fr[12]+((-2.449489742783178*nuVtSqSum[0]*dxvr3R3)+2.449489742783178*nuVtSqSum[0]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[0]*dxvl3R2*dxvr[3])*fl[12]+(12.24744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[1]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[1]*dxvl3R3)*fr[9]+((-2.449489742783178*nuVtSqSum[1]*dxvr3R3)+2.449489742783178*nuVtSqSum[1]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[1]*dxvl3R2*dxvr[3])*fl[9]-25.45584412271572*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fr[5]+25.45584412271572*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fl[5]+(25.45584412271572*nuVtSqSum[1]*fl[2]-25.45584412271572*nuVtSqSum[1]*fr[2])*dxvl[3]*dxvr[3]))/(2.0*dxvr3R3+6.0*dxvl[3]*dxvr3R2+6.0*dxvl3R2*dxvr[3]+2.0*dxvl3R3))-0.5*(1.732050807568877*fjump[12]+fjump[5])+0.25*(2.449489742783178*alphaDrag[0]*favg[12]+2.449489742783178*alphaDrag[1]*favg[9]+1.414213562373095*alphaDrag[0]*favg[5]+1.414213562373095*alphaDrag[1]*favg[2]); 
  Ghat[6] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[0]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[0]*dxvl3R3)*fr[13]+((-2.449489742783178*nuVtSqSum[0]*dxvr3R3)+2.449489742783178*nuVtSqSum[0]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[0]*dxvl3R2*dxvr[3])*fl[13]+(12.24744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[1]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[1]*dxvl3R3)*fr[10]+((-2.449489742783178*nuVtSqSum[1]*dxvr3R3)+2.449489742783178*nuVtSqSum[1]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[1]*dxvl3R2*dxvr[3])*fl[10]-25.45584412271572*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fr[6]+25.45584412271572*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fl[6]-25.45584412271572*nuVtSqSum[1]*dxvl[3]*dxvr[3]*fr[3]+25.45584412271572*nuVtSqSum[1]*dxvl[3]*dxvr[3]*fl[3]))/(2.0*dxvr3R3+6.0*dxvl[3]*dxvr3R2+6.0*dxvl3R2*dxvr[3]+2.0*dxvl3R3))-0.5*(1.732050807568877*fjump[13]+fjump[6])+0.25*(2.449489742783178*alphaDrag[0]*favg[13]+2.449489742783178*alphaDrag[1]*favg[10]+1.414213562373095*alphaDrag[0]*favg[6]+1.414213562373095*alphaDrag[1]*favg[3]); 
  Ghat[7] = (-(1.0*((12.24744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[1]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[1]*dxvl3R3)*fr[15]+((-2.449489742783178*nuVtSqSum[1]*dxvr3R3)+2.449489742783178*nuVtSqSum[1]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[1]*dxvl3R2*dxvr[3])*fl[15]+(12.24744871391589*nuVtSqSum[0]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[0]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[0]*dxvl3R3)*fr[14]+((-2.449489742783178*nuVtSqSum[0]*dxvr3R3)+2.449489742783178*nuVtSqSum[0]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[0]*dxvl3R2*dxvr[3])*fl[14]-25.45584412271572*nuVtSqSum[1]*dxvl[3]*dxvr[3]*fr[11]+25.45584412271572*nuVtSqSum[1]*dxvl[3]*dxvr[3]*fl[11]-25.45584412271572*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fr[7]+25.45584412271572*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fl[7]))/(2.0*dxvr3R3+6.0*dxvl[3]*dxvr3R2+6.0*dxvl3R2*dxvr[3]+2.0*dxvl3R3))+0.25*(2.449489742783178*alphaDrag[1]*favg[15]+2.449489742783178*alphaDrag[0]*favg[14]+1.414213562373095*alphaDrag[1]*favg[11]+1.414213562373095*alphaDrag[0]*favg[7])-0.5*(1.732050807568877*fjump[14]+fjump[7]); 
  Ghat[11] = (-(1.0*((12.24744871391589*nuVtSqSum[0]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[0]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[0]*dxvl3R3)*fr[15]+((-2.449489742783178*nuVtSqSum[0]*dxvr3R3)+2.449489742783178*nuVtSqSum[0]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[0]*dxvl3R2*dxvr[3])*fl[15]+(12.24744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R2+2.449489742783178*nuVtSqSum[1]*dxvl3R2*dxvr[3]-2.449489742783178*nuVtSqSum[1]*dxvl3R3)*fr[14]+((-2.449489742783178*nuVtSqSum[1]*dxvr3R3)+2.449489742783178*nuVtSqSum[1]*dxvl[3]*dxvr3R2+12.24744871391589*nuVtSqSum[1]*dxvl3R2*dxvr[3])*fl[14]-25.45584412271572*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fr[11]+25.45584412271572*nuVtSqSum[0]*dxvl[3]*dxvr[3]*fl[11]-25.45584412271572*nuVtSqSum[1]*dxvl[3]*dxvr[3]*fr[7]+25.45584412271572*nuVtSqSum[1]*dxvl[3]*dxvr[3]*fl[7]))/(2.0*dxvr3R3+6.0*dxvl[3]*dxvr3R2+6.0*dxvl3R2*dxvr[3]+2.0*dxvl3R3))-0.5*(1.732050807568877*fjump[15]+fjump[11])+0.25*(2.449489742783178*alphaDrag[0]*favg[15]+2.449489742783178*alphaDrag[1]*favg[14]+1.414213562373095*alphaDrag[0]*favg[11]+1.414213562373095*alphaDrag[1]*favg[7]); 

  double incr1[16]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[0]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = -0.5*Ghat[6]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = 0.8660254037844386*Ghat[1]; 
  incr1[9] = 0.8660254037844386*Ghat[2]; 
  incr1[10] = 0.8660254037844386*Ghat[3]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = 0.8660254037844386*Ghat[5]; 
  incr1[13] = 0.8660254037844386*Ghat[6]; 
  incr1[14] = 0.8660254037844386*Ghat[7]; 
  incr1[15] = 0.8660254037844386*Ghat[11]; 

  double incr2[16]; 

  incr2[4] = -(1.0*((7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[1]*dxvl3R3*dxvr[3])*fr[8]+((-4.242640687119286*nuVtSqSum[1]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2)*fl[8]+(7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[0]*dxvl3R3*dxvr[3])*fr[4]+((-4.242640687119286*nuVtSqSum[0]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2)*fl[4]+((-4.898979485566357*fl[1]*nuVtSqSum[1])-4.898979485566357*fl[0]*nuVtSqSum[0])*dxvr3R3+((-14.69693845669907*fl[1]*nuVtSqSum[1])-14.69693845669907*fl[0]*nuVtSqSum[0])*dxvl[3]*dxvr3R2+((-14.69693845669907*fr[1]*nuVtSqSum[1])-14.69693845669907*fr[0]*nuVtSqSum[0])*dxvl3R2*dxvr[3]+((-4.898979485566357*fr[1]*nuVtSqSum[1])-4.898979485566357*fr[0]*nuVtSqSum[0])*dxvl3R3))/(8.0*dxvr3R3+24.0*dxvl[3]*dxvr3R2+24.0*dxvl3R2*dxvr[3]+8.0*dxvl3R3); 
  incr2[8] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[0]*dxvl3R3*dxvr[3])*fr[8]+((-4.242640687119286*nuVtSqSum[0]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2)*fl[8]+(7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[1]*dxvl3R3*dxvr[3])*fr[4]+((-4.242640687119286*nuVtSqSum[1]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2)*fl[4]+((-4.898979485566357*fl[0]*nuVtSqSum[1])-4.898979485566357*nuVtSqSum[0]*fl[1])*dxvr3R3+((-14.69693845669907*fl[0]*nuVtSqSum[1])-14.69693845669907*nuVtSqSum[0]*fl[1])*dxvl[3]*dxvr3R2+((-14.69693845669907*fr[0]*nuVtSqSum[1])-14.69693845669907*nuVtSqSum[0]*fr[1])*dxvl3R2*dxvr[3]+((-4.898979485566357*fr[0]*nuVtSqSum[1])-4.898979485566357*nuVtSqSum[0]*fr[1])*dxvl3R3))/(8.0*dxvr3R3+24.0*dxvl[3]*dxvr3R2+24.0*dxvl3R2*dxvr[3]+8.0*dxvl3R3); 
  incr2[9] = -(1.0*((7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[1]*dxvl3R3*dxvr[3])*fr[12]+((-4.242640687119286*nuVtSqSum[1]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2)*fl[12]+(7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[0]*dxvl3R3*dxvr[3])*fr[9]+((-4.242640687119286*nuVtSqSum[0]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2)*fl[9]+((-14.69693845669907*nuVtSqSum[1]*dxvl3R2*dxvr[3])-4.898979485566357*nuVtSqSum[1]*dxvl3R3)*fr[5]+((-4.898979485566357*nuVtSqSum[1]*dxvr3R3)-14.69693845669907*nuVtSqSum[1]*dxvl[3]*dxvr3R2)*fl[5]-4.898979485566357*nuVtSqSum[0]*fl[2]*dxvr3R3-14.69693845669907*nuVtSqSum[0]*fl[2]*dxvl[3]*dxvr3R2-14.69693845669907*nuVtSqSum[0]*fr[2]*dxvl3R2*dxvr[3]-4.898979485566357*nuVtSqSum[0]*fr[2]*dxvl3R3))/(8.0*dxvr3R3+24.0*dxvl[3]*dxvr3R2+24.0*dxvl3R2*dxvr[3]+8.0*dxvl3R3); 
  incr2[10] = -(1.0*((7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[1]*dxvl3R3*dxvr[3])*fr[13]+((-4.242640687119286*nuVtSqSum[1]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2)*fl[13]+(7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[0]*dxvl3R3*dxvr[3])*fr[10]+((-4.242640687119286*nuVtSqSum[0]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2)*fl[10]+((-14.69693845669907*nuVtSqSum[1]*dxvl3R2*dxvr[3])-4.898979485566357*nuVtSqSum[1]*dxvl3R3)*fr[6]+((-4.898979485566357*nuVtSqSum[1]*dxvr3R3)-14.69693845669907*nuVtSqSum[1]*dxvl[3]*dxvr3R2)*fl[6]+((-14.69693845669907*nuVtSqSum[0]*dxvl3R2*dxvr[3])-4.898979485566357*nuVtSqSum[0]*dxvl3R3)*fr[3]+((-4.898979485566357*nuVtSqSum[0]*dxvr3R3)-14.69693845669907*nuVtSqSum[0]*dxvl[3]*dxvr3R2)*fl[3]))/(8.0*dxvr3R3+24.0*dxvl[3]*dxvr3R2+24.0*dxvl3R2*dxvr[3]+8.0*dxvl3R3); 
  incr2[12] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[0]*dxvl3R3*dxvr[3])*fr[12]+((-4.242640687119286*nuVtSqSum[0]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2)*fl[12]+(7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[1]*dxvl3R3*dxvr[3])*fr[9]+((-4.242640687119286*nuVtSqSum[1]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2)*fl[9]+((-14.69693845669907*nuVtSqSum[0]*dxvl3R2*dxvr[3])-4.898979485566357*nuVtSqSum[0]*dxvl3R3)*fr[5]+((-4.898979485566357*nuVtSqSum[0]*dxvr3R3)-14.69693845669907*nuVtSqSum[0]*dxvl[3]*dxvr3R2)*fl[5]-4.898979485566357*nuVtSqSum[1]*fl[2]*dxvr3R3-14.69693845669907*nuVtSqSum[1]*fl[2]*dxvl[3]*dxvr3R2-14.69693845669907*nuVtSqSum[1]*fr[2]*dxvl3R2*dxvr[3]-4.898979485566357*nuVtSqSum[1]*fr[2]*dxvl3R3))/(8.0*dxvr3R3+24.0*dxvl[3]*dxvr3R2+24.0*dxvl3R2*dxvr[3]+8.0*dxvl3R3); 
  incr2[13] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[0]*dxvl3R3*dxvr[3])*fr[13]+((-4.242640687119286*nuVtSqSum[0]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2)*fl[13]+(7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[1]*dxvl3R3*dxvr[3])*fr[10]+((-4.242640687119286*nuVtSqSum[1]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2)*fl[10]+((-14.69693845669907*nuVtSqSum[0]*dxvl3R2*dxvr[3])-4.898979485566357*nuVtSqSum[0]*dxvl3R3)*fr[6]+((-4.898979485566357*nuVtSqSum[0]*dxvr3R3)-14.69693845669907*nuVtSqSum[0]*dxvl[3]*dxvr3R2)*fl[6]+((-14.69693845669907*nuVtSqSum[1]*dxvl3R2*dxvr[3])-4.898979485566357*nuVtSqSum[1]*dxvl3R3)*fr[3]+((-4.898979485566357*nuVtSqSum[1]*dxvr3R3)-14.69693845669907*nuVtSqSum[1]*dxvl[3]*dxvr3R2)*fl[3]))/(8.0*dxvr3R3+24.0*dxvl[3]*dxvr3R2+24.0*dxvl3R2*dxvr[3]+8.0*dxvl3R3); 
  incr2[14] = -(1.0*((7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[1]*dxvl3R3*dxvr[3])*fr[15]+((-4.242640687119286*nuVtSqSum[1]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2)*fl[15]+(7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[0]*dxvl3R3*dxvr[3])*fr[14]+((-4.242640687119286*nuVtSqSum[0]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2)*fl[14]+((-14.69693845669907*nuVtSqSum[1]*dxvl3R2*dxvr[3])-4.898979485566357*nuVtSqSum[1]*dxvl3R3)*fr[11]+((-4.898979485566357*nuVtSqSum[1]*dxvr3R3)-14.69693845669907*nuVtSqSum[1]*dxvl[3]*dxvr3R2)*fl[11]+((-14.69693845669907*nuVtSqSum[0]*dxvl3R2*dxvr[3])-4.898979485566357*nuVtSqSum[0]*dxvl3R3)*fr[7]+((-4.898979485566357*nuVtSqSum[0]*dxvr3R3)-14.69693845669907*nuVtSqSum[0]*dxvl[3]*dxvr3R2)*fl[7]))/(8.0*dxvr3R3+24.0*dxvl[3]*dxvr3R2+24.0*dxvl3R2*dxvr[3]+8.0*dxvl3R3); 
  incr2[15] = -(1.0*((7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[0]*dxvl3R3*dxvr[3])*fr[15]+((-4.242640687119286*nuVtSqSum[0]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[0]*dxvl3R2*dxvr3R2)*fl[15]+(7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2+4.242640687119286*nuVtSqSum[1]*dxvl3R3*dxvr[3])*fr[14]+((-4.242640687119286*nuVtSqSum[1]*dxvl[3]*dxvr3R3)-7.071067811865476*nuVtSqSum[1]*dxvl3R2*dxvr3R2)*fl[14]+((-14.69693845669907*nuVtSqSum[0]*dxvl3R2*dxvr[3])-4.898979485566357*nuVtSqSum[0]*dxvl3R3)*fr[11]+((-4.898979485566357*nuVtSqSum[0]*dxvr3R3)-14.69693845669907*nuVtSqSum[0]*dxvl[3]*dxvr3R2)*fl[11]+((-14.69693845669907*nuVtSqSum[1]*dxvl3R2*dxvr[3])-4.898979485566357*nuVtSqSum[1]*dxvl3R3)*fr[7]+((-4.898979485566357*nuVtSqSum[1]*dxvr3R3)-14.69693845669907*nuVtSqSum[1]*dxvl[3]*dxvr3R2)*fl[7]))/(8.0*dxvr3R3+24.0*dxvl[3]*dxvr3R2+24.0*dxvl3R2*dxvr[3]+8.0*dxvl3R3); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr2[8]*rdvSq4R+incr1[8]*rdv2R; 
  outr[9] += incr2[9]*rdvSq4R+incr1[9]*rdv2R; 
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr2[12]*rdvSq4R+incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += incr1[8]*rdv2L-1.0*incr2[8]*rdvSq4L; 
  outl[9] += incr1[9]*rdv2L-1.0*incr2[9]*rdvSq4L; 
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += incr1[12]*rdv2L-1.0*incr2[12]*rdvSq4L; 
  outl[13] += incr1[13]*rdv2L-1.0*incr2[13]*rdvSq4L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += incr1[15]*rdv2L-1.0*incr2[15]*rdvSq4L; 

  return std::abs(wl[3]-(0.7071067811865475*sumNuUz[0])/nuSum); 
} 
