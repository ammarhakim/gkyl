#include <GkLBOModDecl.h> 
double GkLBOconstNuBoundarySurfPositivity1x1vSer_Vpar_P1(const double m_, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const int *idxl, const int *idxr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:         Cell-center coordinates.
  // dxv[2]:       Cell spacing.
  // idx[2]:       current grid index.
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 
  double limFac, fluxFracL, fluxFracR;

  if (idxr[1] == 1) {

  double outrPosP[2], outrPosM[2]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.5 : positivityWeightByDirR[1]/positivityWeightByDirR[0]; 
  outrPosP[0] = -0.04166666666666666*((7.348469228349534*nuVtSqSum[1]-4.242640687119286*nuVtSqSum[0])*fr[3]+(7.348469228349534*nuVtSqSum[0]-4.242640687119286*nuVtSqSum[1])*fr[2]+(2.449489742783178*fr[0]-4.242640687119286*fr[1])*nuVtSqSum[1]+nuVtSqSum[0]*(2.449489742783178*fr[1]-4.242640687119286*fr[0]))*rdvSq4R; 
  outrPosP[1] = -0.04166666666666666*((7.348469228349534*nuVtSqSum[1]+4.242640687119286*nuVtSqSum[0])*fr[3]+(4.242640687119286*nuVtSqSum[1]+7.348469228349534*nuVtSqSum[0])*fr[2]+((-4.242640687119286*fr[1])-2.449489742783178*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-2.449489742783178*fr[1])-4.242640687119286*fr[0]))*rdvSq4R; 
  outrPosM[0] = 0.04166666666666666*((7.348469228349534*nuVtSqSum[1]-4.242640687119286*nuVtSqSum[0])*fr[3]+(7.348469228349534*nuVtSqSum[0]-4.242640687119286*nuVtSqSum[1])*fr[2]+(2.449489742783178*fr[0]-4.242640687119286*fr[1])*nuVtSqSum[1]+nuVtSqSum[0]*(2.449489742783178*fr[1]-4.242640687119286*fr[0]))*rdvSq4R; 
  outrPosM[1] = 0.04166666666666666*((7.348469228349534*nuVtSqSum[1]+4.242640687119286*nuVtSqSum[0])*fr[3]+(4.242640687119286*nuVtSqSum[1]+7.348469228349534*nuVtSqSum[0])*fr[2]+((-4.242640687119286*fr[1])-2.449489742783178*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-2.449489742783178*fr[1])-4.242640687119286*fr[0]))*rdvSq4R; 
  if(outrPosM[0] < 0.) limFac = std::min(1.0, -fluxFracR*(0.1666666666666667*(fr[3]-1.732050807568877*(fr[2]+fr[1])+3.0*fr[0]))/dtApprox/outrPosM[0]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPosP[0] *= limFac; 
  outrPosM[0] *= limFac; 
  if(outrPosM[1] < 0.) limFac = std::min(1.0, -fluxFracR*(-0.1666666666666667*(fr[3]+1.732050807568877*fr[2]-1.732050807568877*fr[1]-3.0*fr[0]))/dtApprox/outrPosM[1]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPosP[1] *= limFac; 
  outrPosM[1] *= limFac; 
    outr[0] += 0.5*(outrPosP[1]+outrPosM[1]+outrPosP[0]+outrPosM[0]); 
    outr[1] += 0.8660254037844386*(outrPosP[1]+outrPosM[1])-0.8660254037844386*(outrPosP[0]+outrPosM[0]); 
    outr[2] += 0.8660254037844386*outrPosP[1]-0.8660254037844386*outrPosM[1]+0.8660254037844386*outrPosP[0]-0.8660254037844386*outrPosM[0]; 
    outr[3] += 1.5*outrPosP[1]-1.5*(outrPosM[1]+outrPosP[0])+1.5*outrPosM[0]; 

  } else {

  double outlPosP[2], outlPosM[2]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.5 : positivityWeightByDirL[1]/positivityWeightByDirL[0]; 
  outlPosP[0] = -0.04166666666666666*((7.348469228349534*nuVtSqSum[1]-4.242640687119286*nuVtSqSum[0])*fl[3]+(7.348469228349534*nuVtSqSum[0]-4.242640687119286*nuVtSqSum[1])*fl[2]+(4.242640687119286*fl[1]-2.449489742783178*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(4.242640687119286*fl[0]-2.449489742783178*fl[1]))*rdvSq4L; 
  outlPosP[1] = -0.04166666666666666*((7.348469228349534*nuVtSqSum[1]+4.242640687119286*nuVtSqSum[0])*fl[3]+(4.242640687119286*nuVtSqSum[1]+7.348469228349534*nuVtSqSum[0])*fl[2]+(4.242640687119286*fl[1]+2.449489742783178*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(2.449489742783178*fl[1]+4.242640687119286*fl[0]))*rdvSq4L; 
  outlPosM[0] = 0.04166666666666666*((7.348469228349534*nuVtSqSum[1]-4.242640687119286*nuVtSqSum[0])*fl[3]+(7.348469228349534*nuVtSqSum[0]-4.242640687119286*nuVtSqSum[1])*fl[2]+(4.242640687119286*fl[1]-2.449489742783178*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(4.242640687119286*fl[0]-2.449489742783178*fl[1]))*rdvSq4L; 
  outlPosM[1] = 0.04166666666666666*((7.348469228349534*nuVtSqSum[1]+4.242640687119286*nuVtSqSum[0])*fl[3]+(4.242640687119286*nuVtSqSum[1]+7.348469228349534*nuVtSqSum[0])*fl[2]+(4.242640687119286*fl[1]+2.449489742783178*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(2.449489742783178*fl[1]+4.242640687119286*fl[0]))*rdvSq4L; 
  if(outlPosP[0] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.1666666666666667*(fl[3]-1.732050807568877*fl[2]+1.732050807568877*fl[1]-3.0*fl[0]))/dtApprox/outlPosP[0]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPosP[0] *= limFac; 
  outlPosM[0] *= limFac; 
  if(outlPosP[1] < 0.) limFac = std::min(1.0, -fluxFracL*(0.1666666666666667*(fl[3]+1.732050807568877*(fl[2]+fl[1])+3.0*fl[0]))/dtApprox/outlPosP[1]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPosP[1] *= limFac; 
  outlPosM[1] *= limFac; 
    outl[0] += 0.5*(outlPosP[1]+outlPosM[1]+outlPosP[0]+outlPosM[0]); 
    outl[1] += 0.8660254037844386*(outlPosP[1]+outlPosM[1])-0.8660254037844386*(outlPosP[0]+outlPosM[0]); 
    outl[2] += 0.8660254037844386*outlPosP[1]-0.8660254037844386*outlPosM[1]+0.8660254037844386*outlPosP[0]-0.8660254037844386*outlPosM[0]; 
    outl[3] += 1.5*outlPosP[1]-1.5*(outlPosM[1]+outlPosP[0])+1.5*outlPosM[0]; 

  }
  return 0.0; 
} 
