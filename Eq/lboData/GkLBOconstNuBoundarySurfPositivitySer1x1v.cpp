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

  double outrPos[4]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.5 : positivityWeightByDirR[2]/positivityWeightByDirR[0]; 
  outrPos[0] = 0.04166666666666666*((7.348469228349534*nuVtSqSum[1]-4.242640687119286*nuVtSqSum[0])*fr[3]+(7.348469228349534*nuVtSqSum[0]-4.242640687119286*nuVtSqSum[1])*fr[2]+(2.449489742783178*fr[0]-4.242640687119286*fr[1])*nuVtSqSum[1]+nuVtSqSum[0]*(2.449489742783178*fr[1]-4.242640687119286*fr[0]))*rdvSq4R; 
  outrPos[1] = 0.04166666666666666*((7.348469228349534*nuVtSqSum[1]+4.242640687119286*nuVtSqSum[0])*fr[3]+(4.242640687119286*nuVtSqSum[1]+7.348469228349534*nuVtSqSum[0])*fr[2]+((-4.242640687119286*fr[1])-2.449489742783178*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-2.449489742783178*fr[1])-4.242640687119286*fr[0]))*rdvSq4R; 
  outrPos[2] = -0.04166666666666666*((7.348469228349534*nuVtSqSum[1]-4.242640687119286*nuVtSqSum[0])*fr[3]+(7.348469228349534*nuVtSqSum[0]-4.242640687119286*nuVtSqSum[1])*fr[2]+(2.449489742783178*fr[0]-4.242640687119286*fr[1])*nuVtSqSum[1]+nuVtSqSum[0]*(2.449489742783178*fr[1]-4.242640687119286*fr[0]))*rdvSq4R; 
  outrPos[3] = -0.04166666666666666*((7.348469228349534*nuVtSqSum[1]+4.242640687119286*nuVtSqSum[0])*fr[3]+(4.242640687119286*nuVtSqSum[1]+7.348469228349534*nuVtSqSum[0])*fr[2]+((-4.242640687119286*fr[1])-2.449489742783178*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-2.449489742783178*fr[1])-4.242640687119286*fr[0]))*rdvSq4R; 
  if(outrPos[0] < 0.) limFac = std::min(1.0, -fluxFracR*(0.1666666666666667*(fr[3]-1.732050807568877*(fr[2]+fr[1])+3.0*fr[0]))/dtApprox/outrPos[0]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPos[2] *= limFac; 
  outrPos[0] *= limFac; 
  if(outrPos[1] < 0.) limFac = std::min(1.0, -fluxFracR*(-0.1666666666666667*(fr[3]+1.732050807568877*fr[2]-1.732050807568877*fr[1]-3.0*fr[0]))/dtApprox/outrPos[1]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPos[3] *= limFac; 
  outrPos[1] *= limFac; 
    outr[0] += 0.5*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
    outr[1] += 0.8660254037844386*outrPos[3]-0.8660254037844386*outrPos[2]+0.8660254037844386*outrPos[1]-0.8660254037844386*outrPos[0]; 
    outr[2] += 0.8660254037844386*(outrPos[3]+outrPos[2])-0.8660254037844386*(outrPos[1]+outrPos[0]); 
    outr[3] += 1.5*outrPos[3]-1.5*(outrPos[2]+outrPos[1])+1.5*outrPos[0]; 

  } else {

  double outlPos[4]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.5 : positivityWeightByDirL[2]/positivityWeightByDirL[0]; 
  outlPos[0] = 0.04166666666666666*((7.348469228349534*nuVtSqSum[1]-4.242640687119286*nuVtSqSum[0])*fl[3]+(7.348469228349534*nuVtSqSum[0]-4.242640687119286*nuVtSqSum[1])*fl[2]+(4.242640687119286*fl[1]-2.449489742783178*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(4.242640687119286*fl[0]-2.449489742783178*fl[1]))*rdvSq4L; 
  outlPos[1] = 0.04166666666666666*((7.348469228349534*nuVtSqSum[1]+4.242640687119286*nuVtSqSum[0])*fl[3]+(4.242640687119286*nuVtSqSum[1]+7.348469228349534*nuVtSqSum[0])*fl[2]+(4.242640687119286*fl[1]+2.449489742783178*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(2.449489742783178*fl[1]+4.242640687119286*fl[0]))*rdvSq4L; 
  outlPos[2] = -0.04166666666666666*((7.348469228349534*nuVtSqSum[1]-4.242640687119286*nuVtSqSum[0])*fl[3]+(7.348469228349534*nuVtSqSum[0]-4.242640687119286*nuVtSqSum[1])*fl[2]+(4.242640687119286*fl[1]-2.449489742783178*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(4.242640687119286*fl[0]-2.449489742783178*fl[1]))*rdvSq4L; 
  outlPos[3] = -0.04166666666666666*((7.348469228349534*nuVtSqSum[1]+4.242640687119286*nuVtSqSum[0])*fl[3]+(4.242640687119286*nuVtSqSum[1]+7.348469228349534*nuVtSqSum[0])*fl[2]+(4.242640687119286*fl[1]+2.449489742783178*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(2.449489742783178*fl[1]+4.242640687119286*fl[0]))*rdvSq4L; 
  if(outlPos[2] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.1666666666666667*(fl[3]-1.732050807568877*fl[2]+1.732050807568877*fl[1]-3.0*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outlPos[0] *= limFac; 
  if(outlPos[3] < 0.) limFac = std::min(1.0, -fluxFracL*(0.1666666666666667*(fl[3]+1.732050807568877*(fl[2]+fl[1])+3.0*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[1] *= limFac; 
    outl[0] += 0.5*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
    outl[1] += 0.8660254037844386*outlPos[3]-0.8660254037844386*outlPos[2]+0.8660254037844386*outlPos[1]-0.8660254037844386*outlPos[0]; 
    outl[2] += 0.8660254037844386*(outlPos[3]+outlPos[2])-0.8660254037844386*(outlPos[1]+outlPos[0]); 
    outl[3] += 1.5*outlPos[3]-1.5*(outlPos[2]+outlPos[1])+1.5*outlPos[0]; 

  }
  return 0.0; 
} 
