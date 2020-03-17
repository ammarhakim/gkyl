#include <GkLBOModDecl.h> 
double GkLBOconstNuBoundarySurfPositivity1x2vSer_Vpar_P1(const double m_, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const int *idxl, const int *idxr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:         Cell-center coordinates.
  // dxv[3]:       Cell spacing.
  // idx[3]:       current grid index.
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

  double outrPosP[4], outrPosM[4]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.3333333333333333 : positivityWeightByDirR[1]/positivityWeightByDirR[0]; 
  outrPosP[0] = 0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[6]+(nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[5]+(3.0*nuVtSqSum[0]-5.196152422706631*nuVtSqSum[1])*fr[4]+(nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[3]+(3.0*nuVtSqSum[1]-5.196152422706631*nuVtSqSum[0])*fr[2]+(3.0*fr[1]-1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*(3.0*fr[0]-1.732050807568877*fr[1]))*rdvSq4R; 
  outrPosP[1] = 0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fr[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[6]+((-1.732050807568877*nuVtSqSum[1])-1.0*nuVtSqSum[0])*fr[5]+((-5.196152422706631*nuVtSqSum[1])-3.0*nuVtSqSum[0])*fr[4]+((-1.0*nuVtSqSum[1])-1.732050807568877*nuVtSqSum[0])*fr[3]+((-3.0*nuVtSqSum[1])-5.196152422706631*nuVtSqSum[0])*fr[2]+(3.0*fr[1]+1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fr[1]+3.0*fr[0]))*rdvSq4R; 
  outrPosP[2] = -0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[6]+(nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[5]+(5.196152422706631*nuVtSqSum[1]-3.0*nuVtSqSum[0])*fr[4]+(nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[3]+(5.196152422706631*nuVtSqSum[0]-3.0*nuVtSqSum[1])*fr[2]+(1.732050807568877*fr[0]-3.0*fr[1])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fr[1]-3.0*fr[0]))*rdvSq4R; 
  outrPosP[3] = -0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fr[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[6]+((-1.732050807568877*nuVtSqSum[1])-1.0*nuVtSqSum[0])*fr[5]+(5.196152422706631*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[4]+((-1.0*nuVtSqSum[1])-1.732050807568877*nuVtSqSum[0])*fr[3]+(3.0*nuVtSqSum[1]+5.196152422706631*nuVtSqSum[0])*fr[2]+((-3.0*fr[1])-1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-1.732050807568877*fr[1])-3.0*fr[0]))*rdvSq4R; 
  outrPosM[0] = -0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[6]+(nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[5]+(3.0*nuVtSqSum[0]-5.196152422706631*nuVtSqSum[1])*fr[4]+(nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[3]+(3.0*nuVtSqSum[1]-5.196152422706631*nuVtSqSum[0])*fr[2]+(3.0*fr[1]-1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*(3.0*fr[0]-1.732050807568877*fr[1]))*rdvSq4R; 
  outrPosM[1] = -0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fr[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[6]+((-1.732050807568877*nuVtSqSum[1])-1.0*nuVtSqSum[0])*fr[5]+((-5.196152422706631*nuVtSqSum[1])-3.0*nuVtSqSum[0])*fr[4]+((-1.0*nuVtSqSum[1])-1.732050807568877*nuVtSqSum[0])*fr[3]+((-3.0*nuVtSqSum[1])-5.196152422706631*nuVtSqSum[0])*fr[2]+(3.0*fr[1]+1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fr[1]+3.0*fr[0]))*rdvSq4R; 
  outrPosM[2] = 0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[6]+(nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[5]+(5.196152422706631*nuVtSqSum[1]-3.0*nuVtSqSum[0])*fr[4]+(nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[3]+(5.196152422706631*nuVtSqSum[0]-3.0*nuVtSqSum[1])*fr[2]+(1.732050807568877*fr[0]-3.0*fr[1])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fr[1]-3.0*fr[0]))*rdvSq4R; 
  outrPosM[3] = 0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fr[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[6]+((-1.732050807568877*nuVtSqSum[1])-1.0*nuVtSqSum[0])*fr[5]+(5.196152422706631*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[4]+((-1.0*nuVtSqSum[1])-1.732050807568877*nuVtSqSum[0])*fr[3]+(3.0*nuVtSqSum[1]+5.196152422706631*nuVtSqSum[0])*fr[2]+((-3.0*fr[1])-1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-1.732050807568877*fr[1])-3.0*fr[0]))*rdvSq4R; 
  if(outrPosM[0] < 0.) limFac = std::min(1.0, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPosM[0]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPosP[0] *= limFac; 
  outrPosM[0] *= limFac; 
  if(outrPosM[1] < 0.) limFac = std::min(1.0, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPosM[1]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPosP[1] *= limFac; 
  outrPosM[1] *= limFac; 
  if(outrPosM[2] < 0.) limFac = std::min(1.0, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0]))/dtApprox/outrPosM[2]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPosP[2] *= limFac; 
  outrPosM[2] *= limFac; 
  if(outrPosM[3] < 0.) limFac = std::min(1.0, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*fr[5]+4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]-12.72792206135786*fr[0]))/dtApprox/outrPosM[3]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPosP[3] *= limFac; 
  outrPosM[3] *= limFac; 
    outr[0] += 0.3535533905932737*(outrPosP[3]+outrPosM[3]+outrPosP[2]+outrPosM[2]+outrPosP[1]+outrPosM[1]+outrPosP[0]+outrPosM[0]); 
    outr[1] += 0.6123724356957944*(outrPosP[3]+outrPosM[3])-0.6123724356957944*(outrPosP[2]+outrPosM[2])+0.6123724356957944*(outrPosP[1]+outrPosM[1])-0.6123724356957944*(outrPosP[0]+outrPosM[0]); 
    outr[2] += 0.6123724356957944*outrPosP[3]-0.6123724356957944*outrPosM[3]+0.6123724356957944*outrPosP[2]-0.6123724356957944*outrPosM[2]+0.6123724356957944*outrPosP[1]-0.6123724356957944*outrPosM[1]+0.6123724356957944*outrPosP[0]-0.6123724356957944*outrPosM[0]; 
    outr[3] += 0.6123724356957944*(outrPosP[3]+outrPosM[3]+outrPosP[2]+outrPosM[2])-0.6123724356957944*(outrPosP[1]+outrPosM[1]+outrPosP[0]+outrPosM[0]); 
    outr[4] += 1.060660171779821*outrPosP[3]-1.060660171779821*(outrPosM[3]+outrPosP[2])+1.060660171779821*(outrPosM[2]+outrPosP[1])-1.060660171779821*(outrPosM[1]+outrPosP[0])+1.060660171779821*outrPosM[0]; 
    outr[5] += 1.060660171779821*(outrPosP[3]+outrPosM[3])-1.060660171779821*(outrPosP[2]+outrPosM[2]+outrPosP[1]+outrPosM[1])+1.060660171779821*(outrPosP[0]+outrPosM[0]); 
    outr[6] += 1.060660171779821*outrPosP[3]-1.060660171779821*outrPosM[3]+1.060660171779821*outrPosP[2]-1.060660171779821*(outrPosM[2]+outrPosP[1])+1.060660171779821*outrPosM[1]-1.060660171779821*outrPosP[0]+1.060660171779821*outrPosM[0]; 
    outr[7] += 1.837117307087383*outrPosP[3]-1.837117307087383*(outrPosM[3]+outrPosP[2])+1.837117307087383*outrPosM[2]-1.837117307087383*outrPosP[1]+1.837117307087383*(outrPosM[1]+outrPosP[0])-1.837117307087383*outrPosM[0]; 

  } else {

  double outlPosP[4], outlPosM[4]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.3333333333333333 : positivityWeightByDirL[1]/positivityWeightByDirL[0]; 
  outlPosP[0] = 0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fl[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fl[6]+(1.732050807568877*nuVtSqSum[1]-1.0*nuVtSqSum[0])*fl[5]+(3.0*nuVtSqSum[0]-5.196152422706631*nuVtSqSum[1])*fl[4]+(1.732050807568877*nuVtSqSum[0]-1.0*nuVtSqSum[1])*fl[3]+(3.0*nuVtSqSum[1]-5.196152422706631*nuVtSqSum[0])*fl[2]+(1.732050807568877*fl[0]-3.0*fl[1])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fl[1]-3.0*fl[0]))*rdvSq4L; 
  outlPosP[1] = 0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[6]+(1.732050807568877*nuVtSqSum[1]+nuVtSqSum[0])*fl[5]+((-5.196152422706631*nuVtSqSum[1])-3.0*nuVtSqSum[0])*fl[4]+(nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[3]+((-3.0*nuVtSqSum[1])-5.196152422706631*nuVtSqSum[0])*fl[2]+((-3.0*fl[1])-1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-1.732050807568877*fl[1])-3.0*fl[0]))*rdvSq4L; 
  outlPosP[2] = -0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fl[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fl[6]+(1.732050807568877*nuVtSqSum[1]-1.0*nuVtSqSum[0])*fl[5]+(5.196152422706631*nuVtSqSum[1]-3.0*nuVtSqSum[0])*fl[4]+(1.732050807568877*nuVtSqSum[0]-1.0*nuVtSqSum[1])*fl[3]+(5.196152422706631*nuVtSqSum[0]-3.0*nuVtSqSum[1])*fl[2]+(3.0*fl[1]-1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(3.0*fl[0]-1.732050807568877*fl[1]))*rdvSq4L; 
  outlPosP[3] = -0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[6]+(1.732050807568877*nuVtSqSum[1]+nuVtSqSum[0])*fl[5]+(5.196152422706631*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[4]+(nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[3]+(3.0*nuVtSqSum[1]+5.196152422706631*nuVtSqSum[0])*fl[2]+(3.0*fl[1]+1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fl[1]+3.0*fl[0]))*rdvSq4L; 
  outlPosM[0] = -0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fl[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fl[6]+(1.732050807568877*nuVtSqSum[1]-1.0*nuVtSqSum[0])*fl[5]+(3.0*nuVtSqSum[0]-5.196152422706631*nuVtSqSum[1])*fl[4]+(1.732050807568877*nuVtSqSum[0]-1.0*nuVtSqSum[1])*fl[3]+(3.0*nuVtSqSum[1]-5.196152422706631*nuVtSqSum[0])*fl[2]+(1.732050807568877*fl[0]-3.0*fl[1])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fl[1]-3.0*fl[0]))*rdvSq4L; 
  outlPosM[1] = -0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[6]+(1.732050807568877*nuVtSqSum[1]+nuVtSqSum[0])*fl[5]+((-5.196152422706631*nuVtSqSum[1])-3.0*nuVtSqSum[0])*fl[4]+(nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[3]+((-3.0*nuVtSqSum[1])-5.196152422706631*nuVtSqSum[0])*fl[2]+((-3.0*fl[1])-1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-1.732050807568877*fl[1])-3.0*fl[0]))*rdvSq4L; 
  outlPosM[2] = 0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fl[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fl[6]+(1.732050807568877*nuVtSqSum[1]-1.0*nuVtSqSum[0])*fl[5]+(5.196152422706631*nuVtSqSum[1]-3.0*nuVtSqSum[0])*fl[4]+(1.732050807568877*nuVtSqSum[0]-1.0*nuVtSqSum[1])*fl[3]+(5.196152422706631*nuVtSqSum[0]-3.0*nuVtSqSum[1])*fl[2]+(3.0*fl[1]-1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(3.0*fl[0]-1.732050807568877*fl[1]))*rdvSq4L; 
  outlPosM[3] = 0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[6]+(1.732050807568877*nuVtSqSum[1]+nuVtSqSum[0])*fl[5]+(5.196152422706631*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[4]+(nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[3]+(3.0*nuVtSqSum[1]+5.196152422706631*nuVtSqSum[0])*fl[2]+(3.0*fl[1]+1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fl[1]+3.0*fl[0]))*rdvSq4L; 
  if(outlPosP[0] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*fl[5]-4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]+12.72792206135786*fl[0]))/dtApprox/outlPosP[0]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPosP[0] *= limFac; 
  outlPosM[0] *= limFac; 
  if(outlPosP[1] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0]))/dtApprox/outlPosP[1]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPosP[1] *= limFac; 
  outlPosM[1] *= limFac; 
  if(outlPosP[2] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPosP[2]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPosP[2] *= limFac; 
  outlPosM[2] *= limFac; 
  if(outlPosP[3] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPosP[3]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPosP[3] *= limFac; 
  outlPosM[3] *= limFac; 
    outl[0] += 0.3535533905932737*(outlPosP[3]+outlPosM[3]+outlPosP[2]+outlPosM[2]+outlPosP[1]+outlPosM[1]+outlPosP[0]+outlPosM[0]); 
    outl[1] += 0.6123724356957944*(outlPosP[3]+outlPosM[3])-0.6123724356957944*(outlPosP[2]+outlPosM[2])+0.6123724356957944*(outlPosP[1]+outlPosM[1])-0.6123724356957944*(outlPosP[0]+outlPosM[0]); 
    outl[2] += 0.6123724356957944*outlPosP[3]-0.6123724356957944*outlPosM[3]+0.6123724356957944*outlPosP[2]-0.6123724356957944*outlPosM[2]+0.6123724356957944*outlPosP[1]-0.6123724356957944*outlPosM[1]+0.6123724356957944*outlPosP[0]-0.6123724356957944*outlPosM[0]; 
    outl[3] += 0.6123724356957944*(outlPosP[3]+outlPosM[3]+outlPosP[2]+outlPosM[2])-0.6123724356957944*(outlPosP[1]+outlPosM[1]+outlPosP[0]+outlPosM[0]); 
    outl[4] += 1.060660171779821*outlPosP[3]-1.060660171779821*(outlPosM[3]+outlPosP[2])+1.060660171779821*(outlPosM[2]+outlPosP[1])-1.060660171779821*(outlPosM[1]+outlPosP[0])+1.060660171779821*outlPosM[0]; 
    outl[5] += 1.060660171779821*(outlPosP[3]+outlPosM[3])-1.060660171779821*(outlPosP[2]+outlPosM[2]+outlPosP[1]+outlPosM[1])+1.060660171779821*(outlPosP[0]+outlPosM[0]); 
    outl[6] += 1.060660171779821*outlPosP[3]-1.060660171779821*outlPosM[3]+1.060660171779821*outlPosP[2]-1.060660171779821*(outlPosM[2]+outlPosP[1])+1.060660171779821*outlPosM[1]-1.060660171779821*outlPosP[0]+1.060660171779821*outlPosM[0]; 
    outl[7] += 1.837117307087383*outlPosP[3]-1.837117307087383*(outlPosM[3]+outlPosP[2])+1.837117307087383*outlPosM[2]-1.837117307087383*outlPosP[1]+1.837117307087383*(outlPosM[1]+outlPosP[0])-1.837117307087383*outlPosM[0]; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurfPositivity1x2vSer_Mu_P1(const double m_, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const int *idxl, const int *idxr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:         Cell-center coordinates.
  // dxv[3]:       Cell spacing.
  // idx[3]:       current grid index.
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 
  double limFac, fluxFracL, fluxFracR;

  double diffFac[2]; 
  diffFac[0] = 0.7071067811865475*(BmagInv[1]*nuVtSqSum[1]+BmagInv[0]*nuVtSqSum[0])*m_; 
  diffFac[1] = 0.7071067811865475*(BmagInv[0]*nuVtSqSum[1]+nuVtSqSum[0]*BmagInv[1])*m_; 

  if (idxr[2] == 1) {

  double outrPosP[4], outrPosM[4]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.3333333333333333 : positivityWeightByDirR[2]/positivityWeightByDirR[0]; 
  outrPosP[0] = 0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[6]+((6.0*diffFac[0]-10.39230484541326*diffFac[1])*wr[2]+(5.196152422706631*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[5]+((2.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvr[2])*fr[4]+((6.0*diffFac[1]-10.39230484541326*diffFac[0])*wr[2]+(5.196152422706631*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[3]+((2.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[2]+(6.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[1]+fr[0]*(6.0*diffFac[0]-3.464101615137754*diffFac[1]))*wr[2]+dxvr[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*fr[1]+fr[0]*(1.732050807568877*diffFac[1]-3.0*diffFac[0])))*rdvSq4R; 
  outrPosP[1] = 0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*dxvr[2])*fr[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wr[2]+((-1.732050807568877*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[6]+(((-10.39230484541326*diffFac[1])-6.0*diffFac[0])*wr[2]+(5.196152422706631*diffFac[1]+3.0*diffFac[0])*dxvr[2])*fr[5]+(((-3.464101615137754*diffFac[1])-2.0*diffFac[0])*wr[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvr[2])*fr[4]+(((-6.0*diffFac[1])-10.39230484541326*diffFac[0])*wr[2]+(3.0*diffFac[1]+5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+(((-2.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[2]+(6.0*diffFac[1]+3.464101615137754*diffFac[0])*fr[1]+fr[0]*(3.464101615137754*diffFac[1]+6.0*diffFac[0]))*wr[2]+dxvr[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*fr[1]+fr[0]*((-1.732050807568877*diffFac[1])-3.0*diffFac[0])))*rdvSq4R; 
  outrPosP[2] = -0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[6]+((10.39230484541326*diffFac[1]-6.0*diffFac[0])*wr[2]+(3.0*diffFac[0]-5.196152422706631*diffFac[1])*dxvr[2])*fr[5]+((2.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvr[2])*fr[4]+((10.39230484541326*diffFac[0]-6.0*diffFac[1])*wr[2]+(3.0*diffFac[1]-5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+((2.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[2]+(3.464101615137754*diffFac[0]-6.0*diffFac[1])*fr[1]+fr[0]*(3.464101615137754*diffFac[1]-6.0*diffFac[0]))*wr[2]+dxvr[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fr[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*fr[1]+fr[0]*(3.0*diffFac[0]-1.732050807568877*diffFac[1])))*rdvSq4R; 
  outrPosP[3] = -0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*dxvr[2])*fr[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wr[2]+((-1.732050807568877*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[6]+((10.39230484541326*diffFac[1]+6.0*diffFac[0])*wr[2]+((-5.196152422706631*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[5]+(((-3.464101615137754*diffFac[1])-2.0*diffFac[0])*wr[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvr[2])*fr[4]+((6.0*diffFac[1]+10.39230484541326*diffFac[0])*wr[2]+((-3.0*diffFac[1])-5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+(((-2.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[2]+((-6.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[1]+fr[0]*((-3.464101615137754*diffFac[1])-6.0*diffFac[0]))*wr[2]+dxvr[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fr[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*fr[1]+fr[0]*(1.732050807568877*diffFac[1]+3.0*diffFac[0])))*rdvSq4R; 
  outrPosM[0] = -0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[6]+((6.0*diffFac[0]-10.39230484541326*diffFac[1])*wr[2]+(5.196152422706631*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[5]+((2.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvr[2])*fr[4]+((6.0*diffFac[1]-10.39230484541326*diffFac[0])*wr[2]+(5.196152422706631*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[3]+((2.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[2]+(6.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[1]+fr[0]*(6.0*diffFac[0]-3.464101615137754*diffFac[1]))*wr[2]+dxvr[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*fr[1]+fr[0]*(1.732050807568877*diffFac[1]-3.0*diffFac[0])))*rdvSq4R; 
  outrPosM[1] = -0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*dxvr[2])*fr[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wr[2]+((-1.732050807568877*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[6]+(((-10.39230484541326*diffFac[1])-6.0*diffFac[0])*wr[2]+(5.196152422706631*diffFac[1]+3.0*diffFac[0])*dxvr[2])*fr[5]+(((-3.464101615137754*diffFac[1])-2.0*diffFac[0])*wr[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvr[2])*fr[4]+(((-6.0*diffFac[1])-10.39230484541326*diffFac[0])*wr[2]+(3.0*diffFac[1]+5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+(((-2.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[2]+(6.0*diffFac[1]+3.464101615137754*diffFac[0])*fr[1]+fr[0]*(3.464101615137754*diffFac[1]+6.0*diffFac[0]))*wr[2]+dxvr[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*fr[1]+fr[0]*((-1.732050807568877*diffFac[1])-3.0*diffFac[0])))*rdvSq4R; 
  outrPosM[2] = 0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[6]+((10.39230484541326*diffFac[1]-6.0*diffFac[0])*wr[2]+(3.0*diffFac[0]-5.196152422706631*diffFac[1])*dxvr[2])*fr[5]+((2.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvr[2])*fr[4]+((10.39230484541326*diffFac[0]-6.0*diffFac[1])*wr[2]+(3.0*diffFac[1]-5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+((2.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[2]+(3.464101615137754*diffFac[0]-6.0*diffFac[1])*fr[1]+fr[0]*(3.464101615137754*diffFac[1]-6.0*diffFac[0]))*wr[2]+dxvr[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fr[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*fr[1]+fr[0]*(3.0*diffFac[0]-1.732050807568877*diffFac[1])))*rdvSq4R; 
  outrPosM[3] = 0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*dxvr[2])*fr[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wr[2]+((-1.732050807568877*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[6]+((10.39230484541326*diffFac[1]+6.0*diffFac[0])*wr[2]+((-5.196152422706631*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[5]+(((-3.464101615137754*diffFac[1])-2.0*diffFac[0])*wr[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvr[2])*fr[4]+((6.0*diffFac[1]+10.39230484541326*diffFac[0])*wr[2]+((-3.0*diffFac[1])-5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+(((-2.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[2]+((-6.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[1]+fr[0]*((-3.464101615137754*diffFac[1])-6.0*diffFac[0]))*wr[2]+dxvr[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fr[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*fr[1]+fr[0]*(1.732050807568877*diffFac[1]+3.0*diffFac[0])))*rdvSq4R; 
  if(outrPosM[0] < 0.) limFac = std::min(1.0, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPosM[0]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPosP[0] *= limFac; 
  outrPosM[0] *= limFac; 
  if(outrPosM[1] < 0.) limFac = std::min(1.0, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPosM[1]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPosP[1] *= limFac; 
  outrPosM[1] *= limFac; 
  if(outrPosM[2] < 0.) limFac = std::min(1.0, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*fr[5]-4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPosM[2]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPosP[2] *= limFac; 
  outrPosM[2] *= limFac; 
  if(outrPosM[3] < 0.) limFac = std::min(1.0, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*(fr[6]+fr[5])-4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPosM[3]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPosP[3] *= limFac; 
  outrPosM[3] *= limFac; 
    outr[0] += 0.3535533905932737*(outrPosP[3]+outrPosM[3]+outrPosP[2]+outrPosM[2]+outrPosP[1]+outrPosM[1]+outrPosP[0]+outrPosM[0]); 
    outr[1] += 0.6123724356957944*(outrPosP[3]+outrPosM[3])-0.6123724356957944*(outrPosP[2]+outrPosM[2])+0.6123724356957944*(outrPosP[1]+outrPosM[1])-0.6123724356957944*(outrPosP[0]+outrPosM[0]); 
    outr[2] += 0.6123724356957944*(outrPosP[3]+outrPosM[3]+outrPosP[2]+outrPosM[2])-0.6123724356957944*(outrPosP[1]+outrPosM[1]+outrPosP[0]+outrPosM[0]); 
    outr[3] += 0.6123724356957944*outrPosP[3]-0.6123724356957944*outrPosM[3]+0.6123724356957944*outrPosP[2]-0.6123724356957944*outrPosM[2]+0.6123724356957944*outrPosP[1]-0.6123724356957944*outrPosM[1]+0.6123724356957944*outrPosP[0]-0.6123724356957944*outrPosM[0]; 
    outr[4] += 1.060660171779821*(outrPosP[3]+outrPosM[3])-1.060660171779821*(outrPosP[2]+outrPosM[2]+outrPosP[1]+outrPosM[1])+1.060660171779821*(outrPosP[0]+outrPosM[0]); 
    outr[5] += 1.060660171779821*outrPosP[3]-1.060660171779821*(outrPosM[3]+outrPosP[2])+1.060660171779821*(outrPosM[2]+outrPosP[1])-1.060660171779821*(outrPosM[1]+outrPosP[0])+1.060660171779821*outrPosM[0]; 
    outr[6] += 1.060660171779821*outrPosP[3]-1.060660171779821*outrPosM[3]+1.060660171779821*outrPosP[2]-1.060660171779821*(outrPosM[2]+outrPosP[1])+1.060660171779821*outrPosM[1]-1.060660171779821*outrPosP[0]+1.060660171779821*outrPosM[0]; 
    outr[7] += 1.837117307087383*outrPosP[3]-1.837117307087383*(outrPosM[3]+outrPosP[2])+1.837117307087383*outrPosM[2]-1.837117307087383*outrPosP[1]+1.837117307087383*(outrPosM[1]+outrPosP[0])-1.837117307087383*outrPosM[0]; 

  } else {

  double outlPosP[4], outlPosM[4]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.3333333333333333 : positivityWeightByDirL[2]/positivityWeightByDirL[0]; 
  outlPosP[0] = 0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wl[2]+(3.0*diffFac[0]-1.732050807568877*diffFac[1])*dxvl[2])*fl[6]+((6.0*diffFac[0]-10.39230484541326*diffFac[1])*wl[2]+(3.0*diffFac[0]-5.196152422706631*diffFac[1])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]-2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvl[2])*fl[4]+((6.0*diffFac[1]-10.39230484541326*diffFac[0])*wl[2]+(3.0*diffFac[1]-5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((3.464101615137754*diffFac[0]-2.0*diffFac[1])*fl[2]+(3.464101615137754*diffFac[0]-6.0*diffFac[1])*fl[1]+fl[0]*(3.464101615137754*diffFac[1]-6.0*diffFac[0]))*wl[2]+dxvl[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fl[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*fl[1]+fl[0]*(1.732050807568877*diffFac[1]-3.0*diffFac[0])))*rdvSq4L; 
  outlPosP[1] = 0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[6]+(((-10.39230484541326*diffFac[1])-6.0*diffFac[0])*wl[2]+((-5.196152422706631*diffFac[1])-3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]+2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvl[2])*fl[4]+(((-6.0*diffFac[1])-10.39230484541326*diffFac[0])*wl[2]+((-3.0*diffFac[1])-5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((2.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[2]+((-6.0*diffFac[1])-3.464101615137754*diffFac[0])*fl[1]+fl[0]*((-3.464101615137754*diffFac[1])-6.0*diffFac[0]))*wl[2]+dxvl[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fl[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*fl[1]+fl[0]*((-1.732050807568877*diffFac[1])-3.0*diffFac[0])))*rdvSq4L; 
  outlPosP[2] = -0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wl[2]+(3.0*diffFac[0]-1.732050807568877*diffFac[1])*dxvl[2])*fl[6]+((10.39230484541326*diffFac[1]-6.0*diffFac[0])*wl[2]+(5.196152422706631*diffFac[1]-3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]-2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvl[2])*fl[4]+((10.39230484541326*diffFac[0]-6.0*diffFac[1])*wl[2]+(5.196152422706631*diffFac[0]-3.0*diffFac[1])*dxvl[2])*fl[3]+((3.464101615137754*diffFac[0]-2.0*diffFac[1])*fl[2]+(6.0*diffFac[1]-3.464101615137754*diffFac[0])*fl[1]+fl[0]*(6.0*diffFac[0]-3.464101615137754*diffFac[1]))*wl[2]+dxvl[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*fl[1]+fl[0]*(3.0*diffFac[0]-1.732050807568877*diffFac[1])))*rdvSq4L; 
  outlPosP[3] = -0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[6]+((10.39230484541326*diffFac[1]+6.0*diffFac[0])*wl[2]+(5.196152422706631*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]+2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvl[2])*fl[4]+((6.0*diffFac[1]+10.39230484541326*diffFac[0])*wl[2]+(3.0*diffFac[1]+5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((2.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[2]+(6.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[1]+fl[0]*(3.464101615137754*diffFac[1]+6.0*diffFac[0]))*wl[2]+dxvl[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*fl[1]+fl[0]*(1.732050807568877*diffFac[1]+3.0*diffFac[0])))*rdvSq4L; 
  outlPosM[0] = -0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wl[2]+(3.0*diffFac[0]-1.732050807568877*diffFac[1])*dxvl[2])*fl[6]+((6.0*diffFac[0]-10.39230484541326*diffFac[1])*wl[2]+(3.0*diffFac[0]-5.196152422706631*diffFac[1])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]-2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvl[2])*fl[4]+((6.0*diffFac[1]-10.39230484541326*diffFac[0])*wl[2]+(3.0*diffFac[1]-5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((3.464101615137754*diffFac[0]-2.0*diffFac[1])*fl[2]+(3.464101615137754*diffFac[0]-6.0*diffFac[1])*fl[1]+fl[0]*(3.464101615137754*diffFac[1]-6.0*diffFac[0]))*wl[2]+dxvl[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fl[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*fl[1]+fl[0]*(1.732050807568877*diffFac[1]-3.0*diffFac[0])))*rdvSq4L; 
  outlPosM[1] = -0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[6]+(((-10.39230484541326*diffFac[1])-6.0*diffFac[0])*wl[2]+((-5.196152422706631*diffFac[1])-3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]+2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvl[2])*fl[4]+(((-6.0*diffFac[1])-10.39230484541326*diffFac[0])*wl[2]+((-3.0*diffFac[1])-5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((2.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[2]+((-6.0*diffFac[1])-3.464101615137754*diffFac[0])*fl[1]+fl[0]*((-3.464101615137754*diffFac[1])-6.0*diffFac[0]))*wl[2]+dxvl[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fl[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*fl[1]+fl[0]*((-1.732050807568877*diffFac[1])-3.0*diffFac[0])))*rdvSq4L; 
  outlPosM[2] = 0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wl[2]+(3.0*diffFac[0]-1.732050807568877*diffFac[1])*dxvl[2])*fl[6]+((10.39230484541326*diffFac[1]-6.0*diffFac[0])*wl[2]+(5.196152422706631*diffFac[1]-3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]-2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvl[2])*fl[4]+((10.39230484541326*diffFac[0]-6.0*diffFac[1])*wl[2]+(5.196152422706631*diffFac[0]-3.0*diffFac[1])*dxvl[2])*fl[3]+((3.464101615137754*diffFac[0]-2.0*diffFac[1])*fl[2]+(6.0*diffFac[1]-3.464101615137754*diffFac[0])*fl[1]+fl[0]*(6.0*diffFac[0]-3.464101615137754*diffFac[1]))*wl[2]+dxvl[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*fl[1]+fl[0]*(3.0*diffFac[0]-1.732050807568877*diffFac[1])))*rdvSq4L; 
  outlPosM[3] = 0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[6]+((10.39230484541326*diffFac[1]+6.0*diffFac[0])*wl[2]+(5.196152422706631*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]+2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvl[2])*fl[4]+((6.0*diffFac[1]+10.39230484541326*diffFac[0])*wl[2]+(3.0*diffFac[1]+5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((2.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[2]+(6.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[1]+fl[0]*(3.464101615137754*diffFac[1]+6.0*diffFac[0]))*wl[2]+dxvl[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*fl[1]+fl[0]*(1.732050807568877*diffFac[1]+3.0*diffFac[0])))*rdvSq4L; 
  if(outlPosP[0] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*(fl[6]+fl[5])+4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPosP[0]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPosP[0] *= limFac; 
  outlPosM[0] *= limFac; 
  if(outlPosP[1] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*fl[5]+4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPosP[1]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPosP[1] *= limFac; 
  outlPosM[1] *= limFac; 
  if(outlPosP[2] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPosP[2]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPosP[2] *= limFac; 
  outlPosM[2] *= limFac; 
  if(outlPosP[3] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPosP[3]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPosP[3] *= limFac; 
  outlPosM[3] *= limFac; 
    outl[0] += 0.3535533905932737*(outlPosP[3]+outlPosM[3]+outlPosP[2]+outlPosM[2]+outlPosP[1]+outlPosM[1]+outlPosP[0]+outlPosM[0]); 
    outl[1] += 0.6123724356957944*(outlPosP[3]+outlPosM[3])-0.6123724356957944*(outlPosP[2]+outlPosM[2])+0.6123724356957944*(outlPosP[1]+outlPosM[1])-0.6123724356957944*(outlPosP[0]+outlPosM[0]); 
    outl[2] += 0.6123724356957944*(outlPosP[3]+outlPosM[3]+outlPosP[2]+outlPosM[2])-0.6123724356957944*(outlPosP[1]+outlPosM[1]+outlPosP[0]+outlPosM[0]); 
    outl[3] += 0.6123724356957944*outlPosP[3]-0.6123724356957944*outlPosM[3]+0.6123724356957944*outlPosP[2]-0.6123724356957944*outlPosM[2]+0.6123724356957944*outlPosP[1]-0.6123724356957944*outlPosM[1]+0.6123724356957944*outlPosP[0]-0.6123724356957944*outlPosM[0]; 
    outl[4] += 1.060660171779821*(outlPosP[3]+outlPosM[3])-1.060660171779821*(outlPosP[2]+outlPosM[2]+outlPosP[1]+outlPosM[1])+1.060660171779821*(outlPosP[0]+outlPosM[0]); 
    outl[5] += 1.060660171779821*outlPosP[3]-1.060660171779821*(outlPosM[3]+outlPosP[2])+1.060660171779821*(outlPosM[2]+outlPosP[1])-1.060660171779821*(outlPosM[1]+outlPosP[0])+1.060660171779821*outlPosM[0]; 
    outl[6] += 1.060660171779821*outlPosP[3]-1.060660171779821*outlPosM[3]+1.060660171779821*outlPosP[2]-1.060660171779821*(outlPosM[2]+outlPosP[1])+1.060660171779821*outlPosM[1]-1.060660171779821*outlPosP[0]+1.060660171779821*outlPosM[0]; 
    outl[7] += 1.837117307087383*outlPosP[3]-1.837117307087383*(outlPosM[3]+outlPosP[2])+1.837117307087383*outlPosM[2]-1.837117307087383*outlPosP[1]+1.837117307087383*(outlPosM[1]+outlPosP[0])-1.837117307087383*outlPosM[0]; 

  }
  return 0.0; 
} 
