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

  double outrPos[8]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.3333333333333333 : positivityWeightByDirR[2]/positivityWeightByDirR[0]; 
  outrPos[0] = -0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[6]+(nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[5]+(3.0*nuVtSqSum[0]-5.196152422706631*nuVtSqSum[1])*fr[4]+(nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[3]+(3.0*nuVtSqSum[1]-5.196152422706631*nuVtSqSum[0])*fr[2]+(3.0*fr[1]-1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*(3.0*fr[0]-1.732050807568877*fr[1]))*rdvSq4R; 
  outrPos[1] = -0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fr[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[6]+((-1.732050807568877*nuVtSqSum[1])-1.0*nuVtSqSum[0])*fr[5]+((-5.196152422706631*nuVtSqSum[1])-3.0*nuVtSqSum[0])*fr[4]+((-1.0*nuVtSqSum[1])-1.732050807568877*nuVtSqSum[0])*fr[3]+((-3.0*nuVtSqSum[1])-5.196152422706631*nuVtSqSum[0])*fr[2]+(3.0*fr[1]+1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fr[1]+3.0*fr[0]))*rdvSq4R; 
  outrPos[2] = 0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[6]+(nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[5]+(3.0*nuVtSqSum[0]-5.196152422706631*nuVtSqSum[1])*fr[4]+(nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[3]+(3.0*nuVtSqSum[1]-5.196152422706631*nuVtSqSum[0])*fr[2]+(3.0*fr[1]-1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*(3.0*fr[0]-1.732050807568877*fr[1]))*rdvSq4R; 
  outrPos[3] = 0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fr[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[6]+((-1.732050807568877*nuVtSqSum[1])-1.0*nuVtSqSum[0])*fr[5]+((-5.196152422706631*nuVtSqSum[1])-3.0*nuVtSqSum[0])*fr[4]+((-1.0*nuVtSqSum[1])-1.732050807568877*nuVtSqSum[0])*fr[3]+((-3.0*nuVtSqSum[1])-5.196152422706631*nuVtSqSum[0])*fr[2]+(3.0*fr[1]+1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fr[1]+3.0*fr[0]))*rdvSq4R; 
  outrPos[4] = 0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[6]+(nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[5]+(5.196152422706631*nuVtSqSum[1]-3.0*nuVtSqSum[0])*fr[4]+(nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[3]+(5.196152422706631*nuVtSqSum[0]-3.0*nuVtSqSum[1])*fr[2]+(1.732050807568877*fr[0]-3.0*fr[1])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fr[1]-3.0*fr[0]))*rdvSq4R; 
  outrPos[5] = 0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fr[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[6]+((-1.732050807568877*nuVtSqSum[1])-1.0*nuVtSqSum[0])*fr[5]+(5.196152422706631*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[4]+((-1.0*nuVtSqSum[1])-1.732050807568877*nuVtSqSum[0])*fr[3]+(3.0*nuVtSqSum[1]+5.196152422706631*nuVtSqSum[0])*fr[2]+((-3.0*fr[1])-1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-1.732050807568877*fr[1])-3.0*fr[0]))*rdvSq4R; 
  outrPos[6] = -0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[6]+(nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fr[5]+(5.196152422706631*nuVtSqSum[1]-3.0*nuVtSqSum[0])*fr[4]+(nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fr[3]+(5.196152422706631*nuVtSqSum[0]-3.0*nuVtSqSum[1])*fr[2]+(1.732050807568877*fr[0]-3.0*fr[1])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fr[1]-3.0*fr[0]))*rdvSq4R; 
  outrPos[7] = -0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fr[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[6]+((-1.732050807568877*nuVtSqSum[1])-1.0*nuVtSqSum[0])*fr[5]+(5.196152422706631*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fr[4]+((-1.0*nuVtSqSum[1])-1.732050807568877*nuVtSqSum[0])*fr[3]+(3.0*nuVtSqSum[1]+5.196152422706631*nuVtSqSum[0])*fr[2]+((-3.0*fr[1])-1.732050807568877*fr[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-1.732050807568877*fr[1])-3.0*fr[0]))*rdvSq4R; 
  if(outrPos[0] < 0.) limFac = std::min(1.0, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPos[0]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPos[2] *= limFac; 
  outrPos[0] *= limFac; 
  if(outrPos[1] < 0.) limFac = std::min(1.0, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPos[1]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPos[3] *= limFac; 
  outrPos[1] *= limFac; 
  if(outrPos[4] < 0.) limFac = std::min(1.0, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0]))/dtApprox/outrPos[4]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPos[6] *= limFac; 
  outrPos[4] *= limFac; 
  if(outrPos[5] < 0.) limFac = std::min(1.0, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*fr[5]+4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]-12.72792206135786*fr[0]))/dtApprox/outrPos[5]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPos[7] *= limFac; 
  outrPos[5] *= limFac; 
    outr[0] += 0.3535533905932737*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
    outr[1] += 0.6123724356957944*outrPos[7]-0.6123724356957944*outrPos[6]+0.6123724356957944*outrPos[5]-0.6123724356957944*outrPos[4]+0.6123724356957944*outrPos[3]-0.6123724356957944*outrPos[2]+0.6123724356957944*outrPos[1]-0.6123724356957944*outrPos[0]; 
    outr[2] += 0.6123724356957944*(outrPos[7]+outrPos[6])-0.6123724356957944*(outrPos[5]+outrPos[4])+0.6123724356957944*(outrPos[3]+outrPos[2])-0.6123724356957944*(outrPos[1]+outrPos[0]); 
    outr[3] += 0.6123724356957944*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])-0.6123724356957944*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
    outr[4] += 1.060660171779821*outrPos[7]-1.060660171779821*(outrPos[6]+outrPos[5])+1.060660171779821*(outrPos[4]+outrPos[3])-1.060660171779821*(outrPos[2]+outrPos[1])+1.060660171779821*outrPos[0]; 
    outr[5] += 1.060660171779821*outrPos[7]-1.060660171779821*outrPos[6]+1.060660171779821*outrPos[5]-1.060660171779821*(outrPos[4]+outrPos[3])+1.060660171779821*outrPos[2]-1.060660171779821*outrPos[1]+1.060660171779821*outrPos[0]; 
    outr[6] += 1.060660171779821*(outrPos[7]+outrPos[6])-1.060660171779821*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+1.060660171779821*(outrPos[1]+outrPos[0]); 
    outr[7] += 1.837117307087383*outrPos[7]-1.837117307087383*(outrPos[6]+outrPos[5])+1.837117307087383*outrPos[4]-1.837117307087383*outrPos[3]+1.837117307087383*(outrPos[2]+outrPos[1])-1.837117307087383*outrPos[0]; 

  } else {

  double outlPos[8]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.3333333333333333 : positivityWeightByDirL[2]/positivityWeightByDirL[0]; 
  outlPos[0] = -0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fl[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fl[6]+(1.732050807568877*nuVtSqSum[1]-1.0*nuVtSqSum[0])*fl[5]+(3.0*nuVtSqSum[0]-5.196152422706631*nuVtSqSum[1])*fl[4]+(1.732050807568877*nuVtSqSum[0]-1.0*nuVtSqSum[1])*fl[3]+(3.0*nuVtSqSum[1]-5.196152422706631*nuVtSqSum[0])*fl[2]+(1.732050807568877*fl[0]-3.0*fl[1])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fl[1]-3.0*fl[0]))*rdvSq4L; 
  outlPos[1] = -0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[6]+(1.732050807568877*nuVtSqSum[1]+nuVtSqSum[0])*fl[5]+((-5.196152422706631*nuVtSqSum[1])-3.0*nuVtSqSum[0])*fl[4]+(nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[3]+((-3.0*nuVtSqSum[1])-5.196152422706631*nuVtSqSum[0])*fl[2]+((-3.0*fl[1])-1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-1.732050807568877*fl[1])-3.0*fl[0]))*rdvSq4L; 
  outlPos[2] = 0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fl[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fl[6]+(1.732050807568877*nuVtSqSum[1]-1.0*nuVtSqSum[0])*fl[5]+(3.0*nuVtSqSum[0]-5.196152422706631*nuVtSqSum[1])*fl[4]+(1.732050807568877*nuVtSqSum[0]-1.0*nuVtSqSum[1])*fl[3]+(3.0*nuVtSqSum[1]-5.196152422706631*nuVtSqSum[0])*fl[2]+(1.732050807568877*fl[0]-3.0*fl[1])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fl[1]-3.0*fl[0]))*rdvSq4L; 
  outlPos[3] = 0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[6]+(1.732050807568877*nuVtSqSum[1]+nuVtSqSum[0])*fl[5]+((-5.196152422706631*nuVtSqSum[1])-3.0*nuVtSqSum[0])*fl[4]+(nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[3]+((-3.0*nuVtSqSum[1])-5.196152422706631*nuVtSqSum[0])*fl[2]+((-3.0*fl[1])-1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*((-1.732050807568877*fl[1])-3.0*fl[0]))*rdvSq4L; 
  outlPos[4] = 0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fl[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fl[6]+(1.732050807568877*nuVtSqSum[1]-1.0*nuVtSqSum[0])*fl[5]+(5.196152422706631*nuVtSqSum[1]-3.0*nuVtSqSum[0])*fl[4]+(1.732050807568877*nuVtSqSum[0]-1.0*nuVtSqSum[1])*fl[3]+(5.196152422706631*nuVtSqSum[0]-3.0*nuVtSqSum[1])*fl[2]+(3.0*fl[1]-1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(3.0*fl[0]-1.732050807568877*fl[1]))*rdvSq4L; 
  outlPos[5] = 0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[6]+(1.732050807568877*nuVtSqSum[1]+nuVtSqSum[0])*fl[5]+(5.196152422706631*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[4]+(nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[3]+(3.0*nuVtSqSum[1]+5.196152422706631*nuVtSqSum[0])*fl[2]+(3.0*fl[1]+1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fl[1]+3.0*fl[0]))*rdvSq4L; 
  outlPos[6] = -0.04166666666666666*((3.0*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0])*fl[7]+(3.0*nuVtSqSum[0]-1.732050807568877*nuVtSqSum[1])*fl[6]+(1.732050807568877*nuVtSqSum[1]-1.0*nuVtSqSum[0])*fl[5]+(5.196152422706631*nuVtSqSum[1]-3.0*nuVtSqSum[0])*fl[4]+(1.732050807568877*nuVtSqSum[0]-1.0*nuVtSqSum[1])*fl[3]+(5.196152422706631*nuVtSqSum[0]-3.0*nuVtSqSum[1])*fl[2]+(3.0*fl[1]-1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(3.0*fl[0]-1.732050807568877*fl[1]))*rdvSq4L; 
  outlPos[7] = -0.04166666666666666*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[7]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[6]+(1.732050807568877*nuVtSqSum[1]+nuVtSqSum[0])*fl[5]+(5.196152422706631*nuVtSqSum[1]+3.0*nuVtSqSum[0])*fl[4]+(nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*fl[3]+(3.0*nuVtSqSum[1]+5.196152422706631*nuVtSqSum[0])*fl[2]+(3.0*fl[1]+1.732050807568877*fl[0])*nuVtSqSum[1]+nuVtSqSum[0]*(1.732050807568877*fl[1]+3.0*fl[0]))*rdvSq4L; 
  if(outlPos[2] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*fl[5]-4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]+12.72792206135786*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outlPos[0] *= limFac; 
  if(outlPos[3] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[1] *= limFac; 
  if(outlPos[6] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPos[6]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPos[6] *= limFac; 
  outlPos[4] *= limFac; 
  if(outlPos[7] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[5] *= limFac; 
    outl[0] += 0.3535533905932737*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
    outl[1] += 0.6123724356957944*outlPos[7]-0.6123724356957944*outlPos[6]+0.6123724356957944*outlPos[5]-0.6123724356957944*outlPos[4]+0.6123724356957944*outlPos[3]-0.6123724356957944*outlPos[2]+0.6123724356957944*outlPos[1]-0.6123724356957944*outlPos[0]; 
    outl[2] += 0.6123724356957944*(outlPos[7]+outlPos[6])-0.6123724356957944*(outlPos[5]+outlPos[4])+0.6123724356957944*(outlPos[3]+outlPos[2])-0.6123724356957944*(outlPos[1]+outlPos[0]); 
    outl[3] += 0.6123724356957944*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])-0.6123724356957944*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
    outl[4] += 1.060660171779821*outlPos[7]-1.060660171779821*(outlPos[6]+outlPos[5])+1.060660171779821*(outlPos[4]+outlPos[3])-1.060660171779821*(outlPos[2]+outlPos[1])+1.060660171779821*outlPos[0]; 
    outl[5] += 1.060660171779821*outlPos[7]-1.060660171779821*outlPos[6]+1.060660171779821*outlPos[5]-1.060660171779821*(outlPos[4]+outlPos[3])+1.060660171779821*outlPos[2]-1.060660171779821*outlPos[1]+1.060660171779821*outlPos[0]; 
    outl[6] += 1.060660171779821*(outlPos[7]+outlPos[6])-1.060660171779821*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+1.060660171779821*(outlPos[1]+outlPos[0]); 
    outl[7] += 1.837117307087383*outlPos[7]-1.837117307087383*(outlPos[6]+outlPos[5])+1.837117307087383*outlPos[4]-1.837117307087383*outlPos[3]+1.837117307087383*(outlPos[2]+outlPos[1])-1.837117307087383*outlPos[0]; 

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

  double outrPos[8]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.3333333333333333 : positivityWeightByDirR[3]/positivityWeightByDirR[0]; 
  outrPos[0] = -0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[6]+((6.0*diffFac[0]-10.39230484541326*diffFac[1])*wr[2]+(5.196152422706631*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[5]+((2.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvr[2])*fr[4]+((6.0*diffFac[1]-10.39230484541326*diffFac[0])*wr[2]+(5.196152422706631*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[3]+((2.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[2]+(6.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[1]+fr[0]*(6.0*diffFac[0]-3.464101615137754*diffFac[1]))*wr[2]+dxvr[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*fr[1]+fr[0]*(1.732050807568877*diffFac[1]-3.0*diffFac[0])))*rdvSq4R; 
  outrPos[1] = -0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*dxvr[2])*fr[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wr[2]+((-1.732050807568877*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[6]+(((-10.39230484541326*diffFac[1])-6.0*diffFac[0])*wr[2]+(5.196152422706631*diffFac[1]+3.0*diffFac[0])*dxvr[2])*fr[5]+(((-3.464101615137754*diffFac[1])-2.0*diffFac[0])*wr[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvr[2])*fr[4]+(((-6.0*diffFac[1])-10.39230484541326*diffFac[0])*wr[2]+(3.0*diffFac[1]+5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+(((-2.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[2]+(6.0*diffFac[1]+3.464101615137754*diffFac[0])*fr[1]+fr[0]*(3.464101615137754*diffFac[1]+6.0*diffFac[0]))*wr[2]+dxvr[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*fr[1]+fr[0]*((-1.732050807568877*diffFac[1])-3.0*diffFac[0])))*rdvSq4R; 
  outrPos[2] = 0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[6]+((10.39230484541326*diffFac[1]-6.0*diffFac[0])*wr[2]+(3.0*diffFac[0]-5.196152422706631*diffFac[1])*dxvr[2])*fr[5]+((2.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvr[2])*fr[4]+((10.39230484541326*diffFac[0]-6.0*diffFac[1])*wr[2]+(3.0*diffFac[1]-5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+((2.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[2]+(3.464101615137754*diffFac[0]-6.0*diffFac[1])*fr[1]+fr[0]*(3.464101615137754*diffFac[1]-6.0*diffFac[0]))*wr[2]+dxvr[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fr[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*fr[1]+fr[0]*(3.0*diffFac[0]-1.732050807568877*diffFac[1])))*rdvSq4R; 
  outrPos[3] = 0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*dxvr[2])*fr[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wr[2]+((-1.732050807568877*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[6]+((10.39230484541326*diffFac[1]+6.0*diffFac[0])*wr[2]+((-5.196152422706631*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[5]+(((-3.464101615137754*diffFac[1])-2.0*diffFac[0])*wr[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvr[2])*fr[4]+((6.0*diffFac[1]+10.39230484541326*diffFac[0])*wr[2]+((-3.0*diffFac[1])-5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+(((-2.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[2]+((-6.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[1]+fr[0]*((-3.464101615137754*diffFac[1])-6.0*diffFac[0]))*wr[2]+dxvr[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fr[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*fr[1]+fr[0]*(1.732050807568877*diffFac[1]+3.0*diffFac[0])))*rdvSq4R; 
  outrPos[4] = 0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[6]+((6.0*diffFac[0]-10.39230484541326*diffFac[1])*wr[2]+(5.196152422706631*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[5]+((2.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvr[2])*fr[4]+((6.0*diffFac[1]-10.39230484541326*diffFac[0])*wr[2]+(5.196152422706631*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[3]+((2.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[2]+(6.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[1]+fr[0]*(6.0*diffFac[0]-3.464101615137754*diffFac[1]))*wr[2]+dxvr[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*fr[1]+fr[0]*(1.732050807568877*diffFac[1]-3.0*diffFac[0])))*rdvSq4R; 
  outrPos[5] = 0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*dxvr[2])*fr[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wr[2]+((-1.732050807568877*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[6]+(((-10.39230484541326*diffFac[1])-6.0*diffFac[0])*wr[2]+(5.196152422706631*diffFac[1]+3.0*diffFac[0])*dxvr[2])*fr[5]+(((-3.464101615137754*diffFac[1])-2.0*diffFac[0])*wr[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvr[2])*fr[4]+(((-6.0*diffFac[1])-10.39230484541326*diffFac[0])*wr[2]+(3.0*diffFac[1]+5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+(((-2.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[2]+(6.0*diffFac[1]+3.464101615137754*diffFac[0])*fr[1]+fr[0]*(3.464101615137754*diffFac[1]+6.0*diffFac[0]))*wr[2]+dxvr[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*fr[1]+fr[0]*((-1.732050807568877*diffFac[1])-3.0*diffFac[0])))*rdvSq4R; 
  outrPos[6] = -0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wr[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*dxvr[2])*fr[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-3.0*diffFac[0])*dxvr[2])*fr[6]+((10.39230484541326*diffFac[1]-6.0*diffFac[0])*wr[2]+(3.0*diffFac[0]-5.196152422706631*diffFac[1])*dxvr[2])*fr[5]+((2.0*diffFac[0]-3.464101615137754*diffFac[1])*wr[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvr[2])*fr[4]+((10.39230484541326*diffFac[0]-6.0*diffFac[1])*wr[2]+(3.0*diffFac[1]-5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+((2.0*diffFac[1]-3.464101615137754*diffFac[0])*fr[2]+(3.464101615137754*diffFac[0]-6.0*diffFac[1])*fr[1]+fr[0]*(3.464101615137754*diffFac[1]-6.0*diffFac[0]))*wr[2]+dxvr[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fr[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*fr[1]+fr[0]*(3.0*diffFac[0]-1.732050807568877*diffFac[1])))*rdvSq4R; 
  outrPos[7] = -0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wr[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*dxvr[2])*fr[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wr[2]+((-1.732050807568877*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[6]+((10.39230484541326*diffFac[1]+6.0*diffFac[0])*wr[2]+((-5.196152422706631*diffFac[1])-3.0*diffFac[0])*dxvr[2])*fr[5]+(((-3.464101615137754*diffFac[1])-2.0*diffFac[0])*wr[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvr[2])*fr[4]+((6.0*diffFac[1]+10.39230484541326*diffFac[0])*wr[2]+((-3.0*diffFac[1])-5.196152422706631*diffFac[0])*dxvr[2])*fr[3]+(((-2.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[2]+((-6.0*diffFac[1])-3.464101615137754*diffFac[0])*fr[1]+fr[0]*((-3.464101615137754*diffFac[1])-6.0*diffFac[0]))*wr[2]+dxvr[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fr[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*fr[1]+fr[0]*(1.732050807568877*diffFac[1]+3.0*diffFac[0])))*rdvSq4R; 
  if(outrPos[0] < 0.) limFac = std::min(1.0, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPos[0]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPos[4] *= limFac; 
  outrPos[0] *= limFac; 
  if(outrPos[1] < 0.) limFac = std::min(1.0, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPos[1]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPos[5] *= limFac; 
  outrPos[1] *= limFac; 
  if(outrPos[2] < 0.) limFac = std::min(1.0, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*fr[5]-4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPos[2]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPos[6] *= limFac; 
  outrPos[2] *= limFac; 
  if(outrPos[3] < 0.) limFac = std::min(1.0, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*(fr[6]+fr[5])-4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPos[3]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outrPos[7] *= limFac; 
  outrPos[3] *= limFac; 
    outr[0] += 0.3535533905932737*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
    outr[1] += 0.6123724356957944*outrPos[7]-0.6123724356957944*outrPos[6]+0.6123724356957944*outrPos[5]-0.6123724356957944*outrPos[4]+0.6123724356957944*outrPos[3]-0.6123724356957944*outrPos[2]+0.6123724356957944*outrPos[1]-0.6123724356957944*outrPos[0]; 
    outr[2] += 0.6123724356957944*(outrPos[7]+outrPos[6])-0.6123724356957944*(outrPos[5]+outrPos[4])+0.6123724356957944*(outrPos[3]+outrPos[2])-0.6123724356957944*(outrPos[1]+outrPos[0]); 
    outr[3] += 0.6123724356957944*(outrPos[7]+outrPos[6]+outrPos[5]+outrPos[4])-0.6123724356957944*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
    outr[4] += 1.060660171779821*outrPos[7]-1.060660171779821*(outrPos[6]+outrPos[5])+1.060660171779821*(outrPos[4]+outrPos[3])-1.060660171779821*(outrPos[2]+outrPos[1])+1.060660171779821*outrPos[0]; 
    outr[5] += 1.060660171779821*outrPos[7]-1.060660171779821*outrPos[6]+1.060660171779821*outrPos[5]-1.060660171779821*(outrPos[4]+outrPos[3])+1.060660171779821*outrPos[2]-1.060660171779821*outrPos[1]+1.060660171779821*outrPos[0]; 
    outr[6] += 1.060660171779821*(outrPos[7]+outrPos[6])-1.060660171779821*(outrPos[5]+outrPos[4]+outrPos[3]+outrPos[2])+1.060660171779821*(outrPos[1]+outrPos[0]); 
    outr[7] += 1.837117307087383*outrPos[7]-1.837117307087383*(outrPos[6]+outrPos[5])+1.837117307087383*outrPos[4]-1.837117307087383*outrPos[3]+1.837117307087383*(outrPos[2]+outrPos[1])-1.837117307087383*outrPos[0]; 

  } else {

  double outlPos[8]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.3333333333333333 : positivityWeightByDirL[3]/positivityWeightByDirL[0]; 
  outlPos[0] = -0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wl[2]+(3.0*diffFac[0]-1.732050807568877*diffFac[1])*dxvl[2])*fl[6]+((6.0*diffFac[0]-10.39230484541326*diffFac[1])*wl[2]+(3.0*diffFac[0]-5.196152422706631*diffFac[1])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]-2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvl[2])*fl[4]+((6.0*diffFac[1]-10.39230484541326*diffFac[0])*wl[2]+(3.0*diffFac[1]-5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((3.464101615137754*diffFac[0]-2.0*diffFac[1])*fl[2]+(3.464101615137754*diffFac[0]-6.0*diffFac[1])*fl[1]+fl[0]*(3.464101615137754*diffFac[1]-6.0*diffFac[0]))*wl[2]+dxvl[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fl[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*fl[1]+fl[0]*(1.732050807568877*diffFac[1]-3.0*diffFac[0])))*rdvSq4L; 
  outlPos[1] = -0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[6]+(((-10.39230484541326*diffFac[1])-6.0*diffFac[0])*wl[2]+((-5.196152422706631*diffFac[1])-3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]+2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvl[2])*fl[4]+(((-6.0*diffFac[1])-10.39230484541326*diffFac[0])*wl[2]+((-3.0*diffFac[1])-5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((2.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[2]+((-6.0*diffFac[1])-3.464101615137754*diffFac[0])*fl[1]+fl[0]*((-3.464101615137754*diffFac[1])-6.0*diffFac[0]))*wl[2]+dxvl[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fl[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*fl[1]+fl[0]*((-1.732050807568877*diffFac[1])-3.0*diffFac[0])))*rdvSq4L; 
  outlPos[2] = 0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wl[2]+(3.0*diffFac[0]-1.732050807568877*diffFac[1])*dxvl[2])*fl[6]+((10.39230484541326*diffFac[1]-6.0*diffFac[0])*wl[2]+(5.196152422706631*diffFac[1]-3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]-2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvl[2])*fl[4]+((10.39230484541326*diffFac[0]-6.0*diffFac[1])*wl[2]+(5.196152422706631*diffFac[0]-3.0*diffFac[1])*dxvl[2])*fl[3]+((3.464101615137754*diffFac[0]-2.0*diffFac[1])*fl[2]+(6.0*diffFac[1]-3.464101615137754*diffFac[0])*fl[1]+fl[0]*(6.0*diffFac[0]-3.464101615137754*diffFac[1]))*wl[2]+dxvl[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*fl[1]+fl[0]*(3.0*diffFac[0]-1.732050807568877*diffFac[1])))*rdvSq4L; 
  outlPos[3] = 0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[6]+((10.39230484541326*diffFac[1]+6.0*diffFac[0])*wl[2]+(5.196152422706631*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]+2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvl[2])*fl[4]+((6.0*diffFac[1]+10.39230484541326*diffFac[0])*wl[2]+(3.0*diffFac[1]+5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((2.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[2]+(6.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[1]+fl[0]*(3.464101615137754*diffFac[1]+6.0*diffFac[0]))*wl[2]+dxvl[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*fl[1]+fl[0]*(1.732050807568877*diffFac[1]+3.0*diffFac[0])))*rdvSq4L; 
  outlPos[4] = 0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wl[2]+(3.0*diffFac[0]-1.732050807568877*diffFac[1])*dxvl[2])*fl[6]+((6.0*diffFac[0]-10.39230484541326*diffFac[1])*wl[2]+(3.0*diffFac[0]-5.196152422706631*diffFac[1])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]-2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvl[2])*fl[4]+((6.0*diffFac[1]-10.39230484541326*diffFac[0])*wl[2]+(3.0*diffFac[1]-5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((3.464101615137754*diffFac[0]-2.0*diffFac[1])*fl[2]+(3.464101615137754*diffFac[0]-6.0*diffFac[1])*fl[1]+fl[0]*(3.464101615137754*diffFac[1]-6.0*diffFac[0]))*wl[2]+dxvl[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fl[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*fl[1]+fl[0]*(1.732050807568877*diffFac[1]-3.0*diffFac[0])))*rdvSq4L; 
  outlPos[5] = 0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[6]+(((-10.39230484541326*diffFac[1])-6.0*diffFac[0])*wl[2]+((-5.196152422706631*diffFac[1])-3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]+2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvl[2])*fl[4]+(((-6.0*diffFac[1])-10.39230484541326*diffFac[0])*wl[2]+((-3.0*diffFac[1])-5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((2.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[2]+((-6.0*diffFac[1])-3.464101615137754*diffFac[0])*fl[1]+fl[0]*((-3.464101615137754*diffFac[1])-6.0*diffFac[0]))*wl[2]+dxvl[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fl[2]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*fl[1]+fl[0]*((-1.732050807568877*diffFac[1])-3.0*diffFac[0])))*rdvSq4L; 
  outlPos[6] = -0.04166666666666666*(((6.0*diffFac[1]-3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((6.0*diffFac[0]-3.464101615137754*diffFac[1])*wl[2]+(3.0*diffFac[0]-1.732050807568877*diffFac[1])*dxvl[2])*fl[6]+((10.39230484541326*diffFac[1]-6.0*diffFac[0])*wl[2]+(5.196152422706631*diffFac[1]-3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]-2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]-1.0*diffFac[0])*dxvl[2])*fl[4]+((10.39230484541326*diffFac[0]-6.0*diffFac[1])*wl[2]+(5.196152422706631*diffFac[0]-3.0*diffFac[1])*dxvl[2])*fl[3]+((3.464101615137754*diffFac[0]-2.0*diffFac[1])*fl[2]+(6.0*diffFac[1]-3.464101615137754*diffFac[0])*fl[1]+fl[0]*(6.0*diffFac[0]-3.464101615137754*diffFac[1]))*wl[2]+dxvl[2]*((1.732050807568877*diffFac[0]-1.0*diffFac[1])*fl[2]+(3.0*diffFac[1]-1.732050807568877*diffFac[0])*fl[1]+fl[0]*(3.0*diffFac[0]-1.732050807568877*diffFac[1])))*rdvSq4L; 
  outlPos[7] = -0.04166666666666666*(((6.0*diffFac[1]+3.464101615137754*diffFac[0])*wl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*dxvl[2])*fl[7]+((3.464101615137754*diffFac[1]+6.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[6]+((10.39230484541326*diffFac[1]+6.0*diffFac[0])*wl[2]+(5.196152422706631*diffFac[1]+3.0*diffFac[0])*dxvl[2])*fl[5]+((3.464101615137754*diffFac[1]+2.0*diffFac[0])*wl[2]+(1.732050807568877*diffFac[1]+diffFac[0])*dxvl[2])*fl[4]+((6.0*diffFac[1]+10.39230484541326*diffFac[0])*wl[2]+(3.0*diffFac[1]+5.196152422706631*diffFac[0])*dxvl[2])*fl[3]+((2.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[2]+(6.0*diffFac[1]+3.464101615137754*diffFac[0])*fl[1]+fl[0]*(3.464101615137754*diffFac[1]+6.0*diffFac[0]))*wl[2]+dxvl[2]*((diffFac[1]+1.732050807568877*diffFac[0])*fl[2]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*fl[1]+fl[0]*(1.732050807568877*diffFac[1]+3.0*diffFac[0])))*rdvSq4L; 
  if(outlPos[4] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*(fl[6]+fl[5])+4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPos[4]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPos[4] *= limFac; 
  outlPos[0] *= limFac; 
  if(outlPos[5] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*fl[5]+4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPos[5]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPos[5] *= limFac; 
  outlPos[1] *= limFac; 
  if(outlPos[6] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPos[6]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPos[6] *= limFac; 
  outlPos[2] *= limFac; 
  if(outlPos[7] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 0.;
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outlPos[3] *= limFac; 
    outl[0] += 0.3535533905932737*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
    outl[1] += 0.6123724356957944*outlPos[7]-0.6123724356957944*outlPos[6]+0.6123724356957944*outlPos[5]-0.6123724356957944*outlPos[4]+0.6123724356957944*outlPos[3]-0.6123724356957944*outlPos[2]+0.6123724356957944*outlPos[1]-0.6123724356957944*outlPos[0]; 
    outl[2] += 0.6123724356957944*(outlPos[7]+outlPos[6])-0.6123724356957944*(outlPos[5]+outlPos[4])+0.6123724356957944*(outlPos[3]+outlPos[2])-0.6123724356957944*(outlPos[1]+outlPos[0]); 
    outl[3] += 0.6123724356957944*(outlPos[7]+outlPos[6]+outlPos[5]+outlPos[4])-0.6123724356957944*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
    outl[4] += 1.060660171779821*outlPos[7]-1.060660171779821*(outlPos[6]+outlPos[5])+1.060660171779821*(outlPos[4]+outlPos[3])-1.060660171779821*(outlPos[2]+outlPos[1])+1.060660171779821*outlPos[0]; 
    outl[5] += 1.060660171779821*outlPos[7]-1.060660171779821*outlPos[6]+1.060660171779821*outlPos[5]-1.060660171779821*(outlPos[4]+outlPos[3])+1.060660171779821*outlPos[2]-1.060660171779821*outlPos[1]+1.060660171779821*outlPos[0]; 
    outl[6] += 1.060660171779821*(outlPos[7]+outlPos[6])-1.060660171779821*(outlPos[5]+outlPos[4]+outlPos[3]+outlPos[2])+1.060660171779821*(outlPos[1]+outlPos[0]); 
    outl[7] += 1.837117307087383*outlPos[7]-1.837117307087383*(outlPos[6]+outlPos[5])+1.837117307087383*outlPos[4]-1.837117307087383*outlPos[3]+1.837117307087383*(outlPos[2]+outlPos[1])-1.837117307087383*outlPos[0]; 

  }
  return 0.0; 
} 
