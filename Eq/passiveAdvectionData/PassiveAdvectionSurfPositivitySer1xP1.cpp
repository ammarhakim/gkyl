#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionSurfPositivity1xSer_X1_P1(const double *w, const double *dxv, const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  const double *v1 = &fr[2]; 
  double incr[2]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.3535533905932737*(1.732050807568877*v1[1]-1.0*v1[0]); 

  double alpha[1]; 
  alpha[0] = -0.5*(2.449489742783178*v1[1]-1.414213562373095*v1[0]); 
  if (alpha0>0) { 
  double rVal = (2.449489742783178*fl[1])/(2.0*EPSILON+1.414213562373095*fl[0]);  // rVal=f1/f0 
  incr[0] = 0.5*alpha[0]*fl[0]*dfac1*limTheta(rVal,1.0); 
  incr[1] = -0.8660254037844386*alpha[0]*fl[0]*dfac1*limTheta(rVal,1.0); 
  } else { 
  double rVal = (2.449489742783178*fr[1])/(2.0*EPSILON+1.414213562373095*fr[0]);  // rVal=f1/f0 
  incr[0] = 0.5*alpha[0]*fr[0]*dfac1*limTheta(rVal,-1.0); 
  incr[1] = -0.8660254037844386*alpha[0]*fr[0]*dfac1*limTheta(rVal,-1.0); 
  }
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[2], outrPos[2]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 1.0 : positivityWeightByDirL[1]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 1.0 : positivityWeightByDirR[1]/positivityWeightByDirR[0]; 
  outlPos[0] = -0.1666666666666667*(2.449489742783178*incr[1]+4.242640687119286*incr[0]); 
  outlPos[1] = 0.1666666666666667*(2.449489742783178*incr[1]-4.242640687119286*incr[0]); 
  outrPos[0] = -0.1666666666666667*(2.449489742783178*incr[1]-4.242640687119286*incr[0]); 
  outrPos[1] = 0.1666666666666667*(2.449489742783178*incr[1]+4.242640687119286*incr[0]); 
  if(outlPos[1] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.1666666666666667*(2.449489742783178*fl[1]+4.242640687119286*fl[0]))/dtApprox/outlPos[1]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.1666666666666667*(2.449489742783178*fr[1]-4.242640687119286*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[1] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[1] *= limFac; 
  outrPos[0] *= limFac; 
  outr[0] += 0.7071067811865475*(outrPos[1]+outrPos[0]); 
  outr[1] += 0.7071067811865475*(1.732050807568877*outrPos[1]-1.732050807568877*outrPos[0]); 

  outl[0] += 0.7071067811865475*(outlPos[1]+outlPos[0]); 
  outl[1] += 0.7071067811865475*(1.732050807568877*outlPos[1]-1.732050807568877*outlPos[0]); 
  return std::abs(alpha0); 
} 
