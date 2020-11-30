#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionSurfPositivity1xSer_X1_P2(const double *w, const double *dxv, const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  const double *v1 = &fr[3]; 
  double incr[3]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*(2.23606797749979*v1[2]-1.732050807568877*v1[1]+v1[0]); 

  double alpha[1]; 
  alpha[0] = 0.7071067811865475*(2.23606797749979*v1[2]-1.732050807568877*v1[1]+v1[0]); 
  if (alpha0>0) { 
  incr[0] = 0.5*alpha[0]*(2.23606797749979*fl[2]+1.732050807568877*fl[1]+fl[0])*dfac1; 
  incr[1] = -0.5*alpha[0]*(3.872983346207417*fl[2]+3.0*fl[1]+1.732050807568877*fl[0])*dfac1; 
  incr[2] = 0.5*alpha[0]*(5.0*fl[2]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*dfac1; 
  } else { 
  incr[0] = 0.5*alpha[0]*(2.23606797749979*fr[2]-1.732050807568877*fr[1]+fr[0])*dfac1; 
  incr[1] = -0.5*alpha[0]*(3.872983346207417*fr[2]-3.0*fr[1]+1.732050807568877*fr[0])*dfac1; 
  incr[2] = 0.5*alpha[0]*(5.0*fr[2]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*dfac1; 
  }
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[3], outrPos[3]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 1.0 : positivityWeightByDirL[1]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 1.0 : positivityWeightByDirR[1]/positivityWeightByDirR[0]; 
  outlPos[0] = -0.03333333333333333*(6.324555320336761*incr[2]+12.24744871391589*incr[1]+7.071067811865476*incr[0]); 
  outlPos[1] = 0.06666666666666667*(6.324555320336761*incr[2]-14.14213562373095*incr[0]); 
  outlPos[2] = -0.03333333333333333*(6.324555320336761*incr[2]-12.24744871391589*incr[1]+7.071067811865476*incr[0]); 
  outrPos[0] = 0.03333333333333333*(6.324555320336761*incr[2]-12.24744871391589*incr[1]+7.071067811865476*incr[0]); 
  outrPos[1] = -0.06666666666666667*(6.324555320336761*incr[2]-14.14213562373095*incr[0]); 
  outrPos[2] = 0.03333333333333333*(6.324555320336761*incr[2]+12.24744871391589*incr[1]+7.071067811865476*incr[0]); 
  if(outlPos[2] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.03333333333333333*(6.324555320336761*fl[2]+12.24744871391589*fl[1]+7.071067811865476*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.03333333333333333*(6.324555320336761*fr[2]-12.24744871391589*fr[1]+7.071067811865476*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[2] *= limFac; 
  outrPos[0] *= limFac; 
  outr[0] += 0.7071067811865475*(outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.7071067811865475*(1.732050807568877*outrPos[2]-1.732050807568877*outrPos[0]); 
  outr[2] += 0.3535533905932737*(4.47213595499958*outrPos[2]-2.23606797749979*outrPos[1]+4.47213595499958*outrPos[0]); 

  outl[0] += 0.7071067811865475*(outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.7071067811865475*(1.732050807568877*outlPos[2]-1.732050807568877*outlPos[0]); 
  outl[2] += 0.3535533905932737*(4.47213595499958*outlPos[2]-2.23606797749979*outlPos[1]+4.47213595499958*outlPos[0]); 
  return std::abs(alpha0); 
} 
