#include <IncompEulerModDecl.h> 
double IncompEulerSurfPositivity2xSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // H/f: Input Hamiltonian/distribution function.
  // out: Incremented output.

  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // Surface-averaged phase velocity in this direction.
  double alpha0 = -0.25*(3.0*phi[3]-1.732050807568877*phi[2])*dfac_y; 

  double alpha[2]; 
  alpha[0] = -0.5*(4.242640687119286*phi[3]-2.449489742783178*phi[2])*dfac_y; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.3535533905932737*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*(1.732050807568877*fl[3]+fl[2])*dfac_x; 
  incr[3] = -0.3535533905932737*alpha[0]*(3.0*fl[3]+1.732050807568877*fl[2])*dfac_x; 
  } else { 
  incr[0] = -0.3535533905932737*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.3535533905932737*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.3535533905932737*alpha[0]*(1.732050807568877*fr[3]-1.0*fr[2])*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*(3.0*fr[3]-1.732050807568877*fr[2])*dfac_x; 
  }
#elif upwindType == QUAD 
  double fupwind[2];
  double fupwindQuad[2];
  double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[0] = 0.5*(((-0.8660254037844386*(fr[3]+fl[3]))+0.5*fr[2]-0.5*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)+0.8660254037844386*fr[3]-0.8660254037844386*fl[3]-0.5*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[1] = 0.5*((0.8660254037844386*(fr[3]+fl[3])-0.5*fr[2]+0.5*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.8660254037844386*fr[3]+0.8660254037844386*fl[3]+0.5*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  fupwind[0] = 0.7071067811865476*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.7071067811865476*fupwindQuad[1]-0.7071067811865476*fupwindQuad[0]; 
  incr[0] = 0.5*alpha[0]*fupwind[0]*dfac_x; 
  incr[1] = -0.8660254037844386*alpha[0]*fupwind[0]*dfac_x; 
  incr[2] = 0.5*alpha[0]*fupwind[1]*dfac_x; 
  incr[3] = -0.8660254037844386*alpha[0]*fupwind[1]*dfac_x; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[4], outrPos[4]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.5 : positivityWeightByDirL[1]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.5 : positivityWeightByDirR[1]/positivityWeightByDirR[0]; 
  outlPos[0] = 0.1666666666666667*(incr[3]+1.732050807568877*incr[2]-1.732050807568877*incr[1]-3.0*incr[0]); 
  outlPos[1] = -0.1666666666666667*(incr[3]-1.732050807568877*(incr[2]+incr[1])+3.0*incr[0]); 
  outlPos[2] = -0.1666666666666667*(incr[3]+1.732050807568877*(incr[2]+incr[1])+3.0*incr[0]); 
  outlPos[3] = 0.1666666666666667*(incr[3]-1.732050807568877*incr[2]+1.732050807568877*incr[1]-3.0*incr[0]); 
  outrPos[0] = 0.1666666666666667*(incr[3]-1.732050807568877*(incr[2]+incr[1])+3.0*incr[0]); 
  outrPos[1] = -0.1666666666666667*(incr[3]+1.732050807568877*incr[2]-1.732050807568877*incr[1]-3.0*incr[0]); 
  outrPos[2] = -0.1666666666666667*(incr[3]-1.732050807568877*incr[2]+1.732050807568877*incr[1]-3.0*incr[0]); 
  outrPos[3] = 0.1666666666666667*(incr[3]+1.732050807568877*(incr[2]+incr[1])+3.0*incr[0]); 
  if(outlPos[1] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.1666666666666667*(fl[3]+1.732050807568877*fl[2]-1.732050807568877*fl[1]-3.0*fl[0]))/dtApprox/outlPos[1]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.1666666666666667*(fr[3]-1.732050807568877*(fr[2]+fr[1])+3.0*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[1] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[1] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[3] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.1666666666666667*(fl[3]+1.732050807568877*(fl[2]+fl[1])+3.0*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[2] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.1666666666666667*(fr[3]-1.732050807568877*fr[2]+1.732050807568877*fr[1]-3.0*fr[0]))/dtApprox/outrPos[2]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[2] *= limFac; 
  outrPos[3] *= limFac; 
  outrPos[2] *= limFac; 
  outr[0] += 0.5*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.5*(1.732050807568877*outrPos[3]-1.732050807568877*outrPos[2]+1.732050807568877*outrPos[1]-1.732050807568877*outrPos[0]); 
  outr[2] += 0.5*(1.732050807568877*(outrPos[3]+outrPos[2])-1.732050807568877*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.5*(3.0*outrPos[3]-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 

  outl[0] += 0.5*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.5*(1.732050807568877*outlPos[3]-1.732050807568877*outlPos[2]+1.732050807568877*outlPos[1]-1.732050807568877*outlPos[0]); 
  outl[2] += 0.5*(1.732050807568877*(outlPos[3]+outlPos[2])-1.732050807568877*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.5*(3.0*outlPos[3]-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  return std::abs(alpha0); 
} 
double IncompEulerSurfPositivity2xSer_Y_P1(const double q_, const double m_, const double *w, const double *dxv, const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // H/f: Input Hamiltonian/distribution function.
  // out: Incremented output.

  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // Surface-averaged phase velocity in this direction.
  double alpha0 = 0.25*(3.0*phi[3]-1.732050807568877*phi[1])*dfac_x; 

  double alpha[2]; 
  alpha[0] = 0.5*(4.242640687119286*phi[3]-2.449489742783178*phi[1])*dfac_x; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.3535533905932737*alpha[0]*(1.732050807568877*fl[2]+fl[0])*dfac_y; 
  incr[1] = 0.3535533905932737*alpha[0]*(1.732050807568877*fl[3]+fl[1])*dfac_y; 
  incr[2] = -0.3535533905932737*alpha[0]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_y; 
  incr[3] = -0.3535533905932737*alpha[0]*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_y; 
  } else { 
  incr[0] = -0.3535533905932737*alpha[0]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_y; 
  incr[1] = -0.3535533905932737*alpha[0]*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_y; 
  incr[2] = 0.3535533905932737*alpha[0]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_y; 
  incr[3] = 0.3535533905932737*alpha[0]*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_y; 
  }
#elif upwindType == QUAD 
  double fupwind[2];
  double fupwindQuad[2];
  double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[0] = 0.5*(((-0.8660254037844386*(fr[3]+fl[3]))+0.8660254037844386*(fr[2]+fl[2])+0.5*fr[1]-0.5*(fl[1]+fr[0])+0.5*fl[0])*sgn(alphaQuad)+0.8660254037844386*fr[3]-0.8660254037844386*(fl[3]+fr[2])+0.8660254037844386*fl[2]-0.5*(fr[1]+fl[1])+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[1] = 0.5*((0.8660254037844386*(fr[3]+fl[3]+fr[2]+fl[2])-0.5*(fr[1]+fr[0])+0.5*(fl[1]+fl[0]))*sgn(alphaQuad)-0.8660254037844386*(fr[3]+fr[2])+0.8660254037844386*(fl[3]+fl[2])+0.5*(fr[1]+fl[1]+fr[0]+fl[0])); 
  fupwind[0] = 0.7071067811865476*(fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.7071067811865476*fupwindQuad[1]-0.7071067811865476*fupwindQuad[0]; 
  incr[0] = 0.5*alpha[0]*fupwind[0]*dfac_y; 
  incr[1] = 0.5*alpha[0]*fupwind[1]*dfac_y; 
  incr[2] = -0.8660254037844386*alpha[0]*fupwind[0]*dfac_y; 
  incr[3] = -0.8660254037844386*alpha[0]*fupwind[1]*dfac_y; 

#endif 
  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[4], outrPos[4]; 
  fluxFracL = positivityWeightByDirL[0] == 0. ? 0.5 : positivityWeightByDirL[2]/positivityWeightByDirL[0]; 
  fluxFracR = positivityWeightByDirR[0] == 0. ? 0.5 : positivityWeightByDirR[2]/positivityWeightByDirR[0]; 
  outlPos[0] = 0.1666666666666667*(incr[3]-1.732050807568877*incr[2]+1.732050807568877*incr[1]-3.0*incr[0]); 
  outlPos[1] = -0.1666666666666667*(incr[3]+1.732050807568877*(incr[2]+incr[1])+3.0*incr[0]); 
  outlPos[2] = -0.1666666666666667*(incr[3]-1.732050807568877*(incr[2]+incr[1])+3.0*incr[0]); 
  outlPos[3] = 0.1666666666666667*(incr[3]+1.732050807568877*incr[2]-1.732050807568877*incr[1]-3.0*incr[0]); 
  outrPos[0] = 0.1666666666666667*(incr[3]-1.732050807568877*(incr[2]+incr[1])+3.0*incr[0]); 
  outrPos[1] = -0.1666666666666667*(incr[3]+1.732050807568877*incr[2]-1.732050807568877*incr[1]-3.0*incr[0]); 
  outrPos[2] = -0.1666666666666667*(incr[3]-1.732050807568877*incr[2]+1.732050807568877*incr[1]-3.0*incr[0]); 
  outrPos[3] = 0.1666666666666667*(incr[3]+1.732050807568877*(incr[2]+incr[1])+3.0*incr[0]); 
  if(outlPos[2] < EPSILON) limFac = std::min(1.0, -fluxFracL*(-0.1666666666666667*(fl[3]-1.732050807568877*fl[2]+1.732050807568877*fl[1]-3.0*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 1.0; 
  if(outrPos[0] < EPSILON) limFac = std::min(limFac, -fluxFracR*(0.1666666666666667*(fr[3]-1.732050807568877*(fr[2]+fr[1])+3.0*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outlPos[0] *= limFac; 
  outrPos[2] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[3] < EPSILON) limFac = std::min(1.0, -fluxFracL*(0.1666666666666667*(fl[3]+1.732050807568877*(fl[2]+fl[1])+3.0*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[1] < EPSILON) limFac = std::min(limFac, -fluxFracR*(-0.1666666666666667*(fr[3]+1.732050807568877*fr[2]-1.732050807568877*fr[1]-3.0*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outlPos[1] *= limFac; 
  outrPos[3] *= limFac; 
  outrPos[1] *= limFac; 
  outr[0] += 0.5*(outrPos[3]+outrPos[2]+outrPos[1]+outrPos[0]); 
  outr[1] += 0.5*(1.732050807568877*outrPos[3]-1.732050807568877*outrPos[2]+1.732050807568877*outrPos[1]-1.732050807568877*outrPos[0]); 
  outr[2] += 0.5*(1.732050807568877*(outrPos[3]+outrPos[2])-1.732050807568877*(outrPos[1]+outrPos[0])); 
  outr[3] += 0.5*(3.0*outrPos[3]-3.0*(outrPos[2]+outrPos[1])+3.0*outrPos[0]); 

  outl[0] += 0.5*(outlPos[3]+outlPos[2]+outlPos[1]+outlPos[0]); 
  outl[1] += 0.5*(1.732050807568877*outlPos[3]-1.732050807568877*outlPos[2]+1.732050807568877*outlPos[1]-1.732050807568877*outlPos[0]); 
  outl[2] += 0.5*(1.732050807568877*(outlPos[3]+outlPos[2])-1.732050807568877*(outlPos[1]+outlPos[0])); 
  outl[3] += 0.5*(3.0*outlPos[3]-3.0*(outlPos[2]+outlPos[1])+3.0*outlPos[0]); 
  return std::abs(alpha0); 
} 
