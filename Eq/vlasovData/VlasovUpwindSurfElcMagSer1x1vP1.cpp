#include <VlasovModDecl.h> 
__host__ __device__ double VlasovUpwindSurfElcMag1x1vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double *E0 = &EM[0]; 

  double alpha[2]; 
  double incr[4]; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 

  const double amid = 0.7071067811865475*alpha[0]; 

  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha[0]-0.7071067811865475*alpha[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fr[3]+fl[3]))+0.8660254037844386*(fr[2]+fl[2])+0.5*fr[1]-0.5*(fl[1]+fr[0])+0.5*fl[0])*sgn(alphaOrdR)+0.8660254037844386*fr[3]-0.8660254037844386*(fl[3]+fr[2])+0.8660254037844386*fl[2]-0.5*(fr[1]+fl[1])+0.5*(fr[0]+fl[0])); 
  alphaOrdR = 0.7071067811865475*(alpha[1]+alpha[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fr[3]+fl[3]+fr[2]+fl[2])-0.5*(fr[1]+fr[0])+0.5*(fl[1]+fl[0]))*sgn(alphaOrdR)-0.8660254037844386*(fr[3]+fr[2])+0.8660254037844386*(fl[3]+fl[2])+0.5*(fr[1]+fl[1]+fr[0]+fl[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[1] = 0.5*(alpha[0]*fUp[1]+fUp[0]*alpha[1]); 
  incr[2] = -0.8660254037844386*(alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alpha[0]*fUp[1]+fUp[0]*alpha[1]); 

  outr[0] += incr[0]*dv10r; 
  outr[1] += incr[1]*dv10r; 
  outr[2] += incr[2]*dv10r; 
  outr[3] += incr[3]*dv10r; 

  outl[0] += -1.0*incr[0]*dv10l; 
  outl[1] += -1.0*incr[1]*dv10l; 
  outl[2] += incr[2]*dv10l; 
  outl[3] += incr[3]*dv10l; 


  return std::abs(amid); 
} 
