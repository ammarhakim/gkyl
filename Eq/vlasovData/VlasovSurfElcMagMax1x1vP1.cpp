#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag1x1vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 

  double Ghat[2]; 
  double favg[2]; 
  double alpha[2]; 

  favg[0] = (-1.224744871391589*fr[2])+1.224744871391589*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = 0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 

  const double amid = 0.7071067811865475*alpha[0]; 

  Ghat[0] = 0.3535533905932737*((1.732050807568877*(fr[2]+fl[2])-1.0*fr[0]+fl[0])*amax+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(alpha[0]*favg[1]+favg[0]*alpha[1])-0.3535533905932737*(fr[1]-1.0*fl[1])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += -1.224744871391589*Ghat[0]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -1.224744871391589*Ghat[0]*dv10l; 

  return std::abs(amid); 
} 
