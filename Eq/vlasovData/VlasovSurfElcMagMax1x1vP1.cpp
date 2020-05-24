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

  double Ghat[3]; 
  double alpha[3]; 

  double favg[3]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 

  double fjump[3]; 
  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 

  alpha[0] = 1.414213562373095*E0[0]; 
  alpha[1] = 1.414213562373095*E0[1]; 
  const double amid = 0.5*alpha[0]; 

  Ghat[0] = (-0.8660254037844386*fjump[2])+alpha[0]*(0.4330127018922193*favg[2]+0.25*favg[0])+0.25*alpha[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.4330127018922193*favg[2]+0.25*favg[0])-0.5*fjump[1]+0.25*alpha[0]*favg[1]; 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10l; 

  return std::abs(amid); 
} 
