#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag2x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[2]; 
  double dv10r = 2/dxvr[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double *B2 = &EM[15]; 

  double Ghat[5]; 
  double alpha[5]; 

  double favg[5]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 

  double fjump[5]; 
  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0]); 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1]); 
  alpha[2] = 2.0*(B2[2]*wv2+E0[2]); 
  alpha[4] = 0.5773502691896258*B2[0]*dv2; 
  const double amid = 0.25*alpha[0]; 

  Ghat[0] = 0.125*alpha[4]*favg[4]-0.8660254037844386*fjump[3]+alpha[0]*(0.2165063509461096*favg[3]+0.125*favg[0])+0.125*(alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.2165063509461096*favg[3]+0.125*favg[0])-0.5*fjump[1]+0.125*alpha[0]*favg[1]; 
  Ghat[2] = alpha[2]*(0.2165063509461096*favg[3]+0.125*favg[0])-0.5*fjump[2]+0.125*alpha[0]*favg[2]; 
  Ghat[4] = (-0.5*fjump[4])+0.125*alpha[0]*favg[4]+(0.2165063509461096*favg[3]+0.125*favg[0])*alpha[4]; 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += 0.5*Ghat[2]*dv10r; 
  outr[3] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[4] += 0.5*Ghat[4]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.5*Ghat[2]*dv10l; 
  outl[3] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[4] += -0.5*Ghat[4]*dv10l; 

  return std::abs(amid); 
} 
__host__ __device__ double VlasovSurfElcMag2x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[3]; 
  double dv11r = 2/dxvr[3]; 
  const double *E1 = &EM[3]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double *B2 = &EM[15]; 

  double Ghat[5]; 
  double alpha[5]; 

  double favg[5]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 

  double fjump[5]; 
  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(-1*fr[4]-fl[4]); 

  alpha[0] = 2.0*E1[0]-2.0*B2[0]*wv1; 
  alpha[1] = 2.0*E1[1]-2.0*B2[1]*wv1; 
  alpha[2] = 2.0*E1[2]-2.0*B2[2]*wv1; 
  alpha[3] = -0.5773502691896258*B2[0]*dv1; 
  const double amid = 0.25*alpha[0]; 

  Ghat[0] = (-0.8660254037844386*fjump[4])+alpha[0]*(0.2165063509461096*favg[4]+0.125*favg[0])+0.125*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.2165063509461096*favg[4]+0.125*favg[0])-0.5*fjump[1]+0.125*alpha[0]*favg[1]; 
  Ghat[2] = alpha[2]*(0.2165063509461096*favg[4]+0.125*favg[0])-0.5*fjump[2]+0.125*alpha[0]*favg[2]; 
  Ghat[3] = alpha[3]*(0.2165063509461096*favg[4]+0.125*favg[0])-0.5*fjump[3]+0.125*alpha[0]*favg[3]; 

  outr[0] += 0.5*Ghat[0]*dv11r; 
  outr[1] += 0.5*Ghat[1]*dv11r; 
  outr[2] += 0.5*Ghat[2]*dv11r; 
  outr[3] += 0.5*Ghat[3]*dv11r; 
  outr[4] += -0.8660254037844386*Ghat[0]*dv11r; 

  outl[0] += -0.5*Ghat[0]*dv11l; 
  outl[1] += -0.5*Ghat[1]*dv11l; 
  outl[2] += -0.5*Ghat[2]*dv11l; 
  outl[3] += -0.5*Ghat[3]*dv11l; 
  outl[4] += -0.8660254037844386*Ghat[0]*dv11l; 

  return std::abs(amid); 
} 
