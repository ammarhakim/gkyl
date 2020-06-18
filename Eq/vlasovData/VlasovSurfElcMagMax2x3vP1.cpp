#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag2x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[2]; 
  double dv10r = 2/dxvr[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 
  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double Ghat[5]; 
  double favg[5]; 
  double alpha[5]; 

  favg[0] = (-1.224744871391589*fr[3])+1.224744871391589*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = 0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = 0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = 0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = 0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 2.0*(B2[2]*wv2+E0[2])-2.0*B1[2]*wv3; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[4] = -0.5773502691896258*B1[0]*dv3; 

  const double amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[3]+fl[3])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.125*(alpha[0]*favg[1]+favg[0]*alpha[1])-0.3535533905932737*(fr[1]-1.0*fl[1])*amax; 
  Ghat[2] = 0.125*(alpha[0]*favg[2]+favg[0]*alpha[2])-0.3535533905932737*(fr[2]-1.0*fl[2])*amax; 
  Ghat[3] = 0.125*(alpha[0]*favg[3]+favg[0]*alpha[3])-0.3535533905932737*(fr[4]-1.0*fl[4])*amax; 
  Ghat[4] = 0.125*(alpha[0]*favg[4]+favg[0]*alpha[4])-0.3535533905932737*(fr[5]-1.0*fl[5])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[3] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv10r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[3] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv10l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv10l; 

  return std::abs(amid); 
} 
__host__ __device__ double VlasovSurfElcMag2x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[3]; 
  double dv11r = 2/dxvr[3]; 
  const double *E1 = &EM[3]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 
  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double Ghat[5]; 
  double favg[5]; 
  double alpha[5]; 

  favg[0] = (-1.224744871391589*fr[4])+1.224744871391589*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = 0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = 0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = 0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = 0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 

  alpha[0] = 2.0*B0[0]*wv3-2.0*B2[0]*wv1+2.0*E1[0]; 
  alpha[1] = 2.0*B0[1]*wv3-2.0*B2[1]*wv1+2.0*E1[1]; 
  alpha[2] = 2.0*B0[2]*wv3-2.0*B2[2]*wv1+2.0*E1[2]; 
  alpha[3] = -0.5773502691896258*B2[0]*dv1; 
  alpha[4] = 0.5773502691896258*B0[0]*dv3; 

  const double amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[4]+fl[4])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.125*(alpha[0]*favg[1]+favg[0]*alpha[1])-0.3535533905932737*(fr[1]-1.0*fl[1])*amax; 
  Ghat[2] = 0.125*(alpha[0]*favg[2]+favg[0]*alpha[2])-0.3535533905932737*(fr[2]-1.0*fl[2])*amax; 
  Ghat[3] = 0.125*(alpha[0]*favg[3]+favg[0]*alpha[3])-0.3535533905932737*(fr[3]-1.0*fl[3])*amax; 
  Ghat[4] = 0.125*(alpha[0]*favg[4]+favg[0]*alpha[4])-0.3535533905932737*(fr[5]-1.0*fl[5])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv11r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv11r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv11r; 
  outr[3] += 0.7071067811865475*Ghat[3]*dv11r; 
  outr[4] += -1.224744871391589*Ghat[0]*dv11r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv11r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv11l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv11l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv11l; 
  outl[3] += -0.7071067811865475*Ghat[3]*dv11l; 
  outl[4] += -1.224744871391589*Ghat[0]*dv11l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv11l; 

  return std::abs(amid); 
} 
__host__ __device__ double VlasovSurfElcMag2x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12l = 2/dxvl[4]; 
  double dv12r = 2/dxvr[4]; 
  const double *E2 = &EM[6]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 
  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double Ghat[5]; 
  double favg[5]; 
  double alpha[5]; 

  favg[0] = (-1.224744871391589*fr[5])+1.224744871391589*fl[5]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = 0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = 0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = 0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = 0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 

  alpha[0] = 2.0*(B1[0]*wv1+E2[0])-2.0*B0[0]*wv2; 
  alpha[1] = 2.0*(B1[1]*wv1+E2[1])-2.0*B0[1]*wv2; 
  alpha[2] = 2.0*(B1[2]*wv1+E2[2])-2.0*B0[2]*wv2; 
  alpha[3] = 0.5773502691896258*B1[0]*dv1; 
  alpha[4] = -0.5773502691896258*B0[0]*dv2; 

  const double amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[5]+fl[5])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.125*(alpha[0]*favg[1]+favg[0]*alpha[1])-0.3535533905932737*(fr[1]-1.0*fl[1])*amax; 
  Ghat[2] = 0.125*(alpha[0]*favg[2]+favg[0]*alpha[2])-0.3535533905932737*(fr[2]-1.0*fl[2])*amax; 
  Ghat[3] = 0.125*(alpha[0]*favg[3]+favg[0]*alpha[3])-0.3535533905932737*(fr[3]-1.0*fl[3])*amax; 
  Ghat[4] = 0.125*(alpha[0]*favg[4]+favg[0]*alpha[4])-0.3535533905932737*(fr[4]-1.0*fl[4])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv12r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv12r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv12r; 
  outr[3] += 0.7071067811865475*Ghat[3]*dv12r; 
  outr[4] += 0.7071067811865475*Ghat[4]*dv12r; 
  outr[5] += -1.224744871391589*Ghat[0]*dv12r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv12l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv12l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv12l; 
  outl[3] += -0.7071067811865475*Ghat[3]*dv12l; 
  outl[4] += -0.7071067811865475*Ghat[4]*dv12l; 
  outl[5] += -1.224744871391589*Ghat[0]*dv12l; 

  return std::abs(amid); 
} 
