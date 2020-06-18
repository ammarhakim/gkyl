#include <VlasovModDecl.h> 
__host__ __device__ double VlasovPhiSurf1x1vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *phi, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // phi:       electrostatic potential.
  // EM:        external EM field vectors.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2./dxvl[1]; 
  double dv10r = 2./dxvr[1]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 

  const double *E0 = &EM[0]; 

  double Ghat[3]; 
  double favg[3]; 
  double alpha[3]; 

  favg[0] = 1.58113883008419*fr[5]+1.58113883008419*fl[5]-1.224744871391589*fr[2]+1.224744871391589*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = 1.58113883008419*fr[7]+1.58113883008419*fl[7]-1.224744871391589*fr[3]+1.224744871391589*fl[3]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[6])+1.224744871391589*fl[6]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 

  alpha[0] = E0[0]-1.732050807568877*phi[1]; 
  alpha[1] = E0[1]-3.872983346207417*phi[2]; 
  alpha[2] = E0[2]; 

  const double amid = 0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2]; 

  Ghat[0] = 0.3535533905932737*(alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[5]-1.0*(2.23606797749979*fl[5]+1.732050807568877*(fr[2]+fl[2]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.07071067811865475*(4.47213595499958*(alpha[1]*favg[2]+favg[1]*alpha[2])+5.0*(alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[7]-1.0*(3.872983346207417*fl[7]+3.0*(fr[3]+fl[3])))+3.0*(fr[1]-1.0*fl[1]))*amax; 
  Ghat[2] = 0.07071067811865475*(8.660254037844387*(fr[6]+fl[6])+5.0*(fl[4]-1.0*fr[4]))*amax+0.01010152544552211*((22.3606797749979*alpha[2]+35.0*alpha[0])*favg[2]+35.0*favg[0]*alpha[2]+31.30495168499706*alpha[1]*favg[1]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[3] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[4] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[5] += 1.58113883008419*Ghat[0]*dv10r; 
  outr[6] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[7] += 1.58113883008419*Ghat[1]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[3] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[4] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[5] += -1.58113883008419*Ghat[0]*dv10l; 
  outl[6] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[7] += -1.58113883008419*Ghat[1]*dv10l; 

  return std::abs(amid); 
} 
