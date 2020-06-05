#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag1x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double *B2 = &EM[15]; 

  double Ghat[6]; 
  double favg[6]; 
  double alpha[6]; 

  favg[0] = 1.58113883008419*fr[8]+1.58113883008419*fl[8]-1.224744871391589*fr[2]+1.224744871391589*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[4])+1.224744871391589*fl[4]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[6])+1.224744871391589*fl[6]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[3] = 0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[4] = 0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[5] = 0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 0.408248290463863*B2[0]*dv2; 
  alpha[3] = 0.408248290463863*B2[1]*dv2; 
  alpha[4] = 1.414213562373095*(B2[2]*wv2+E0[2]); 

  const double amid = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 

  Ghat[0] = 0.25*(alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[8]-1.0*(2.23606797749979*fl[8]+1.732050807568877*(fr[2]+fl[2]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[4]+fl[4])-1.0*fr[1]+fl[1])*amax+0.05*(4.47213595499958*(alpha[1]*favg[4]+favg[1]*alpha[4])+5.0*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[6]+fl[6])-1.0*fr[3]+fl[3])*amax+0.05*(4.47213595499958*alpha[2]*favg[5]+5.0*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = 0.05*(4.47213595499958*(alpha[3]*(favg[5]+favg[4])+favg[3]*alpha[4])+5.0*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2]))-0.3535533905932737*(fr[5]-1.0*fl[5])*amax; 
  Ghat[4] = 0.007142857142857143*((22.3606797749979*alpha[4]+35.0*alpha[0])*favg[4]+35.0*favg[0]*alpha[4]+31.30495168499706*(alpha[3]*favg[3]+alpha[1]*favg[1]))-0.3535533905932737*(fr[7]-1.0*fl[7])*amax; 
  Ghat[5] = 0.05*(5.0*alpha[0]*favg[5]+4.47213595499958*(alpha[3]*favg[3]+alpha[2]*favg[2]))-0.3535533905932737*(fr[9]-1.0*fl[9])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[3] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[4] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[5] += 0.7071067811865475*Ghat[3]*dv10r; 
  outr[6] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[7] += 0.7071067811865475*Ghat[4]*dv10r; 
  outr[8] += 1.58113883008419*Ghat[0]*dv10r; 
  outr[9] += 0.7071067811865475*Ghat[5]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[3] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[4] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[5] += -0.7071067811865475*Ghat[3]*dv10l; 
  outl[6] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[7] += -0.7071067811865475*Ghat[4]*dv10l; 
  outl[8] += -1.58113883008419*Ghat[0]*dv10l; 
  outl[9] += -0.7071067811865475*Ghat[5]*dv10l; 

  return std::abs(amid); 
} 
__host__ __device__ double VlasovSurfElcMag1x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[2]; 
  double dv11r = 2/dxvr[2]; 
  const double *E1 = &EM[3]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double *B2 = &EM[15]; 

  double Ghat[6]; 
  double favg[6]; 
  double alpha[6]; 

  favg[0] = 1.58113883008419*fr[9]+1.58113883008419*fl[9]-1.224744871391589*fr[3]+1.224744871391589*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[5])+1.224744871391589*fl[5]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[6])+1.224744871391589*fl[6]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = 0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = 0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[5] = 0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = -0.408248290463863*B2[0]*dv1; 
  alpha[3] = -0.408248290463863*B2[1]*dv1; 
  alpha[4] = 1.414213562373095*E1[2]-1.414213562373095*B2[2]*wv1; 

  const double amid = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 

  Ghat[0] = 0.25*(alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[9]-1.0*(2.23606797749979*fl[9]+1.732050807568877*(fr[3]+fl[3]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[5]+fl[5])-1.0*fr[1]+fl[1])*amax+0.05*(4.47213595499958*(alpha[1]*favg[4]+favg[1]*alpha[4])+5.0*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[6]+fl[6])-1.0*fr[2]+fl[2])*amax+0.05*(4.47213595499958*alpha[2]*favg[5]+5.0*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = 0.05*(4.47213595499958*(alpha[3]*(favg[5]+favg[4])+favg[3]*alpha[4])+5.0*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2]))-0.3535533905932737*(fr[4]-1.0*fl[4])*amax; 
  Ghat[4] = 0.007142857142857143*((22.3606797749979*alpha[4]+35.0*alpha[0])*favg[4]+35.0*favg[0]*alpha[4]+31.30495168499706*(alpha[3]*favg[3]+alpha[1]*favg[1]))-0.3535533905932737*(fr[7]-1.0*fl[7])*amax; 
  Ghat[5] = 0.05*(5.0*alpha[0]*favg[5]+4.47213595499958*(alpha[3]*favg[3]+alpha[2]*favg[2]))-0.3535533905932737*(fr[8]-1.0*fl[8])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv11r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv11r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv11r; 
  outr[3] += -1.224744871391589*Ghat[0]*dv11r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv11r; 
  outr[5] += -1.224744871391589*Ghat[1]*dv11r; 
  outr[6] += -1.224744871391589*Ghat[2]*dv11r; 
  outr[7] += 0.7071067811865475*Ghat[4]*dv11r; 
  outr[8] += 0.7071067811865475*Ghat[5]*dv11r; 
  outr[9] += 1.58113883008419*Ghat[0]*dv11r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv11l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv11l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv11l; 
  outl[3] += -1.224744871391589*Ghat[0]*dv11l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv11l; 
  outl[5] += -1.224744871391589*Ghat[1]*dv11l; 
  outl[6] += -1.224744871391589*Ghat[2]*dv11l; 
  outl[7] += -0.7071067811865475*Ghat[4]*dv11l; 
  outl[8] += -0.7071067811865475*Ghat[5]*dv11l; 
  outl[9] += -1.58113883008419*Ghat[0]*dv11l; 

  return std::abs(amid); 
} 
