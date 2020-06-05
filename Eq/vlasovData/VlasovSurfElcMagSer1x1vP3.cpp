#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag1x1vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 

  double Ghat[4]; 
  double favg[4]; 
  double alpha[4]; 

  favg[0] = (-1.870828693386971*fr[9])+1.870828693386971*fl[9]+1.58113883008419*fr[5]+1.58113883008419*fl[5]-1.224744871391589*fr[2]+1.224744871391589*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.870828693386971*fr[11])+1.870828693386971*fl[11]+1.58113883008419*fr[7]+1.58113883008419*fl[7]-1.224744871391589*fr[3]+1.224744871391589*fl[3]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[6])+1.224744871391589*fl[6]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[3] = (-1.224744871391589*fr[10])+1.224744871391589*fl[10]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 
  alpha[2] = E0[2]; 
  alpha[3] = E0[3]; 

  const double amid = 0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2]; 

  Ghat[0] = 0.3535533905932737*((2.645751311064591*(fr[9]+fl[9])+2.23606797749979*(fl[5]-1.0*fr[5])+1.732050807568877*(fr[2]+fl[2])-1.0*fr[0]+fl[0])*amax+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.1178511301977579*(1.732050807568877*(4.58257569495584*(fr[11]+fl[11])+3.872983346207417*(fl[7]-1.0*fr[7])+3.0*(fr[3]+fl[3]))+3.0*(fl[1]-1.0*fr[1]))*amax+0.01010152544552211*(30.7408522978788*(alpha[2]*favg[3]+favg[2]*alpha[3])+31.30495168499706*(alpha[1]*favg[2]+favg[1]*alpha[2])+35.0*(alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.07071067811865475*(8.660254037844387*(fr[6]+fl[6])+5.0*(fl[4]-1.0*fr[4]))*amax+0.003367175148507369*((62.60990336999411*alpha[3]+92.22255689363637*alpha[1])*favg[3]+92.22255689363637*favg[1]*alpha[3]+(67.0820393249937*alpha[2]+105.0*alpha[0])*favg[2]+105.0*favg[0]*alpha[2]+93.91485505499116*alpha[1]*favg[1]); 
  Ghat[3] = 0.05050762722761053*(12.12435565298214*(fr[10]+fl[10])+7.0*(fl[8]-1.0*fr[8]))*amax+0.003367175148507369*((62.60990336999411*alpha[2]+105.0*alpha[0])*favg[3]+(62.60990336999411*favg[2]+105.0*favg[0])*alpha[3]+92.22255689363637*(alpha[1]*favg[2]+favg[1]*alpha[2])); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[3] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[4] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[5] += 1.58113883008419*Ghat[0]*dv10r; 
  outr[6] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[7] += 1.58113883008419*Ghat[1]*dv10r; 
  outr[8] += 0.7071067811865475*Ghat[3]*dv10r; 
  outr[9] += -1.870828693386971*Ghat[0]*dv10r; 
  outr[10] += -1.224744871391589*Ghat[3]*dv10r; 
  outr[11] += -1.870828693386971*Ghat[1]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[3] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[4] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[5] += -1.58113883008419*Ghat[0]*dv10l; 
  outl[6] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[7] += -1.58113883008419*Ghat[1]*dv10l; 
  outl[8] += -0.7071067811865475*Ghat[3]*dv10l; 
  outl[9] += -1.870828693386971*Ghat[0]*dv10l; 
  outl[10] += -1.224744871391589*Ghat[3]*dv10l; 
  outl[11] += -1.870828693386971*Ghat[1]*dv10l; 

  return std::abs(amid); 
} 
