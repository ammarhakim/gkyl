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

  double Ghat[10]; 
  double alpha[10]; 

  double favg[10]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 

  double fjump[10]; 
  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(-1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(1*fr[9]-fl[9]); 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0]); 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1]); 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[5] = 0.5773502691896258*B2[1]*dv2; 
  alpha[7] = 2.0*(B2[2]*wv2+E0[2]); 
  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 

  Ghat[0] = (-1.118033988749895*fjump[8])+alpha[0]*(0.3952847075210473*favg[8]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+0.1767766952966368*alpha[7]*favg[7]+alpha[3]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+0.1767766952966368*alpha[5]*favg[5]+alpha[1]*(0.3061862178478971*favg[4]+0.1767766952966368*favg[1])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.3952847075210473*favg[8]+0.1581138830084189*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+(0.273861278752583*favg[4]+0.1581138830084189*favg[1])*alpha[7]+alpha[5]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+0.1767766952966368*alpha[3]*favg[5]-0.8660254037844386*fjump[4]+alpha[0]*(0.3061862178478971*favg[4]+0.1767766952966368*favg[1])-0.5*fjump[1]; 
  Ghat[3] = alpha[3]*(0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[6]+alpha[0]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+0.1767766952966368*alpha[1]*favg[5]+(0.3061862178478971*favg[4]+0.1767766952966368*favg[1])*alpha[5]-0.5*fjump[3]; 
  Ghat[5] = alpha[5]*(0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.1581138830084189*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+0.1581138830084189*favg[5]*alpha[7]+alpha[1]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])-0.5*fjump[5]+0.1767766952966368*alpha[0]*favg[5]+alpha[3]*(0.3061862178478971*favg[4]+0.1767766952966368*favg[1]); 
  Ghat[7] = alpha[7]*(0.3952847075210473*favg[8]+0.1129384878631564*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-0.5*fjump[7]+0.1767766952966368*alpha[0]*favg[7]+0.1581138830084189*alpha[5]*favg[5]+alpha[1]*(0.273861278752583*favg[4]+0.1581138830084189*favg[1]); 
  Ghat[9] = (-0.5*fjump[9])+0.1767766952966368*alpha[0]*favg[9]+alpha[3]*(0.273861278752583*favg[6]+0.1581138830084189*favg[3])+0.1581138830084189*alpha[5]*favg[5]; 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[3] += 0.5*Ghat[3]*dv10r; 
  outr[4] += -0.8660254037844386*Ghat[1]*dv10r; 
  outr[5] += 0.5*Ghat[5]*dv10r; 
  outr[6] += -0.8660254037844386*Ghat[3]*dv10r; 
  outr[7] += 0.5*Ghat[7]*dv10r; 
  outr[8] += 1.118033988749895*Ghat[0]*dv10r; 
  outr[9] += 0.5*Ghat[9]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[3] += -0.5*Ghat[3]*dv10l; 
  outl[4] += -0.8660254037844386*Ghat[1]*dv10l; 
  outl[5] += -0.5*Ghat[5]*dv10l; 
  outl[6] += -0.8660254037844386*Ghat[3]*dv10l; 
  outl[7] += -0.5*Ghat[7]*dv10l; 
  outl[8] += -1.118033988749895*Ghat[0]*dv10l; 
  outl[9] += -0.5*Ghat[9]*dv10l; 

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

  double Ghat[10]; 
  double alpha[10]; 

  double favg[10]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 

  double fjump[10]; 
  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(-1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(1*fr[9]-fl[9]); 

  alpha[0] = 2.0*E1[0]-2.0*B2[0]*wv1; 
  alpha[1] = 2.0*E1[1]-2.0*B2[1]*wv1; 
  alpha[2] = -0.5773502691896258*B2[0]*dv1; 
  alpha[4] = -0.5773502691896258*B2[1]*dv1; 
  alpha[7] = 2.0*E1[2]-2.0*B2[2]*wv1; 
  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 

  Ghat[0] = (-1.118033988749895*fjump[9])+alpha[0]*(0.3952847075210473*favg[9]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+0.1767766952966368*alpha[7]*favg[7]+alpha[2]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])+alpha[1]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[1])+0.1767766952966368*alpha[4]*favg[4]-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.3952847075210473*favg[9]+0.1581138830084189*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+(0.273861278752583*favg[5]+0.1581138830084189*favg[1])*alpha[7]+alpha[4]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])-0.8660254037844386*fjump[5]+alpha[0]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[1])+0.1767766952966368*alpha[2]*favg[4]-0.5*fjump[1]; 
  Ghat[2] = alpha[2]*(0.3952847075210473*favg[9]+0.1581138830084189*favg[8]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[6]+alpha[0]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])+alpha[4]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[1])+0.1767766952966368*alpha[1]*favg[4]-0.5*fjump[2]; 
  Ghat[4] = alpha[4]*(0.3952847075210473*favg[9]+0.1581138830084189*(favg[8]+favg[7])+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+0.1581138830084189*favg[4]*alpha[7]+alpha[1]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])+alpha[2]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[1])-0.5*fjump[4]+0.1767766952966368*alpha[0]*favg[4]; 
  Ghat[7] = alpha[7]*(0.3952847075210473*favg[9]+0.1129384878631564*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])-0.5*fjump[7]+0.1767766952966368*alpha[0]*favg[7]+alpha[1]*(0.273861278752583*favg[5]+0.1581138830084189*favg[1])+0.1581138830084189*alpha[4]*favg[4]; 
  Ghat[8] = (-0.5*fjump[8])+0.1767766952966368*alpha[0]*favg[8]+alpha[2]*(0.273861278752583*favg[6]+0.1581138830084189*favg[2])+0.1581138830084189*alpha[4]*favg[4]; 

  outr[0] += 0.5*Ghat[0]*dv11r; 
  outr[1] += 0.5*Ghat[1]*dv11r; 
  outr[2] += 0.5*Ghat[2]*dv11r; 
  outr[3] += -0.8660254037844386*Ghat[0]*dv11r; 
  outr[4] += 0.5*Ghat[4]*dv11r; 
  outr[5] += -0.8660254037844386*Ghat[1]*dv11r; 
  outr[6] += -0.8660254037844386*Ghat[2]*dv11r; 
  outr[7] += 0.5*Ghat[7]*dv11r; 
  outr[8] += 0.5*Ghat[8]*dv11r; 
  outr[9] += 1.118033988749895*Ghat[0]*dv11r; 

  outl[0] += -0.5*Ghat[0]*dv11l; 
  outl[1] += -0.5*Ghat[1]*dv11l; 
  outl[2] += -0.5*Ghat[2]*dv11l; 
  outl[3] += -0.8660254037844386*Ghat[0]*dv11l; 
  outl[4] += -0.5*Ghat[4]*dv11l; 
  outl[5] += -0.8660254037844386*Ghat[1]*dv11l; 
  outl[6] += -0.8660254037844386*Ghat[2]*dv11l; 
  outl[7] += -0.5*Ghat[7]*dv11l; 
  outl[8] += -0.5*Ghat[8]*dv11l; 
  outl[9] += -1.118033988749895*Ghat[0]*dv11l; 

  return std::abs(amid); 
} 
