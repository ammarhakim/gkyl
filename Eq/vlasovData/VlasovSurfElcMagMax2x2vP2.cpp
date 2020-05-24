#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag2x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[2]; 
  double dv10r = 2/dxvr[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double *B2 = &EM[30]; 

  double Ghat[15]; 
  double alpha[15]; 

  double favg[15]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 

  double fjump[15]; 
  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(-1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(1*fr[9]-fl[9]); 
  fjump[10] = amax*(-1*fr[10]-fl[10]); 
  fjump[11] = amax*(1*fr[11]-fl[11]); 
  fjump[12] = amax*(1*fr[12]-fl[12]); 
  fjump[13] = amax*(1*fr[13]-fl[13]); 
  fjump[14] = amax*(1*fr[14]-fl[14]); 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0]); 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1]); 
  alpha[2] = 2.0*(B2[2]*wv2+E0[2]); 
  alpha[4] = 0.5773502691896258*B2[0]*dv2; 
  alpha[5] = 2.0*(B2[3]*wv2+E0[3]); 
  alpha[8] = 0.5773502691896258*B2[1]*dv2; 
  alpha[9] = 0.5773502691896258*B2[2]*dv2; 
  alpha[11] = 2.0*(B2[4]*wv2+E0[4]); 
  alpha[12] = 2.0*(B2[5]*wv2+E0[5]); 
  const double amid = (-0.2795084971874737*alpha[12])-0.2795084971874737*alpha[11]+0.25*alpha[0]; 

  Ghat[0] = (-1.118033988749895*fjump[13])+alpha[0]*(0.2795084971874737*favg[13]+0.2165063509461096*favg[3]+0.125*favg[0])+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11])+alpha[4]*(0.2165063509461096*favg[10]+0.125*favg[4])+0.125*(alpha[9]*favg[9]+alpha[8]*favg[8])+alpha[2]*(0.2165063509461096*favg[7]+0.125*favg[2])+alpha[1]*(0.2165063509461096*favg[6]+0.125*favg[1])+0.125*alpha[5]*favg[5]-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.2795084971874737*favg[13]+0.1118033988749895*favg[11]+0.2165063509461096*favg[3]+0.125*favg[0])+(0.1936491673103708*favg[6]+0.1118033988749895*favg[1])*alpha[11]+alpha[8]*(0.2165063509461096*favg[10]+0.125*favg[4])+0.125*alpha[4]*favg[8]+alpha[5]*(0.2165063509461096*favg[7]+0.125*favg[2])-0.8660254037844386*fjump[6]+alpha[0]*(0.2165063509461096*favg[6]+0.125*favg[1])+0.125*alpha[2]*favg[5]-0.5*fjump[1]; 
  Ghat[2] = alpha[2]*(0.2795084971874737*favg[13]+0.1118033988749895*favg[12]+0.2165063509461096*favg[3]+0.125*favg[0])+(0.1936491673103708*favg[7]+0.1118033988749895*favg[2])*alpha[12]+alpha[9]*(0.2165063509461096*favg[10]+0.125*favg[4])+0.125*alpha[4]*favg[9]-0.8660254037844386*fjump[7]+alpha[0]*(0.2165063509461096*favg[7]+0.125*favg[2])+alpha[5]*(0.2165063509461096*favg[6]+0.125*favg[1])+0.125*alpha[1]*favg[5]-0.5*fjump[2]; 
  Ghat[4] = alpha[4]*(0.1118033988749895*favg[14]+0.2795084971874737*favg[13]+0.2165063509461096*favg[3]+0.125*favg[0])-0.8660254037844386*fjump[10]+alpha[0]*(0.2165063509461096*favg[10]+0.125*favg[4])+0.125*alpha[2]*favg[9]+(0.2165063509461096*favg[7]+0.125*favg[2])*alpha[9]+0.125*alpha[1]*favg[8]+(0.2165063509461096*favg[6]+0.125*favg[1])*alpha[8]-0.5*fjump[4]; 
  Ghat[5] = alpha[5]*(0.2795084971874737*favg[13]+0.1118033988749895*(favg[12]+favg[11])+0.2165063509461096*favg[3]+0.125*favg[0])+0.1118033988749895*favg[5]*(alpha[12]+alpha[11])+0.125*(alpha[8]*favg[9]+favg[8]*alpha[9])+alpha[1]*(0.2165063509461096*favg[7]+0.125*favg[2])+alpha[2]*(0.2165063509461096*favg[6]+0.125*favg[1])-0.5*fjump[5]+0.125*alpha[0]*favg[5]; 
  Ghat[8] = alpha[8]*(0.1118033988749895*favg[14]+0.2795084971874737*favg[13]+0.1118033988749895*favg[11]+0.2165063509461096*favg[3]+0.125*favg[0])+0.1118033988749895*favg[8]*alpha[11]+alpha[1]*(0.2165063509461096*favg[10]+0.125*favg[4])+0.125*(alpha[5]*favg[9]+favg[5]*alpha[9])-0.5*fjump[8]+0.125*alpha[0]*favg[8]+alpha[4]*(0.2165063509461096*favg[6]+0.125*favg[1]); 
  Ghat[9] = alpha[9]*(0.1118033988749895*favg[14]+0.2795084971874737*favg[13]+0.1118033988749895*favg[12]+0.2165063509461096*favg[3]+0.125*favg[0])+0.1118033988749895*favg[9]*alpha[12]+alpha[2]*(0.2165063509461096*favg[10]+0.125*favg[4])-0.5*fjump[9]+0.125*(alpha[0]*favg[9]+alpha[5]*favg[8]+favg[5]*alpha[8])+alpha[4]*(0.2165063509461096*favg[7]+0.125*favg[2]); 
  Ghat[11] = alpha[11]*(0.2795084971874737*favg[13]+0.07985957062499249*favg[11]+0.2165063509461096*favg[3]+0.125*favg[0])-0.5*fjump[11]+0.125*alpha[0]*favg[11]+0.1118033988749895*alpha[8]*favg[8]+alpha[1]*(0.1936491673103708*favg[6]+0.1118033988749895*favg[1])+0.1118033988749895*alpha[5]*favg[5]; 
  Ghat[12] = alpha[12]*(0.2795084971874737*favg[13]+0.07985957062499249*favg[12]+0.2165063509461096*favg[3]+0.125*favg[0])-0.5*fjump[12]+0.125*alpha[0]*favg[12]+0.1118033988749895*alpha[9]*favg[9]+alpha[2]*(0.1936491673103708*favg[7]+0.1118033988749895*favg[2])+0.1118033988749895*alpha[5]*favg[5]; 
  Ghat[14] = (-0.5*fjump[14])+0.125*alpha[0]*favg[14]+alpha[4]*(0.1936491673103708*favg[10]+0.1118033988749895*favg[4])+0.1118033988749895*(alpha[9]*favg[9]+alpha[8]*favg[8]); 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += 0.5*Ghat[2]*dv10r; 
  outr[3] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[4] += 0.5*Ghat[4]*dv10r; 
  outr[5] += 0.5*Ghat[5]*dv10r; 
  outr[6] += -0.8660254037844386*Ghat[1]*dv10r; 
  outr[7] += -0.8660254037844386*Ghat[2]*dv10r; 
  outr[8] += 0.5*Ghat[8]*dv10r; 
  outr[9] += 0.5*Ghat[9]*dv10r; 
  outr[10] += -0.8660254037844386*Ghat[4]*dv10r; 
  outr[11] += 0.5*Ghat[11]*dv10r; 
  outr[12] += 0.5*Ghat[12]*dv10r; 
  outr[13] += 1.118033988749895*Ghat[0]*dv10r; 
  outr[14] += 0.5*Ghat[14]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.5*Ghat[2]*dv10l; 
  outl[3] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[4] += -0.5*Ghat[4]*dv10l; 
  outl[5] += -0.5*Ghat[5]*dv10l; 
  outl[6] += -0.8660254037844386*Ghat[1]*dv10l; 
  outl[7] += -0.8660254037844386*Ghat[2]*dv10l; 
  outl[8] += -0.5*Ghat[8]*dv10l; 
  outl[9] += -0.5*Ghat[9]*dv10l; 
  outl[10] += -0.8660254037844386*Ghat[4]*dv10l; 
  outl[11] += -0.5*Ghat[11]*dv10l; 
  outl[12] += -0.5*Ghat[12]*dv10l; 
  outl[13] += -1.118033988749895*Ghat[0]*dv10l; 
  outl[14] += -0.5*Ghat[14]*dv10l; 

  return std::abs(amid); 
} 
__host__ __device__ double VlasovSurfElcMag2x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[3]; 
  double dv11r = 2/dxvr[3]; 
  const double *E1 = &EM[6]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double *B2 = &EM[30]; 

  double Ghat[15]; 
  double alpha[15]; 

  double favg[15]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = -1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 

  double fjump[15]; 
  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(-1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(-1*fr[8]-fl[8]); 
  fjump[9] = amax*(-1*fr[9]-fl[9]); 
  fjump[10] = amax*(-1*fr[10]-fl[10]); 
  fjump[11] = amax*(1*fr[11]-fl[11]); 
  fjump[12] = amax*(1*fr[12]-fl[12]); 
  fjump[13] = amax*(1*fr[13]-fl[13]); 
  fjump[14] = amax*(1*fr[14]-fl[14]); 

  alpha[0] = 2.0*E1[0]-2.0*B2[0]*wv1; 
  alpha[1] = 2.0*E1[1]-2.0*B2[1]*wv1; 
  alpha[2] = 2.0*E1[2]-2.0*B2[2]*wv1; 
  alpha[3] = -0.5773502691896258*B2[0]*dv1; 
  alpha[5] = 2.0*E1[3]-2.0*B2[3]*wv1; 
  alpha[6] = -0.5773502691896258*B2[1]*dv1; 
  alpha[7] = -0.5773502691896258*B2[2]*dv1; 
  alpha[11] = 2.0*E1[4]-2.0*B2[4]*wv1; 
  alpha[12] = 2.0*E1[5]-2.0*B2[5]*wv1; 
  const double amid = (-0.2795084971874737*alpha[12])-0.2795084971874737*alpha[11]+0.25*alpha[0]; 

  Ghat[0] = (-1.118033988749895*fjump[14])+alpha[0]*(0.2795084971874737*favg[14]+0.2165063509461096*favg[4]+0.125*favg[0])+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11])+alpha[3]*(0.2165063509461096*favg[10]+0.125*favg[3])+alpha[2]*(0.2165063509461096*favg[9]+0.125*favg[2])+alpha[1]*(0.2165063509461096*favg[8]+0.125*favg[1])+0.125*(alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5])-0.8660254037844386*fjump[4]-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.2795084971874737*favg[14]+0.1118033988749895*favg[11]+0.2165063509461096*favg[4]+0.125*favg[0])+(0.1936491673103708*favg[8]+0.1118033988749895*favg[1])*alpha[11]+alpha[6]*(0.2165063509461096*favg[10]+0.125*favg[3])+alpha[5]*(0.2165063509461096*favg[9]+0.125*favg[2])-0.8660254037844386*fjump[8]+alpha[0]*(0.2165063509461096*favg[8]+0.125*favg[1])+0.125*(alpha[3]*favg[6]+alpha[2]*favg[5])-0.5*fjump[1]; 
  Ghat[2] = alpha[2]*(0.2795084971874737*favg[14]+0.1118033988749895*favg[12]+0.2165063509461096*favg[4]+0.125*favg[0])+(0.1936491673103708*favg[9]+0.1118033988749895*favg[2])*alpha[12]+alpha[7]*(0.2165063509461096*favg[10]+0.125*favg[3])-0.8660254037844386*fjump[9]+alpha[0]*(0.2165063509461096*favg[9]+0.125*favg[2])+alpha[5]*(0.2165063509461096*favg[8]+0.125*favg[1])+0.125*(alpha[3]*favg[7]+alpha[1]*favg[5])-0.5*fjump[2]; 
  Ghat[3] = alpha[3]*(0.2795084971874737*favg[14]+0.1118033988749895*favg[13]+0.2165063509461096*favg[4]+0.125*favg[0])-0.8660254037844386*fjump[10]+alpha[0]*(0.2165063509461096*favg[10]+0.125*favg[3])+alpha[7]*(0.2165063509461096*favg[9]+0.125*favg[2])+alpha[6]*(0.2165063509461096*favg[8]+0.125*favg[1])+0.125*(alpha[2]*favg[7]+alpha[1]*favg[6])-0.5*fjump[3]; 
  Ghat[5] = alpha[5]*(0.2795084971874737*favg[14]+0.1118033988749895*(favg[12]+favg[11])+0.2165063509461096*favg[4]+0.125*favg[0])+0.1118033988749895*favg[5]*(alpha[12]+alpha[11])+alpha[1]*(0.2165063509461096*favg[9]+0.125*favg[2])+alpha[2]*(0.2165063509461096*favg[8]+0.125*favg[1])+0.125*(alpha[6]*favg[7]+favg[6]*alpha[7])-0.5*fjump[5]+0.125*alpha[0]*favg[5]; 
  Ghat[6] = alpha[6]*(0.2795084971874737*favg[14]+0.1118033988749895*(favg[13]+favg[11])+0.2165063509461096*favg[4]+0.125*favg[0])+0.1118033988749895*favg[6]*alpha[11]+alpha[1]*(0.2165063509461096*favg[10]+0.125*favg[3])+alpha[3]*(0.2165063509461096*favg[8]+0.125*favg[1])+0.125*(alpha[5]*favg[7]+favg[5]*alpha[7])-0.5*fjump[6]+0.125*alpha[0]*favg[6]; 
  Ghat[7] = alpha[7]*(0.2795084971874737*favg[14]+0.1118033988749895*(favg[13]+favg[12])+0.2165063509461096*favg[4]+0.125*favg[0])+0.1118033988749895*favg[7]*alpha[12]+alpha[2]*(0.2165063509461096*favg[10]+0.125*favg[3])+alpha[3]*(0.2165063509461096*favg[9]+0.125*favg[2])-0.5*fjump[7]+0.125*(alpha[0]*favg[7]+alpha[5]*favg[6]+favg[5]*alpha[6]); 
  Ghat[11] = alpha[11]*(0.2795084971874737*favg[14]+0.07985957062499249*favg[11]+0.2165063509461096*favg[4]+0.125*favg[0])-0.5*fjump[11]+0.125*alpha[0]*favg[11]+alpha[1]*(0.1936491673103708*favg[8]+0.1118033988749895*favg[1])+0.1118033988749895*(alpha[6]*favg[6]+alpha[5]*favg[5]); 
  Ghat[12] = alpha[12]*(0.2795084971874737*favg[14]+0.07985957062499249*favg[12]+0.2165063509461096*favg[4]+0.125*favg[0])-0.5*fjump[12]+0.125*alpha[0]*favg[12]+alpha[2]*(0.1936491673103708*favg[9]+0.1118033988749895*favg[2])+0.1118033988749895*(alpha[7]*favg[7]+alpha[5]*favg[5]); 
  Ghat[13] = (-0.5*fjump[13])+0.125*alpha[0]*favg[13]+alpha[3]*(0.1936491673103708*favg[10]+0.1118033988749895*favg[3])+0.1118033988749895*(alpha[7]*favg[7]+alpha[6]*favg[6]); 

  outr[0] += 0.5*Ghat[0]*dv11r; 
  outr[1] += 0.5*Ghat[1]*dv11r; 
  outr[2] += 0.5*Ghat[2]*dv11r; 
  outr[3] += 0.5*Ghat[3]*dv11r; 
  outr[4] += -0.8660254037844386*Ghat[0]*dv11r; 
  outr[5] += 0.5*Ghat[5]*dv11r; 
  outr[6] += 0.5*Ghat[6]*dv11r; 
  outr[7] += 0.5*Ghat[7]*dv11r; 
  outr[8] += -0.8660254037844386*Ghat[1]*dv11r; 
  outr[9] += -0.8660254037844386*Ghat[2]*dv11r; 
  outr[10] += -0.8660254037844386*Ghat[3]*dv11r; 
  outr[11] += 0.5*Ghat[11]*dv11r; 
  outr[12] += 0.5*Ghat[12]*dv11r; 
  outr[13] += 0.5*Ghat[13]*dv11r; 
  outr[14] += 1.118033988749895*Ghat[0]*dv11r; 

  outl[0] += -0.5*Ghat[0]*dv11l; 
  outl[1] += -0.5*Ghat[1]*dv11l; 
  outl[2] += -0.5*Ghat[2]*dv11l; 
  outl[3] += -0.5*Ghat[3]*dv11l; 
  outl[4] += -0.8660254037844386*Ghat[0]*dv11l; 
  outl[5] += -0.5*Ghat[5]*dv11l; 
  outl[6] += -0.5*Ghat[6]*dv11l; 
  outl[7] += -0.5*Ghat[7]*dv11l; 
  outl[8] += -0.8660254037844386*Ghat[1]*dv11l; 
  outl[9] += -0.8660254037844386*Ghat[2]*dv11l; 
  outl[10] += -0.8660254037844386*Ghat[3]*dv11l; 
  outl[11] += -0.5*Ghat[11]*dv11l; 
  outl[12] += -0.5*Ghat[12]*dv11l; 
  outl[13] += -0.5*Ghat[13]*dv11l; 
  outl[14] += -1.118033988749895*Ghat[0]*dv11l; 

  return std::abs(amid); 
} 
