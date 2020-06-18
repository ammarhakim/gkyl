#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag1x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double dv3 = dxvr[3], wv3 = wr[3]; 
  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double Ghat[10]; 
  double favg[10]; 
  double alpha[10]; 

  favg[0] = 1.58113883008419*fr[12]+1.58113883008419*fl[12]-1.224744871391589*fr[2]+1.224744871391589*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[5])+1.224744871391589*fl[5]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[7])+1.224744871391589*fl[7]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[3] = (-1.224744871391589*fr[9])+1.224744871391589*fl[9]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = 0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[5] = 0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[6] = 0.7071067811865475*fr[10]+0.7071067811865475*fl[10]; 
  favg[7] = 0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 
  favg[8] = 0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 
  favg[9] = 0.7071067811865475*fr[14]+0.7071067811865475*fl[14]; 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 0.5773502691896258*B2[0]*dv2; 
  alpha[3] = -0.5773502691896258*B1[0]*dv3; 
  alpha[4] = 0.5773502691896258*B2[1]*dv2; 
  alpha[5] = -0.5773502691896258*B1[1]*dv3; 
  alpha[7] = 2.0*(B2[2]*wv2+E0[2])-2.0*B1[2]*wv3; 

  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 

  Ghat[0] = 0.1767766952966368*(alpha[7]*favg[7]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[12]-1.0*(2.23606797749979*fl[12]+1.732050807568877*(fr[2]+fl[2]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[5]+fl[5])-1.0*fr[1]+fl[1])*amax+0.03535533905932737*(4.47213595499958*(alpha[1]*favg[7]+favg[1]*alpha[7])+5.0*(alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[7]+fl[7])-1.0*fr[3]+fl[3])*amax+0.03535533905932737*(4.47213595499958*alpha[2]*favg[8]+5.0*(alpha[3]*favg[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[9]+fl[9])-1.0*fr[4]+fl[4])*amax+0.03535533905932737*(4.47213595499958*alpha[3]*favg[9]+5.0*(alpha[2]*favg[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3])); 
  Ghat[4] = 0.03535533905932737*(4.47213595499958*(alpha[4]*(favg[8]+favg[7])+favg[4]*alpha[7])+5.0*(alpha[5]*favg[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]))-0.3535533905932737*(fr[6]-1.0*fl[6])*amax; 
  Ghat[5] = 0.03535533905932737*(4.47213595499958*(alpha[5]*(favg[9]+favg[7])+favg[5]*alpha[7])+5.0*(alpha[4]*favg[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3]))-0.3535533905932737*(fr[8]-1.0*fl[8])*amax; 
  Ghat[6] = 0.1767766952966368*(alpha[0]*favg[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3])-0.3535533905932737*(fr[10]-1.0*fl[10])*amax; 
  Ghat[7] = 0.005050762722761053*((22.3606797749979*alpha[7]+35.0*alpha[0])*favg[7]+35.0*favg[0]*alpha[7]+31.30495168499706*(alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[1]*favg[1]))-0.3535533905932737*(fr[11]-1.0*fl[11])*amax; 
  Ghat[8] = 0.03535533905932737*(5.0*alpha[0]*favg[8]+4.47213595499958*(alpha[4]*favg[4]+alpha[2]*favg[2]))-0.3535533905932737*(fr[13]-1.0*fl[13])*amax; 
  Ghat[9] = 0.03535533905932737*(5.0*alpha[0]*favg[9]+4.47213595499958*(alpha[5]*favg[5]+alpha[3]*favg[3]))-0.3535533905932737*(fr[14]-1.0*fl[14])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[3] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv10r; 
  outr[5] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[6] += 0.7071067811865475*Ghat[4]*dv10r; 
  outr[7] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[8] += 0.7071067811865475*Ghat[5]*dv10r; 
  outr[9] += -1.224744871391589*Ghat[3]*dv10r; 
  outr[10] += 0.7071067811865475*Ghat[6]*dv10r; 
  outr[11] += 0.7071067811865475*Ghat[7]*dv10r; 
  outr[12] += 1.58113883008419*Ghat[0]*dv10r; 
  outr[13] += 0.7071067811865475*Ghat[8]*dv10r; 
  outr[14] += 0.7071067811865475*Ghat[9]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[3] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv10l; 
  outl[5] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[6] += -0.7071067811865475*Ghat[4]*dv10l; 
  outl[7] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[8] += -0.7071067811865475*Ghat[5]*dv10l; 
  outl[9] += -1.224744871391589*Ghat[3]*dv10l; 
  outl[10] += -0.7071067811865475*Ghat[6]*dv10l; 
  outl[11] += -0.7071067811865475*Ghat[7]*dv10l; 
  outl[12] += -1.58113883008419*Ghat[0]*dv10l; 
  outl[13] += -0.7071067811865475*Ghat[8]*dv10l; 
  outl[14] += -0.7071067811865475*Ghat[9]*dv10l; 

  return std::abs(amid); 
} 
__host__ __device__ double VlasovSurfElcMag1x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[2]; 
  double dv11r = 2/dxvr[2]; 
  const double *E1 = &EM[3]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double dv3 = dxvr[3], wv3 = wr[3]; 
  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double Ghat[10]; 
  double favg[10]; 
  double alpha[10]; 

  favg[0] = 1.58113883008419*fr[13]+1.58113883008419*fl[13]-1.224744871391589*fr[3]+1.224744871391589*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[6])+1.224744871391589*fl[6]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[7])+1.224744871391589*fl[7]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[10])+1.224744871391589*fl[10]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = 0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = 0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[6] = 0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[7] = 0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 
  favg[8] = 0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[9] = 0.7071067811865475*fr[14]+0.7071067811865475*fl[14]; 

  alpha[0] = 2.0*B0[0]*wv3-2.0*B2[0]*wv1+2.0*E1[0]; 
  alpha[1] = 2.0*B0[1]*wv3-2.0*B2[1]*wv1+2.0*E1[1]; 
  alpha[2] = -0.5773502691896258*B2[0]*dv1; 
  alpha[3] = 0.5773502691896258*B0[0]*dv3; 
  alpha[4] = -0.5773502691896258*B2[1]*dv1; 
  alpha[5] = 0.5773502691896258*B0[1]*dv3; 
  alpha[7] = 2.0*B0[2]*wv3-2.0*B2[2]*wv1+2.0*E1[2]; 

  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 

  Ghat[0] = 0.1767766952966368*(alpha[7]*favg[7]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[13]-1.0*(2.23606797749979*fl[13]+1.732050807568877*(fr[3]+fl[3]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[6]+fl[6])-1.0*fr[1]+fl[1])*amax+0.03535533905932737*(4.47213595499958*(alpha[1]*favg[7]+favg[1]*alpha[7])+5.0*(alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[7]+fl[7])-1.0*fr[2]+fl[2])*amax+0.03535533905932737*(4.47213595499958*alpha[2]*favg[8]+5.0*(alpha[3]*favg[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[10]+fl[10])-1.0*fr[4]+fl[4])*amax+0.03535533905932737*(4.47213595499958*alpha[3]*favg[9]+5.0*(alpha[2]*favg[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3])); 
  Ghat[4] = 0.03535533905932737*(4.47213595499958*(alpha[4]*(favg[8]+favg[7])+favg[4]*alpha[7])+5.0*(alpha[5]*favg[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]))-0.3535533905932737*(fr[5]-1.0*fl[5])*amax; 
  Ghat[5] = 0.03535533905932737*(4.47213595499958*(alpha[5]*(favg[9]+favg[7])+favg[5]*alpha[7])+5.0*(alpha[4]*favg[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3]))-0.3535533905932737*(fr[8]-1.0*fl[8])*amax; 
  Ghat[6] = 0.1767766952966368*(alpha[0]*favg[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3])-0.3535533905932737*(fr[9]-1.0*fl[9])*amax; 
  Ghat[7] = 0.005050762722761053*((22.3606797749979*alpha[7]+35.0*alpha[0])*favg[7]+35.0*favg[0]*alpha[7]+31.30495168499706*(alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[1]*favg[1]))-0.3535533905932737*(fr[11]-1.0*fl[11])*amax; 
  Ghat[8] = 0.03535533905932737*(5.0*alpha[0]*favg[8]+4.47213595499958*(alpha[4]*favg[4]+alpha[2]*favg[2]))-0.3535533905932737*(fr[12]-1.0*fl[12])*amax; 
  Ghat[9] = 0.03535533905932737*(5.0*alpha[0]*favg[9]+4.47213595499958*(alpha[5]*favg[5]+alpha[3]*favg[3]))-0.3535533905932737*(fr[14]-1.0*fl[14])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv11r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv11r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv11r; 
  outr[3] += -1.224744871391589*Ghat[0]*dv11r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv11r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv11r; 
  outr[6] += -1.224744871391589*Ghat[1]*dv11r; 
  outr[7] += -1.224744871391589*Ghat[2]*dv11r; 
  outr[8] += 0.7071067811865475*Ghat[5]*dv11r; 
  outr[9] += 0.7071067811865475*Ghat[6]*dv11r; 
  outr[10] += -1.224744871391589*Ghat[3]*dv11r; 
  outr[11] += 0.7071067811865475*Ghat[7]*dv11r; 
  outr[12] += 0.7071067811865475*Ghat[8]*dv11r; 
  outr[13] += 1.58113883008419*Ghat[0]*dv11r; 
  outr[14] += 0.7071067811865475*Ghat[9]*dv11r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv11l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv11l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv11l; 
  outl[3] += -1.224744871391589*Ghat[0]*dv11l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv11l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv11l; 
  outl[6] += -1.224744871391589*Ghat[1]*dv11l; 
  outl[7] += -1.224744871391589*Ghat[2]*dv11l; 
  outl[8] += -0.7071067811865475*Ghat[5]*dv11l; 
  outl[9] += -0.7071067811865475*Ghat[6]*dv11l; 
  outl[10] += -1.224744871391589*Ghat[3]*dv11l; 
  outl[11] += -0.7071067811865475*Ghat[7]*dv11l; 
  outl[12] += -0.7071067811865475*Ghat[8]*dv11l; 
  outl[13] += -1.58113883008419*Ghat[0]*dv11l; 
  outl[14] += -0.7071067811865475*Ghat[9]*dv11l; 

  return std::abs(amid); 
} 
__host__ __device__ double VlasovSurfElcMag1x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12l = 2/dxvl[3]; 
  double dv12r = 2/dxvr[3]; 
  const double *E2 = &EM[6]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double dv3 = dxvr[3], wv3 = wr[3]; 
  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double Ghat[10]; 
  double favg[10]; 
  double alpha[10]; 

  favg[0] = 1.58113883008419*fr[14]+1.58113883008419*fl[14]-1.224744871391589*fr[4]+1.224744871391589*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[8])+1.224744871391589*fl[8]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[9])+1.224744871391589*fl[9]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[10])+1.224744871391589*fl[10]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = 0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = 0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = 0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = 0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 
  favg[8] = 0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[9] = 0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 

  alpha[0] = 2.0*(B1[0]*wv1+E2[0])-2.0*B0[0]*wv2; 
  alpha[1] = 2.0*(B1[1]*wv1+E2[1])-2.0*B0[1]*wv2; 
  alpha[2] = 0.5773502691896258*B1[0]*dv1; 
  alpha[3] = -0.5773502691896258*B0[0]*dv2; 
  alpha[4] = 0.5773502691896258*B1[1]*dv1; 
  alpha[5] = -0.5773502691896258*B0[1]*dv2; 
  alpha[7] = 2.0*(B1[2]*wv1+E2[2])-2.0*B0[2]*wv2; 

  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 

  Ghat[0] = 0.1767766952966368*(alpha[7]*favg[7]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[14]-1.0*(2.23606797749979*fl[14]+1.732050807568877*(fr[4]+fl[4]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[8]+fl[8])-1.0*fr[1]+fl[1])*amax+0.03535533905932737*(4.47213595499958*(alpha[1]*favg[7]+favg[1]*alpha[7])+5.0*(alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[9]+fl[9])-1.0*fr[2]+fl[2])*amax+0.03535533905932737*(4.47213595499958*alpha[2]*favg[8]+5.0*(alpha[3]*favg[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[10]+fl[10])-1.0*fr[3]+fl[3])*amax+0.03535533905932737*(4.47213595499958*alpha[3]*favg[9]+5.0*(alpha[2]*favg[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3])); 
  Ghat[4] = 0.03535533905932737*(4.47213595499958*(alpha[4]*(favg[8]+favg[7])+favg[4]*alpha[7])+5.0*(alpha[5]*favg[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]))-0.3535533905932737*(fr[5]-1.0*fl[5])*amax; 
  Ghat[5] = 0.03535533905932737*(4.47213595499958*(alpha[5]*(favg[9]+favg[7])+favg[5]*alpha[7])+5.0*(alpha[4]*favg[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3]))-0.3535533905932737*(fr[6]-1.0*fl[6])*amax; 
  Ghat[6] = 0.1767766952966368*(alpha[0]*favg[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3])-0.3535533905932737*(fr[7]-1.0*fl[7])*amax; 
  Ghat[7] = 0.005050762722761053*((22.3606797749979*alpha[7]+35.0*alpha[0])*favg[7]+35.0*favg[0]*alpha[7]+31.30495168499706*(alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[1]*favg[1]))-0.3535533905932737*(fr[11]-1.0*fl[11])*amax; 
  Ghat[8] = 0.03535533905932737*(5.0*alpha[0]*favg[8]+4.47213595499958*(alpha[4]*favg[4]+alpha[2]*favg[2]))-0.3535533905932737*(fr[12]-1.0*fl[12])*amax; 
  Ghat[9] = 0.03535533905932737*(5.0*alpha[0]*favg[9]+4.47213595499958*(alpha[5]*favg[5]+alpha[3]*favg[3]))-0.3535533905932737*(fr[13]-1.0*fl[13])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv12r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv12r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv12r; 
  outr[3] += 0.7071067811865475*Ghat[3]*dv12r; 
  outr[4] += -1.224744871391589*Ghat[0]*dv12r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv12r; 
  outr[6] += 0.7071067811865475*Ghat[5]*dv12r; 
  outr[7] += 0.7071067811865475*Ghat[6]*dv12r; 
  outr[8] += -1.224744871391589*Ghat[1]*dv12r; 
  outr[9] += -1.224744871391589*Ghat[2]*dv12r; 
  outr[10] += -1.224744871391589*Ghat[3]*dv12r; 
  outr[11] += 0.7071067811865475*Ghat[7]*dv12r; 
  outr[12] += 0.7071067811865475*Ghat[8]*dv12r; 
  outr[13] += 0.7071067811865475*Ghat[9]*dv12r; 
  outr[14] += 1.58113883008419*Ghat[0]*dv12r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv12l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv12l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv12l; 
  outl[3] += -0.7071067811865475*Ghat[3]*dv12l; 
  outl[4] += -1.224744871391589*Ghat[0]*dv12l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv12l; 
  outl[6] += -0.7071067811865475*Ghat[5]*dv12l; 
  outl[7] += -0.7071067811865475*Ghat[6]*dv12l; 
  outl[8] += -1.224744871391589*Ghat[1]*dv12l; 
  outl[9] += -1.224744871391589*Ghat[2]*dv12l; 
  outl[10] += -1.224744871391589*Ghat[3]*dv12l; 
  outl[11] += -0.7071067811865475*Ghat[7]*dv12l; 
  outl[12] += -0.7071067811865475*Ghat[8]*dv12l; 
  outl[13] += -0.7071067811865475*Ghat[9]*dv12l; 
  outl[14] += -1.58113883008419*Ghat[0]*dv12l; 

  return std::abs(amid); 
} 
