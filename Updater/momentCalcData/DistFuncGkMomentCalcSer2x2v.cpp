#include <math.h> 
#include <DistFuncMomentCalcModDecl.h> 
__host__ __device__ void GkMomentCalc2x2vSer_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
  out[6] += 2.0*f[19]*volFact; 
  out[7] += 2.0*f[20]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M0_P3(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
  out[6] += 2.0*f[19]*volFact; 
  out[7] += 2.0*f[20]*volFact; 
  out[8] += 2.0*f[31]*volFact; 
  out[9] += 2.0*f[32]*volFact; 
  out[10] += 2.0*f[48]*volFact; 
  out[11] += 2.0*f[49]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1); 
  out[2] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1); 
  out[3] += volFact*(2.0*f[5]*wx1+0.5773502691896258*f[11]*dv1); 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1); 
  out[2] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1); 
  out[3] += volFact*(2.0*f[5]*wx1+0.5773502691896258*f[15]*dv1); 
  out[4] += volFact*(2.0*f[11]*wx1+0.5773502691896257*f[21]*dv1); 
  out[5] += volFact*(2.0*f[12]*wx1+0.5773502691896257*f[22]*dv1); 
  out[6] += volFact*(2.0*f[19]*wx1+0.5773502691896257*f[32]*dv1); 
  out[7] += volFact*(2.0*f[20]*wx1+0.5773502691896257*f[33]*dv1); 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M1_P3(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1); 
  out[2] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1); 
  out[3] += volFact*(2.0*f[5]*wx1+0.5773502691896258*f[15]*dv1); 
  out[4] += volFact*(2.0*f[11]*wx1+0.5773502691896257*f[21]*dv1); 
  out[5] += volFact*(2.0*f[12]*wx1+0.5773502691896257*f[22]*dv1); 
  out[6] += volFact*(2.0*f[19]*wx1+0.5773502691896257*f[36]*dv1); 
  out[7] += volFact*(2.0*f[20]*wx1+0.5773502691896257*f[37]*dv1); 
  out[8] += volFact*(2.0*f[31]*wx1+0.5773502691896256*f[50]*dv1); 
  out[9] += volFact*(2.0*f[32]*wx1+0.5773502691896256*f[51]*dv1); 
  out[10] += volFact*(2.0*f[48]*wx1+0.5773502691896256*f[64]*dv1); 
  out[11] += volFact*(2.0*f[49]*wx1+0.5773502691896256*f[65]*dv1); 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M1proj_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += 2.0*f[0]*volFact*wx1; 
  out[1] += 2.0*f[1]*volFact*wx1; 
  out[2] += 2.0*f[2]*volFact*wx1; 
  out[3] += 2.0*f[5]*volFact*wx1; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M1proj_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1); 
  out[2] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1); 
  out[3] += volFact*(2.0*f[5]*wx1+0.5773502691896258*f[15]*dv1); 
  out[4] += volFact*(2.0*f[11]*wx1+0.5773502691896257*f[21]*dv1); 
  out[5] += volFact*(2.0*f[12]*wx1+0.5773502691896257*f[22]*dv1); 
  out[6] += volFact*(2.0*f[19]*wx1+0.5773502691896257*f[32]*dv1); 
  out[7] += volFact*(2.0*f[20]*wx1+0.5773502691896257*f[33]*dv1); 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M1proj_P3(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1); 
  out[2] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1); 
  out[3] += volFact*(2.0*f[5]*wx1+0.5773502691896258*f[15]*dv1); 
  out[4] += volFact*(2.0*f[11]*wx1+0.5773502691896257*f[21]*dv1); 
  out[5] += volFact*(2.0*f[12]*wx1+0.5773502691896257*f[22]*dv1); 
  out[6] += volFact*(2.0*f[19]*wx1+0.5773502691896257*f[36]*dv1); 
  out[7] += volFact*(2.0*f[20]*wx1+0.5773502691896257*f[37]*dv1); 
  out[8] += volFact*(2.0*f[31]*wx1+0.5773502691896256*f[50]*dv1); 
  out[9] += volFact*(2.0*f[32]*wx1+0.5773502691896256*f[51]*dv1); 
  out[10] += volFact*(2.0*f[48]*wx1+0.5773502691896256*f[64]*dv1); 
  out[11] += volFact*(2.0*f[49]*wx1+0.5773502691896256*f[65]*dv1); 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.1666666666666667*f[5]*dv1_sq); 
  double tmp[4]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[12]*dv2; 
  out[0] += (2.0*(0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[23]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[24]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[15]*dv1*wx1+0.149071198499986*f[34]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
  out[4] += volFact*(2.0*f[11]*wx1_sq+1.154700538379251*f[21]*dv1*wx1+0.1666666666666667*f[11]*dv1_sq); 
  out[5] += volFact*(2.0*f[12]*wx1_sq+1.154700538379251*f[22]*dv1*wx1+0.1666666666666667*f[12]*dv1_sq); 
  out[6] += volFact*(2.0*f[19]*wx1_sq+1.154700538379251*f[32]*dv1*wx1+0.1666666666666667*f[19]*dv1_sq); 
  out[7] += volFact*(2.0*f[20]*wx1_sq+1.154700538379251*f[33]*dv1*wx1+0.1666666666666667*f[20]*dv1_sq); 
  double tmp[8]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[16]*dv2; 
  tmp[4] = 2.0*f[11]*wx2+0.5773502691896257*f[25]*dv2; 
  tmp[5] = 2.0*f[12]*wx2+0.5773502691896257*f[26]*dv2; 
  tmp[6] = 2.0*f[19]*wx2+0.5773502691896257*f[35]*dv2; 
  tmp[7] = 2.0*f[20]*wx2+0.5773502691896257*f[36]*dv2; 
  out[0] += (2.0*(0.5*Bmag[7]*tmp[7]+0.5*Bmag[6]*tmp[6]+0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.5000000000000001*Bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*Bmag[7]+0.447213595499958*Bmag[3]*tmp[6]+0.447213595499958*tmp[3]*Bmag[6]+0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.447213595499958*Bmag[3]*tmp[7]+0.447213595499958*tmp[3]*Bmag[7]+0.5000000000000001*Bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*Bmag[6]+0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.4*Bmag[6]*tmp[7]+0.447213595499958*Bmag[2]*tmp[7]+0.4*tmp[6]*Bmag[7]+0.447213595499958*tmp[2]*Bmag[7]+0.447213595499958*Bmag[1]*tmp[6]+0.447213595499958*tmp[1]*Bmag[6]+0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  out[4] += (2.0*(0.4472135954999579*Bmag[7]*tmp[7]+0.31943828249997*Bmag[6]*tmp[6]+0.5000000000000001*Bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*Bmag[6]+0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  out[5] += (2.0*(0.31943828249997*Bmag[7]*tmp[7]+0.5000000000000001*Bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*Bmag[7]+0.4472135954999579*Bmag[6]*tmp[6]+0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
  out[6] += (2.0*(0.4*Bmag[3]*tmp[7]+0.4*tmp[3]*Bmag[7]+0.4472135954999579*Bmag[5]*tmp[6]+0.31943828249997*Bmag[4]*tmp[6]+0.5*Bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*Bmag[6]+0.31943828249997*tmp[4]*Bmag[6]+0.5*tmp[0]*Bmag[6]+0.5000000000000001*Bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*Bmag[4]+0.447213595499958*Bmag[1]*tmp[3]+0.447213595499958*tmp[1]*Bmag[3])*volFact)/m_; 
  out[7] += (2.0*(0.31943828249997*Bmag[5]*tmp[7]+0.4472135954999579*Bmag[4]*tmp[7]+0.5*Bmag[0]*tmp[7]+0.31943828249997*tmp[5]*Bmag[7]+0.4472135954999579*tmp[4]*Bmag[7]+0.5*tmp[0]*Bmag[7]+0.4*Bmag[3]*tmp[6]+0.4*tmp[3]*Bmag[6]+0.5000000000000001*Bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*Bmag[5]+0.447213595499958*Bmag[2]*tmp[3]+0.447213595499958*tmp[2]*Bmag[3])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M2_P3(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[23]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[24]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[15]*dv1*wx1+0.149071198499986*f[38]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
  out[4] += volFact*(2.0*f[11]*wx1_sq+1.154700538379251*f[21]*dv1*wx1+0.1666666666666667*f[11]*dv1_sq); 
  out[5] += volFact*(2.0*f[12]*wx1_sq+1.154700538379251*f[22]*dv1*wx1+0.1666666666666667*f[12]*dv1_sq); 
  out[6] += volFact*(2.0*f[19]*wx1_sq+1.154700538379251*f[36]*dv1*wx1+0.1666666666666667*f[19]*dv1_sq); 
  out[7] += volFact*(2.0*f[20]*wx1_sq+1.154700538379251*f[37]*dv1*wx1+0.1666666666666667*f[20]*dv1_sq); 
  out[8] += volFact*(2.0*f[31]*wx1_sq+1.154700538379251*f[50]*dv1*wx1+0.1666666666666667*f[31]*dv1_sq); 
  out[9] += volFact*(2.0*f[32]*wx1_sq+1.154700538379251*f[51]*dv1*wx1+0.1666666666666667*f[32]*dv1_sq); 
  out[10] += volFact*(2.0*f[48]*wx1_sq+1.154700538379251*f[64]*dv1*wx1+0.1666666666666667*f[48]*dv1_sq); 
  out[11] += volFact*(2.0*f[49]*wx1_sq+1.154700538379251*f[65]*dv1*wx1+0.1666666666666667*f[49]*dv1_sq); 
  double tmp[12]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[16]*dv2; 
  tmp[4] = 2.0*f[11]*wx2+0.5773502691896257*f[25]*dv2; 
  tmp[5] = 2.0*f[12]*wx2+0.5773502691896257*f[26]*dv2; 
  tmp[6] = 2.0*f[19]*wx2+0.5773502691896257*f[39]*dv2; 
  tmp[7] = 2.0*f[20]*wx2+0.5773502691896257*f[40]*dv2; 
  tmp[8] = 2.0*f[31]*wx2+0.5773502691896256*f[54]*dv2; 
  tmp[9] = 2.0*f[32]*wx2+0.5773502691896256*f[55]*dv2; 
  tmp[10] = 2.0*f[48]*wx2+0.5773502691896256*f[67]*dv2; 
  tmp[11] = 2.0*f[49]*wx2+0.5773502691896256*f[68]*dv2; 
  out[0] += (2.0*(0.5*Bmag[11]*tmp[11]+0.5*Bmag[10]*tmp[10]+0.5*Bmag[9]*tmp[9]+0.5*Bmag[8]*tmp[8]+0.5*Bmag[7]*tmp[7]+0.5*Bmag[6]*tmp[6]+0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.5*Bmag[9]*tmp[11]+0.5*tmp[9]*Bmag[11]+0.4391550328268399*Bmag[6]*tmp[10]+0.4391550328268399*tmp[6]*Bmag[10]+0.4391550328268398*Bmag[4]*tmp[8]+0.4391550328268398*tmp[4]*Bmag[8]+0.5000000000000001*Bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*Bmag[7]+0.447213595499958*Bmag[3]*tmp[6]+0.447213595499958*tmp[3]*Bmag[6]+0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.4391550328268399*Bmag[7]*tmp[11]+0.4391550328268399*tmp[7]*Bmag[11]+0.5*Bmag[8]*tmp[10]+0.5*tmp[8]*Bmag[10]+0.4391550328268398*Bmag[5]*tmp[9]+0.4391550328268398*tmp[5]*Bmag[9]+0.447213595499958*Bmag[3]*tmp[7]+0.447213595499958*tmp[3]*Bmag[7]+0.5000000000000001*Bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*Bmag[6]+0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.4391550328268399*Bmag[5]*tmp[11]+0.4391550328268399*tmp[5]*Bmag[11]+0.4391550328268399*Bmag[4]*tmp[10]+0.4391550328268399*tmp[4]*Bmag[10]+0.4391550328268399*Bmag[7]*tmp[9]+0.4391550328268399*tmp[7]*Bmag[9]+0.4391550328268399*Bmag[6]*tmp[8]+0.4391550328268399*tmp[6]*Bmag[8]+0.4*Bmag[6]*tmp[7]+0.447213595499958*Bmag[2]*tmp[7]+0.4*tmp[6]*Bmag[7]+0.447213595499958*tmp[2]*Bmag[7]+0.447213595499958*Bmag[1]*tmp[6]+0.447213595499958*tmp[1]*Bmag[6]+0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  out[4] += (2.0*(0.4472135954999579*Bmag[11]*tmp[11]+0.2981423969999719*Bmag[10]*tmp[10]+0.4391550328268399*Bmag[3]*tmp[10]+0.4391550328268399*tmp[3]*Bmag[10]+0.2981423969999719*Bmag[8]*tmp[8]+0.4391550328268398*Bmag[1]*tmp[8]+0.4391550328268398*tmp[1]*Bmag[8]+0.4472135954999579*Bmag[7]*tmp[7]+0.31943828249997*Bmag[6]*tmp[6]+0.5000000000000001*Bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*Bmag[6]+0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  out[5] += (2.0*(0.2981423969999719*Bmag[11]*tmp[11]+0.4391550328268399*Bmag[3]*tmp[11]+0.4391550328268399*tmp[3]*Bmag[11]+0.4472135954999579*Bmag[10]*tmp[10]+0.2981423969999719*Bmag[9]*tmp[9]+0.4391550328268398*Bmag[2]*tmp[9]+0.4391550328268398*tmp[2]*Bmag[9]+0.31943828249997*Bmag[7]*tmp[7]+0.5000000000000001*Bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*Bmag[7]+0.4472135954999579*Bmag[6]*tmp[6]+0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
  out[6] += (2.0*(0.3927922024247863*Bmag[7]*tmp[11]+0.3927922024247863*tmp[7]*Bmag[11]+0.2981423969999719*Bmag[8]*tmp[10]+0.3927922024247863*Bmag[7]*tmp[10]+0.4391550328268399*Bmag[1]*tmp[10]+0.2981423969999719*tmp[8]*Bmag[10]+0.3927922024247863*tmp[7]*Bmag[10]+0.4391550328268399*tmp[1]*Bmag[10]+0.4391550328268399*Bmag[3]*tmp[8]+0.4391550328268399*tmp[3]*Bmag[8]+0.4*Bmag[3]*tmp[7]+0.4*tmp[3]*Bmag[7]+0.4472135954999579*Bmag[5]*tmp[6]+0.31943828249997*Bmag[4]*tmp[6]+0.5*Bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*Bmag[6]+0.31943828249997*tmp[4]*Bmag[6]+0.5*tmp[0]*Bmag[6]+0.5000000000000001*Bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*Bmag[4]+0.447213595499958*Bmag[1]*tmp[3]+0.447213595499958*tmp[1]*Bmag[3])*volFact)/m_; 
  out[7] += (2.0*(0.2981423969999719*Bmag[9]*tmp[11]+0.3927922024247863*Bmag[6]*tmp[11]+0.4391550328268399*Bmag[2]*tmp[11]+0.2981423969999719*tmp[9]*Bmag[11]+0.3927922024247863*tmp[6]*Bmag[11]+0.4391550328268399*tmp[2]*Bmag[11]+0.3927922024247863*Bmag[6]*tmp[10]+0.3927922024247863*tmp[6]*Bmag[10]+0.4391550328268399*Bmag[3]*tmp[9]+0.4391550328268399*tmp[3]*Bmag[9]+0.31943828249997*Bmag[5]*tmp[7]+0.4472135954999579*Bmag[4]*tmp[7]+0.5*Bmag[0]*tmp[7]+0.31943828249997*tmp[5]*Bmag[7]+0.4472135954999579*tmp[4]*Bmag[7]+0.5*tmp[0]*Bmag[7]+0.4*Bmag[3]*tmp[6]+0.4*tmp[3]*Bmag[6]+0.5000000000000001*Bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*Bmag[5]+0.447213595499958*Bmag[2]*tmp[3]+0.447213595499958*tmp[2]*Bmag[3])*volFact)/m_; 
  out[8] += (2.0*(0.2981423969999719*Bmag[6]*tmp[10]+0.5*Bmag[2]*tmp[10]+0.2981423969999719*tmp[6]*Bmag[10]+0.5*tmp[2]*Bmag[10]+0.2981423969999719*Bmag[4]*tmp[8]+0.5*Bmag[0]*tmp[8]+0.2981423969999719*tmp[4]*Bmag[8]+0.5*tmp[0]*Bmag[8]+0.4391550328268399*Bmag[3]*tmp[6]+0.4391550328268399*tmp[3]*Bmag[6]+0.4391550328268398*Bmag[1]*tmp[4]+0.4391550328268398*tmp[1]*Bmag[4])*volFact)/m_; 
  out[9] += (2.0*(0.2981423969999719*Bmag[7]*tmp[11]+0.5*Bmag[1]*tmp[11]+0.2981423969999719*tmp[7]*Bmag[11]+0.5*tmp[1]*Bmag[11]+0.2981423969999719*Bmag[5]*tmp[9]+0.5*Bmag[0]*tmp[9]+0.2981423969999719*tmp[5]*Bmag[9]+0.5*tmp[0]*Bmag[9]+0.4391550328268399*Bmag[3]*tmp[7]+0.4391550328268399*tmp[3]*Bmag[7]+0.4391550328268398*Bmag[2]*tmp[5]+0.4391550328268398*tmp[2]*Bmag[5])*volFact)/m_; 
  out[10] += (2.0*(0.4472135954999579*Bmag[5]*tmp[10]+0.2981423969999719*Bmag[4]*tmp[10]+0.5*Bmag[0]*tmp[10]+0.4472135954999579*tmp[5]*Bmag[10]+0.2981423969999719*tmp[4]*Bmag[10]+0.5*tmp[0]*Bmag[10]+0.2981423969999719*Bmag[6]*tmp[8]+0.5*Bmag[2]*tmp[8]+0.2981423969999719*tmp[6]*Bmag[8]+0.5*tmp[2]*Bmag[8]+0.3927922024247863*Bmag[6]*tmp[7]+0.3927922024247863*tmp[6]*Bmag[7]+0.4391550328268399*Bmag[1]*tmp[6]+0.4391550328268399*tmp[1]*Bmag[6]+0.4391550328268399*Bmag[3]*tmp[4]+0.4391550328268399*tmp[3]*Bmag[4])*volFact)/m_; 
  out[11] += (2.0*(0.2981423969999719*Bmag[5]*tmp[11]+0.4472135954999579*Bmag[4]*tmp[11]+0.5*Bmag[0]*tmp[11]+0.2981423969999719*tmp[5]*Bmag[11]+0.4472135954999579*tmp[4]*Bmag[11]+0.5*tmp[0]*Bmag[11]+0.2981423969999719*Bmag[7]*tmp[9]+0.5*Bmag[1]*tmp[9]+0.2981423969999719*tmp[7]*Bmag[9]+0.5*tmp[1]*Bmag[9]+0.3927922024247863*Bmag[6]*tmp[7]+0.4391550328268399*Bmag[2]*tmp[7]+0.3927922024247863*tmp[6]*Bmag[7]+0.4391550328268399*tmp[2]*Bmag[7]+0.4391550328268399*Bmag[3]*tmp[5]+0.4391550328268399*tmp[3]*Bmag[5])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.1666666666666667*f[5]*dv1_sq); 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[23]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[24]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[15]*dv1*wx1+0.149071198499986*f[34]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
  out[4] += volFact*(2.0*f[11]*wx1_sq+1.154700538379251*f[21]*dv1*wx1+0.1666666666666667*f[11]*dv1_sq); 
  out[5] += volFact*(2.0*f[12]*wx1_sq+1.154700538379251*f[22]*dv1*wx1+0.1666666666666667*f[12]*dv1_sq); 
  out[6] += volFact*(2.0*f[19]*wx1_sq+1.154700538379251*f[32]*dv1*wx1+0.1666666666666667*f[19]*dv1_sq); 
  out[7] += volFact*(2.0*f[20]*wx1_sq+1.154700538379251*f[33]*dv1*wx1+0.1666666666666667*f[20]*dv1_sq); 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M2par_P3(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[23]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[24]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[15]*dv1*wx1+0.149071198499986*f[38]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
  out[4] += volFact*(2.0*f[11]*wx1_sq+1.154700538379251*f[21]*dv1*wx1+0.1666666666666667*f[11]*dv1_sq); 
  out[5] += volFact*(2.0*f[12]*wx1_sq+1.154700538379251*f[22]*dv1*wx1+0.1666666666666667*f[12]*dv1_sq); 
  out[6] += volFact*(2.0*f[19]*wx1_sq+1.154700538379251*f[36]*dv1*wx1+0.1666666666666667*f[19]*dv1_sq); 
  out[7] += volFact*(2.0*f[20]*wx1_sq+1.154700538379251*f[37]*dv1*wx1+0.1666666666666667*f[20]*dv1_sq); 
  out[8] += volFact*(2.0*f[31]*wx1_sq+1.154700538379251*f[50]*dv1*wx1+0.1666666666666667*f[31]*dv1_sq); 
  out[9] += volFact*(2.0*f[32]*wx1_sq+1.154700538379251*f[51]*dv1*wx1+0.1666666666666667*f[32]*dv1_sq); 
  out[10] += volFact*(2.0*f[48]*wx1_sq+1.154700538379251*f[64]*dv1*wx1+0.1666666666666667*f[48]*dv1_sq); 
  out[11] += volFact*(2.0*f[49]*wx1_sq+1.154700538379251*f[65]*dv1*wx1+0.1666666666666667*f[49]*dv1_sq); 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  double tmp[4]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[12]*dv2; 
  out[0] += ((0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += ((0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += ((0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  double tmp[8]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[16]*dv2; 
  tmp[4] = 2.0*f[11]*wx2+0.5773502691896257*f[25]*dv2; 
  tmp[5] = 2.0*f[12]*wx2+0.5773502691896257*f[26]*dv2; 
  tmp[6] = 2.0*f[19]*wx2+0.5773502691896257*f[35]*dv2; 
  tmp[7] = 2.0*f[20]*wx2+0.5773502691896257*f[36]*dv2; 
  out[0] += ((0.5*Bmag[7]*tmp[7]+0.5*Bmag[6]*tmp[6]+0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.5000000000000001*Bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*Bmag[7]+0.447213595499958*Bmag[3]*tmp[6]+0.447213595499958*tmp[3]*Bmag[6]+0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += ((0.447213595499958*Bmag[3]*tmp[7]+0.447213595499958*tmp[3]*Bmag[7]+0.5000000000000001*Bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*Bmag[6]+0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += ((0.4*Bmag[6]*tmp[7]+0.447213595499958*Bmag[2]*tmp[7]+0.4*tmp[6]*Bmag[7]+0.447213595499958*tmp[2]*Bmag[7]+0.447213595499958*Bmag[1]*tmp[6]+0.447213595499958*tmp[1]*Bmag[6]+0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  out[4] += ((0.4472135954999579*Bmag[7]*tmp[7]+0.31943828249997*Bmag[6]*tmp[6]+0.5000000000000001*Bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*Bmag[6]+0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  out[5] += ((0.31943828249997*Bmag[7]*tmp[7]+0.5000000000000001*Bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*Bmag[7]+0.4472135954999579*Bmag[6]*tmp[6]+0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
  out[6] += ((0.4*Bmag[3]*tmp[7]+0.4*tmp[3]*Bmag[7]+0.4472135954999579*Bmag[5]*tmp[6]+0.31943828249997*Bmag[4]*tmp[6]+0.5*Bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*Bmag[6]+0.31943828249997*tmp[4]*Bmag[6]+0.5*tmp[0]*Bmag[6]+0.5000000000000001*Bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*Bmag[4]+0.447213595499958*Bmag[1]*tmp[3]+0.447213595499958*tmp[1]*Bmag[3])*volFact)/m_; 
  out[7] += ((0.31943828249997*Bmag[5]*tmp[7]+0.4472135954999579*Bmag[4]*tmp[7]+0.5*Bmag[0]*tmp[7]+0.31943828249997*tmp[5]*Bmag[7]+0.4472135954999579*tmp[4]*Bmag[7]+0.5*tmp[0]*Bmag[7]+0.4*Bmag[3]*tmp[6]+0.4*tmp[3]*Bmag[6]+0.5000000000000001*Bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*Bmag[5]+0.447213595499958*Bmag[2]*tmp[3]+0.447213595499958*tmp[2]*Bmag[3])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M2perp_P3(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  double tmp[12]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[16]*dv2; 
  tmp[4] = 2.0*f[11]*wx2+0.5773502691896257*f[25]*dv2; 
  tmp[5] = 2.0*f[12]*wx2+0.5773502691896257*f[26]*dv2; 
  tmp[6] = 2.0*f[19]*wx2+0.5773502691896257*f[39]*dv2; 
  tmp[7] = 2.0*f[20]*wx2+0.5773502691896257*f[40]*dv2; 
  tmp[8] = 2.0*f[31]*wx2+0.5773502691896256*f[54]*dv2; 
  tmp[9] = 2.0*f[32]*wx2+0.5773502691896256*f[55]*dv2; 
  tmp[10] = 2.0*f[48]*wx2+0.5773502691896256*f[67]*dv2; 
  tmp[11] = 2.0*f[49]*wx2+0.5773502691896256*f[68]*dv2; 
  out[0] += ((0.5*Bmag[11]*tmp[11]+0.5*Bmag[10]*tmp[10]+0.5*Bmag[9]*tmp[9]+0.5*Bmag[8]*tmp[8]+0.5*Bmag[7]*tmp[7]+0.5*Bmag[6]*tmp[6]+0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.5*Bmag[9]*tmp[11]+0.5*tmp[9]*Bmag[11]+0.4391550328268399*Bmag[6]*tmp[10]+0.4391550328268399*tmp[6]*Bmag[10]+0.4391550328268398*Bmag[4]*tmp[8]+0.4391550328268398*tmp[4]*Bmag[8]+0.5000000000000001*Bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*Bmag[7]+0.447213595499958*Bmag[3]*tmp[6]+0.447213595499958*tmp[3]*Bmag[6]+0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += ((0.4391550328268399*Bmag[7]*tmp[11]+0.4391550328268399*tmp[7]*Bmag[11]+0.5*Bmag[8]*tmp[10]+0.5*tmp[8]*Bmag[10]+0.4391550328268398*Bmag[5]*tmp[9]+0.4391550328268398*tmp[5]*Bmag[9]+0.447213595499958*Bmag[3]*tmp[7]+0.447213595499958*tmp[3]*Bmag[7]+0.5000000000000001*Bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*Bmag[6]+0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += ((0.4391550328268399*Bmag[5]*tmp[11]+0.4391550328268399*tmp[5]*Bmag[11]+0.4391550328268399*Bmag[4]*tmp[10]+0.4391550328268399*tmp[4]*Bmag[10]+0.4391550328268399*Bmag[7]*tmp[9]+0.4391550328268399*tmp[7]*Bmag[9]+0.4391550328268399*Bmag[6]*tmp[8]+0.4391550328268399*tmp[6]*Bmag[8]+0.4*Bmag[6]*tmp[7]+0.447213595499958*Bmag[2]*tmp[7]+0.4*tmp[6]*Bmag[7]+0.447213595499958*tmp[2]*Bmag[7]+0.447213595499958*Bmag[1]*tmp[6]+0.447213595499958*tmp[1]*Bmag[6]+0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  out[4] += ((0.4472135954999579*Bmag[11]*tmp[11]+0.2981423969999719*Bmag[10]*tmp[10]+0.4391550328268399*Bmag[3]*tmp[10]+0.4391550328268399*tmp[3]*Bmag[10]+0.2981423969999719*Bmag[8]*tmp[8]+0.4391550328268398*Bmag[1]*tmp[8]+0.4391550328268398*tmp[1]*Bmag[8]+0.4472135954999579*Bmag[7]*tmp[7]+0.31943828249997*Bmag[6]*tmp[6]+0.5000000000000001*Bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*Bmag[6]+0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  out[5] += ((0.2981423969999719*Bmag[11]*tmp[11]+0.4391550328268399*Bmag[3]*tmp[11]+0.4391550328268399*tmp[3]*Bmag[11]+0.4472135954999579*Bmag[10]*tmp[10]+0.2981423969999719*Bmag[9]*tmp[9]+0.4391550328268398*Bmag[2]*tmp[9]+0.4391550328268398*tmp[2]*Bmag[9]+0.31943828249997*Bmag[7]*tmp[7]+0.5000000000000001*Bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*Bmag[7]+0.4472135954999579*Bmag[6]*tmp[6]+0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
  out[6] += ((0.3927922024247863*Bmag[7]*tmp[11]+0.3927922024247863*tmp[7]*Bmag[11]+0.2981423969999719*Bmag[8]*tmp[10]+0.3927922024247863*Bmag[7]*tmp[10]+0.4391550328268399*Bmag[1]*tmp[10]+0.2981423969999719*tmp[8]*Bmag[10]+0.3927922024247863*tmp[7]*Bmag[10]+0.4391550328268399*tmp[1]*Bmag[10]+0.4391550328268399*Bmag[3]*tmp[8]+0.4391550328268399*tmp[3]*Bmag[8]+0.4*Bmag[3]*tmp[7]+0.4*tmp[3]*Bmag[7]+0.4472135954999579*Bmag[5]*tmp[6]+0.31943828249997*Bmag[4]*tmp[6]+0.5*Bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*Bmag[6]+0.31943828249997*tmp[4]*Bmag[6]+0.5*tmp[0]*Bmag[6]+0.5000000000000001*Bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*Bmag[4]+0.447213595499958*Bmag[1]*tmp[3]+0.447213595499958*tmp[1]*Bmag[3])*volFact)/m_; 
  out[7] += ((0.2981423969999719*Bmag[9]*tmp[11]+0.3927922024247863*Bmag[6]*tmp[11]+0.4391550328268399*Bmag[2]*tmp[11]+0.2981423969999719*tmp[9]*Bmag[11]+0.3927922024247863*tmp[6]*Bmag[11]+0.4391550328268399*tmp[2]*Bmag[11]+0.3927922024247863*Bmag[6]*tmp[10]+0.3927922024247863*tmp[6]*Bmag[10]+0.4391550328268399*Bmag[3]*tmp[9]+0.4391550328268399*tmp[3]*Bmag[9]+0.31943828249997*Bmag[5]*tmp[7]+0.4472135954999579*Bmag[4]*tmp[7]+0.5*Bmag[0]*tmp[7]+0.31943828249997*tmp[5]*Bmag[7]+0.4472135954999579*tmp[4]*Bmag[7]+0.5*tmp[0]*Bmag[7]+0.4*Bmag[3]*tmp[6]+0.4*tmp[3]*Bmag[6]+0.5000000000000001*Bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*Bmag[5]+0.447213595499958*Bmag[2]*tmp[3]+0.447213595499958*tmp[2]*Bmag[3])*volFact)/m_; 
  out[8] += ((0.2981423969999719*Bmag[6]*tmp[10]+0.5*Bmag[2]*tmp[10]+0.2981423969999719*tmp[6]*Bmag[10]+0.5*tmp[2]*Bmag[10]+0.2981423969999719*Bmag[4]*tmp[8]+0.5*Bmag[0]*tmp[8]+0.2981423969999719*tmp[4]*Bmag[8]+0.5*tmp[0]*Bmag[8]+0.4391550328268399*Bmag[3]*tmp[6]+0.4391550328268399*tmp[3]*Bmag[6]+0.4391550328268398*Bmag[1]*tmp[4]+0.4391550328268398*tmp[1]*Bmag[4])*volFact)/m_; 
  out[9] += ((0.2981423969999719*Bmag[7]*tmp[11]+0.5*Bmag[1]*tmp[11]+0.2981423969999719*tmp[7]*Bmag[11]+0.5*tmp[1]*Bmag[11]+0.2981423969999719*Bmag[5]*tmp[9]+0.5*Bmag[0]*tmp[9]+0.2981423969999719*tmp[5]*Bmag[9]+0.5*tmp[0]*Bmag[9]+0.4391550328268399*Bmag[3]*tmp[7]+0.4391550328268399*tmp[3]*Bmag[7]+0.4391550328268398*Bmag[2]*tmp[5]+0.4391550328268398*tmp[2]*Bmag[5])*volFact)/m_; 
  out[10] += ((0.4472135954999579*Bmag[5]*tmp[10]+0.2981423969999719*Bmag[4]*tmp[10]+0.5*Bmag[0]*tmp[10]+0.4472135954999579*tmp[5]*Bmag[10]+0.2981423969999719*tmp[4]*Bmag[10]+0.5*tmp[0]*Bmag[10]+0.2981423969999719*Bmag[6]*tmp[8]+0.5*Bmag[2]*tmp[8]+0.2981423969999719*tmp[6]*Bmag[8]+0.5*tmp[2]*Bmag[8]+0.3927922024247863*Bmag[6]*tmp[7]+0.3927922024247863*tmp[6]*Bmag[7]+0.4391550328268399*Bmag[1]*tmp[6]+0.4391550328268399*tmp[1]*Bmag[6]+0.4391550328268399*Bmag[3]*tmp[4]+0.4391550328268399*tmp[3]*Bmag[4])*volFact)/m_; 
  out[11] += ((0.2981423969999719*Bmag[5]*tmp[11]+0.4472135954999579*Bmag[4]*tmp[11]+0.5*Bmag[0]*tmp[11]+0.2981423969999719*tmp[5]*Bmag[11]+0.4472135954999579*tmp[4]*Bmag[11]+0.5*tmp[0]*Bmag[11]+0.2981423969999719*Bmag[7]*tmp[9]+0.5*Bmag[1]*tmp[9]+0.2981423969999719*tmp[7]*Bmag[9]+0.5*tmp[1]*Bmag[9]+0.3927922024247863*Bmag[6]*tmp[7]+0.4391550328268399*Bmag[2]*tmp[7]+0.3927922024247863*tmp[6]*Bmag[7]+0.4391550328268399*tmp[2]*Bmag[7]+0.4391550328268399*Bmag[3]*tmp[5]+0.4391550328268399*tmp[3]*Bmag[5])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[3]*dv1*wx1_sq+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[3]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[6]*dv1*wx1_sq+0.5*f[1]*dv1_sq*wx1+0.08660254037844387*f[6]*dv1*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1*wx1_sq+1.732050807568877*f[7]*dv1*wx1_sq+0.5*f[2]*dv1_sq*wx1+0.08660254037844387*f[7]*dv1*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1*wx1_sq+1.732050807568877*f[11]*dv1*wx1_sq+0.5*f[5]*dv1_sq*wx1+0.08660254037844387*f[11]*dv1*dv1_sq); 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[3]*dv1*wx1_sq+0.4472135954999579*f[13]*dv1_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[3]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[6]*dv1*wx1_sq+0.447213595499958*f[23]*dv1_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.08660254037844387*f[6]*dv1*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1*wx1_sq+1.732050807568877*f[7]*dv1*wx1_sq+0.447213595499958*f[24]*dv1_sq*wx1+0.5*f[2]*dv1_sq*wx1+0.08660254037844387*f[7]*dv1*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1*wx1_sq+1.732050807568877*f[15]*dv1*wx1_sq+0.4472135954999579*f[34]*dv1_sq*wx1+0.5*f[5]*dv1_sq*wx1+0.08660254037844387*f[15]*dv1*dv1_sq); 
  out[4] += volFact*(2.0*f[11]*wx1*wx1_sq+1.732050807568877*f[21]*dv1*wx1_sq+0.5*f[11]*dv1_sq*wx1+0.08660254037844385*f[21]*dv1*dv1_sq); 
  out[5] += volFact*(2.0*f[12]*wx1*wx1_sq+1.732050807568877*f[22]*dv1*wx1_sq+0.5*f[12]*dv1_sq*wx1+0.08660254037844385*f[22]*dv1*dv1_sq); 
  out[6] += volFact*(2.0*f[19]*wx1*wx1_sq+1.732050807568877*f[32]*dv1*wx1_sq+0.5*f[19]*dv1_sq*wx1+0.08660254037844385*f[32]*dv1*dv1_sq); 
  out[7] += volFact*(2.0*f[20]*wx1*wx1_sq+1.732050807568877*f[33]*dv1*wx1_sq+0.5*f[20]*dv1_sq*wx1+0.08660254037844385*f[33]*dv1*dv1_sq); 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M3par_P3(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[3]*dv1*wx1_sq+0.4472135954999579*f[13]*dv1_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.03779644730092272*f[33]*dv1*dv1_sq+0.08660254037844387*f[3]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[6]*dv1*wx1_sq+0.447213595499958*f[23]*dv1_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.03779644730092273*f[52]*dv1*dv1_sq+0.08660254037844387*f[6]*dv1*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1*wx1_sq+1.732050807568877*f[7]*dv1*wx1_sq+0.447213595499958*f[24]*dv1_sq*wx1+0.5*f[2]*dv1_sq*wx1+0.03779644730092273*f[53]*dv1*dv1_sq+0.08660254037844387*f[7]*dv1*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1*wx1_sq+1.732050807568877*f[15]*dv1*wx1_sq+0.4472135954999579*f[38]*dv1_sq*wx1+0.5*f[5]*dv1_sq*wx1+0.03779644730092272*f[66]*dv1*dv1_sq+0.08660254037844387*f[15]*dv1*dv1_sq); 
  out[4] += volFact*(2.0*f[11]*wx1*wx1_sq+1.732050807568877*f[21]*dv1*wx1_sq+0.5*f[11]*dv1_sq*wx1+0.08660254037844385*f[21]*dv1*dv1_sq); 
  out[5] += volFact*(2.0*f[12]*wx1*wx1_sq+1.732050807568877*f[22]*dv1*wx1_sq+0.5*f[12]*dv1_sq*wx1+0.08660254037844385*f[22]*dv1*dv1_sq); 
  out[6] += volFact*(2.0*f[19]*wx1*wx1_sq+1.732050807568877*f[36]*dv1*wx1_sq+0.5*f[19]*dv1_sq*wx1+0.08660254037844385*f[36]*dv1*dv1_sq); 
  out[7] += volFact*(2.0*f[20]*wx1*wx1_sq+1.732050807568877*f[37]*dv1*wx1_sq+0.5*f[20]*dv1_sq*wx1+0.08660254037844385*f[37]*dv1*dv1_sq); 
  out[8] += volFact*(2.0*f[31]*wx1*wx1_sq+1.732050807568877*f[50]*dv1*wx1_sq+0.5*f[31]*dv1_sq*wx1+0.08660254037844385*f[50]*dv1*dv1_sq); 
  out[9] += volFact*(2.0*f[32]*wx1*wx1_sq+1.732050807568877*f[51]*dv1*wx1_sq+0.5*f[32]*dv1_sq*wx1+0.08660254037844385*f[51]*dv1*dv1_sq); 
  out[10] += volFact*(2.0*f[48]*wx1*wx1_sq+1.732050807568877*f[64]*dv1*wx1_sq+0.5*f[48]*dv1_sq*wx1+0.08660254037844385*f[64]*dv1*dv1_sq); 
  out[11] += volFact*(2.0*f[49]*wx1*wx1_sq+1.732050807568877*f[65]*dv1*wx1_sq+0.5*f[49]*dv1_sq*wx1+0.08660254037844385*f[65]*dv1*dv1_sq); 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (volFact*(Bmag[3]*f[5]*wx1*wx2+Bmag[2]*f[2]*wx1*wx2+Bmag[1]*f[1]*wx1*wx2+Bmag[0]*f[0]*wx1*wx2+0.2886751345948129*Bmag[3]*f[11]*dv1*wx2+0.2886751345948129*Bmag[2]*f[7]*dv1*wx2+0.2886751345948129*Bmag[1]*f[6]*dv1*wx2+0.2886751345948129*Bmag[0]*f[3]*dv1*wx2+0.2886751345948129*Bmag[3]*f[12]*dv2*wx1+0.2886751345948129*Bmag[2]*f[9]*dv2*wx1+0.2886751345948129*Bmag[1]*f[8]*dv2*wx1+0.2886751345948129*Bmag[0]*f[4]*dv2*wx1+0.08333333333333333*Bmag[3]*f[15]*dv1*dv2+0.08333333333333333*Bmag[2]*f[14]*dv1*dv2+0.08333333333333333*Bmag[1]*f[13]*dv1*dv2+0.08333333333333333*Bmag[0]*f[10]*dv1*dv2))/m_; 
  out[1] += (volFact*(Bmag[2]*f[5]*wx1*wx2+f[2]*Bmag[3]*wx1*wx2+Bmag[0]*f[1]*wx1*wx2+f[0]*Bmag[1]*wx1*wx2+0.2886751345948129*Bmag[2]*f[11]*dv1*wx2+0.2886751345948129*Bmag[3]*f[7]*dv1*wx2+0.2886751345948129*Bmag[0]*f[6]*dv1*wx2+0.2886751345948129*Bmag[1]*f[3]*dv1*wx2+0.2886751345948129*Bmag[2]*f[12]*dv2*wx1+0.2886751345948129*Bmag[3]*f[9]*dv2*wx1+0.2886751345948129*Bmag[0]*f[8]*dv2*wx1+0.2886751345948129*Bmag[1]*f[4]*dv2*wx1+0.08333333333333333*Bmag[2]*f[15]*dv1*dv2+0.08333333333333333*Bmag[3]*f[14]*dv1*dv2+0.08333333333333333*Bmag[0]*f[13]*dv1*dv2+0.08333333333333333*Bmag[1]*f[10]*dv1*dv2))/m_; 
  out[2] += (volFact*(Bmag[1]*f[5]*wx1*wx2+f[1]*Bmag[3]*wx1*wx2+Bmag[0]*f[2]*wx1*wx2+f[0]*Bmag[2]*wx1*wx2+0.2886751345948129*Bmag[1]*f[11]*dv1*wx2+0.2886751345948129*Bmag[0]*f[7]*dv1*wx2+0.2886751345948129*Bmag[3]*f[6]*dv1*wx2+0.2886751345948129*Bmag[2]*f[3]*dv1*wx2+0.2886751345948129*Bmag[1]*f[12]*dv2*wx1+0.2886751345948129*Bmag[0]*f[9]*dv2*wx1+0.2886751345948129*Bmag[3]*f[8]*dv2*wx1+0.2886751345948129*Bmag[2]*f[4]*dv2*wx1+0.08333333333333333*Bmag[1]*f[15]*dv1*dv2+0.08333333333333333*Bmag[0]*f[14]*dv1*dv2+0.08333333333333333*Bmag[3]*f[13]*dv1*dv2+0.08333333333333333*Bmag[2]*f[10]*dv1*dv2))/m_; 
  out[3] += (volFact*(Bmag[0]*f[5]*wx1*wx2+f[0]*Bmag[3]*wx1*wx2+Bmag[1]*f[2]*wx1*wx2+f[1]*Bmag[2]*wx1*wx2+0.2886751345948129*Bmag[0]*f[11]*dv1*wx2+0.2886751345948129*Bmag[1]*f[7]*dv1*wx2+0.2886751345948129*Bmag[2]*f[6]*dv1*wx2+0.2886751345948129*Bmag[3]*f[3]*dv1*wx2+0.2886751345948129*Bmag[0]*f[12]*dv2*wx1+0.2886751345948129*Bmag[1]*f[9]*dv2*wx1+0.2886751345948129*Bmag[2]*f[8]*dv2*wx1+0.2886751345948129*Bmag[3]*f[4]*dv2*wx1+0.08333333333333333*Bmag[0]*f[15]*dv1*dv2+0.08333333333333333*Bmag[1]*f[14]*dv1*dv2+0.08333333333333333*Bmag[2]*f[13]*dv1*dv2+0.08333333333333333*Bmag[3]*f[10]*dv1*dv2))/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (volFact*(Bmag[7]*f[20]*wx1*wx2+Bmag[6]*f[19]*wx1*wx2+Bmag[5]*f[12]*wx1*wx2+Bmag[4]*f[11]*wx1*wx2+Bmag[3]*f[5]*wx1*wx2+Bmag[2]*f[2]*wx1*wx2+Bmag[1]*f[1]*wx1*wx2+Bmag[0]*f[0]*wx1*wx2+0.2886751345948129*Bmag[7]*f[33]*dv1*wx2+0.2886751345948129*Bmag[6]*f[32]*dv1*wx2+0.2886751345948129*Bmag[5]*f[22]*dv1*wx2+0.2886751345948129*Bmag[4]*f[21]*dv1*wx2+0.2886751345948129*Bmag[3]*f[15]*dv1*wx2+0.2886751345948129*Bmag[2]*f[7]*dv1*wx2+0.2886751345948129*Bmag[1]*f[6]*dv1*wx2+0.2886751345948129*Bmag[0]*f[3]*dv1*wx2+0.2886751345948129*Bmag[7]*f[36]*dv2*wx1+0.2886751345948129*Bmag[6]*f[35]*dv2*wx1+0.2886751345948129*Bmag[5]*f[26]*dv2*wx1+0.2886751345948129*Bmag[4]*f[25]*dv2*wx1+0.2886751345948129*Bmag[3]*f[16]*dv2*wx1+0.2886751345948129*Bmag[2]*f[9]*dv2*wx1+0.2886751345948129*Bmag[1]*f[8]*dv2*wx1+0.2886751345948129*Bmag[0]*f[4]*dv2*wx1+0.08333333333333333*Bmag[7]*f[45]*dv1*dv2+0.08333333333333333*Bmag[6]*f[44]*dv1*dv2+0.08333333333333333*Bmag[5]*f[38]*dv1*dv2+0.08333333333333333*Bmag[4]*f[37]*dv1*dv2+0.08333333333333333*Bmag[3]*f[31]*dv1*dv2+0.08333333333333333*Bmag[2]*f[18]*dv1*dv2+0.08333333333333333*Bmag[1]*f[17]*dv1*dv2+0.08333333333333333*Bmag[0]*f[10]*dv1*dv2))/m_; 
  out[1] += (volFact*(1.0*Bmag[5]*f[20]*wx1*wx2+0.8944271909999161*Bmag[3]*f[19]*wx1*wx2+1.0*Bmag[7]*f[12]*wx1*wx2+0.8944271909999159*Bmag[1]*f[11]*wx1*wx2+0.8944271909999161*f[5]*Bmag[6]*wx1*wx2+Bmag[2]*f[5]*wx1*wx2+0.8944271909999159*f[1]*Bmag[4]*wx1*wx2+f[2]*Bmag[3]*wx1*wx2+Bmag[0]*f[1]*wx1*wx2+f[0]*Bmag[1]*wx1*wx2+0.2886751345948129*Bmag[5]*f[33]*dv1*wx2+0.2581988897471612*Bmag[3]*f[32]*dv1*wx2+0.2886751345948129*Bmag[7]*f[22]*dv1*wx2+0.2581988897471611*Bmag[1]*f[21]*dv1*wx2+0.2581988897471611*Bmag[6]*f[15]*dv1*wx2+0.2886751345948129*Bmag[2]*f[15]*dv1*wx2+0.2886751345948129*Bmag[3]*f[7]*dv1*wx2+0.2581988897471612*Bmag[4]*f[6]*dv1*wx2+0.2886751345948129*Bmag[0]*f[6]*dv1*wx2+0.2886751345948129*Bmag[1]*f[3]*dv1*wx2+0.2886751345948129*Bmag[5]*f[36]*dv2*wx1+0.2581988897471612*Bmag[3]*f[35]*dv2*wx1+0.2886751345948129*Bmag[7]*f[26]*dv2*wx1+0.2581988897471611*Bmag[1]*f[25]*dv2*wx1+0.2581988897471611*Bmag[6]*f[16]*dv2*wx1+0.2886751345948129*Bmag[2]*f[16]*dv2*wx1+0.2886751345948129*Bmag[3]*f[9]*dv2*wx1+0.2581988897471612*Bmag[4]*f[8]*dv2*wx1+0.2886751345948129*Bmag[0]*f[8]*dv2*wx1+0.2886751345948129*Bmag[1]*f[4]*dv2*wx1+0.08333333333333336*Bmag[5]*f[45]*dv1*dv2+0.07453559924999302*Bmag[3]*f[44]*dv1*dv2+0.08333333333333336*Bmag[7]*f[38]*dv1*dv2+0.07453559924999298*Bmag[1]*f[37]*dv1*dv2+0.07453559924999302*Bmag[6]*f[31]*dv1*dv2+0.08333333333333333*Bmag[2]*f[31]*dv1*dv2+0.08333333333333333*Bmag[3]*f[18]*dv1*dv2+0.07453559924999298*Bmag[4]*f[17]*dv1*dv2+0.08333333333333333*Bmag[0]*f[17]*dv1*dv2+0.08333333333333333*Bmag[1]*f[10]*dv1*dv2))/m_; 
  out[2] += (volFact*(0.8944271909999161*Bmag[3]*f[20]*wx1*wx2+1.0*Bmag[4]*f[19]*wx1*wx2+0.8944271909999159*Bmag[2]*f[12]*wx1*wx2+1.0*Bmag[6]*f[11]*wx1*wx2+0.8944271909999161*f[5]*Bmag[7]*wx1*wx2+Bmag[1]*f[5]*wx1*wx2+0.8944271909999159*f[2]*Bmag[5]*wx1*wx2+f[1]*Bmag[3]*wx1*wx2+Bmag[0]*f[2]*wx1*wx2+f[0]*Bmag[2]*wx1*wx2+0.2581988897471612*Bmag[3]*f[33]*dv1*wx2+0.2886751345948129*Bmag[4]*f[32]*dv1*wx2+0.2581988897471611*Bmag[2]*f[22]*dv1*wx2+0.2886751345948129*Bmag[6]*f[21]*dv1*wx2+0.2581988897471611*Bmag[7]*f[15]*dv1*wx2+0.2886751345948129*Bmag[1]*f[15]*dv1*wx2+0.2581988897471612*Bmag[5]*f[7]*dv1*wx2+0.2886751345948129*Bmag[0]*f[7]*dv1*wx2+0.2886751345948129*Bmag[3]*f[6]*dv1*wx2+0.2886751345948129*Bmag[2]*f[3]*dv1*wx2+0.2581988897471612*Bmag[3]*f[36]*dv2*wx1+0.2886751345948129*Bmag[4]*f[35]*dv2*wx1+0.2581988897471611*Bmag[2]*f[26]*dv2*wx1+0.2886751345948129*Bmag[6]*f[25]*dv2*wx1+0.2581988897471611*Bmag[7]*f[16]*dv2*wx1+0.2886751345948129*Bmag[1]*f[16]*dv2*wx1+0.2581988897471612*Bmag[5]*f[9]*dv2*wx1+0.2886751345948129*Bmag[0]*f[9]*dv2*wx1+0.2886751345948129*Bmag[3]*f[8]*dv2*wx1+0.2886751345948129*Bmag[2]*f[4]*dv2*wx1+0.07453559924999302*Bmag[3]*f[45]*dv1*dv2+0.08333333333333336*Bmag[4]*f[44]*dv1*dv2+0.07453559924999298*Bmag[2]*f[38]*dv1*dv2+0.08333333333333336*Bmag[6]*f[37]*dv1*dv2+0.07453559924999302*Bmag[7]*f[31]*dv1*dv2+0.08333333333333333*Bmag[1]*f[31]*dv1*dv2+0.07453559924999298*Bmag[5]*f[18]*dv1*dv2+0.08333333333333333*Bmag[0]*f[18]*dv1*dv2+0.08333333333333333*Bmag[3]*f[17]*dv1*dv2+0.08333333333333333*Bmag[2]*f[10]*dv1*dv2))/m_; 
  out[3] += (volFact*(0.8*Bmag[6]*f[20]*wx1*wx2+0.8944271909999161*Bmag[2]*f[20]*wx1*wx2+0.8*Bmag[7]*f[19]*wx1*wx2+0.8944271909999161*Bmag[1]*f[19]*wx1*wx2+0.8944271909999159*Bmag[3]*f[12]*wx1*wx2+0.8944271909999159*Bmag[3]*f[11]*wx1*wx2+0.8944271909999161*f[2]*Bmag[7]*wx1*wx2+0.8944271909999161*f[1]*Bmag[6]*wx1*wx2+0.8944271909999159*Bmag[5]*f[5]*wx1*wx2+0.8944271909999159*Bmag[4]*f[5]*wx1*wx2+Bmag[0]*f[5]*wx1*wx2+f[0]*Bmag[3]*wx1*wx2+Bmag[1]*f[2]*wx1*wx2+f[1]*Bmag[2]*wx1*wx2+0.2309401076758503*Bmag[6]*f[33]*dv1*wx2+0.2581988897471612*Bmag[2]*f[33]*dv1*wx2+0.2309401076758503*Bmag[7]*f[32]*dv1*wx2+0.2581988897471612*Bmag[1]*f[32]*dv1*wx2+0.2581988897471611*Bmag[3]*f[22]*dv1*wx2+0.2581988897471611*Bmag[3]*f[21]*dv1*wx2+0.2581988897471612*Bmag[5]*f[15]*dv1*wx2+0.2581988897471612*Bmag[4]*f[15]*dv1*wx2+0.2886751345948129*Bmag[0]*f[15]*dv1*wx2+0.2581988897471611*Bmag[7]*f[7]*dv1*wx2+0.2886751345948129*Bmag[1]*f[7]*dv1*wx2+0.2581988897471611*Bmag[6]*f[6]*dv1*wx2+0.2886751345948129*Bmag[2]*f[6]*dv1*wx2+0.2886751345948129*Bmag[3]*f[3]*dv1*wx2+0.2309401076758503*Bmag[6]*f[36]*dv2*wx1+0.2581988897471612*Bmag[2]*f[36]*dv2*wx1+0.2309401076758503*Bmag[7]*f[35]*dv2*wx1+0.2581988897471612*Bmag[1]*f[35]*dv2*wx1+0.2581988897471611*Bmag[3]*f[26]*dv2*wx1+0.2581988897471611*Bmag[3]*f[25]*dv2*wx1+0.2581988897471612*Bmag[5]*f[16]*dv2*wx1+0.2581988897471612*Bmag[4]*f[16]*dv2*wx1+0.2886751345948129*Bmag[0]*f[16]*dv2*wx1+0.2581988897471611*Bmag[7]*f[9]*dv2*wx1+0.2886751345948129*Bmag[1]*f[9]*dv2*wx1+0.2581988897471611*Bmag[6]*f[8]*dv2*wx1+0.2886751345948129*Bmag[2]*f[8]*dv2*wx1+0.2886751345948129*Bmag[3]*f[4]*dv2*wx1+0.06666666666666667*Bmag[6]*f[45]*dv1*dv2+0.07453559924999302*Bmag[2]*f[45]*dv1*dv2+0.06666666666666667*Bmag[7]*f[44]*dv1*dv2+0.07453559924999302*Bmag[1]*f[44]*dv1*dv2+0.07453559924999298*Bmag[3]*f[38]*dv1*dv2+0.07453559924999298*Bmag[3]*f[37]*dv1*dv2+0.07453559924999298*Bmag[5]*f[31]*dv1*dv2+0.07453559924999298*Bmag[4]*f[31]*dv1*dv2+0.08333333333333333*Bmag[0]*f[31]*dv1*dv2+0.07453559924999302*Bmag[7]*f[18]*dv1*dv2+0.08333333333333333*Bmag[1]*f[18]*dv1*dv2+0.07453559924999302*Bmag[6]*f[17]*dv1*dv2+0.08333333333333333*Bmag[2]*f[17]*dv1*dv2+0.08333333333333333*Bmag[3]*f[10]*dv1*dv2))/m_; 
  out[4] += (volFact*(0.8944271909999159*Bmag[7]*f[20]*wx1*wx2+0.6388765649999399*Bmag[6]*f[19]*wx1*wx2+1.0*Bmag[2]*f[19]*wx1*wx2+0.6388765649999399*Bmag[4]*f[11]*wx1*wx2+Bmag[0]*f[11]*wx1*wx2+1.0*f[2]*Bmag[6]*wx1*wx2+0.8944271909999159*Bmag[3]*f[5]*wx1*wx2+f[0]*Bmag[4]*wx1*wx2+0.8944271909999159*Bmag[1]*f[1]*wx1*wx2+0.2581988897471611*Bmag[7]*f[33]*dv1*wx2+0.1844277783908294*Bmag[6]*f[32]*dv1*wx2+0.2886751345948129*Bmag[2]*f[32]*dv1*wx2+0.1844277783908294*Bmag[4]*f[21]*dv1*wx2+0.2886751345948129*Bmag[0]*f[21]*dv1*wx2+0.2581988897471612*Bmag[3]*f[15]*dv1*wx2+0.2886751345948129*Bmag[6]*f[7]*dv1*wx2+0.2581988897471612*Bmag[1]*f[6]*dv1*wx2+0.2886751345948129*f[3]*Bmag[4]*dv1*wx2+0.2581988897471611*Bmag[7]*f[36]*dv2*wx1+0.1844277783908294*Bmag[6]*f[35]*dv2*wx1+0.2886751345948129*Bmag[2]*f[35]*dv2*wx1+0.1844277783908294*Bmag[4]*f[25]*dv2*wx1+0.2886751345948129*Bmag[0]*f[25]*dv2*wx1+0.2581988897471612*Bmag[3]*f[16]*dv2*wx1+0.2886751345948129*Bmag[6]*f[9]*dv2*wx1+0.2581988897471612*Bmag[1]*f[8]*dv2*wx1+0.2886751345948129*Bmag[4]*f[4]*dv2*wx1+0.07453559924999298*Bmag[7]*f[45]*dv1*dv2+0.05323971374999499*Bmag[6]*f[44]*dv1*dv2+0.08333333333333336*Bmag[2]*f[44]*dv1*dv2+0.05323971374999499*Bmag[4]*f[37]*dv1*dv2+0.08333333333333333*Bmag[0]*f[37]*dv1*dv2+0.07453559924999298*Bmag[3]*f[31]*dv1*dv2+0.08333333333333336*Bmag[6]*f[18]*dv1*dv2+0.07453559924999298*Bmag[1]*f[17]*dv1*dv2+0.08333333333333333*Bmag[4]*f[10]*dv1*dv2))/m_; 
  out[5] += (volFact*(0.6388765649999399*Bmag[7]*f[20]*wx1*wx2+1.0*Bmag[1]*f[20]*wx1*wx2+0.8944271909999159*Bmag[6]*f[19]*wx1*wx2+0.6388765649999399*Bmag[5]*f[12]*wx1*wx2+Bmag[0]*f[12]*wx1*wx2+1.0*f[1]*Bmag[7]*wx1*wx2+0.8944271909999159*Bmag[3]*f[5]*wx1*wx2+f[0]*Bmag[5]*wx1*wx2+0.8944271909999159*Bmag[2]*f[2]*wx1*wx2+0.1844277783908294*Bmag[7]*f[33]*dv1*wx2+0.2886751345948129*Bmag[1]*f[33]*dv1*wx2+0.2581988897471611*Bmag[6]*f[32]*dv1*wx2+0.1844277783908294*Bmag[5]*f[22]*dv1*wx2+0.2886751345948129*Bmag[0]*f[22]*dv1*wx2+0.2581988897471612*Bmag[3]*f[15]*dv1*wx2+0.2581988897471612*Bmag[2]*f[7]*dv1*wx2+0.2886751345948129*f[6]*Bmag[7]*dv1*wx2+0.2886751345948129*f[3]*Bmag[5]*dv1*wx2+0.1844277783908294*Bmag[7]*f[36]*dv2*wx1+0.2886751345948129*Bmag[1]*f[36]*dv2*wx1+0.2581988897471611*Bmag[6]*f[35]*dv2*wx1+0.1844277783908294*Bmag[5]*f[26]*dv2*wx1+0.2886751345948129*Bmag[0]*f[26]*dv2*wx1+0.2581988897471612*Bmag[3]*f[16]*dv2*wx1+0.2581988897471612*Bmag[2]*f[9]*dv2*wx1+0.2886751345948129*Bmag[7]*f[8]*dv2*wx1+0.2886751345948129*f[4]*Bmag[5]*dv2*wx1+0.05323971374999499*Bmag[7]*f[45]*dv1*dv2+0.08333333333333336*Bmag[1]*f[45]*dv1*dv2+0.07453559924999298*Bmag[6]*f[44]*dv1*dv2+0.05323971374999499*Bmag[5]*f[38]*dv1*dv2+0.08333333333333333*Bmag[0]*f[38]*dv1*dv2+0.07453559924999298*Bmag[3]*f[31]*dv1*dv2+0.07453559924999298*Bmag[2]*f[18]*dv1*dv2+0.08333333333333336*Bmag[7]*f[17]*dv1*dv2+0.08333333333333333*Bmag[5]*f[10]*dv1*dv2))/m_; 
  out[6] += (volFact*(0.8*Bmag[3]*f[20]*wx1*wx2+0.8944271909999159*Bmag[5]*f[19]*wx1*wx2+0.6388765649999399*Bmag[4]*f[19]*wx1*wx2+Bmag[0]*f[19]*wx1*wx2+0.8944271909999159*Bmag[6]*f[12]*wx1*wx2+0.6388765649999399*Bmag[6]*f[11]*wx1*wx2+1.0*Bmag[2]*f[11]*wx1*wx2+0.8*f[5]*Bmag[7]*wx1*wx2+f[0]*Bmag[6]*wx1*wx2+0.8944271909999161*Bmag[1]*f[5]*wx1*wx2+1.0*f[2]*Bmag[4]*wx1*wx2+0.8944271909999161*f[1]*Bmag[3]*wx1*wx2+0.2309401076758503*Bmag[3]*f[33]*dv1*wx2+0.2581988897471611*Bmag[5]*f[32]*dv1*wx2+0.1844277783908294*Bmag[4]*f[32]*dv1*wx2+0.2886751345948129*Bmag[0]*f[32]*dv1*wx2+0.2581988897471611*Bmag[6]*f[22]*dv1*wx2+0.1844277783908294*Bmag[6]*f[21]*dv1*wx2+0.2886751345948129*Bmag[2]*f[21]*dv1*wx2+0.2309401076758504*Bmag[7]*f[15]*dv1*wx2+0.2581988897471611*Bmag[1]*f[15]*dv1*wx2+0.2886751345948129*Bmag[4]*f[7]*dv1*wx2+0.2581988897471611*Bmag[3]*f[6]*dv1*wx2+0.2886751345948129*f[3]*Bmag[6]*dv1*wx2+0.2309401076758503*Bmag[3]*f[36]*dv2*wx1+0.2581988897471611*Bmag[5]*f[35]*dv2*wx1+0.1844277783908294*Bmag[4]*f[35]*dv2*wx1+0.2886751345948129*Bmag[0]*f[35]*dv2*wx1+0.2581988897471611*Bmag[6]*f[26]*dv2*wx1+0.1844277783908294*Bmag[6]*f[25]*dv2*wx1+0.2886751345948129*Bmag[2]*f[25]*dv2*wx1+0.2309401076758504*Bmag[7]*f[16]*dv2*wx1+0.2581988897471611*Bmag[1]*f[16]*dv2*wx1+0.2886751345948129*Bmag[4]*f[9]*dv2*wx1+0.2581988897471611*Bmag[3]*f[8]*dv2*wx1+0.2886751345948129*f[4]*Bmag[6]*dv2*wx1+0.06666666666666667*Bmag[3]*f[45]*dv1*dv2+0.07453559924999298*Bmag[5]*f[44]*dv1*dv2+0.05323971374999499*Bmag[4]*f[44]*dv1*dv2+0.08333333333333333*Bmag[0]*f[44]*dv1*dv2+0.07453559924999298*Bmag[6]*f[38]*dv1*dv2+0.05323971374999499*Bmag[6]*f[37]*dv1*dv2+0.08333333333333336*Bmag[2]*f[37]*dv1*dv2+0.06666666666666667*Bmag[7]*f[31]*dv1*dv2+0.07453559924999302*Bmag[1]*f[31]*dv1*dv2+0.08333333333333336*Bmag[4]*f[18]*dv1*dv2+0.07453559924999302*Bmag[3]*f[17]*dv1*dv2+0.08333333333333333*Bmag[6]*f[10]*dv1*dv2))/m_; 
  out[7] += (volFact*(0.6388765649999399*Bmag[5]*f[20]*wx1*wx2+0.8944271909999159*Bmag[4]*f[20]*wx1*wx2+Bmag[0]*f[20]*wx1*wx2+0.8*Bmag[3]*f[19]*wx1*wx2+0.6388765649999399*Bmag[7]*f[12]*wx1*wx2+1.0*Bmag[1]*f[12]*wx1*wx2+0.8944271909999159*Bmag[7]*f[11]*wx1*wx2+f[0]*Bmag[7]*wx1*wx2+0.8*f[5]*Bmag[6]*wx1*wx2+0.8944271909999161*Bmag[2]*f[5]*wx1*wx2+1.0*f[1]*Bmag[5]*wx1*wx2+0.8944271909999161*f[2]*Bmag[3]*wx1*wx2+0.1844277783908294*Bmag[5]*f[33]*dv1*wx2+0.2581988897471611*Bmag[4]*f[33]*dv1*wx2+0.2886751345948129*Bmag[0]*f[33]*dv1*wx2+0.2309401076758503*Bmag[3]*f[32]*dv1*wx2+0.1844277783908294*Bmag[7]*f[22]*dv1*wx2+0.2886751345948129*Bmag[1]*f[22]*dv1*wx2+0.2581988897471611*Bmag[7]*f[21]*dv1*wx2+0.2309401076758504*Bmag[6]*f[15]*dv1*wx2+0.2581988897471611*Bmag[2]*f[15]*dv1*wx2+0.2581988897471611*Bmag[3]*f[7]*dv1*wx2+0.2886751345948129*f[3]*Bmag[7]*dv1*wx2+0.2886751345948129*Bmag[5]*f[6]*dv1*wx2+0.1844277783908294*Bmag[5]*f[36]*dv2*wx1+0.2581988897471611*Bmag[4]*f[36]*dv2*wx1+0.2886751345948129*Bmag[0]*f[36]*dv2*wx1+0.2309401076758503*Bmag[3]*f[35]*dv2*wx1+0.1844277783908294*Bmag[7]*f[26]*dv2*wx1+0.2886751345948129*Bmag[1]*f[26]*dv2*wx1+0.2581988897471611*Bmag[7]*f[25]*dv2*wx1+0.2309401076758504*Bmag[6]*f[16]*dv2*wx1+0.2581988897471611*Bmag[2]*f[16]*dv2*wx1+0.2581988897471611*Bmag[3]*f[9]*dv2*wx1+0.2886751345948129*Bmag[5]*f[8]*dv2*wx1+0.2886751345948129*f[4]*Bmag[7]*dv2*wx1+0.05323971374999499*Bmag[5]*f[45]*dv1*dv2+0.07453559924999298*Bmag[4]*f[45]*dv1*dv2+0.08333333333333333*Bmag[0]*f[45]*dv1*dv2+0.06666666666666667*Bmag[3]*f[44]*dv1*dv2+0.05323971374999499*Bmag[7]*f[38]*dv1*dv2+0.08333333333333336*Bmag[1]*f[38]*dv1*dv2+0.07453559924999298*Bmag[7]*f[37]*dv1*dv2+0.06666666666666667*Bmag[6]*f[31]*dv1*dv2+0.07453559924999302*Bmag[2]*f[31]*dv1*dv2+0.07453559924999302*Bmag[3]*f[18]*dv1*dv2+0.08333333333333336*Bmag[5]*f[17]*dv1*dv2+0.08333333333333333*Bmag[7]*f[10]*dv1*dv2))/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M3perp_P3(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (volFact*(Bmag[11]*f[49]*wx1*wx2+Bmag[10]*f[48]*wx1*wx2+Bmag[9]*f[32]*wx1*wx2+Bmag[8]*f[31]*wx1*wx2+Bmag[7]*f[20]*wx1*wx2+Bmag[6]*f[19]*wx1*wx2+Bmag[5]*f[12]*wx1*wx2+Bmag[4]*f[11]*wx1*wx2+Bmag[3]*f[5]*wx1*wx2+Bmag[2]*f[2]*wx1*wx2+Bmag[1]*f[1]*wx1*wx2+Bmag[0]*f[0]*wx1*wx2+0.2886751345948128*Bmag[11]*f[65]*dv1*wx2+0.2886751345948128*Bmag[10]*f[64]*dv1*wx2+0.2886751345948128*Bmag[9]*f[51]*dv1*wx2+0.2886751345948128*Bmag[8]*f[50]*dv1*wx2+0.2886751345948129*Bmag[7]*f[37]*dv1*wx2+0.2886751345948129*Bmag[6]*f[36]*dv1*wx2+0.2886751345948129*Bmag[5]*f[22]*dv1*wx2+0.2886751345948129*Bmag[4]*f[21]*dv1*wx2+0.2886751345948129*Bmag[3]*f[15]*dv1*wx2+0.2886751345948129*Bmag[2]*f[7]*dv1*wx2+0.2886751345948129*Bmag[1]*f[6]*dv1*wx2+0.2886751345948129*Bmag[0]*f[3]*dv1*wx2+0.2886751345948128*Bmag[11]*f[68]*dv2*wx1+0.2886751345948128*Bmag[10]*f[67]*dv2*wx1+0.2886751345948128*Bmag[9]*f[55]*dv2*wx1+0.2886751345948128*Bmag[8]*f[54]*dv2*wx1+0.2886751345948129*Bmag[7]*f[40]*dv2*wx1+0.2886751345948129*Bmag[6]*f[39]*dv2*wx1+0.2886751345948129*Bmag[5]*f[26]*dv2*wx1+0.2886751345948129*Bmag[4]*f[25]*dv2*wx1+0.2886751345948129*Bmag[3]*f[16]*dv2*wx1+0.2886751345948129*Bmag[2]*f[9]*dv2*wx1+0.2886751345948129*Bmag[1]*f[8]*dv2*wx1+0.2886751345948129*Bmag[0]*f[4]*dv2*wx1+0.08333333333333333*Bmag[11]*f[77]*dv1*dv2+0.08333333333333333*Bmag[10]*f[76]*dv1*dv2+0.08333333333333333*Bmag[9]*f[70]*dv1*dv2+0.08333333333333333*Bmag[8]*f[69]*dv1*dv2+0.08333333333333333*Bmag[7]*f[61]*dv1*dv2+0.08333333333333333*Bmag[6]*f[60]*dv1*dv2+0.08333333333333333*Bmag[5]*f[42]*dv1*dv2+0.08333333333333333*Bmag[4]*f[41]*dv1*dv2+0.08333333333333333*Bmag[3]*f[35]*dv1*dv2+0.08333333333333333*Bmag[2]*f[18]*dv1*dv2+0.08333333333333333*Bmag[1]*f[17]*dv1*dv2+0.08333333333333333*Bmag[0]*f[10]*dv1*dv2))/m_; 
  out[1] += (volFact*(1.0*Bmag[9]*f[49]*wx1*wx2+0.8783100656536798*Bmag[6]*f[48]*wx1*wx2+1.0*Bmag[11]*f[32]*wx1*wx2+0.8783100656536796*Bmag[4]*f[31]*wx1*wx2+1.0*Bmag[5]*f[20]*wx1*wx2+0.8783100656536798*Bmag[10]*f[19]*wx1*wx2+0.8944271909999161*Bmag[3]*f[19]*wx1*wx2+1.0*Bmag[7]*f[12]*wx1*wx2+0.8783100656536796*Bmag[8]*f[11]*wx1*wx2+0.8944271909999159*Bmag[1]*f[11]*wx1*wx2+0.8944271909999161*f[5]*Bmag[6]*wx1*wx2+Bmag[2]*f[5]*wx1*wx2+0.8944271909999159*f[1]*Bmag[4]*wx1*wx2+f[2]*Bmag[3]*wx1*wx2+Bmag[0]*f[1]*wx1*wx2+f[0]*Bmag[1]*wx1*wx2+0.2886751345948129*Bmag[9]*f[65]*dv1*wx2+0.253546276418555*Bmag[6]*f[64]*dv1*wx2+0.2886751345948129*Bmag[11]*f[51]*dv1*wx2+0.2535462764185549*Bmag[4]*f[50]*dv1*wx2+0.2886751345948129*Bmag[5]*f[37]*dv1*wx2+0.2535462764185549*Bmag[10]*f[36]*dv1*wx2+0.2581988897471612*Bmag[3]*f[36]*dv1*wx2+0.2886751345948129*Bmag[7]*f[22]*dv1*wx2+0.253546276418555*Bmag[8]*f[21]*dv1*wx2+0.2581988897471611*Bmag[1]*f[21]*dv1*wx2+0.2581988897471611*Bmag[6]*f[15]*dv1*wx2+0.2886751345948129*Bmag[2]*f[15]*dv1*wx2+0.2886751345948129*Bmag[3]*f[7]*dv1*wx2+0.2581988897471612*Bmag[4]*f[6]*dv1*wx2+0.2886751345948129*Bmag[0]*f[6]*dv1*wx2+0.2886751345948129*Bmag[1]*f[3]*dv1*wx2+0.2886751345948129*Bmag[9]*f[68]*dv2*wx1+0.253546276418555*Bmag[6]*f[67]*dv2*wx1+0.2886751345948129*Bmag[11]*f[55]*dv2*wx1+0.2535462764185549*Bmag[4]*f[54]*dv2*wx1+0.2886751345948129*Bmag[5]*f[40]*dv2*wx1+0.2535462764185549*Bmag[10]*f[39]*dv2*wx1+0.2581988897471612*Bmag[3]*f[39]*dv2*wx1+0.2886751345948129*Bmag[7]*f[26]*dv2*wx1+0.253546276418555*Bmag[8]*f[25]*dv2*wx1+0.2581988897471611*Bmag[1]*f[25]*dv2*wx1+0.2581988897471611*Bmag[6]*f[16]*dv2*wx1+0.2886751345948129*Bmag[2]*f[16]*dv2*wx1+0.2886751345948129*Bmag[3]*f[9]*dv2*wx1+0.2581988897471612*Bmag[4]*f[8]*dv2*wx1+0.2886751345948129*Bmag[0]*f[8]*dv2*wx1+0.2886751345948129*Bmag[1]*f[4]*dv2*wx1+0.08333333333333334*Bmag[9]*f[77]*dv1*dv2+0.07319250547114*Bmag[6]*f[76]*dv1*dv2+0.08333333333333334*Bmag[11]*f[70]*dv1*dv2+0.07319250547113998*Bmag[4]*f[69]*dv1*dv2+0.08333333333333336*Bmag[5]*f[61]*dv1*dv2+0.07319250547114*Bmag[10]*f[60]*dv1*dv2+0.07453559924999302*Bmag[3]*f[60]*dv1*dv2+0.08333333333333336*Bmag[7]*f[42]*dv1*dv2+0.07319250547113998*Bmag[8]*f[41]*dv1*dv2+0.07453559924999298*Bmag[1]*f[41]*dv1*dv2+0.07453559924999302*Bmag[6]*f[35]*dv1*dv2+0.08333333333333333*Bmag[2]*f[35]*dv1*dv2+0.08333333333333333*Bmag[3]*f[18]*dv1*dv2+0.07453559924999298*Bmag[4]*f[17]*dv1*dv2+0.08333333333333333*Bmag[0]*f[17]*dv1*dv2+0.08333333333333333*Bmag[1]*f[10]*dv1*dv2))/m_; 
  out[2] += (volFact*(0.8783100656536798*Bmag[7]*f[49]*wx1*wx2+1.0*Bmag[8]*f[48]*wx1*wx2+0.8783100656536796*Bmag[5]*f[32]*wx1*wx2+1.0*Bmag[10]*f[31]*wx1*wx2+0.8783100656536798*Bmag[11]*f[20]*wx1*wx2+0.8944271909999161*Bmag[3]*f[20]*wx1*wx2+1.0*Bmag[4]*f[19]*wx1*wx2+0.8783100656536796*Bmag[9]*f[12]*wx1*wx2+0.8944271909999159*Bmag[2]*f[12]*wx1*wx2+1.0*Bmag[6]*f[11]*wx1*wx2+0.8944271909999161*f[5]*Bmag[7]*wx1*wx2+Bmag[1]*f[5]*wx1*wx2+0.8944271909999159*f[2]*Bmag[5]*wx1*wx2+f[1]*Bmag[3]*wx1*wx2+Bmag[0]*f[2]*wx1*wx2+f[0]*Bmag[2]*wx1*wx2+0.253546276418555*Bmag[7]*f[65]*dv1*wx2+0.2886751345948129*Bmag[8]*f[64]*dv1*wx2+0.2535462764185549*Bmag[5]*f[51]*dv1*wx2+0.2886751345948129*Bmag[10]*f[50]*dv1*wx2+0.2535462764185549*Bmag[11]*f[37]*dv1*wx2+0.2581988897471612*Bmag[3]*f[37]*dv1*wx2+0.2886751345948129*Bmag[4]*f[36]*dv1*wx2+0.253546276418555*Bmag[9]*f[22]*dv1*wx2+0.2581988897471611*Bmag[2]*f[22]*dv1*wx2+0.2886751345948129*Bmag[6]*f[21]*dv1*wx2+0.2581988897471611*Bmag[7]*f[15]*dv1*wx2+0.2886751345948129*Bmag[1]*f[15]*dv1*wx2+0.2581988897471612*Bmag[5]*f[7]*dv1*wx2+0.2886751345948129*Bmag[0]*f[7]*dv1*wx2+0.2886751345948129*Bmag[3]*f[6]*dv1*wx2+0.2886751345948129*Bmag[2]*f[3]*dv1*wx2+0.253546276418555*Bmag[7]*f[68]*dv2*wx1+0.2886751345948129*Bmag[8]*f[67]*dv2*wx1+0.2535462764185549*Bmag[5]*f[55]*dv2*wx1+0.2886751345948129*Bmag[10]*f[54]*dv2*wx1+0.2535462764185549*Bmag[11]*f[40]*dv2*wx1+0.2581988897471612*Bmag[3]*f[40]*dv2*wx1+0.2886751345948129*Bmag[4]*f[39]*dv2*wx1+0.253546276418555*Bmag[9]*f[26]*dv2*wx1+0.2581988897471611*Bmag[2]*f[26]*dv2*wx1+0.2886751345948129*Bmag[6]*f[25]*dv2*wx1+0.2581988897471611*Bmag[7]*f[16]*dv2*wx1+0.2886751345948129*Bmag[1]*f[16]*dv2*wx1+0.2581988897471612*Bmag[5]*f[9]*dv2*wx1+0.2886751345948129*Bmag[0]*f[9]*dv2*wx1+0.2886751345948129*Bmag[3]*f[8]*dv2*wx1+0.2886751345948129*Bmag[2]*f[4]*dv2*wx1+0.07319250547114*Bmag[7]*f[77]*dv1*dv2+0.08333333333333334*Bmag[8]*f[76]*dv1*dv2+0.07319250547113998*Bmag[5]*f[70]*dv1*dv2+0.08333333333333334*Bmag[10]*f[69]*dv1*dv2+0.07319250547114*Bmag[11]*f[61]*dv1*dv2+0.07453559924999302*Bmag[3]*f[61]*dv1*dv2+0.08333333333333336*Bmag[4]*f[60]*dv1*dv2+0.07319250547113998*Bmag[9]*f[42]*dv1*dv2+0.07453559924999298*Bmag[2]*f[42]*dv1*dv2+0.08333333333333336*Bmag[6]*f[41]*dv1*dv2+0.07453559924999302*Bmag[7]*f[35]*dv1*dv2+0.08333333333333333*Bmag[1]*f[35]*dv1*dv2+0.07453559924999298*Bmag[5]*f[18]*dv1*dv2+0.08333333333333333*Bmag[0]*f[18]*dv1*dv2+0.08333333333333333*Bmag[3]*f[17]*dv1*dv2+0.08333333333333333*Bmag[2]*f[10]*dv1*dv2))/m_; 
  out[3] += (volFact*(0.8783100656536798*Bmag[5]*f[49]*wx1*wx2+0.8783100656536798*Bmag[4]*f[48]*wx1*wx2+0.87831006565368*Bmag[7]*f[32]*wx1*wx2+0.87831006565368*Bmag[6]*f[31]*wx1*wx2+0.87831006565368*Bmag[9]*f[20]*wx1*wx2+0.8*Bmag[6]*f[20]*wx1*wx2+0.8944271909999161*Bmag[2]*f[20]*wx1*wx2+0.87831006565368*Bmag[8]*f[19]*wx1*wx2+0.8*Bmag[7]*f[19]*wx1*wx2+0.8944271909999161*Bmag[1]*f[19]*wx1*wx2+0.8783100656536798*Bmag[11]*f[12]*wx1*wx2+0.8944271909999159*Bmag[3]*f[12]*wx1*wx2+0.8783100656536798*Bmag[10]*f[11]*wx1*wx2+0.8944271909999159*Bmag[3]*f[11]*wx1*wx2+0.8944271909999161*f[2]*Bmag[7]*wx1*wx2+0.8944271909999161*f[1]*Bmag[6]*wx1*wx2+0.8944271909999159*Bmag[5]*f[5]*wx1*wx2+0.8944271909999159*Bmag[4]*f[5]*wx1*wx2+Bmag[0]*f[5]*wx1*wx2+f[0]*Bmag[3]*wx1*wx2+Bmag[1]*f[2]*wx1*wx2+f[1]*Bmag[2]*wx1*wx2+0.253546276418555*Bmag[5]*f[65]*dv1*wx2+0.253546276418555*Bmag[4]*f[64]*dv1*wx2+0.253546276418555*Bmag[7]*f[51]*dv1*wx2+0.253546276418555*Bmag[6]*f[50]*dv1*wx2+0.253546276418555*Bmag[9]*f[37]*dv1*wx2+0.2309401076758503*Bmag[6]*f[37]*dv1*wx2+0.2581988897471612*Bmag[2]*f[37]*dv1*wx2+0.253546276418555*Bmag[8]*f[36]*dv1*wx2+0.2309401076758503*Bmag[7]*f[36]*dv1*wx2+0.2581988897471612*Bmag[1]*f[36]*dv1*wx2+0.253546276418555*Bmag[11]*f[22]*dv1*wx2+0.2581988897471611*Bmag[3]*f[22]*dv1*wx2+0.253546276418555*Bmag[10]*f[21]*dv1*wx2+0.2581988897471611*Bmag[3]*f[21]*dv1*wx2+0.2581988897471612*Bmag[5]*f[15]*dv1*wx2+0.2581988897471612*Bmag[4]*f[15]*dv1*wx2+0.2886751345948129*Bmag[0]*f[15]*dv1*wx2+0.2581988897471611*Bmag[7]*f[7]*dv1*wx2+0.2886751345948129*Bmag[1]*f[7]*dv1*wx2+0.2581988897471611*Bmag[6]*f[6]*dv1*wx2+0.2886751345948129*Bmag[2]*f[6]*dv1*wx2+0.2886751345948129*Bmag[3]*f[3]*dv1*wx2+0.253546276418555*Bmag[5]*f[68]*dv2*wx1+0.253546276418555*Bmag[4]*f[67]*dv2*wx1+0.253546276418555*Bmag[7]*f[55]*dv2*wx1+0.253546276418555*Bmag[6]*f[54]*dv2*wx1+0.253546276418555*Bmag[9]*f[40]*dv2*wx1+0.2309401076758503*Bmag[6]*f[40]*dv2*wx1+0.2581988897471612*Bmag[2]*f[40]*dv2*wx1+0.253546276418555*Bmag[8]*f[39]*dv2*wx1+0.2309401076758503*Bmag[7]*f[39]*dv2*wx1+0.2581988897471612*Bmag[1]*f[39]*dv2*wx1+0.253546276418555*Bmag[11]*f[26]*dv2*wx1+0.2581988897471611*Bmag[3]*f[26]*dv2*wx1+0.253546276418555*Bmag[10]*f[25]*dv2*wx1+0.2581988897471611*Bmag[3]*f[25]*dv2*wx1+0.2581988897471612*Bmag[5]*f[16]*dv2*wx1+0.2581988897471612*Bmag[4]*f[16]*dv2*wx1+0.2886751345948129*Bmag[0]*f[16]*dv2*wx1+0.2581988897471611*Bmag[7]*f[9]*dv2*wx1+0.2886751345948129*Bmag[1]*f[9]*dv2*wx1+0.2581988897471611*Bmag[6]*f[8]*dv2*wx1+0.2886751345948129*Bmag[2]*f[8]*dv2*wx1+0.2886751345948129*Bmag[3]*f[4]*dv2*wx1+0.07319250547113998*Bmag[5]*f[77]*dv1*dv2+0.07319250547113998*Bmag[4]*f[76]*dv1*dv2+0.07319250547113999*Bmag[7]*f[70]*dv1*dv2+0.07319250547113999*Bmag[6]*f[69]*dv1*dv2+0.07319250547113999*Bmag[9]*f[61]*dv1*dv2+0.06666666666666667*Bmag[6]*f[61]*dv1*dv2+0.07453559924999302*Bmag[2]*f[61]*dv1*dv2+0.07319250547113999*Bmag[8]*f[60]*dv1*dv2+0.06666666666666667*Bmag[7]*f[60]*dv1*dv2+0.07453559924999302*Bmag[1]*f[60]*dv1*dv2+0.07319250547113998*Bmag[11]*f[42]*dv1*dv2+0.07453559924999298*Bmag[3]*f[42]*dv1*dv2+0.07319250547113998*Bmag[10]*f[41]*dv1*dv2+0.07453559924999298*Bmag[3]*f[41]*dv1*dv2+0.07453559924999298*Bmag[5]*f[35]*dv1*dv2+0.07453559924999298*Bmag[4]*f[35]*dv1*dv2+0.08333333333333333*Bmag[0]*f[35]*dv1*dv2+0.07453559924999302*Bmag[7]*f[18]*dv1*dv2+0.08333333333333333*Bmag[1]*f[18]*dv1*dv2+0.07453559924999302*Bmag[6]*f[17]*dv1*dv2+0.08333333333333333*Bmag[2]*f[17]*dv1*dv2+0.08333333333333333*Bmag[3]*f[10]*dv1*dv2))/m_; 
  out[4] += (volFact*(0.8944271909999159*Bmag[11]*f[49]*wx1*wx2+0.5962847939999438*Bmag[10]*f[48]*wx1*wx2+0.8783100656536798*Bmag[3]*f[48]*wx1*wx2+0.5962847939999438*Bmag[8]*f[31]*wx1*wx2+0.8783100656536796*Bmag[1]*f[31]*wx1*wx2+0.8944271909999159*Bmag[7]*f[20]*wx1*wx2+0.6388765649999399*Bmag[6]*f[19]*wx1*wx2+1.0*Bmag[2]*f[19]*wx1*wx2+0.6388765649999399*Bmag[4]*f[11]*wx1*wx2+Bmag[0]*f[11]*wx1*wx2+0.8783100656536798*f[5]*Bmag[10]*wx1*wx2+0.8783100656536796*f[1]*Bmag[8]*wx1*wx2+1.0*f[2]*Bmag[6]*wx1*wx2+0.8944271909999159*Bmag[3]*f[5]*wx1*wx2+f[0]*Bmag[4]*wx1*wx2+0.8944271909999159*Bmag[1]*f[1]*wx1*wx2+0.258198889747161*Bmag[11]*f[65]*dv1*wx2+0.172132593164774*Bmag[10]*f[64]*dv1*wx2+0.253546276418555*Bmag[3]*f[64]*dv1*wx2+0.172132593164774*Bmag[8]*f[50]*dv1*wx2+0.2535462764185549*Bmag[1]*f[50]*dv1*wx2+0.2581988897471611*Bmag[7]*f[37]*dv1*wx2+0.1844277783908294*Bmag[6]*f[36]*dv1*wx2+0.2886751345948129*Bmag[2]*f[36]*dv1*wx2+0.1844277783908294*Bmag[4]*f[21]*dv1*wx2+0.2886751345948129*Bmag[0]*f[21]*dv1*wx2+0.2535462764185549*Bmag[10]*f[15]*dv1*wx2+0.2581988897471612*Bmag[3]*f[15]*dv1*wx2+0.253546276418555*f[6]*Bmag[8]*dv1*wx2+0.2886751345948129*Bmag[6]*f[7]*dv1*wx2+0.2581988897471612*Bmag[1]*f[6]*dv1*wx2+0.2886751345948129*f[3]*Bmag[4]*dv1*wx2+0.258198889747161*Bmag[11]*f[68]*dv2*wx1+0.172132593164774*Bmag[10]*f[67]*dv2*wx1+0.253546276418555*Bmag[3]*f[67]*dv2*wx1+0.172132593164774*Bmag[8]*f[54]*dv2*wx1+0.2535462764185549*Bmag[1]*f[54]*dv2*wx1+0.2581988897471611*Bmag[7]*f[40]*dv2*wx1+0.1844277783908294*Bmag[6]*f[39]*dv2*wx1+0.2886751345948129*Bmag[2]*f[39]*dv2*wx1+0.1844277783908294*Bmag[4]*f[25]*dv2*wx1+0.2886751345948129*Bmag[0]*f[25]*dv2*wx1+0.2535462764185549*Bmag[10]*f[16]*dv2*wx1+0.2581988897471612*Bmag[3]*f[16]*dv2*wx1+0.2886751345948129*Bmag[6]*f[9]*dv2*wx1+0.253546276418555*Bmag[8]*f[8]*dv2*wx1+0.2581988897471612*Bmag[1]*f[8]*dv2*wx1+0.2886751345948129*Bmag[4]*f[4]*dv2*wx1+0.07453559924999298*Bmag[11]*f[77]*dv1*dv2+0.04969039949999532*Bmag[10]*f[76]*dv1*dv2+0.07319250547113998*Bmag[3]*f[76]*dv1*dv2+0.04969039949999532*Bmag[8]*f[69]*dv1*dv2+0.07319250547113998*Bmag[1]*f[69]*dv1*dv2+0.07453559924999298*Bmag[7]*f[61]*dv1*dv2+0.05323971374999499*Bmag[6]*f[60]*dv1*dv2+0.08333333333333336*Bmag[2]*f[60]*dv1*dv2+0.05323971374999499*Bmag[4]*f[41]*dv1*dv2+0.08333333333333333*Bmag[0]*f[41]*dv1*dv2+0.07319250547113998*Bmag[10]*f[35]*dv1*dv2+0.07453559924999298*Bmag[3]*f[35]*dv1*dv2+0.08333333333333336*Bmag[6]*f[18]*dv1*dv2+0.07319250547113998*Bmag[8]*f[17]*dv1*dv2+0.07453559924999298*Bmag[1]*f[17]*dv1*dv2+0.08333333333333333*Bmag[4]*f[10]*dv1*dv2))/m_; 
  out[5] += (volFact*(0.5962847939999438*Bmag[11]*f[49]*wx1*wx2+0.8783100656536798*Bmag[3]*f[49]*wx1*wx2+0.8944271909999159*Bmag[10]*f[48]*wx1*wx2+0.5962847939999438*Bmag[9]*f[32]*wx1*wx2+0.8783100656536796*Bmag[2]*f[32]*wx1*wx2+0.6388765649999399*Bmag[7]*f[20]*wx1*wx2+1.0*Bmag[1]*f[20]*wx1*wx2+0.8944271909999159*Bmag[6]*f[19]*wx1*wx2+0.6388765649999399*Bmag[5]*f[12]*wx1*wx2+Bmag[0]*f[12]*wx1*wx2+0.8783100656536798*f[5]*Bmag[11]*wx1*wx2+0.8783100656536796*f[2]*Bmag[9]*wx1*wx2+1.0*f[1]*Bmag[7]*wx1*wx2+0.8944271909999159*Bmag[3]*f[5]*wx1*wx2+f[0]*Bmag[5]*wx1*wx2+0.8944271909999159*Bmag[2]*f[2]*wx1*wx2+0.172132593164774*Bmag[11]*f[65]*dv1*wx2+0.253546276418555*Bmag[3]*f[65]*dv1*wx2+0.258198889747161*Bmag[10]*f[64]*dv1*wx2+0.172132593164774*Bmag[9]*f[51]*dv1*wx2+0.2535462764185549*Bmag[2]*f[51]*dv1*wx2+0.1844277783908294*Bmag[7]*f[37]*dv1*wx2+0.2886751345948129*Bmag[1]*f[37]*dv1*wx2+0.2581988897471611*Bmag[6]*f[36]*dv1*wx2+0.1844277783908294*Bmag[5]*f[22]*dv1*wx2+0.2886751345948129*Bmag[0]*f[22]*dv1*wx2+0.2535462764185549*Bmag[11]*f[15]*dv1*wx2+0.2581988897471612*Bmag[3]*f[15]*dv1*wx2+0.253546276418555*f[7]*Bmag[9]*dv1*wx2+0.2581988897471612*Bmag[2]*f[7]*dv1*wx2+0.2886751345948129*f[6]*Bmag[7]*dv1*wx2+0.2886751345948129*f[3]*Bmag[5]*dv1*wx2+0.172132593164774*Bmag[11]*f[68]*dv2*wx1+0.253546276418555*Bmag[3]*f[68]*dv2*wx1+0.258198889747161*Bmag[10]*f[67]*dv2*wx1+0.172132593164774*Bmag[9]*f[55]*dv2*wx1+0.2535462764185549*Bmag[2]*f[55]*dv2*wx1+0.1844277783908294*Bmag[7]*f[40]*dv2*wx1+0.2886751345948129*Bmag[1]*f[40]*dv2*wx1+0.2581988897471611*Bmag[6]*f[39]*dv2*wx1+0.1844277783908294*Bmag[5]*f[26]*dv2*wx1+0.2886751345948129*Bmag[0]*f[26]*dv2*wx1+0.2535462764185549*Bmag[11]*f[16]*dv2*wx1+0.2581988897471612*Bmag[3]*f[16]*dv2*wx1+0.253546276418555*Bmag[9]*f[9]*dv2*wx1+0.2581988897471612*Bmag[2]*f[9]*dv2*wx1+0.2886751345948129*Bmag[7]*f[8]*dv2*wx1+0.2886751345948129*f[4]*Bmag[5]*dv2*wx1+0.04969039949999532*Bmag[11]*f[77]*dv1*dv2+0.07319250547113998*Bmag[3]*f[77]*dv1*dv2+0.07453559924999298*Bmag[10]*f[76]*dv1*dv2+0.04969039949999532*Bmag[9]*f[70]*dv1*dv2+0.07319250547113998*Bmag[2]*f[70]*dv1*dv2+0.05323971374999499*Bmag[7]*f[61]*dv1*dv2+0.08333333333333336*Bmag[1]*f[61]*dv1*dv2+0.07453559924999298*Bmag[6]*f[60]*dv1*dv2+0.05323971374999499*Bmag[5]*f[42]*dv1*dv2+0.08333333333333333*Bmag[0]*f[42]*dv1*dv2+0.07319250547113998*Bmag[11]*f[35]*dv1*dv2+0.07453559924999298*Bmag[3]*f[35]*dv1*dv2+0.07319250547113998*Bmag[9]*f[18]*dv1*dv2+0.07453559924999298*Bmag[2]*f[18]*dv1*dv2+0.08333333333333336*Bmag[7]*f[17]*dv1*dv2+0.08333333333333333*Bmag[5]*f[10]*dv1*dv2))/m_; 
  out[6] += (volFact*(0.7855844048495726*Bmag[7]*f[49]*wx1*wx2+0.5962847939999437*Bmag[8]*f[48]*wx1*wx2+0.7855844048495726*Bmag[7]*f[48]*wx1*wx2+0.8783100656536798*Bmag[1]*f[48]*wx1*wx2+0.5962847939999437*Bmag[10]*f[31]*wx1*wx2+0.87831006565368*Bmag[3]*f[31]*wx1*wx2+0.7855844048495726*Bmag[11]*f[20]*wx1*wx2+0.7855844048495726*Bmag[10]*f[20]*wx1*wx2+0.8*Bmag[3]*f[20]*wx1*wx2+0.8944271909999159*Bmag[5]*f[19]*wx1*wx2+0.6388765649999399*Bmag[4]*f[19]*wx1*wx2+Bmag[0]*f[19]*wx1*wx2+0.8944271909999159*Bmag[6]*f[12]*wx1*wx2+0.6388765649999399*Bmag[6]*f[11]*wx1*wx2+1.0*Bmag[2]*f[11]*wx1*wx2+0.8783100656536798*f[1]*Bmag[10]*wx1*wx2+0.87831006565368*f[5]*Bmag[8]*wx1*wx2+0.8*f[5]*Bmag[7]*wx1*wx2+f[0]*Bmag[6]*wx1*wx2+0.8944271909999161*Bmag[1]*f[5]*wx1*wx2+1.0*f[2]*Bmag[4]*wx1*wx2+0.8944271909999161*f[1]*Bmag[3]*wx1*wx2+0.2267786838055363*Bmag[7]*f[65]*dv1*wx2+0.1721325931647741*Bmag[8]*f[64]*dv1*wx2+0.2267786838055363*Bmag[7]*f[64]*dv1*wx2+0.253546276418555*Bmag[1]*f[64]*dv1*wx2+0.1721325931647741*Bmag[10]*f[50]*dv1*wx2+0.253546276418555*Bmag[3]*f[50]*dv1*wx2+0.2267786838055363*Bmag[11]*f[37]*dv1*wx2+0.2267786838055363*Bmag[10]*f[37]*dv1*wx2+0.2309401076758503*Bmag[3]*f[37]*dv1*wx2+0.2581988897471611*Bmag[5]*f[36]*dv1*wx2+0.1844277783908294*Bmag[4]*f[36]*dv1*wx2+0.2886751345948129*Bmag[0]*f[36]*dv1*wx2+0.2581988897471611*Bmag[6]*f[22]*dv1*wx2+0.1844277783908294*Bmag[6]*f[21]*dv1*wx2+0.2886751345948129*Bmag[2]*f[21]*dv1*wx2+0.253546276418555*Bmag[8]*f[15]*dv1*wx2+0.2309401076758504*Bmag[7]*f[15]*dv1*wx2+0.2581988897471611*Bmag[1]*f[15]*dv1*wx2+0.253546276418555*f[6]*Bmag[10]*dv1*wx2+0.2886751345948129*Bmag[4]*f[7]*dv1*wx2+0.2581988897471611*Bmag[3]*f[6]*dv1*wx2+0.2886751345948129*f[3]*Bmag[6]*dv1*wx2+0.2267786838055363*Bmag[7]*f[68]*dv2*wx1+0.1721325931647741*Bmag[8]*f[67]*dv2*wx1+0.2267786838055363*Bmag[7]*f[67]*dv2*wx1+0.253546276418555*Bmag[1]*f[67]*dv2*wx1+0.1721325931647741*Bmag[10]*f[54]*dv2*wx1+0.253546276418555*Bmag[3]*f[54]*dv2*wx1+0.2267786838055363*Bmag[11]*f[40]*dv2*wx1+0.2267786838055363*Bmag[10]*f[40]*dv2*wx1+0.2309401076758503*Bmag[3]*f[40]*dv2*wx1+0.2581988897471611*Bmag[5]*f[39]*dv2*wx1+0.1844277783908294*Bmag[4]*f[39]*dv2*wx1+0.2886751345948129*Bmag[0]*f[39]*dv2*wx1+0.2581988897471611*Bmag[6]*f[26]*dv2*wx1+0.1844277783908294*Bmag[6]*f[25]*dv2*wx1+0.2886751345948129*Bmag[2]*f[25]*dv2*wx1+0.253546276418555*Bmag[8]*f[16]*dv2*wx1+0.2309401076758504*Bmag[7]*f[16]*dv2*wx1+0.2581988897471611*Bmag[1]*f[16]*dv2*wx1+0.253546276418555*f[8]*Bmag[10]*dv2*wx1+0.2886751345948129*Bmag[4]*f[9]*dv2*wx1+0.2581988897471611*Bmag[3]*f[8]*dv2*wx1+0.2886751345948129*f[4]*Bmag[6]*dv2*wx1+0.06546536707079771*Bmag[7]*f[77]*dv1*dv2+0.04969039949999531*Bmag[8]*f[76]*dv1*dv2+0.06546536707079771*Bmag[7]*f[76]*dv1*dv2+0.07319250547114*Bmag[1]*f[76]*dv1*dv2+0.04969039949999531*Bmag[10]*f[69]*dv1*dv2+0.07319250547113999*Bmag[3]*f[69]*dv1*dv2+0.06546536707079771*Bmag[11]*f[61]*dv1*dv2+0.06546536707079771*Bmag[10]*f[61]*dv1*dv2+0.06666666666666667*Bmag[3]*f[61]*dv1*dv2+0.07453559924999298*Bmag[5]*f[60]*dv1*dv2+0.05323971374999499*Bmag[4]*f[60]*dv1*dv2+0.08333333333333333*Bmag[0]*f[60]*dv1*dv2+0.07453559924999298*Bmag[6]*f[42]*dv1*dv2+0.05323971374999499*Bmag[6]*f[41]*dv1*dv2+0.08333333333333336*Bmag[2]*f[41]*dv1*dv2+0.07319250547113999*Bmag[8]*f[35]*dv1*dv2+0.06666666666666667*Bmag[7]*f[35]*dv1*dv2+0.07453559924999302*Bmag[1]*f[35]*dv1*dv2+0.08333333333333336*Bmag[4]*f[18]*dv1*dv2+0.07319250547114*Bmag[10]*f[17]*dv1*dv2+0.07453559924999302*Bmag[3]*f[17]*dv1*dv2+0.08333333333333333*Bmag[6]*f[10]*dv1*dv2))/m_; 
  out[7] += (volFact*(0.5962847939999437*Bmag[9]*f[49]*wx1*wx2+0.7855844048495726*Bmag[6]*f[49]*wx1*wx2+0.8783100656536798*Bmag[2]*f[49]*wx1*wx2+0.7855844048495726*Bmag[6]*f[48]*wx1*wx2+0.5962847939999437*Bmag[11]*f[32]*wx1*wx2+0.87831006565368*Bmag[3]*f[32]*wx1*wx2+0.6388765649999399*Bmag[5]*f[20]*wx1*wx2+0.8944271909999159*Bmag[4]*f[20]*wx1*wx2+Bmag[0]*f[20]*wx1*wx2+0.7855844048495726*Bmag[11]*f[19]*wx1*wx2+0.7855844048495726*Bmag[10]*f[19]*wx1*wx2+0.8*Bmag[3]*f[19]*wx1*wx2+0.6388765649999399*Bmag[7]*f[12]*wx1*wx2+1.0*Bmag[1]*f[12]*wx1*wx2+0.8944271909999159*Bmag[7]*f[11]*wx1*wx2+0.8783100656536798*f[2]*Bmag[11]*wx1*wx2+0.87831006565368*f[5]*Bmag[9]*wx1*wx2+f[0]*Bmag[7]*wx1*wx2+0.8*f[5]*Bmag[6]*wx1*wx2+0.8944271909999161*Bmag[2]*f[5]*wx1*wx2+1.0*f[1]*Bmag[5]*wx1*wx2+0.8944271909999161*f[2]*Bmag[3]*wx1*wx2+0.1721325931647741*Bmag[9]*f[65]*dv1*wx2+0.2267786838055363*Bmag[6]*f[65]*dv1*wx2+0.253546276418555*Bmag[2]*f[65]*dv1*wx2+0.2267786838055363*Bmag[6]*f[64]*dv1*wx2+0.1721325931647741*Bmag[11]*f[51]*dv1*wx2+0.253546276418555*Bmag[3]*f[51]*dv1*wx2+0.1844277783908294*Bmag[5]*f[37]*dv1*wx2+0.2581988897471611*Bmag[4]*f[37]*dv1*wx2+0.2886751345948129*Bmag[0]*f[37]*dv1*wx2+0.2267786838055363*Bmag[11]*f[36]*dv1*wx2+0.2267786838055363*Bmag[10]*f[36]*dv1*wx2+0.2309401076758503*Bmag[3]*f[36]*dv1*wx2+0.1844277783908294*Bmag[7]*f[22]*dv1*wx2+0.2886751345948129*Bmag[1]*f[22]*dv1*wx2+0.2581988897471611*Bmag[7]*f[21]*dv1*wx2+0.253546276418555*Bmag[9]*f[15]*dv1*wx2+0.2309401076758504*Bmag[6]*f[15]*dv1*wx2+0.2581988897471611*Bmag[2]*f[15]*dv1*wx2+0.253546276418555*f[7]*Bmag[11]*dv1*wx2+0.2581988897471611*Bmag[3]*f[7]*dv1*wx2+0.2886751345948129*f[3]*Bmag[7]*dv1*wx2+0.2886751345948129*Bmag[5]*f[6]*dv1*wx2+0.1721325931647741*Bmag[9]*f[68]*dv2*wx1+0.2267786838055363*Bmag[6]*f[68]*dv2*wx1+0.253546276418555*Bmag[2]*f[68]*dv2*wx1+0.2267786838055363*Bmag[6]*f[67]*dv2*wx1+0.1721325931647741*Bmag[11]*f[55]*dv2*wx1+0.253546276418555*Bmag[3]*f[55]*dv2*wx1+0.1844277783908294*Bmag[5]*f[40]*dv2*wx1+0.2581988897471611*Bmag[4]*f[40]*dv2*wx1+0.2886751345948129*Bmag[0]*f[40]*dv2*wx1+0.2267786838055363*Bmag[11]*f[39]*dv2*wx1+0.2267786838055363*Bmag[10]*f[39]*dv2*wx1+0.2309401076758503*Bmag[3]*f[39]*dv2*wx1+0.1844277783908294*Bmag[7]*f[26]*dv2*wx1+0.2886751345948129*Bmag[1]*f[26]*dv2*wx1+0.2581988897471611*Bmag[7]*f[25]*dv2*wx1+0.253546276418555*Bmag[9]*f[16]*dv2*wx1+0.2309401076758504*Bmag[6]*f[16]*dv2*wx1+0.2581988897471611*Bmag[2]*f[16]*dv2*wx1+0.253546276418555*f[9]*Bmag[11]*dv2*wx1+0.2581988897471611*Bmag[3]*f[9]*dv2*wx1+0.2886751345948129*Bmag[5]*f[8]*dv2*wx1+0.2886751345948129*f[4]*Bmag[7]*dv2*wx1+0.04969039949999531*Bmag[9]*f[77]*dv1*dv2+0.06546536707079771*Bmag[6]*f[77]*dv1*dv2+0.07319250547114*Bmag[2]*f[77]*dv1*dv2+0.06546536707079771*Bmag[6]*f[76]*dv1*dv2+0.04969039949999531*Bmag[11]*f[70]*dv1*dv2+0.07319250547113999*Bmag[3]*f[70]*dv1*dv2+0.05323971374999499*Bmag[5]*f[61]*dv1*dv2+0.07453559924999298*Bmag[4]*f[61]*dv1*dv2+0.08333333333333333*Bmag[0]*f[61]*dv1*dv2+0.06546536707079771*Bmag[11]*f[60]*dv1*dv2+0.06546536707079771*Bmag[10]*f[60]*dv1*dv2+0.06666666666666667*Bmag[3]*f[60]*dv1*dv2+0.05323971374999499*Bmag[7]*f[42]*dv1*dv2+0.08333333333333336*Bmag[1]*f[42]*dv1*dv2+0.07453559924999298*Bmag[7]*f[41]*dv1*dv2+0.07319250547113999*Bmag[9]*f[35]*dv1*dv2+0.06666666666666667*Bmag[6]*f[35]*dv1*dv2+0.07453559924999302*Bmag[2]*f[35]*dv1*dv2+0.07319250547114*Bmag[11]*f[18]*dv1*dv2+0.07453559924999302*Bmag[3]*f[18]*dv1*dv2+0.08333333333333336*Bmag[5]*f[17]*dv1*dv2+0.08333333333333333*Bmag[7]*f[10]*dv1*dv2))/m_; 
  out[8] += (volFact*(0.5962847939999437*Bmag[6]*f[48]*wx1*wx2+1.0*Bmag[2]*f[48]*wx1*wx2+0.5962847939999438*Bmag[4]*f[31]*wx1*wx2+Bmag[0]*f[31]*wx1*wx2+0.5962847939999437*Bmag[10]*f[19]*wx1*wx2+0.8783100656536798*Bmag[3]*f[19]*wx1*wx2+0.5962847939999438*Bmag[8]*f[11]*wx1*wx2+0.8783100656536796*Bmag[1]*f[11]*wx1*wx2+1.0*f[2]*Bmag[10]*wx1*wx2+f[0]*Bmag[8]*wx1*wx2+0.8783100656536798*f[5]*Bmag[6]*wx1*wx2+0.8783100656536796*f[1]*Bmag[4]*wx1*wx2+0.1721325931647741*Bmag[6]*f[64]*dv1*wx2+0.2886751345948129*Bmag[2]*f[64]*dv1*wx2+0.172132593164774*Bmag[4]*f[50]*dv1*wx2+0.2886751345948128*Bmag[0]*f[50]*dv1*wx2+0.172132593164774*Bmag[10]*f[36]*dv1*wx2+0.253546276418555*Bmag[3]*f[36]*dv1*wx2+0.1721325931647741*Bmag[8]*f[21]*dv1*wx2+0.253546276418555*Bmag[1]*f[21]*dv1*wx2+0.253546276418555*Bmag[6]*f[15]*dv1*wx2+0.2886751345948128*f[7]*Bmag[10]*dv1*wx2+0.2886751345948129*f[3]*Bmag[8]*dv1*wx2+0.253546276418555*Bmag[4]*f[6]*dv1*wx2+0.1721325931647741*Bmag[6]*f[67]*dv2*wx1+0.2886751345948129*Bmag[2]*f[67]*dv2*wx1+0.172132593164774*Bmag[4]*f[54]*dv2*wx1+0.2886751345948128*Bmag[0]*f[54]*dv2*wx1+0.172132593164774*Bmag[10]*f[39]*dv2*wx1+0.253546276418555*Bmag[3]*f[39]*dv2*wx1+0.1721325931647741*Bmag[8]*f[25]*dv2*wx1+0.253546276418555*Bmag[1]*f[25]*dv2*wx1+0.253546276418555*Bmag[6]*f[16]*dv2*wx1+0.2886751345948128*f[9]*Bmag[10]*dv2*wx1+0.253546276418555*Bmag[4]*f[8]*dv2*wx1+0.2886751345948129*f[4]*Bmag[8]*dv2*wx1+0.04969039949999531*Bmag[6]*f[76]*dv1*dv2+0.08333333333333334*Bmag[2]*f[76]*dv1*dv2+0.04969039949999532*Bmag[4]*f[69]*dv1*dv2+0.08333333333333333*Bmag[0]*f[69]*dv1*dv2+0.04969039949999531*Bmag[10]*f[60]*dv1*dv2+0.07319250547113999*Bmag[3]*f[60]*dv1*dv2+0.04969039949999532*Bmag[8]*f[41]*dv1*dv2+0.07319250547113998*Bmag[1]*f[41]*dv1*dv2+0.07319250547113999*Bmag[6]*f[35]*dv1*dv2+0.08333333333333334*Bmag[10]*f[18]*dv1*dv2+0.07319250547113998*Bmag[4]*f[17]*dv1*dv2+0.08333333333333333*Bmag[8]*f[10]*dv1*dv2))/m_; 
  out[9] += (volFact*(0.5962847939999437*Bmag[7]*f[49]*wx1*wx2+1.0*Bmag[1]*f[49]*wx1*wx2+0.5962847939999438*Bmag[5]*f[32]*wx1*wx2+Bmag[0]*f[32]*wx1*wx2+0.5962847939999437*Bmag[11]*f[20]*wx1*wx2+0.8783100656536798*Bmag[3]*f[20]*wx1*wx2+0.5962847939999438*Bmag[9]*f[12]*wx1*wx2+0.8783100656536796*Bmag[2]*f[12]*wx1*wx2+1.0*f[1]*Bmag[11]*wx1*wx2+f[0]*Bmag[9]*wx1*wx2+0.8783100656536798*f[5]*Bmag[7]*wx1*wx2+0.8783100656536796*f[2]*Bmag[5]*wx1*wx2+0.1721325931647741*Bmag[7]*f[65]*dv1*wx2+0.2886751345948129*Bmag[1]*f[65]*dv1*wx2+0.172132593164774*Bmag[5]*f[51]*dv1*wx2+0.2886751345948128*Bmag[0]*f[51]*dv1*wx2+0.172132593164774*Bmag[11]*f[37]*dv1*wx2+0.253546276418555*Bmag[3]*f[37]*dv1*wx2+0.1721325931647741*Bmag[9]*f[22]*dv1*wx2+0.253546276418555*Bmag[2]*f[22]*dv1*wx2+0.253546276418555*Bmag[7]*f[15]*dv1*wx2+0.2886751345948128*f[6]*Bmag[11]*dv1*wx2+0.2886751345948129*f[3]*Bmag[9]*dv1*wx2+0.253546276418555*Bmag[5]*f[7]*dv1*wx2+0.1721325931647741*Bmag[7]*f[68]*dv2*wx1+0.2886751345948129*Bmag[1]*f[68]*dv2*wx1+0.172132593164774*Bmag[5]*f[55]*dv2*wx1+0.2886751345948128*Bmag[0]*f[55]*dv2*wx1+0.172132593164774*Bmag[11]*f[40]*dv2*wx1+0.253546276418555*Bmag[3]*f[40]*dv2*wx1+0.1721325931647741*Bmag[9]*f[26]*dv2*wx1+0.253546276418555*Bmag[2]*f[26]*dv2*wx1+0.253546276418555*Bmag[7]*f[16]*dv2*wx1+0.2886751345948128*f[8]*Bmag[11]*dv2*wx1+0.253546276418555*Bmag[5]*f[9]*dv2*wx1+0.2886751345948129*f[4]*Bmag[9]*dv2*wx1+0.04969039949999531*Bmag[7]*f[77]*dv1*dv2+0.08333333333333334*Bmag[1]*f[77]*dv1*dv2+0.04969039949999532*Bmag[5]*f[70]*dv1*dv2+0.08333333333333333*Bmag[0]*f[70]*dv1*dv2+0.04969039949999531*Bmag[11]*f[61]*dv1*dv2+0.07319250547113999*Bmag[3]*f[61]*dv1*dv2+0.04969039949999532*Bmag[9]*f[42]*dv1*dv2+0.07319250547113998*Bmag[2]*f[42]*dv1*dv2+0.07319250547113999*Bmag[7]*f[35]*dv1*dv2+0.07319250547113998*Bmag[5]*f[18]*dv1*dv2+0.08333333333333334*Bmag[11]*f[17]*dv1*dv2+0.08333333333333333*Bmag[9]*f[10]*dv1*dv2))/m_; 
  out[10] += (volFact*(0.8944271909999159*Bmag[5]*f[48]*wx1*wx2+0.5962847939999438*Bmag[4]*f[48]*wx1*wx2+Bmag[0]*f[48]*wx1*wx2+0.5962847939999437*Bmag[6]*f[31]*wx1*wx2+1.0*Bmag[2]*f[31]*wx1*wx2+0.7855844048495726*Bmag[6]*f[20]*wx1*wx2+0.5962847939999437*Bmag[8]*f[19]*wx1*wx2+0.7855844048495726*Bmag[7]*f[19]*wx1*wx2+0.8783100656536798*Bmag[1]*f[19]*wx1*wx2+0.8944271909999159*Bmag[10]*f[12]*wx1*wx2+0.5962847939999438*Bmag[10]*f[11]*wx1*wx2+0.8783100656536798*Bmag[3]*f[11]*wx1*wx2+f[0]*Bmag[10]*wx1*wx2+1.0*f[2]*Bmag[8]*wx1*wx2+0.8783100656536798*f[1]*Bmag[6]*wx1*wx2+0.8783100656536798*Bmag[4]*f[5]*wx1*wx2+0.258198889747161*Bmag[5]*f[64]*dv1*wx2+0.172132593164774*Bmag[4]*f[64]*dv1*wx2+0.2886751345948128*Bmag[0]*f[64]*dv1*wx2+0.1721325931647741*Bmag[6]*f[50]*dv1*wx2+0.2886751345948129*Bmag[2]*f[50]*dv1*wx2+0.2267786838055363*Bmag[6]*f[37]*dv1*wx2+0.172132593164774*Bmag[8]*f[36]*dv1*wx2+0.2267786838055363*Bmag[7]*f[36]*dv1*wx2+0.2535462764185549*Bmag[1]*f[36]*dv1*wx2+0.2581988897471611*Bmag[10]*f[22]*dv1*wx2+0.1721325931647741*Bmag[10]*f[21]*dv1*wx2+0.253546276418555*Bmag[3]*f[21]*dv1*wx2+0.2535462764185549*Bmag[4]*f[15]*dv1*wx2+0.2886751345948129*f[3]*Bmag[10]*dv1*wx2+0.2886751345948128*f[7]*Bmag[8]*dv1*wx2+0.253546276418555*Bmag[6]*f[6]*dv1*wx2+0.258198889747161*Bmag[5]*f[67]*dv2*wx1+0.172132593164774*Bmag[4]*f[67]*dv2*wx1+0.2886751345948128*Bmag[0]*f[67]*dv2*wx1+0.1721325931647741*Bmag[6]*f[54]*dv2*wx1+0.2886751345948129*Bmag[2]*f[54]*dv2*wx1+0.2267786838055363*Bmag[6]*f[40]*dv2*wx1+0.172132593164774*Bmag[8]*f[39]*dv2*wx1+0.2267786838055363*Bmag[7]*f[39]*dv2*wx1+0.2535462764185549*Bmag[1]*f[39]*dv2*wx1+0.2581988897471611*Bmag[10]*f[26]*dv2*wx1+0.1721325931647741*Bmag[10]*f[25]*dv2*wx1+0.253546276418555*Bmag[3]*f[25]*dv2*wx1+0.2535462764185549*Bmag[4]*f[16]*dv2*wx1+0.2886751345948129*f[4]*Bmag[10]*dv2*wx1+0.2886751345948128*Bmag[8]*f[9]*dv2*wx1+0.253546276418555*Bmag[6]*f[8]*dv2*wx1+0.07453559924999298*Bmag[5]*f[76]*dv1*dv2+0.04969039949999532*Bmag[4]*f[76]*dv1*dv2+0.08333333333333333*Bmag[0]*f[76]*dv1*dv2+0.04969039949999531*Bmag[6]*f[69]*dv1*dv2+0.08333333333333334*Bmag[2]*f[69]*dv1*dv2+0.06546536707079771*Bmag[6]*f[61]*dv1*dv2+0.04969039949999531*Bmag[8]*f[60]*dv1*dv2+0.06546536707079771*Bmag[7]*f[60]*dv1*dv2+0.07319250547114*Bmag[1]*f[60]*dv1*dv2+0.07453559924999298*Bmag[10]*f[42]*dv1*dv2+0.04969039949999532*Bmag[10]*f[41]*dv1*dv2+0.07319250547113998*Bmag[3]*f[41]*dv1*dv2+0.07319250547113998*Bmag[4]*f[35]*dv1*dv2+0.08333333333333334*Bmag[8]*f[18]*dv1*dv2+0.07319250547114*Bmag[6]*f[17]*dv1*dv2+0.08333333333333333*Bmag[10]*f[10]*dv1*dv2))/m_; 
  out[11] += (volFact*(0.5962847939999438*Bmag[5]*f[49]*wx1*wx2+0.8944271909999159*Bmag[4]*f[49]*wx1*wx2+Bmag[0]*f[49]*wx1*wx2+0.5962847939999437*Bmag[7]*f[32]*wx1*wx2+1.0*Bmag[1]*f[32]*wx1*wx2+0.5962847939999437*Bmag[9]*f[20]*wx1*wx2+0.7855844048495726*Bmag[6]*f[20]*wx1*wx2+0.8783100656536798*Bmag[2]*f[20]*wx1*wx2+0.7855844048495726*Bmag[7]*f[19]*wx1*wx2+0.5962847939999438*Bmag[11]*f[12]*wx1*wx2+0.8783100656536798*Bmag[3]*f[12]*wx1*wx2+0.8944271909999159*Bmag[11]*f[11]*wx1*wx2+f[0]*Bmag[11]*wx1*wx2+1.0*f[1]*Bmag[9]*wx1*wx2+0.8783100656536798*f[2]*Bmag[7]*wx1*wx2+0.8783100656536798*Bmag[5]*f[5]*wx1*wx2+0.172132593164774*Bmag[5]*f[65]*dv1*wx2+0.258198889747161*Bmag[4]*f[65]*dv1*wx2+0.2886751345948128*Bmag[0]*f[65]*dv1*wx2+0.1721325931647741*Bmag[7]*f[51]*dv1*wx2+0.2886751345948129*Bmag[1]*f[51]*dv1*wx2+0.172132593164774*Bmag[9]*f[37]*dv1*wx2+0.2267786838055363*Bmag[6]*f[37]*dv1*wx2+0.2535462764185549*Bmag[2]*f[37]*dv1*wx2+0.2267786838055363*Bmag[7]*f[36]*dv1*wx2+0.1721325931647741*Bmag[11]*f[22]*dv1*wx2+0.253546276418555*Bmag[3]*f[22]*dv1*wx2+0.2581988897471611*Bmag[11]*f[21]*dv1*wx2+0.2535462764185549*Bmag[5]*f[15]*dv1*wx2+0.2886751345948129*f[3]*Bmag[11]*dv1*wx2+0.2886751345948128*f[6]*Bmag[9]*dv1*wx2+0.253546276418555*Bmag[7]*f[7]*dv1*wx2+0.172132593164774*Bmag[5]*f[68]*dv2*wx1+0.258198889747161*Bmag[4]*f[68]*dv2*wx1+0.2886751345948128*Bmag[0]*f[68]*dv2*wx1+0.1721325931647741*Bmag[7]*f[55]*dv2*wx1+0.2886751345948129*Bmag[1]*f[55]*dv2*wx1+0.172132593164774*Bmag[9]*f[40]*dv2*wx1+0.2267786838055363*Bmag[6]*f[40]*dv2*wx1+0.2535462764185549*Bmag[2]*f[40]*dv2*wx1+0.2267786838055363*Bmag[7]*f[39]*dv2*wx1+0.1721325931647741*Bmag[11]*f[26]*dv2*wx1+0.253546276418555*Bmag[3]*f[26]*dv2*wx1+0.2581988897471611*Bmag[11]*f[25]*dv2*wx1+0.2535462764185549*Bmag[5]*f[16]*dv2*wx1+0.2886751345948129*f[4]*Bmag[11]*dv2*wx1+0.253546276418555*Bmag[7]*f[9]*dv2*wx1+0.2886751345948128*f[8]*Bmag[9]*dv2*wx1+0.04969039949999532*Bmag[5]*f[77]*dv1*dv2+0.07453559924999298*Bmag[4]*f[77]*dv1*dv2+0.08333333333333333*Bmag[0]*f[77]*dv1*dv2+0.04969039949999531*Bmag[7]*f[70]*dv1*dv2+0.08333333333333334*Bmag[1]*f[70]*dv1*dv2+0.04969039949999531*Bmag[9]*f[61]*dv1*dv2+0.06546536707079771*Bmag[6]*f[61]*dv1*dv2+0.07319250547114*Bmag[2]*f[61]*dv1*dv2+0.06546536707079771*Bmag[7]*f[60]*dv1*dv2+0.04969039949999532*Bmag[11]*f[42]*dv1*dv2+0.07319250547113998*Bmag[3]*f[42]*dv1*dv2+0.07453559924999298*Bmag[11]*f[41]*dv1*dv2+0.07319250547113998*Bmag[5]*f[35]*dv1*dv2+0.07319250547114*Bmag[7]*f[18]*dv1*dv2+0.08333333333333334*Bmag[9]*f[17]*dv1*dv2+0.08333333333333333*f[10]*Bmag[11]*dv1*dv2))/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM0[2] += 2.0*f[2]*volFact; 
  outM0[3] += 2.0*f[5]*volFact; 
  outM1[0] += 0.3333333333333333*volFact*(6.0*f[0]*wx1+1.732050807568877*f[3]*dv1); 
  outM1[1] += 0.3333333333333333*volFact*(6.0*f[1]*wx1+1.732050807568877*f[6]*dv1); 
  outM1[2] += 0.3333333333333333*volFact*(6.0*f[2]*wx1+1.732050807568877*f[7]*dv1); 
  outM1[3] += 0.3333333333333333*volFact*(6.0*f[5]*wx1+1.732050807568877*f[11]*dv1); 
  outM2[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq); 
  outM2[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq); 
  outM2[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq); 
  outM2[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.1666666666666667*f[5]*dv1_sq); 
  double tmp[4]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[12]*dv2; 
  outM2[0] += (2.0*(0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  outM2[2] += (2.0*(0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  outM2[3] += (2.0*(0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM0[2] += 2.0*f[2]*volFact; 
  outM0[3] += 2.0*f[5]*volFact; 
  outM0[4] += 2.0*f[11]*volFact; 
  outM0[5] += 2.0*f[12]*volFact; 
  outM0[6] += 2.0*f[19]*volFact; 
  outM0[7] += 2.0*f[20]*volFact; 
  outM1[0] += 0.3333333333333333*volFact*(6.0*f[0]*wx1+1.732050807568877*f[3]*dv1); 
  outM1[1] += 0.3333333333333333*volFact*(6.0*f[1]*wx1+1.732050807568877*f[6]*dv1); 
  outM1[2] += 0.3333333333333333*volFact*(6.0*f[2]*wx1+1.732050807568877*f[7]*dv1); 
  outM1[3] += 0.3333333333333333*volFact*(6.0*f[5]*wx1+1.732050807568877*f[15]*dv1); 
  outM1[4] += 0.06666666666666667*volFact*(30.0*f[11]*wx1+8.660254037844387*f[21]*dv1); 
  outM1[5] += 0.06666666666666667*volFact*(30.0*f[12]*wx1+8.660254037844387*f[22]*dv1); 
  outM1[6] += 0.06666666666666667*volFact*(30.0*f[19]*wx1+8.660254037844387*f[32]*dv1); 
  outM1[7] += 0.06666666666666667*volFact*(30.0*f[20]*wx1+8.660254037844387*f[33]*dv1); 
  outM2[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  outM2[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[23]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  outM2[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[24]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  outM2[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[15]*dv1*wx1+0.149071198499986*f[34]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
  outM2[4] += volFact*(2.0*f[11]*wx1_sq+1.154700538379251*f[21]*dv1*wx1+0.1666666666666667*f[11]*dv1_sq); 
  outM2[5] += volFact*(2.0*f[12]*wx1_sq+1.154700538379251*f[22]*dv1*wx1+0.1666666666666667*f[12]*dv1_sq); 
  outM2[6] += volFact*(2.0*f[19]*wx1_sq+1.154700538379251*f[32]*dv1*wx1+0.1666666666666667*f[19]*dv1_sq); 
  outM2[7] += volFact*(2.0*f[20]*wx1_sq+1.154700538379251*f[33]*dv1*wx1+0.1666666666666667*f[20]*dv1_sq); 
  double tmp[8]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[16]*dv2; 
  tmp[4] = 2.0*f[11]*wx2+0.5773502691896257*f[25]*dv2; 
  tmp[5] = 2.0*f[12]*wx2+0.5773502691896257*f[26]*dv2; 
  tmp[6] = 2.0*f[19]*wx2+0.5773502691896257*f[35]*dv2; 
  tmp[7] = 2.0*f[20]*wx2+0.5773502691896257*f[36]*dv2; 
  outM2[0] += (2.0*(0.5*Bmag[7]*tmp[7]+0.5*Bmag[6]*tmp[6]+0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.5000000000000001*Bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*Bmag[7]+0.447213595499958*Bmag[3]*tmp[6]+0.447213595499958*tmp[3]*Bmag[6]+0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  outM2[2] += (2.0*(0.447213595499958*Bmag[3]*tmp[7]+0.447213595499958*tmp[3]*Bmag[7]+0.5000000000000001*Bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*Bmag[6]+0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  outM2[3] += (2.0*(0.4*Bmag[6]*tmp[7]+0.447213595499958*Bmag[2]*tmp[7]+0.4*tmp[6]*Bmag[7]+0.447213595499958*tmp[2]*Bmag[7]+0.447213595499958*Bmag[1]*tmp[6]+0.447213595499958*tmp[1]*Bmag[6]+0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  outM2[4] += (2.0*(0.4472135954999579*Bmag[7]*tmp[7]+0.31943828249997*Bmag[6]*tmp[6]+0.5000000000000001*Bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*Bmag[6]+0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  outM2[5] += (2.0*(0.31943828249997*Bmag[7]*tmp[7]+0.5000000000000001*Bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*Bmag[7]+0.4472135954999579*Bmag[6]*tmp[6]+0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
  outM2[6] += (2.0*(0.4*Bmag[3]*tmp[7]+0.4*tmp[3]*Bmag[7]+0.4472135954999579*Bmag[5]*tmp[6]+0.31943828249997*Bmag[4]*tmp[6]+0.5*Bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*Bmag[6]+0.31943828249997*tmp[4]*Bmag[6]+0.5*tmp[0]*Bmag[6]+0.5000000000000001*Bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*Bmag[4]+0.447213595499958*Bmag[1]*tmp[3]+0.447213595499958*tmp[1]*Bmag[3])*volFact)/m_; 
  outM2[7] += (2.0*(0.31943828249997*Bmag[5]*tmp[7]+0.4472135954999579*Bmag[4]*tmp[7]+0.5*Bmag[0]*tmp[7]+0.31943828249997*tmp[5]*Bmag[7]+0.4472135954999579*tmp[4]*Bmag[7]+0.5*tmp[0]*Bmag[7]+0.4*Bmag[3]*tmp[6]+0.4*tmp[3]*Bmag[6]+0.5000000000000001*Bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*Bmag[5]+0.447213595499958*Bmag[2]*tmp[3]+0.447213595499958*tmp[2]*Bmag[3])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_ThreeMoments_P3(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM0[2] += 2.0*f[2]*volFact; 
  outM0[3] += 2.0*f[5]*volFact; 
  outM0[4] += 2.0*f[11]*volFact; 
  outM0[5] += 2.0*f[12]*volFact; 
  outM0[6] += 2.0*f[19]*volFact; 
  outM0[7] += 2.0*f[20]*volFact; 
  outM0[8] += 2.0*f[31]*volFact; 
  outM0[9] += 2.0*f[32]*volFact; 
  outM0[10] += 2.0*f[48]*volFact; 
  outM0[11] += 2.0*f[49]*volFact; 
  outM1[0] += 0.3333333333333333*volFact*(6.0*f[0]*wx1+1.732050807568877*f[3]*dv1); 
  outM1[1] += 0.3333333333333333*volFact*(6.0*f[1]*wx1+1.732050807568877*f[6]*dv1); 
  outM1[2] += 0.3333333333333333*volFact*(6.0*f[2]*wx1+1.732050807568877*f[7]*dv1); 
  outM1[3] += 0.3333333333333333*volFact*(6.0*f[5]*wx1+1.732050807568877*f[15]*dv1); 
  outM1[4] += 0.06666666666666667*volFact*(30.0*f[11]*wx1+8.660254037844387*f[21]*dv1); 
  outM1[5] += 0.06666666666666667*volFact*(30.0*f[12]*wx1+8.660254037844387*f[22]*dv1); 
  outM1[6] += 0.06666666666666667*volFact*(30.0*f[19]*wx1+8.660254037844387*f[36]*dv1); 
  outM1[7] += 0.06666666666666667*volFact*(30.0*f[20]*wx1+8.660254037844387*f[37]*dv1); 
  outM1[8] += 0.04761904761904762*volFact*(42.0*f[31]*wx1+12.12435565298214*f[50]*dv1); 
  outM1[9] += 0.04761904761904762*volFact*(42.0*f[32]*wx1+12.12435565298214*f[51]*dv1); 
  outM1[10] += 0.04761904761904762*volFact*(42.0*f[48]*wx1+12.12435565298214*f[64]*dv1); 
  outM1[11] += 0.04761904761904762*volFact*(42.0*f[49]*wx1+12.12435565298214*f[65]*dv1); 
  outM2[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  outM2[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[23]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  outM2[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[24]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  outM2[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[15]*dv1*wx1+0.149071198499986*f[38]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
  outM2[4] += volFact*(2.0*f[11]*wx1_sq+1.154700538379251*f[21]*dv1*wx1+0.1666666666666667*f[11]*dv1_sq); 
  outM2[5] += volFact*(2.0*f[12]*wx1_sq+1.154700538379251*f[22]*dv1*wx1+0.1666666666666667*f[12]*dv1_sq); 
  outM2[6] += volFact*(2.0*f[19]*wx1_sq+1.154700538379251*f[36]*dv1*wx1+0.1666666666666667*f[19]*dv1_sq); 
  outM2[7] += volFact*(2.0*f[20]*wx1_sq+1.154700538379251*f[37]*dv1*wx1+0.1666666666666667*f[20]*dv1_sq); 
  outM2[8] += volFact*(2.0*f[31]*wx1_sq+1.154700538379251*f[50]*dv1*wx1+0.1666666666666667*f[31]*dv1_sq); 
  outM2[9] += volFact*(2.0*f[32]*wx1_sq+1.154700538379251*f[51]*dv1*wx1+0.1666666666666667*f[32]*dv1_sq); 
  outM2[10] += volFact*(2.0*f[48]*wx1_sq+1.154700538379251*f[64]*dv1*wx1+0.1666666666666667*f[48]*dv1_sq); 
  outM2[11] += volFact*(2.0*f[49]*wx1_sq+1.154700538379251*f[65]*dv1*wx1+0.1666666666666667*f[49]*dv1_sq); 
  double tmp[12]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[16]*dv2; 
  tmp[4] = 2.0*f[11]*wx2+0.5773502691896257*f[25]*dv2; 
  tmp[5] = 2.0*f[12]*wx2+0.5773502691896257*f[26]*dv2; 
  tmp[6] = 2.0*f[19]*wx2+0.5773502691896257*f[39]*dv2; 
  tmp[7] = 2.0*f[20]*wx2+0.5773502691896257*f[40]*dv2; 
  tmp[8] = 2.0*f[31]*wx2+0.5773502691896256*f[54]*dv2; 
  tmp[9] = 2.0*f[32]*wx2+0.5773502691896256*f[55]*dv2; 
  tmp[10] = 2.0*f[48]*wx2+0.5773502691896256*f[67]*dv2; 
  tmp[11] = 2.0*f[49]*wx2+0.5773502691896256*f[68]*dv2; 
  outM2[0] += (2.0*(0.5*Bmag[11]*tmp[11]+0.5*Bmag[10]*tmp[10]+0.5*Bmag[9]*tmp[9]+0.5*Bmag[8]*tmp[8]+0.5*Bmag[7]*tmp[7]+0.5*Bmag[6]*tmp[6]+0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.5*Bmag[9]*tmp[11]+0.5*tmp[9]*Bmag[11]+0.4391550328268399*Bmag[6]*tmp[10]+0.4391550328268399*tmp[6]*Bmag[10]+0.4391550328268398*Bmag[4]*tmp[8]+0.4391550328268398*tmp[4]*Bmag[8]+0.5000000000000001*Bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*Bmag[7]+0.447213595499958*Bmag[3]*tmp[6]+0.447213595499958*tmp[3]*Bmag[6]+0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  outM2[2] += (2.0*(0.4391550328268399*Bmag[7]*tmp[11]+0.4391550328268399*tmp[7]*Bmag[11]+0.5*Bmag[8]*tmp[10]+0.5*tmp[8]*Bmag[10]+0.4391550328268398*Bmag[5]*tmp[9]+0.4391550328268398*tmp[5]*Bmag[9]+0.447213595499958*Bmag[3]*tmp[7]+0.447213595499958*tmp[3]*Bmag[7]+0.5000000000000001*Bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*Bmag[6]+0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  outM2[3] += (2.0*(0.4391550328268399*Bmag[5]*tmp[11]+0.4391550328268399*tmp[5]*Bmag[11]+0.4391550328268399*Bmag[4]*tmp[10]+0.4391550328268399*tmp[4]*Bmag[10]+0.4391550328268399*Bmag[7]*tmp[9]+0.4391550328268399*tmp[7]*Bmag[9]+0.4391550328268399*Bmag[6]*tmp[8]+0.4391550328268399*tmp[6]*Bmag[8]+0.4*Bmag[6]*tmp[7]+0.447213595499958*Bmag[2]*tmp[7]+0.4*tmp[6]*Bmag[7]+0.447213595499958*tmp[2]*Bmag[7]+0.447213595499958*Bmag[1]*tmp[6]+0.447213595499958*tmp[1]*Bmag[6]+0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  outM2[4] += (2.0*(0.4472135954999579*Bmag[11]*tmp[11]+0.2981423969999719*Bmag[10]*tmp[10]+0.4391550328268399*Bmag[3]*tmp[10]+0.4391550328268399*tmp[3]*Bmag[10]+0.2981423969999719*Bmag[8]*tmp[8]+0.4391550328268398*Bmag[1]*tmp[8]+0.4391550328268398*tmp[1]*Bmag[8]+0.4472135954999579*Bmag[7]*tmp[7]+0.31943828249997*Bmag[6]*tmp[6]+0.5000000000000001*Bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*Bmag[6]+0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  outM2[5] += (2.0*(0.2981423969999719*Bmag[11]*tmp[11]+0.4391550328268399*Bmag[3]*tmp[11]+0.4391550328268399*tmp[3]*Bmag[11]+0.4472135954999579*Bmag[10]*tmp[10]+0.2981423969999719*Bmag[9]*tmp[9]+0.4391550328268398*Bmag[2]*tmp[9]+0.4391550328268398*tmp[2]*Bmag[9]+0.31943828249997*Bmag[7]*tmp[7]+0.5000000000000001*Bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*Bmag[7]+0.4472135954999579*Bmag[6]*tmp[6]+0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
  outM2[6] += (2.0*(0.3927922024247863*Bmag[7]*tmp[11]+0.3927922024247863*tmp[7]*Bmag[11]+0.2981423969999719*Bmag[8]*tmp[10]+0.3927922024247863*Bmag[7]*tmp[10]+0.4391550328268399*Bmag[1]*tmp[10]+0.2981423969999719*tmp[8]*Bmag[10]+0.3927922024247863*tmp[7]*Bmag[10]+0.4391550328268399*tmp[1]*Bmag[10]+0.4391550328268399*Bmag[3]*tmp[8]+0.4391550328268399*tmp[3]*Bmag[8]+0.4*Bmag[3]*tmp[7]+0.4*tmp[3]*Bmag[7]+0.4472135954999579*Bmag[5]*tmp[6]+0.31943828249997*Bmag[4]*tmp[6]+0.5*Bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*Bmag[6]+0.31943828249997*tmp[4]*Bmag[6]+0.5*tmp[0]*Bmag[6]+0.5000000000000001*Bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*Bmag[4]+0.447213595499958*Bmag[1]*tmp[3]+0.447213595499958*tmp[1]*Bmag[3])*volFact)/m_; 
  outM2[7] += (2.0*(0.2981423969999719*Bmag[9]*tmp[11]+0.3927922024247863*Bmag[6]*tmp[11]+0.4391550328268399*Bmag[2]*tmp[11]+0.2981423969999719*tmp[9]*Bmag[11]+0.3927922024247863*tmp[6]*Bmag[11]+0.4391550328268399*tmp[2]*Bmag[11]+0.3927922024247863*Bmag[6]*tmp[10]+0.3927922024247863*tmp[6]*Bmag[10]+0.4391550328268399*Bmag[3]*tmp[9]+0.4391550328268399*tmp[3]*Bmag[9]+0.31943828249997*Bmag[5]*tmp[7]+0.4472135954999579*Bmag[4]*tmp[7]+0.5*Bmag[0]*tmp[7]+0.31943828249997*tmp[5]*Bmag[7]+0.4472135954999579*tmp[4]*Bmag[7]+0.5*tmp[0]*Bmag[7]+0.4*Bmag[3]*tmp[6]+0.4*tmp[3]*Bmag[6]+0.5000000000000001*Bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*Bmag[5]+0.447213595499958*Bmag[2]*tmp[3]+0.447213595499958*tmp[2]*Bmag[3])*volFact)/m_; 
  outM2[8] += (2.0*(0.2981423969999719*Bmag[6]*tmp[10]+0.5*Bmag[2]*tmp[10]+0.2981423969999719*tmp[6]*Bmag[10]+0.5*tmp[2]*Bmag[10]+0.2981423969999719*Bmag[4]*tmp[8]+0.5*Bmag[0]*tmp[8]+0.2981423969999719*tmp[4]*Bmag[8]+0.5*tmp[0]*Bmag[8]+0.4391550328268399*Bmag[3]*tmp[6]+0.4391550328268399*tmp[3]*Bmag[6]+0.4391550328268398*Bmag[1]*tmp[4]+0.4391550328268398*tmp[1]*Bmag[4])*volFact)/m_; 
  outM2[9] += (2.0*(0.2981423969999719*Bmag[7]*tmp[11]+0.5*Bmag[1]*tmp[11]+0.2981423969999719*tmp[7]*Bmag[11]+0.5*tmp[1]*Bmag[11]+0.2981423969999719*Bmag[5]*tmp[9]+0.5*Bmag[0]*tmp[9]+0.2981423969999719*tmp[5]*Bmag[9]+0.5*tmp[0]*Bmag[9]+0.4391550328268399*Bmag[3]*tmp[7]+0.4391550328268399*tmp[3]*Bmag[7]+0.4391550328268398*Bmag[2]*tmp[5]+0.4391550328268398*tmp[2]*Bmag[5])*volFact)/m_; 
  outM2[10] += (2.0*(0.4472135954999579*Bmag[5]*tmp[10]+0.2981423969999719*Bmag[4]*tmp[10]+0.5*Bmag[0]*tmp[10]+0.4472135954999579*tmp[5]*Bmag[10]+0.2981423969999719*tmp[4]*Bmag[10]+0.5*tmp[0]*Bmag[10]+0.2981423969999719*Bmag[6]*tmp[8]+0.5*Bmag[2]*tmp[8]+0.2981423969999719*tmp[6]*Bmag[8]+0.5*tmp[2]*Bmag[8]+0.3927922024247863*Bmag[6]*tmp[7]+0.3927922024247863*tmp[6]*Bmag[7]+0.4391550328268399*Bmag[1]*tmp[6]+0.4391550328268399*tmp[1]*Bmag[6]+0.4391550328268399*Bmag[3]*tmp[4]+0.4391550328268399*tmp[3]*Bmag[4])*volFact)/m_; 
  outM2[11] += (2.0*(0.2981423969999719*Bmag[5]*tmp[11]+0.4472135954999579*Bmag[4]*tmp[11]+0.5*Bmag[0]*tmp[11]+0.2981423969999719*tmp[5]*Bmag[11]+0.4472135954999579*tmp[4]*Bmag[11]+0.5*tmp[0]*Bmag[11]+0.2981423969999719*Bmag[7]*tmp[9]+0.5*Bmag[1]*tmp[9]+0.2981423969999719*tmp[7]*Bmag[9]+0.5*tmp[1]*Bmag[9]+0.3927922024247863*Bmag[6]*tmp[7]+0.4391550328268399*Bmag[2]*tmp[7]+0.3927922024247863*tmp[6]*Bmag[7]+0.4391550328268399*tmp[2]*Bmag[7]+0.4391550328268399*Bmag[3]*tmp[5]+0.4391550328268399*tmp[3]*Bmag[5])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M0_step1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[2]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[2]*volFact; 
  out[3] += 1.414213562373095*f[4]*volFact; 
  out[4] += 1.414213562373095*f[5]*volFact; 
  out[5] += 1.414213562373095*f[8]*volFact; 
  out[6] += 1.414213562373095*f[9]*volFact; 
  out[7] += 1.414213562373095*f[12]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M0_step1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[2]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[2]*volFact; 
  out[3] += 1.414213562373095*f[4]*volFact; 
  out[4] += 1.414213562373095*f[5]*volFact; 
  out[5] += 1.414213562373095*f[8]*volFact; 
  out[6] += 1.414213562373095*f[9]*volFact; 
  out[7] += 1.414213562373095*f[11]*volFact; 
  out[8] += 1.414213562373095*f[12]*volFact; 
  out[9] += 1.414213562373095*f[14]*volFact; 
  out[10] += 1.414213562373095*f[16]*volFact; 
  out[11] += 1.414213562373095*f[19]*volFact; 
  out[12] += 1.414213562373095*f[20]*volFact; 
  out[13] += 1.414213562373095*f[25]*volFact; 
  out[14] += 1.414213562373095*f[26]*volFact; 
  out[15] += 1.414213562373095*f[28]*volFact; 
  out[16] += 1.414213562373095*f[29]*volFact; 
  out[17] += 1.414213562373095*f[35]*volFact; 
  out[18] += 1.414213562373095*f[36]*volFact; 
  out[19] += 1.414213562373095*f[41]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M0_step1_P3(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[2]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[2]*volFact; 
  out[3] += 1.414213562373095*f[4]*volFact; 
  out[4] += 1.414213562373095*f[5]*volFact; 
  out[5] += 1.414213562373095*f[8]*volFact; 
  out[6] += 1.414213562373095*f[9]*volFact; 
  out[7] += 1.414213562373095*f[11]*volFact; 
  out[8] += 1.414213562373095*f[12]*volFact; 
  out[9] += 1.414213562373095*f[14]*volFact; 
  out[10] += 1.414213562373095*f[16]*volFact; 
  out[11] += 1.414213562373095*f[19]*volFact; 
  out[12] += 1.414213562373095*f[20]*volFact; 
  out[13] += 1.414213562373095*f[25]*volFact; 
  out[14] += 1.414213562373095*f[26]*volFact; 
  out[15] += 1.414213562373095*f[28]*volFact; 
  out[16] += 1.414213562373095*f[29]*volFact; 
  out[17] += 1.414213562373095*f[31]*volFact; 
  out[18] += 1.414213562373095*f[32]*volFact; 
  out[19] += 1.414213562373095*f[34]*volFact; 
  out[20] += 1.414213562373095*f[39]*volFact; 
  out[21] += 1.414213562373095*f[40]*volFact; 
  out[22] += 1.414213562373095*f[45]*volFact; 
  out[23] += 1.414213562373095*f[48]*volFact; 
  out[24] += 1.414213562373095*f[49]*volFact; 
  out[25] += 1.414213562373095*f[54]*volFact; 
  out[26] += 1.414213562373095*f[55]*volFact; 
  out[27] += 1.414213562373095*f[57]*volFact; 
  out[28] += 1.414213562373095*f[58]*volFact; 
  out[29] += 1.414213562373095*f[67]*volFact; 
  out[30] += 1.414213562373095*f[68]*volFact; 
  out[31] += 1.414213562373095*f[73]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M0_step2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]/2; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[4]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M0_step2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]/2; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[4]*volFact; 
  out[4] += 2.828427124746191*f[7]*volFact; 
  out[5] += 2.828427124746191*f[8]*volFact; 
  out[6] += 2.828427124746191*f[11]*volFact; 
  out[7] += 2.828427124746191*f[12]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M0_step2_P3(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]/2; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[4]*volFact; 
  out[4] += 2.828427124746191*f[7]*volFact; 
  out[5] += 2.828427124746191*f[8]*volFact; 
  out[6] += 2.828427124746191*f[11]*volFact; 
  out[7] += 2.828427124746191*f[12]*volFact; 
  out[8] += 2.828427124746191*f[17]*volFact; 
  out[9] += 2.828427124746191*f[18]*volFact; 
  out[10] += 2.828427124746191*f[23]*volFact; 
  out[11] += 2.828427124746191*f[24]*volFact; 
} 
