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
__host__ __device__ void GkMomentCalc2x2vSer_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  double tmp[4]; 
  tmp[0] = 4.0*f[0]*wx2+1.154700538379252*f[4]*dv2; 
  tmp[1] = 4.0*f[1]*wx2+1.154700538379252*f[8]*dv2; 
  tmp[2] = 4.0*f[2]*wx2+1.154700538379252*f[9]*dv2; 
  tmp[3] = 4.0*f[5]*wx2+1.154700538379252*f[12]*dv2; 
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
  tmp[0] = 4.0*f[0]*wx2+1.154700538379252*f[4]*dv2; 
  tmp[1] = 4.0*f[1]*wx2+1.154700538379252*f[8]*dv2; 
  tmp[2] = 4.0*f[2]*wx2+1.154700538379252*f[9]*dv2; 
  tmp[3] = 4.0*f[5]*wx2+1.154700538379252*f[16]*dv2; 
  tmp[4] = 4.0*f[11]*wx2+1.154700538379251*f[25]*dv2; 
  tmp[5] = 4.0*f[12]*wx2+1.154700538379251*f[26]*dv2; 
  tmp[6] = 4.0*f[19]*wx2+1.154700538379251*f[35]*dv2; 
  tmp[7] = 4.0*f[20]*wx2+1.154700538379251*f[36]*dv2; 
  out[0] += ((0.5*Bmag[7]*tmp[7]+0.5*Bmag[6]*tmp[6]+0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.5000000000000001*Bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*Bmag[7]+0.447213595499958*Bmag[3]*tmp[6]+0.447213595499958*tmp[3]*Bmag[6]+0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += ((0.447213595499958*Bmag[3]*tmp[7]+0.447213595499958*tmp[3]*Bmag[7]+0.5000000000000001*Bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*Bmag[6]+0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += ((0.4*Bmag[6]*tmp[7]+0.447213595499958*Bmag[2]*tmp[7]+0.4*tmp[6]*Bmag[7]+0.447213595499958*tmp[2]*Bmag[7]+0.447213595499958*Bmag[1]*tmp[6]+0.447213595499958*tmp[1]*Bmag[6]+0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  out[4] += ((0.4472135954999579*Bmag[7]*tmp[7]+0.31943828249997*Bmag[6]*tmp[6]+0.5000000000000001*Bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*Bmag[6]+0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  out[5] += ((0.31943828249997*Bmag[7]*tmp[7]+0.5000000000000001*Bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*Bmag[7]+0.4472135954999579*Bmag[6]*tmp[6]+0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
  out[6] += ((0.4*Bmag[3]*tmp[7]+0.4*tmp[3]*Bmag[7]+0.4472135954999579*Bmag[5]*tmp[6]+0.31943828249997*Bmag[4]*tmp[6]+0.5*Bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*Bmag[6]+0.31943828249997*tmp[4]*Bmag[6]+0.5*tmp[0]*Bmag[6]+0.5000000000000001*Bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*Bmag[4]+0.447213595499958*Bmag[1]*tmp[3]+0.447213595499958*tmp[1]*Bmag[3])*volFact)/m_; 
  out[7] += ((0.31943828249997*Bmag[5]*tmp[7]+0.4472135954999579*Bmag[4]*tmp[7]+0.5*Bmag[0]*tmp[7]+0.31943828249997*tmp[5]*Bmag[7]+0.4472135954999579*tmp[4]*Bmag[7]+0.5*tmp[0]*Bmag[7]+0.4*Bmag[3]*tmp[6]+0.4*tmp[3]*Bmag[6]+0.5000000000000001*Bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*Bmag[5]+0.447213595499958*Bmag[2]*tmp[3]+0.447213595499958*tmp[2]*Bmag[3])*volFact)/m_; 
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
__host__ __device__ void GkMomentCalc2x2vSer_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (volFact*(2.0*Bmag[3]*f[5]*wx1*wx2+2.0*Bmag[2]*f[2]*wx1*wx2+2.0*Bmag[1]*f[1]*wx1*wx2+2.0*Bmag[0]*f[0]*wx1*wx2+0.5773502691896258*Bmag[3]*f[11]*dv1*wx2+0.5773502691896258*Bmag[2]*f[7]*dv1*wx2+0.5773502691896258*Bmag[1]*f[6]*dv1*wx2+0.5773502691896258*Bmag[0]*f[3]*dv1*wx2+0.5773502691896258*Bmag[3]*f[12]*dv2*wx1+0.5773502691896258*Bmag[2]*f[9]*dv2*wx1+0.5773502691896258*Bmag[1]*f[8]*dv2*wx1+0.5773502691896258*Bmag[0]*f[4]*dv2*wx1+0.1666666666666667*Bmag[3]*f[15]*dv1*dv2+0.1666666666666667*Bmag[2]*f[14]*dv1*dv2+0.1666666666666667*Bmag[1]*f[13]*dv1*dv2+0.1666666666666667*Bmag[0]*f[10]*dv1*dv2))/m_; 
  out[1] += (volFact*(2.0*Bmag[2]*f[5]*wx1*wx2+2.0*f[2]*Bmag[3]*wx1*wx2+2.0*Bmag[0]*f[1]*wx1*wx2+2.0*f[0]*Bmag[1]*wx1*wx2+0.5773502691896258*Bmag[2]*f[11]*dv1*wx2+0.5773502691896258*Bmag[3]*f[7]*dv1*wx2+0.5773502691896258*Bmag[0]*f[6]*dv1*wx2+0.5773502691896258*Bmag[1]*f[3]*dv1*wx2+0.5773502691896258*Bmag[2]*f[12]*dv2*wx1+0.5773502691896258*Bmag[3]*f[9]*dv2*wx1+0.5773502691896258*Bmag[0]*f[8]*dv2*wx1+0.5773502691896258*Bmag[1]*f[4]*dv2*wx1+0.1666666666666667*Bmag[2]*f[15]*dv1*dv2+0.1666666666666667*Bmag[3]*f[14]*dv1*dv2+0.1666666666666667*Bmag[0]*f[13]*dv1*dv2+0.1666666666666667*Bmag[1]*f[10]*dv1*dv2))/m_; 
  out[2] += (volFact*(2.0*Bmag[1]*f[5]*wx1*wx2+2.0*f[1]*Bmag[3]*wx1*wx2+2.0*Bmag[0]*f[2]*wx1*wx2+2.0*f[0]*Bmag[2]*wx1*wx2+0.5773502691896258*Bmag[1]*f[11]*dv1*wx2+0.5773502691896258*Bmag[0]*f[7]*dv1*wx2+0.5773502691896258*Bmag[3]*f[6]*dv1*wx2+0.5773502691896258*Bmag[2]*f[3]*dv1*wx2+0.5773502691896258*Bmag[1]*f[12]*dv2*wx1+0.5773502691896258*Bmag[0]*f[9]*dv2*wx1+0.5773502691896258*Bmag[3]*f[8]*dv2*wx1+0.5773502691896258*Bmag[2]*f[4]*dv2*wx1+0.1666666666666667*Bmag[1]*f[15]*dv1*dv2+0.1666666666666667*Bmag[0]*f[14]*dv1*dv2+0.1666666666666667*Bmag[3]*f[13]*dv1*dv2+0.1666666666666667*Bmag[2]*f[10]*dv1*dv2))/m_; 
  out[3] += (volFact*(2.0*Bmag[0]*f[5]*wx1*wx2+2.0*f[0]*Bmag[3]*wx1*wx2+2.0*Bmag[1]*f[2]*wx1*wx2+2.0*f[1]*Bmag[2]*wx1*wx2+0.5773502691896258*Bmag[0]*f[11]*dv1*wx2+0.5773502691896258*Bmag[1]*f[7]*dv1*wx2+0.5773502691896258*Bmag[2]*f[6]*dv1*wx2+0.5773502691896258*Bmag[3]*f[3]*dv1*wx2+0.5773502691896258*Bmag[0]*f[12]*dv2*wx1+0.5773502691896258*Bmag[1]*f[9]*dv2*wx1+0.5773502691896258*Bmag[2]*f[8]*dv2*wx1+0.5773502691896258*Bmag[3]*f[4]*dv2*wx1+0.1666666666666667*Bmag[0]*f[15]*dv1*dv2+0.1666666666666667*Bmag[1]*f[14]*dv1*dv2+0.1666666666666667*Bmag[2]*f[13]*dv1*dv2+0.1666666666666667*Bmag[3]*f[10]*dv1*dv2))/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vSer_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (volFact*(2.0*Bmag[7]*f[20]*wx1*wx2+2.0*Bmag[6]*f[19]*wx1*wx2+2.0*Bmag[5]*f[12]*wx1*wx2+2.0*Bmag[4]*f[11]*wx1*wx2+2.0*Bmag[3]*f[5]*wx1*wx2+2.0*Bmag[2]*f[2]*wx1*wx2+2.0*Bmag[1]*f[1]*wx1*wx2+2.0*Bmag[0]*f[0]*wx1*wx2+0.5773502691896257*Bmag[7]*f[33]*dv1*wx2+0.5773502691896257*Bmag[6]*f[32]*dv1*wx2+0.5773502691896257*Bmag[5]*f[22]*dv1*wx2+0.5773502691896257*Bmag[4]*f[21]*dv1*wx2+0.5773502691896258*Bmag[3]*f[15]*dv1*wx2+0.5773502691896258*Bmag[2]*f[7]*dv1*wx2+0.5773502691896258*Bmag[1]*f[6]*dv1*wx2+0.5773502691896258*Bmag[0]*f[3]*dv1*wx2+0.5773502691896257*Bmag[7]*f[36]*dv2*wx1+0.5773502691896257*Bmag[6]*f[35]*dv2*wx1+0.5773502691896257*Bmag[5]*f[26]*dv2*wx1+0.5773502691896257*Bmag[4]*f[25]*dv2*wx1+0.5773502691896258*Bmag[3]*f[16]*dv2*wx1+0.5773502691896258*Bmag[2]*f[9]*dv2*wx1+0.5773502691896258*Bmag[1]*f[8]*dv2*wx1+0.5773502691896258*Bmag[0]*f[4]*dv2*wx1+0.1666666666666667*Bmag[7]*f[45]*dv1*dv2+0.1666666666666667*Bmag[6]*f[44]*dv1*dv2+0.1666666666666667*Bmag[5]*f[38]*dv1*dv2+0.1666666666666667*Bmag[4]*f[37]*dv1*dv2+0.1666666666666667*Bmag[3]*f[31]*dv1*dv2+0.1666666666666667*Bmag[2]*f[18]*dv1*dv2+0.1666666666666667*Bmag[1]*f[17]*dv1*dv2+0.1666666666666667*Bmag[0]*f[10]*dv1*dv2))/m_; 
  out[1] += (volFact*(2.0*Bmag[5]*f[20]*wx1*wx2+1.788854381999832*Bmag[3]*f[19]*wx1*wx2+2.0*Bmag[7]*f[12]*wx1*wx2+1.788854381999832*Bmag[1]*f[11]*wx1*wx2+1.788854381999832*f[5]*Bmag[6]*wx1*wx2+2.0*Bmag[2]*f[5]*wx1*wx2+1.788854381999832*f[1]*Bmag[4]*wx1*wx2+2.0*f[2]*Bmag[3]*wx1*wx2+2.0*Bmag[0]*f[1]*wx1*wx2+2.0*f[0]*Bmag[1]*wx1*wx2+0.5773502691896258*Bmag[5]*f[33]*dv1*wx2+0.5163977794943223*Bmag[3]*f[32]*dv1*wx2+0.5773502691896258*Bmag[7]*f[22]*dv1*wx2+0.5163977794943222*Bmag[1]*f[21]*dv1*wx2+0.5163977794943222*Bmag[6]*f[15]*dv1*wx2+0.5773502691896258*Bmag[2]*f[15]*dv1*wx2+0.5773502691896258*Bmag[3]*f[7]*dv1*wx2+0.5163977794943223*Bmag[4]*f[6]*dv1*wx2+0.5773502691896258*Bmag[0]*f[6]*dv1*wx2+0.5773502691896258*Bmag[1]*f[3]*dv1*wx2+0.5773502691896258*Bmag[5]*f[36]*dv2*wx1+0.5163977794943223*Bmag[3]*f[35]*dv2*wx1+0.5773502691896258*Bmag[7]*f[26]*dv2*wx1+0.5163977794943222*Bmag[1]*f[25]*dv2*wx1+0.5163977794943222*Bmag[6]*f[16]*dv2*wx1+0.5773502691896258*Bmag[2]*f[16]*dv2*wx1+0.5773502691896258*Bmag[3]*f[9]*dv2*wx1+0.5163977794943223*Bmag[4]*f[8]*dv2*wx1+0.5773502691896258*Bmag[0]*f[8]*dv2*wx1+0.5773502691896258*Bmag[1]*f[4]*dv2*wx1+0.1666666666666667*Bmag[5]*f[45]*dv1*dv2+0.149071198499986*Bmag[3]*f[44]*dv1*dv2+0.1666666666666667*Bmag[7]*f[38]*dv1*dv2+0.149071198499986*Bmag[1]*f[37]*dv1*dv2+0.149071198499986*Bmag[6]*f[31]*dv1*dv2+0.1666666666666667*Bmag[2]*f[31]*dv1*dv2+0.1666666666666667*Bmag[3]*f[18]*dv1*dv2+0.149071198499986*Bmag[4]*f[17]*dv1*dv2+0.1666666666666667*Bmag[0]*f[17]*dv1*dv2+0.1666666666666667*Bmag[1]*f[10]*dv1*dv2))/m_; 
  out[2] += (volFact*(1.788854381999832*Bmag[3]*f[20]*wx1*wx2+2.0*Bmag[4]*f[19]*wx1*wx2+1.788854381999832*Bmag[2]*f[12]*wx1*wx2+2.0*Bmag[6]*f[11]*wx1*wx2+1.788854381999832*f[5]*Bmag[7]*wx1*wx2+2.0*Bmag[1]*f[5]*wx1*wx2+1.788854381999832*f[2]*Bmag[5]*wx1*wx2+2.0*f[1]*Bmag[3]*wx1*wx2+2.0*Bmag[0]*f[2]*wx1*wx2+2.0*f[0]*Bmag[2]*wx1*wx2+0.5163977794943223*Bmag[3]*f[33]*dv1*wx2+0.5773502691896258*Bmag[4]*f[32]*dv1*wx2+0.5163977794943222*Bmag[2]*f[22]*dv1*wx2+0.5773502691896258*Bmag[6]*f[21]*dv1*wx2+0.5163977794943222*Bmag[7]*f[15]*dv1*wx2+0.5773502691896258*Bmag[1]*f[15]*dv1*wx2+0.5163977794943223*Bmag[5]*f[7]*dv1*wx2+0.5773502691896258*Bmag[0]*f[7]*dv1*wx2+0.5773502691896258*Bmag[3]*f[6]*dv1*wx2+0.5773502691896258*Bmag[2]*f[3]*dv1*wx2+0.5163977794943223*Bmag[3]*f[36]*dv2*wx1+0.5773502691896258*Bmag[4]*f[35]*dv2*wx1+0.5163977794943222*Bmag[2]*f[26]*dv2*wx1+0.5773502691896258*Bmag[6]*f[25]*dv2*wx1+0.5163977794943222*Bmag[7]*f[16]*dv2*wx1+0.5773502691896258*Bmag[1]*f[16]*dv2*wx1+0.5163977794943223*Bmag[5]*f[9]*dv2*wx1+0.5773502691896258*Bmag[0]*f[9]*dv2*wx1+0.5773502691896258*Bmag[3]*f[8]*dv2*wx1+0.5773502691896258*Bmag[2]*f[4]*dv2*wx1+0.149071198499986*Bmag[3]*f[45]*dv1*dv2+0.1666666666666667*Bmag[4]*f[44]*dv1*dv2+0.149071198499986*Bmag[2]*f[38]*dv1*dv2+0.1666666666666667*Bmag[6]*f[37]*dv1*dv2+0.149071198499986*Bmag[7]*f[31]*dv1*dv2+0.1666666666666667*Bmag[1]*f[31]*dv1*dv2+0.149071198499986*Bmag[5]*f[18]*dv1*dv2+0.1666666666666667*Bmag[0]*f[18]*dv1*dv2+0.1666666666666667*Bmag[3]*f[17]*dv1*dv2+0.1666666666666667*Bmag[2]*f[10]*dv1*dv2))/m_; 
  out[3] += (volFact*(1.6*Bmag[6]*f[20]*wx1*wx2+1.788854381999832*Bmag[2]*f[20]*wx1*wx2+1.6*Bmag[7]*f[19]*wx1*wx2+1.788854381999832*Bmag[1]*f[19]*wx1*wx2+1.788854381999832*Bmag[3]*f[12]*wx1*wx2+1.788854381999832*Bmag[3]*f[11]*wx1*wx2+1.788854381999832*f[2]*Bmag[7]*wx1*wx2+1.788854381999832*f[1]*Bmag[6]*wx1*wx2+1.788854381999832*Bmag[5]*f[5]*wx1*wx2+1.788854381999832*Bmag[4]*f[5]*wx1*wx2+2.0*Bmag[0]*f[5]*wx1*wx2+2.0*f[0]*Bmag[3]*wx1*wx2+2.0*Bmag[1]*f[2]*wx1*wx2+2.0*f[1]*Bmag[2]*wx1*wx2+0.4618802153517005*Bmag[6]*f[33]*dv1*wx2+0.5163977794943223*Bmag[2]*f[33]*dv1*wx2+0.4618802153517005*Bmag[7]*f[32]*dv1*wx2+0.5163977794943223*Bmag[1]*f[32]*dv1*wx2+0.5163977794943222*Bmag[3]*f[22]*dv1*wx2+0.5163977794943222*Bmag[3]*f[21]*dv1*wx2+0.5163977794943223*Bmag[5]*f[15]*dv1*wx2+0.5163977794943223*Bmag[4]*f[15]*dv1*wx2+0.5773502691896258*Bmag[0]*f[15]*dv1*wx2+0.5163977794943222*Bmag[7]*f[7]*dv1*wx2+0.5773502691896258*Bmag[1]*f[7]*dv1*wx2+0.5163977794943222*Bmag[6]*f[6]*dv1*wx2+0.5773502691896258*Bmag[2]*f[6]*dv1*wx2+0.5773502691896258*Bmag[3]*f[3]*dv1*wx2+0.4618802153517005*Bmag[6]*f[36]*dv2*wx1+0.5163977794943223*Bmag[2]*f[36]*dv2*wx1+0.4618802153517005*Bmag[7]*f[35]*dv2*wx1+0.5163977794943223*Bmag[1]*f[35]*dv2*wx1+0.5163977794943222*Bmag[3]*f[26]*dv2*wx1+0.5163977794943222*Bmag[3]*f[25]*dv2*wx1+0.5163977794943223*Bmag[5]*f[16]*dv2*wx1+0.5163977794943223*Bmag[4]*f[16]*dv2*wx1+0.5773502691896258*Bmag[0]*f[16]*dv2*wx1+0.5163977794943222*Bmag[7]*f[9]*dv2*wx1+0.5773502691896258*Bmag[1]*f[9]*dv2*wx1+0.5163977794943222*Bmag[6]*f[8]*dv2*wx1+0.5773502691896258*Bmag[2]*f[8]*dv2*wx1+0.5773502691896258*Bmag[3]*f[4]*dv2*wx1+0.1333333333333333*Bmag[6]*f[45]*dv1*dv2+0.149071198499986*Bmag[2]*f[45]*dv1*dv2+0.1333333333333333*Bmag[7]*f[44]*dv1*dv2+0.149071198499986*Bmag[1]*f[44]*dv1*dv2+0.149071198499986*Bmag[3]*f[38]*dv1*dv2+0.149071198499986*Bmag[3]*f[37]*dv1*dv2+0.149071198499986*Bmag[5]*f[31]*dv1*dv2+0.149071198499986*Bmag[4]*f[31]*dv1*dv2+0.1666666666666667*Bmag[0]*f[31]*dv1*dv2+0.149071198499986*Bmag[7]*f[18]*dv1*dv2+0.1666666666666667*Bmag[1]*f[18]*dv1*dv2+0.149071198499986*Bmag[6]*f[17]*dv1*dv2+0.1666666666666667*Bmag[2]*f[17]*dv1*dv2+0.1666666666666667*Bmag[3]*f[10]*dv1*dv2))/m_; 
  out[4] += (volFact*(1.788854381999832*Bmag[7]*f[20]*wx1*wx2+1.27775312999988*Bmag[6]*f[19]*wx1*wx2+2.0*Bmag[2]*f[19]*wx1*wx2+1.27775312999988*Bmag[4]*f[11]*wx1*wx2+2.0*Bmag[0]*f[11]*wx1*wx2+2.0*f[2]*Bmag[6]*wx1*wx2+1.788854381999832*Bmag[3]*f[5]*wx1*wx2+2.0*f[0]*Bmag[4]*wx1*wx2+1.788854381999832*Bmag[1]*f[1]*wx1*wx2+0.5163977794943222*Bmag[7]*f[33]*dv1*wx2+0.3688555567816588*Bmag[6]*f[32]*dv1*wx2+0.5773502691896258*Bmag[2]*f[32]*dv1*wx2+0.3688555567816588*Bmag[4]*f[21]*dv1*wx2+0.5773502691896257*Bmag[0]*f[21]*dv1*wx2+0.5163977794943223*Bmag[3]*f[15]*dv1*wx2+0.5773502691896257*Bmag[6]*f[7]*dv1*wx2+0.5163977794943223*Bmag[1]*f[6]*dv1*wx2+0.5773502691896258*f[3]*Bmag[4]*dv1*wx2+0.5163977794943222*Bmag[7]*f[36]*dv2*wx1+0.3688555567816588*Bmag[6]*f[35]*dv2*wx1+0.5773502691896258*Bmag[2]*f[35]*dv2*wx1+0.3688555567816588*Bmag[4]*f[25]*dv2*wx1+0.5773502691896257*Bmag[0]*f[25]*dv2*wx1+0.5163977794943223*Bmag[3]*f[16]*dv2*wx1+0.5773502691896257*Bmag[6]*f[9]*dv2*wx1+0.5163977794943223*Bmag[1]*f[8]*dv2*wx1+0.5773502691896258*Bmag[4]*f[4]*dv2*wx1+0.149071198499986*Bmag[7]*f[45]*dv1*dv2+0.10647942749999*Bmag[6]*f[44]*dv1*dv2+0.1666666666666667*Bmag[2]*f[44]*dv1*dv2+0.10647942749999*Bmag[4]*f[37]*dv1*dv2+0.1666666666666667*Bmag[0]*f[37]*dv1*dv2+0.149071198499986*Bmag[3]*f[31]*dv1*dv2+0.1666666666666667*Bmag[6]*f[18]*dv1*dv2+0.149071198499986*Bmag[1]*f[17]*dv1*dv2+0.1666666666666667*Bmag[4]*f[10]*dv1*dv2))/m_; 
  out[5] += (volFact*(1.27775312999988*Bmag[7]*f[20]*wx1*wx2+2.0*Bmag[1]*f[20]*wx1*wx2+1.788854381999832*Bmag[6]*f[19]*wx1*wx2+1.27775312999988*Bmag[5]*f[12]*wx1*wx2+2.0*Bmag[0]*f[12]*wx1*wx2+2.0*f[1]*Bmag[7]*wx1*wx2+1.788854381999832*Bmag[3]*f[5]*wx1*wx2+2.0*f[0]*Bmag[5]*wx1*wx2+1.788854381999832*Bmag[2]*f[2]*wx1*wx2+0.3688555567816588*Bmag[7]*f[33]*dv1*wx2+0.5773502691896258*Bmag[1]*f[33]*dv1*wx2+0.5163977794943222*Bmag[6]*f[32]*dv1*wx2+0.3688555567816588*Bmag[5]*f[22]*dv1*wx2+0.5773502691896257*Bmag[0]*f[22]*dv1*wx2+0.5163977794943223*Bmag[3]*f[15]*dv1*wx2+0.5163977794943223*Bmag[2]*f[7]*dv1*wx2+0.5773502691896257*f[6]*Bmag[7]*dv1*wx2+0.5773502691896258*f[3]*Bmag[5]*dv1*wx2+0.3688555567816588*Bmag[7]*f[36]*dv2*wx1+0.5773502691896258*Bmag[1]*f[36]*dv2*wx1+0.5163977794943222*Bmag[6]*f[35]*dv2*wx1+0.3688555567816588*Bmag[5]*f[26]*dv2*wx1+0.5773502691896257*Bmag[0]*f[26]*dv2*wx1+0.5163977794943223*Bmag[3]*f[16]*dv2*wx1+0.5163977794943223*Bmag[2]*f[9]*dv2*wx1+0.5773502691896257*Bmag[7]*f[8]*dv2*wx1+0.5773502691896258*f[4]*Bmag[5]*dv2*wx1+0.10647942749999*Bmag[7]*f[45]*dv1*dv2+0.1666666666666667*Bmag[1]*f[45]*dv1*dv2+0.149071198499986*Bmag[6]*f[44]*dv1*dv2+0.10647942749999*Bmag[5]*f[38]*dv1*dv2+0.1666666666666667*Bmag[0]*f[38]*dv1*dv2+0.149071198499986*Bmag[3]*f[31]*dv1*dv2+0.149071198499986*Bmag[2]*f[18]*dv1*dv2+0.1666666666666667*Bmag[7]*f[17]*dv1*dv2+0.1666666666666667*Bmag[5]*f[10]*dv1*dv2))/m_; 
  out[6] += (volFact*(1.6*Bmag[3]*f[20]*wx1*wx2+1.788854381999832*Bmag[5]*f[19]*wx1*wx2+1.27775312999988*Bmag[4]*f[19]*wx1*wx2+2.0*Bmag[0]*f[19]*wx1*wx2+1.788854381999832*Bmag[6]*f[12]*wx1*wx2+1.27775312999988*Bmag[6]*f[11]*wx1*wx2+2.0*Bmag[2]*f[11]*wx1*wx2+1.6*f[5]*Bmag[7]*wx1*wx2+2.0*f[0]*Bmag[6]*wx1*wx2+1.788854381999832*Bmag[1]*f[5]*wx1*wx2+2.0*f[2]*Bmag[4]*wx1*wx2+1.788854381999832*f[1]*Bmag[3]*wx1*wx2+0.4618802153517005*Bmag[3]*f[33]*dv1*wx2+0.5163977794943222*Bmag[5]*f[32]*dv1*wx2+0.3688555567816588*Bmag[4]*f[32]*dv1*wx2+0.5773502691896257*Bmag[0]*f[32]*dv1*wx2+0.5163977794943222*Bmag[6]*f[22]*dv1*wx2+0.3688555567816588*Bmag[6]*f[21]*dv1*wx2+0.5773502691896258*Bmag[2]*f[21]*dv1*wx2+0.4618802153517007*Bmag[7]*f[15]*dv1*wx2+0.5163977794943222*Bmag[1]*f[15]*dv1*wx2+0.5773502691896257*Bmag[4]*f[7]*dv1*wx2+0.5163977794943222*Bmag[3]*f[6]*dv1*wx2+0.5773502691896258*f[3]*Bmag[6]*dv1*wx2+0.4618802153517005*Bmag[3]*f[36]*dv2*wx1+0.5163977794943222*Bmag[5]*f[35]*dv2*wx1+0.3688555567816588*Bmag[4]*f[35]*dv2*wx1+0.5773502691896257*Bmag[0]*f[35]*dv2*wx1+0.5163977794943222*Bmag[6]*f[26]*dv2*wx1+0.3688555567816588*Bmag[6]*f[25]*dv2*wx1+0.5773502691896258*Bmag[2]*f[25]*dv2*wx1+0.4618802153517007*Bmag[7]*f[16]*dv2*wx1+0.5163977794943222*Bmag[1]*f[16]*dv2*wx1+0.5773502691896257*Bmag[4]*f[9]*dv2*wx1+0.5163977794943222*Bmag[3]*f[8]*dv2*wx1+0.5773502691896258*f[4]*Bmag[6]*dv2*wx1+0.1333333333333333*Bmag[3]*f[45]*dv1*dv2+0.149071198499986*Bmag[5]*f[44]*dv1*dv2+0.10647942749999*Bmag[4]*f[44]*dv1*dv2+0.1666666666666667*Bmag[0]*f[44]*dv1*dv2+0.149071198499986*Bmag[6]*f[38]*dv1*dv2+0.10647942749999*Bmag[6]*f[37]*dv1*dv2+0.1666666666666667*Bmag[2]*f[37]*dv1*dv2+0.1333333333333333*Bmag[7]*f[31]*dv1*dv2+0.149071198499986*Bmag[1]*f[31]*dv1*dv2+0.1666666666666667*Bmag[4]*f[18]*dv1*dv2+0.149071198499986*Bmag[3]*f[17]*dv1*dv2+0.1666666666666667*Bmag[6]*f[10]*dv1*dv2))/m_; 
  out[7] += (volFact*(1.27775312999988*Bmag[5]*f[20]*wx1*wx2+1.788854381999832*Bmag[4]*f[20]*wx1*wx2+2.0*Bmag[0]*f[20]*wx1*wx2+1.6*Bmag[3]*f[19]*wx1*wx2+1.27775312999988*Bmag[7]*f[12]*wx1*wx2+2.0*Bmag[1]*f[12]*wx1*wx2+1.788854381999832*Bmag[7]*f[11]*wx1*wx2+2.0*f[0]*Bmag[7]*wx1*wx2+1.6*f[5]*Bmag[6]*wx1*wx2+1.788854381999832*Bmag[2]*f[5]*wx1*wx2+2.0*f[1]*Bmag[5]*wx1*wx2+1.788854381999832*f[2]*Bmag[3]*wx1*wx2+0.3688555567816588*Bmag[5]*f[33]*dv1*wx2+0.5163977794943222*Bmag[4]*f[33]*dv1*wx2+0.5773502691896257*Bmag[0]*f[33]*dv1*wx2+0.4618802153517005*Bmag[3]*f[32]*dv1*wx2+0.3688555567816588*Bmag[7]*f[22]*dv1*wx2+0.5773502691896258*Bmag[1]*f[22]*dv1*wx2+0.5163977794943222*Bmag[7]*f[21]*dv1*wx2+0.4618802153517007*Bmag[6]*f[15]*dv1*wx2+0.5163977794943222*Bmag[2]*f[15]*dv1*wx2+0.5163977794943222*Bmag[3]*f[7]*dv1*wx2+0.5773502691896258*f[3]*Bmag[7]*dv1*wx2+0.5773502691896257*Bmag[5]*f[6]*dv1*wx2+0.3688555567816588*Bmag[5]*f[36]*dv2*wx1+0.5163977794943222*Bmag[4]*f[36]*dv2*wx1+0.5773502691896257*Bmag[0]*f[36]*dv2*wx1+0.4618802153517005*Bmag[3]*f[35]*dv2*wx1+0.3688555567816588*Bmag[7]*f[26]*dv2*wx1+0.5773502691896258*Bmag[1]*f[26]*dv2*wx1+0.5163977794943222*Bmag[7]*f[25]*dv2*wx1+0.4618802153517007*Bmag[6]*f[16]*dv2*wx1+0.5163977794943222*Bmag[2]*f[16]*dv2*wx1+0.5163977794943222*Bmag[3]*f[9]*dv2*wx1+0.5773502691896257*Bmag[5]*f[8]*dv2*wx1+0.5773502691896258*f[4]*Bmag[7]*dv2*wx1+0.10647942749999*Bmag[5]*f[45]*dv1*dv2+0.149071198499986*Bmag[4]*f[45]*dv1*dv2+0.1666666666666667*Bmag[0]*f[45]*dv1*dv2+0.1333333333333333*Bmag[3]*f[44]*dv1*dv2+0.10647942749999*Bmag[7]*f[38]*dv1*dv2+0.1666666666666667*Bmag[1]*f[38]*dv1*dv2+0.149071198499986*Bmag[7]*f[37]*dv1*dv2+0.1333333333333333*Bmag[6]*f[31]*dv1*dv2+0.149071198499986*Bmag[2]*f[31]*dv1*dv2+0.149071198499986*Bmag[3]*f[18]*dv1*dv2+0.1666666666666667*Bmag[5]*f[17]*dv1*dv2+0.1666666666666667*Bmag[7]*f[10]*dv1*dv2))/m_; 
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
