#include <math.h> 
#include <DistFuncMomentCalcModDecl.h> 
__host__ __device__ void GkMomentCalc2x2vTensor_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
  out[8] += 2.0*f[44]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1); 
  out[2] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1); 
  out[3] += volFact*(2.0*f[5]*wx1+0.5773502691896258*f[11]*dv1); 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
  out[8] += volFact*(2.0*f[44]*wx1+0.5773502691896258*f[54]*dv1); 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M1proj_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += 2.0*f[0]*volFact*wx1; 
  out[1] += 2.0*f[1]*volFact*wx1; 
  out[2] += 2.0*f[2]*volFact*wx1; 
  out[3] += 2.0*f[5]*volFact*wx1; 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M1proj_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
  out[8] += volFact*(2.0*f[44]*wx1+0.5773502691896258*f[54]*dv1); 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
__host__ __device__ void GkMomentCalc2x2vTensor_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
  out[4] += volFact*(2.0*f[11]*wx1_sq+1.154700538379251*f[21]*dv1*wx1+0.149071198499986*f[45]*dv1_sq+0.1666666666666667*f[11]*dv1_sq); 
  out[5] += volFact*(2.0*f[12]*wx1_sq+1.154700538379251*f[22]*dv1*wx1+0.149071198499986*f[46]*dv1_sq+0.1666666666666667*f[12]*dv1_sq); 
  out[6] += volFact*(2.0*f[19]*wx1_sq+1.154700538379251*f[32]*dv1*wx1+0.149071198499986*f[55]*dv1_sq+0.1666666666666667*f[19]*dv1_sq); 
  out[7] += volFact*(2.0*f[20]*wx1_sq+1.154700538379251*f[33]*dv1*wx1+0.149071198499986*f[56]*dv1_sq+0.1666666666666667*f[20]*dv1_sq); 
  out[8] += volFact*(2.0*f[44]*wx1_sq+1.154700538379252*f[54]*dv1*wx1+0.149071198499986*f[72]*dv1_sq+0.1666666666666667*f[44]*dv1_sq); 
  double tmp[9]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[16]*dv2; 
  tmp[4] = 2.0*f[11]*wx2+0.5773502691896257*f[25]*dv2; 
  tmp[5] = 2.0*f[12]*wx2+0.5773502691896257*f[26]*dv2; 
  tmp[6] = 2.0*f[19]*wx2+0.5773502691896257*f[35]*dv2; 
  tmp[7] = 2.0*f[20]*wx2+0.5773502691896257*f[36]*dv2; 
  tmp[8] = 2.0*f[44]*wx2+0.5773502691896258*f[57]*dv2; 
  out[0] += (2.0*(0.5*Bmag[8]*tmp[8]+0.5*Bmag[7]*tmp[7]+0.5*Bmag[6]*tmp[6]+0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.447213595499958*Bmag[7]*tmp[8]+0.447213595499958*tmp[7]*Bmag[8]+0.5000000000000001*Bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*Bmag[7]+0.447213595499958*Bmag[3]*tmp[6]+0.447213595499958*tmp[3]*Bmag[6]+0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.447213595499958*Bmag[6]*tmp[8]+0.447213595499958*tmp[6]*Bmag[8]+0.447213595499958*Bmag[3]*tmp[7]+0.447213595499958*tmp[3]*Bmag[7]+0.5000000000000001*Bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*Bmag[6]+0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.4*Bmag[3]*tmp[8]+0.4*tmp[3]*Bmag[8]+0.4*Bmag[6]*tmp[7]+0.447213595499958*Bmag[2]*tmp[7]+0.4*tmp[6]*Bmag[7]+0.447213595499958*tmp[2]*Bmag[7]+0.447213595499958*Bmag[1]*tmp[6]+0.447213595499958*tmp[1]*Bmag[6]+0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  out[4] += (2.0*(0.31943828249997*Bmag[8]*tmp[8]+0.5*Bmag[5]*tmp[8]+0.5*tmp[5]*Bmag[8]+0.4472135954999579*Bmag[7]*tmp[7]+0.31943828249997*Bmag[6]*tmp[6]+0.5000000000000001*Bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*Bmag[6]+0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  out[5] += (2.0*(0.31943828249997*Bmag[8]*tmp[8]+0.5*Bmag[4]*tmp[8]+0.5*tmp[4]*Bmag[8]+0.31943828249997*Bmag[7]*tmp[7]+0.5000000000000001*Bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*Bmag[7]+0.4472135954999579*Bmag[6]*tmp[6]+0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
  out[6] += (2.0*(0.2857142857142857*Bmag[6]*tmp[8]+0.447213595499958*Bmag[2]*tmp[8]+0.2857142857142857*tmp[6]*Bmag[8]+0.447213595499958*tmp[2]*Bmag[8]+0.4*Bmag[3]*tmp[7]+0.4*tmp[3]*Bmag[7]+0.4472135954999579*Bmag[5]*tmp[6]+0.31943828249997*Bmag[4]*tmp[6]+0.5*Bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*Bmag[6]+0.31943828249997*tmp[4]*Bmag[6]+0.5*tmp[0]*Bmag[6]+0.5000000000000001*Bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*Bmag[4]+0.447213595499958*Bmag[1]*tmp[3]+0.447213595499958*tmp[1]*Bmag[3])*volFact)/m_; 
  out[7] += (2.0*(0.2857142857142857*Bmag[7]*tmp[8]+0.447213595499958*Bmag[1]*tmp[8]+0.2857142857142857*tmp[7]*Bmag[8]+0.447213595499958*tmp[1]*Bmag[8]+0.31943828249997*Bmag[5]*tmp[7]+0.4472135954999579*Bmag[4]*tmp[7]+0.5*Bmag[0]*tmp[7]+0.31943828249997*tmp[5]*Bmag[7]+0.4472135954999579*tmp[4]*Bmag[7]+0.5*tmp[0]*Bmag[7]+0.4*Bmag[3]*tmp[6]+0.4*tmp[3]*Bmag[6]+0.5000000000000001*Bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*Bmag[5]+0.447213595499958*Bmag[2]*tmp[3]+0.447213595499958*tmp[2]*Bmag[3])*volFact)/m_; 
  out[8] += (2.0*(0.2040816326530612*Bmag[8]*tmp[8]+0.31943828249997*Bmag[5]*tmp[8]+0.31943828249997*Bmag[4]*tmp[8]+0.5*Bmag[0]*tmp[8]+0.31943828249997*tmp[5]*Bmag[8]+0.31943828249997*tmp[4]*Bmag[8]+0.5*tmp[0]*Bmag[8]+0.2857142857142857*Bmag[7]*tmp[7]+0.447213595499958*Bmag[1]*tmp[7]+0.447213595499958*tmp[1]*Bmag[7]+0.2857142857142857*Bmag[6]*tmp[6]+0.447213595499958*Bmag[2]*tmp[6]+0.447213595499958*tmp[2]*Bmag[6]+0.5*Bmag[4]*tmp[5]+0.5*tmp[4]*Bmag[5]+0.4*Bmag[3]*tmp[3])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
__host__ __device__ void GkMomentCalc2x2vTensor_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
  out[4] += volFact*(2.0*f[11]*wx1_sq+1.154700538379251*f[21]*dv1*wx1+0.149071198499986*f[45]*dv1_sq+0.1666666666666667*f[11]*dv1_sq); 
  out[5] += volFact*(2.0*f[12]*wx1_sq+1.154700538379251*f[22]*dv1*wx1+0.149071198499986*f[46]*dv1_sq+0.1666666666666667*f[12]*dv1_sq); 
  out[6] += volFact*(2.0*f[19]*wx1_sq+1.154700538379251*f[32]*dv1*wx1+0.149071198499986*f[55]*dv1_sq+0.1666666666666667*f[19]*dv1_sq); 
  out[7] += volFact*(2.0*f[20]*wx1_sq+1.154700538379251*f[33]*dv1*wx1+0.149071198499986*f[56]*dv1_sq+0.1666666666666667*f[20]*dv1_sq); 
  out[8] += volFact*(2.0*f[44]*wx1_sq+1.154700538379252*f[54]*dv1*wx1+0.149071198499986*f[72]*dv1_sq+0.1666666666666667*f[44]*dv1_sq); 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
__host__ __device__ void GkMomentCalc2x2vTensor_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  double tmp[9]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[16]*dv2; 
  tmp[4] = 2.0*f[11]*wx2+0.5773502691896257*f[25]*dv2; 
  tmp[5] = 2.0*f[12]*wx2+0.5773502691896257*f[26]*dv2; 
  tmp[6] = 2.0*f[19]*wx2+0.5773502691896257*f[35]*dv2; 
  tmp[7] = 2.0*f[20]*wx2+0.5773502691896257*f[36]*dv2; 
  tmp[8] = 2.0*f[44]*wx2+0.5773502691896258*f[57]*dv2; 
  out[0] += ((0.5*Bmag[8]*tmp[8]+0.5*Bmag[7]*tmp[7]+0.5*Bmag[6]*tmp[6]+0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.447213595499958*Bmag[7]*tmp[8]+0.447213595499958*tmp[7]*Bmag[8]+0.5000000000000001*Bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*Bmag[7]+0.447213595499958*Bmag[3]*tmp[6]+0.447213595499958*tmp[3]*Bmag[6]+0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += ((0.447213595499958*Bmag[6]*tmp[8]+0.447213595499958*tmp[6]*Bmag[8]+0.447213595499958*Bmag[3]*tmp[7]+0.447213595499958*tmp[3]*Bmag[7]+0.5000000000000001*Bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*Bmag[6]+0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += ((0.4*Bmag[3]*tmp[8]+0.4*tmp[3]*Bmag[8]+0.4*Bmag[6]*tmp[7]+0.447213595499958*Bmag[2]*tmp[7]+0.4*tmp[6]*Bmag[7]+0.447213595499958*tmp[2]*Bmag[7]+0.447213595499958*Bmag[1]*tmp[6]+0.447213595499958*tmp[1]*Bmag[6]+0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  out[4] += ((0.31943828249997*Bmag[8]*tmp[8]+0.5*Bmag[5]*tmp[8]+0.5*tmp[5]*Bmag[8]+0.4472135954999579*Bmag[7]*tmp[7]+0.31943828249997*Bmag[6]*tmp[6]+0.5000000000000001*Bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*Bmag[6]+0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  out[5] += ((0.31943828249997*Bmag[8]*tmp[8]+0.5*Bmag[4]*tmp[8]+0.5*tmp[4]*Bmag[8]+0.31943828249997*Bmag[7]*tmp[7]+0.5000000000000001*Bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*Bmag[7]+0.4472135954999579*Bmag[6]*tmp[6]+0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
  out[6] += ((0.2857142857142857*Bmag[6]*tmp[8]+0.447213595499958*Bmag[2]*tmp[8]+0.2857142857142857*tmp[6]*Bmag[8]+0.447213595499958*tmp[2]*Bmag[8]+0.4*Bmag[3]*tmp[7]+0.4*tmp[3]*Bmag[7]+0.4472135954999579*Bmag[5]*tmp[6]+0.31943828249997*Bmag[4]*tmp[6]+0.5*Bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*Bmag[6]+0.31943828249997*tmp[4]*Bmag[6]+0.5*tmp[0]*Bmag[6]+0.5000000000000001*Bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*Bmag[4]+0.447213595499958*Bmag[1]*tmp[3]+0.447213595499958*tmp[1]*Bmag[3])*volFact)/m_; 
  out[7] += ((0.2857142857142857*Bmag[7]*tmp[8]+0.447213595499958*Bmag[1]*tmp[8]+0.2857142857142857*tmp[7]*Bmag[8]+0.447213595499958*tmp[1]*Bmag[8]+0.31943828249997*Bmag[5]*tmp[7]+0.4472135954999579*Bmag[4]*tmp[7]+0.5*Bmag[0]*tmp[7]+0.31943828249997*tmp[5]*Bmag[7]+0.4472135954999579*tmp[4]*Bmag[7]+0.5*tmp[0]*Bmag[7]+0.4*Bmag[3]*tmp[6]+0.4*tmp[3]*Bmag[6]+0.5000000000000001*Bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*Bmag[5]+0.447213595499958*Bmag[2]*tmp[3]+0.447213595499958*tmp[2]*Bmag[3])*volFact)/m_; 
  out[8] += ((0.2040816326530612*Bmag[8]*tmp[8]+0.31943828249997*Bmag[5]*tmp[8]+0.31943828249997*Bmag[4]*tmp[8]+0.5*Bmag[0]*tmp[8]+0.31943828249997*tmp[5]*Bmag[8]+0.31943828249997*tmp[4]*Bmag[8]+0.5*tmp[0]*Bmag[8]+0.2857142857142857*Bmag[7]*tmp[7]+0.447213595499958*Bmag[1]*tmp[7]+0.447213595499958*tmp[1]*Bmag[7]+0.2857142857142857*Bmag[6]*tmp[6]+0.447213595499958*Bmag[2]*tmp[6]+0.447213595499958*tmp[2]*Bmag[6]+0.5*Bmag[4]*tmp[5]+0.5*tmp[4]*Bmag[5]+0.4*Bmag[3]*tmp[3])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
__host__ __device__ void GkMomentCalc2x2vTensor_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
  out[4] += volFact*(2.0*f[11]*wx1*wx1_sq+1.732050807568877*f[21]*dv1*wx1_sq+0.4472135954999579*f[45]*dv1_sq*wx1+0.5*f[11]*dv1_sq*wx1+0.08660254037844385*f[21]*dv1*dv1_sq); 
  out[5] += volFact*(2.0*f[12]*wx1*wx1_sq+1.732050807568877*f[22]*dv1*wx1_sq+0.4472135954999579*f[46]*dv1_sq*wx1+0.5*f[12]*dv1_sq*wx1+0.08660254037844385*f[22]*dv1*dv1_sq); 
  out[6] += volFact*(2.0*f[19]*wx1*wx1_sq+1.732050807568877*f[32]*dv1*wx1_sq+0.447213595499958*f[55]*dv1_sq*wx1+0.5*f[19]*dv1_sq*wx1+0.08660254037844385*f[32]*dv1*dv1_sq); 
  out[7] += volFact*(2.0*f[20]*wx1*wx1_sq+1.732050807568877*f[33]*dv1*wx1_sq+0.447213595499958*f[56]*dv1_sq*wx1+0.5*f[20]*dv1_sq*wx1+0.08660254037844385*f[33]*dv1*dv1_sq); 
  out[8] += volFact*(2.0*f[44]*wx1*wx1_sq+1.732050807568877*f[54]*dv1*wx1_sq+0.4472135954999579*f[72]*dv1_sq*wx1+0.5*f[44]*dv1_sq*wx1+0.08660254037844387*f[54]*dv1*dv1_sq); 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (volFact*(Bmag[3]*f[5]*wx1*wx2+Bmag[2]*f[2]*wx1*wx2+Bmag[1]*f[1]*wx1*wx2+Bmag[0]*f[0]*wx1*wx2+0.2886751345948129*Bmag[3]*f[11]*dv1*wx2+0.2886751345948129*Bmag[2]*f[7]*dv1*wx2+0.2886751345948129*Bmag[1]*f[6]*dv1*wx2+0.2886751345948129*Bmag[0]*f[3]*dv1*wx2+0.2886751345948129*Bmag[3]*f[12]*dv2*wx1+0.2886751345948129*Bmag[2]*f[9]*dv2*wx1+0.2886751345948129*Bmag[1]*f[8]*dv2*wx1+0.2886751345948129*Bmag[0]*f[4]*dv2*wx1+0.08333333333333333*Bmag[3]*f[15]*dv1*dv2+0.08333333333333333*Bmag[2]*f[14]*dv1*dv2+0.08333333333333333*Bmag[1]*f[13]*dv1*dv2+0.08333333333333333*Bmag[0]*f[10]*dv1*dv2))/m_; 
  out[1] += (volFact*(Bmag[2]*f[5]*wx1*wx2+f[2]*Bmag[3]*wx1*wx2+Bmag[0]*f[1]*wx1*wx2+f[0]*Bmag[1]*wx1*wx2+0.2886751345948129*Bmag[2]*f[11]*dv1*wx2+0.2886751345948129*Bmag[3]*f[7]*dv1*wx2+0.2886751345948129*Bmag[0]*f[6]*dv1*wx2+0.2886751345948129*Bmag[1]*f[3]*dv1*wx2+0.2886751345948129*Bmag[2]*f[12]*dv2*wx1+0.2886751345948129*Bmag[3]*f[9]*dv2*wx1+0.2886751345948129*Bmag[0]*f[8]*dv2*wx1+0.2886751345948129*Bmag[1]*f[4]*dv2*wx1+0.08333333333333333*Bmag[2]*f[15]*dv1*dv2+0.08333333333333333*Bmag[3]*f[14]*dv1*dv2+0.08333333333333333*Bmag[0]*f[13]*dv1*dv2+0.08333333333333333*Bmag[1]*f[10]*dv1*dv2))/m_; 
  out[2] += (volFact*(Bmag[1]*f[5]*wx1*wx2+f[1]*Bmag[3]*wx1*wx2+Bmag[0]*f[2]*wx1*wx2+f[0]*Bmag[2]*wx1*wx2+0.2886751345948129*Bmag[1]*f[11]*dv1*wx2+0.2886751345948129*Bmag[0]*f[7]*dv1*wx2+0.2886751345948129*Bmag[3]*f[6]*dv1*wx2+0.2886751345948129*Bmag[2]*f[3]*dv1*wx2+0.2886751345948129*Bmag[1]*f[12]*dv2*wx1+0.2886751345948129*Bmag[0]*f[9]*dv2*wx1+0.2886751345948129*Bmag[3]*f[8]*dv2*wx1+0.2886751345948129*Bmag[2]*f[4]*dv2*wx1+0.08333333333333333*Bmag[1]*f[15]*dv1*dv2+0.08333333333333333*Bmag[0]*f[14]*dv1*dv2+0.08333333333333333*Bmag[3]*f[13]*dv1*dv2+0.08333333333333333*Bmag[2]*f[10]*dv1*dv2))/m_; 
  out[3] += (volFact*(Bmag[0]*f[5]*wx1*wx2+f[0]*Bmag[3]*wx1*wx2+Bmag[1]*f[2]*wx1*wx2+f[1]*Bmag[2]*wx1*wx2+0.2886751345948129*Bmag[0]*f[11]*dv1*wx2+0.2886751345948129*Bmag[1]*f[7]*dv1*wx2+0.2886751345948129*Bmag[2]*f[6]*dv1*wx2+0.2886751345948129*Bmag[3]*f[3]*dv1*wx2+0.2886751345948129*Bmag[0]*f[12]*dv2*wx1+0.2886751345948129*Bmag[1]*f[9]*dv2*wx1+0.2886751345948129*Bmag[2]*f[8]*dv2*wx1+0.2886751345948129*Bmag[3]*f[4]*dv2*wx1+0.08333333333333333*Bmag[0]*f[15]*dv1*dv2+0.08333333333333333*Bmag[1]*f[14]*dv1*dv2+0.08333333333333333*Bmag[2]*f[13]*dv1*dv2+0.08333333333333333*Bmag[3]*f[10]*dv1*dv2))/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (volFact*(Bmag[8]*f[44]*wx1*wx2+Bmag[7]*f[20]*wx1*wx2+Bmag[6]*f[19]*wx1*wx2+Bmag[5]*f[12]*wx1*wx2+Bmag[4]*f[11]*wx1*wx2+Bmag[3]*f[5]*wx1*wx2+Bmag[2]*f[2]*wx1*wx2+Bmag[1]*f[1]*wx1*wx2+Bmag[0]*f[0]*wx1*wx2+0.2886751345948129*Bmag[8]*f[54]*dv1*wx2+0.2886751345948129*Bmag[7]*f[33]*dv1*wx2+0.2886751345948129*Bmag[6]*f[32]*dv1*wx2+0.2886751345948129*Bmag[5]*f[22]*dv1*wx2+0.2886751345948129*Bmag[4]*f[21]*dv1*wx2+0.2886751345948129*Bmag[3]*f[15]*dv1*wx2+0.2886751345948129*Bmag[2]*f[7]*dv1*wx2+0.2886751345948129*Bmag[1]*f[6]*dv1*wx2+0.2886751345948129*Bmag[0]*f[3]*dv1*wx2+0.2886751345948129*Bmag[8]*f[57]*dv2*wx1+0.2886751345948129*Bmag[7]*f[36]*dv2*wx1+0.2886751345948129*Bmag[6]*f[35]*dv2*wx1+0.2886751345948129*Bmag[5]*f[26]*dv2*wx1+0.2886751345948129*Bmag[4]*f[25]*dv2*wx1+0.2886751345948129*Bmag[3]*f[16]*dv2*wx1+0.2886751345948129*Bmag[2]*f[9]*dv2*wx1+0.2886751345948129*Bmag[1]*f[8]*dv2*wx1+0.2886751345948129*Bmag[0]*f[4]*dv2*wx1+0.08333333333333333*Bmag[8]*f[66]*dv1*dv2+0.08333333333333333*Bmag[7]*f[51]*dv1*dv2+0.08333333333333333*Bmag[6]*f[50]*dv1*dv2+0.08333333333333333*Bmag[5]*f[38]*dv1*dv2+0.08333333333333333*Bmag[4]*f[37]*dv1*dv2+0.08333333333333333*Bmag[3]*f[31]*dv1*dv2+0.08333333333333333*Bmag[2]*f[18]*dv1*dv2+0.08333333333333333*Bmag[1]*f[17]*dv1*dv2+0.08333333333333333*Bmag[0]*f[10]*dv1*dv2))/m_; 
  out[1] += (volFact*(0.8944271909999161*Bmag[7]*f[44]*wx1*wx2+0.8944271909999161*Bmag[8]*f[20]*wx1*wx2+1.0*Bmag[5]*f[20]*wx1*wx2+0.8944271909999161*Bmag[3]*f[19]*wx1*wx2+1.0*Bmag[7]*f[12]*wx1*wx2+0.8944271909999159*Bmag[1]*f[11]*wx1*wx2+0.8944271909999161*f[5]*Bmag[6]*wx1*wx2+Bmag[2]*f[5]*wx1*wx2+0.8944271909999159*f[1]*Bmag[4]*wx1*wx2+f[2]*Bmag[3]*wx1*wx2+Bmag[0]*f[1]*wx1*wx2+f[0]*Bmag[1]*wx1*wx2+0.2581988897471611*Bmag[7]*f[54]*dv1*wx2+0.2581988897471612*Bmag[8]*f[33]*dv1*wx2+0.2886751345948129*Bmag[5]*f[33]*dv1*wx2+0.2581988897471612*Bmag[3]*f[32]*dv1*wx2+0.2886751345948129*Bmag[7]*f[22]*dv1*wx2+0.2581988897471611*Bmag[1]*f[21]*dv1*wx2+0.2581988897471611*Bmag[6]*f[15]*dv1*wx2+0.2886751345948129*Bmag[2]*f[15]*dv1*wx2+0.2886751345948129*Bmag[3]*f[7]*dv1*wx2+0.2581988897471612*Bmag[4]*f[6]*dv1*wx2+0.2886751345948129*Bmag[0]*f[6]*dv1*wx2+0.2886751345948129*Bmag[1]*f[3]*dv1*wx2+0.2581988897471611*Bmag[7]*f[57]*dv2*wx1+0.2581988897471612*Bmag[8]*f[36]*dv2*wx1+0.2886751345948129*Bmag[5]*f[36]*dv2*wx1+0.2581988897471612*Bmag[3]*f[35]*dv2*wx1+0.2886751345948129*Bmag[7]*f[26]*dv2*wx1+0.2581988897471611*Bmag[1]*f[25]*dv2*wx1+0.2581988897471611*Bmag[6]*f[16]*dv2*wx1+0.2886751345948129*Bmag[2]*f[16]*dv2*wx1+0.2886751345948129*Bmag[3]*f[9]*dv2*wx1+0.2581988897471612*Bmag[4]*f[8]*dv2*wx1+0.2886751345948129*Bmag[0]*f[8]*dv2*wx1+0.2886751345948129*Bmag[1]*f[4]*dv2*wx1+0.07453559924999302*Bmag[7]*f[66]*dv1*dv2+0.07453559924999302*Bmag[8]*f[51]*dv1*dv2+0.08333333333333336*Bmag[5]*f[51]*dv1*dv2+0.07453559924999302*Bmag[3]*f[50]*dv1*dv2+0.08333333333333336*Bmag[7]*f[38]*dv1*dv2+0.07453559924999298*Bmag[1]*f[37]*dv1*dv2+0.07453559924999302*Bmag[6]*f[31]*dv1*dv2+0.08333333333333333*Bmag[2]*f[31]*dv1*dv2+0.08333333333333333*Bmag[3]*f[18]*dv1*dv2+0.07453559924999298*Bmag[4]*f[17]*dv1*dv2+0.08333333333333333*Bmag[0]*f[17]*dv1*dv2+0.08333333333333333*Bmag[1]*f[10]*dv1*dv2))/m_; 
  out[2] += (volFact*(0.8944271909999161*Bmag[6]*f[44]*wx1*wx2+0.8944271909999161*Bmag[3]*f[20]*wx1*wx2+0.8944271909999161*Bmag[8]*f[19]*wx1*wx2+1.0*Bmag[4]*f[19]*wx1*wx2+0.8944271909999159*Bmag[2]*f[12]*wx1*wx2+1.0*Bmag[6]*f[11]*wx1*wx2+0.8944271909999161*f[5]*Bmag[7]*wx1*wx2+Bmag[1]*f[5]*wx1*wx2+0.8944271909999159*f[2]*Bmag[5]*wx1*wx2+f[1]*Bmag[3]*wx1*wx2+Bmag[0]*f[2]*wx1*wx2+f[0]*Bmag[2]*wx1*wx2+0.2581988897471611*Bmag[6]*f[54]*dv1*wx2+0.2581988897471612*Bmag[3]*f[33]*dv1*wx2+0.2581988897471612*Bmag[8]*f[32]*dv1*wx2+0.2886751345948129*Bmag[4]*f[32]*dv1*wx2+0.2581988897471611*Bmag[2]*f[22]*dv1*wx2+0.2886751345948129*Bmag[6]*f[21]*dv1*wx2+0.2581988897471611*Bmag[7]*f[15]*dv1*wx2+0.2886751345948129*Bmag[1]*f[15]*dv1*wx2+0.2581988897471612*Bmag[5]*f[7]*dv1*wx2+0.2886751345948129*Bmag[0]*f[7]*dv1*wx2+0.2886751345948129*Bmag[3]*f[6]*dv1*wx2+0.2886751345948129*Bmag[2]*f[3]*dv1*wx2+0.2581988897471611*Bmag[6]*f[57]*dv2*wx1+0.2581988897471612*Bmag[3]*f[36]*dv2*wx1+0.2581988897471612*Bmag[8]*f[35]*dv2*wx1+0.2886751345948129*Bmag[4]*f[35]*dv2*wx1+0.2581988897471611*Bmag[2]*f[26]*dv2*wx1+0.2886751345948129*Bmag[6]*f[25]*dv2*wx1+0.2581988897471611*Bmag[7]*f[16]*dv2*wx1+0.2886751345948129*Bmag[1]*f[16]*dv2*wx1+0.2581988897471612*Bmag[5]*f[9]*dv2*wx1+0.2886751345948129*Bmag[0]*f[9]*dv2*wx1+0.2886751345948129*Bmag[3]*f[8]*dv2*wx1+0.2886751345948129*Bmag[2]*f[4]*dv2*wx1+0.07453559924999302*Bmag[6]*f[66]*dv1*dv2+0.07453559924999302*Bmag[3]*f[51]*dv1*dv2+0.07453559924999302*Bmag[8]*f[50]*dv1*dv2+0.08333333333333336*Bmag[4]*f[50]*dv1*dv2+0.07453559924999298*Bmag[2]*f[38]*dv1*dv2+0.08333333333333336*Bmag[6]*f[37]*dv1*dv2+0.07453559924999302*Bmag[7]*f[31]*dv1*dv2+0.08333333333333333*Bmag[1]*f[31]*dv1*dv2+0.07453559924999298*Bmag[5]*f[18]*dv1*dv2+0.08333333333333333*Bmag[0]*f[18]*dv1*dv2+0.08333333333333333*Bmag[3]*f[17]*dv1*dv2+0.08333333333333333*Bmag[2]*f[10]*dv1*dv2))/m_; 
  out[3] += (volFact*(0.8*Bmag[3]*f[44]*wx1*wx2+0.8*Bmag[6]*f[20]*wx1*wx2+0.8944271909999161*Bmag[2]*f[20]*wx1*wx2+0.8*Bmag[7]*f[19]*wx1*wx2+0.8944271909999161*Bmag[1]*f[19]*wx1*wx2+0.8944271909999159*Bmag[3]*f[12]*wx1*wx2+0.8944271909999159*Bmag[3]*f[11]*wx1*wx2+0.8*f[5]*Bmag[8]*wx1*wx2+0.8944271909999161*f[2]*Bmag[7]*wx1*wx2+0.8944271909999161*f[1]*Bmag[6]*wx1*wx2+0.8944271909999159*Bmag[5]*f[5]*wx1*wx2+0.8944271909999159*Bmag[4]*f[5]*wx1*wx2+Bmag[0]*f[5]*wx1*wx2+f[0]*Bmag[3]*wx1*wx2+Bmag[1]*f[2]*wx1*wx2+f[1]*Bmag[2]*wx1*wx2+0.2309401076758504*Bmag[3]*f[54]*dv1*wx2+0.2309401076758503*Bmag[6]*f[33]*dv1*wx2+0.2581988897471612*Bmag[2]*f[33]*dv1*wx2+0.2309401076758503*Bmag[7]*f[32]*dv1*wx2+0.2581988897471612*Bmag[1]*f[32]*dv1*wx2+0.2581988897471611*Bmag[3]*f[22]*dv1*wx2+0.2581988897471611*Bmag[3]*f[21]*dv1*wx2+0.2309401076758504*Bmag[8]*f[15]*dv1*wx2+0.2581988897471612*Bmag[5]*f[15]*dv1*wx2+0.2581988897471612*Bmag[4]*f[15]*dv1*wx2+0.2886751345948129*Bmag[0]*f[15]*dv1*wx2+0.2581988897471611*Bmag[7]*f[7]*dv1*wx2+0.2886751345948129*Bmag[1]*f[7]*dv1*wx2+0.2581988897471611*Bmag[6]*f[6]*dv1*wx2+0.2886751345948129*Bmag[2]*f[6]*dv1*wx2+0.2886751345948129*Bmag[3]*f[3]*dv1*wx2+0.2309401076758504*Bmag[3]*f[57]*dv2*wx1+0.2309401076758503*Bmag[6]*f[36]*dv2*wx1+0.2581988897471612*Bmag[2]*f[36]*dv2*wx1+0.2309401076758503*Bmag[7]*f[35]*dv2*wx1+0.2581988897471612*Bmag[1]*f[35]*dv2*wx1+0.2581988897471611*Bmag[3]*f[26]*dv2*wx1+0.2581988897471611*Bmag[3]*f[25]*dv2*wx1+0.2309401076758504*Bmag[8]*f[16]*dv2*wx1+0.2581988897471612*Bmag[5]*f[16]*dv2*wx1+0.2581988897471612*Bmag[4]*f[16]*dv2*wx1+0.2886751345948129*Bmag[0]*f[16]*dv2*wx1+0.2581988897471611*Bmag[7]*f[9]*dv2*wx1+0.2886751345948129*Bmag[1]*f[9]*dv2*wx1+0.2581988897471611*Bmag[6]*f[8]*dv2*wx1+0.2886751345948129*Bmag[2]*f[8]*dv2*wx1+0.2886751345948129*Bmag[3]*f[4]*dv2*wx1+0.06666666666666667*Bmag[3]*f[66]*dv1*dv2+0.06666666666666667*Bmag[6]*f[51]*dv1*dv2+0.07453559924999302*Bmag[2]*f[51]*dv1*dv2+0.06666666666666667*Bmag[7]*f[50]*dv1*dv2+0.07453559924999302*Bmag[1]*f[50]*dv1*dv2+0.07453559924999298*Bmag[3]*f[38]*dv1*dv2+0.07453559924999298*Bmag[3]*f[37]*dv1*dv2+0.06666666666666667*Bmag[8]*f[31]*dv1*dv2+0.07453559924999298*Bmag[5]*f[31]*dv1*dv2+0.07453559924999298*Bmag[4]*f[31]*dv1*dv2+0.08333333333333333*Bmag[0]*f[31]*dv1*dv2+0.07453559924999302*Bmag[7]*f[18]*dv1*dv2+0.08333333333333333*Bmag[1]*f[18]*dv1*dv2+0.07453559924999302*Bmag[6]*f[17]*dv1*dv2+0.08333333333333333*Bmag[2]*f[17]*dv1*dv2+0.08333333333333333*Bmag[3]*f[10]*dv1*dv2))/m_; 
  out[4] += (volFact*(0.6388765649999399*Bmag[8]*f[44]*wx1*wx2+Bmag[5]*f[44]*wx1*wx2+0.8944271909999159*Bmag[7]*f[20]*wx1*wx2+0.6388765649999399*Bmag[6]*f[19]*wx1*wx2+1.0*Bmag[2]*f[19]*wx1*wx2+Bmag[8]*f[12]*wx1*wx2+0.6388765649999399*Bmag[4]*f[11]*wx1*wx2+Bmag[0]*f[11]*wx1*wx2+1.0*f[2]*Bmag[6]*wx1*wx2+0.8944271909999159*Bmag[3]*f[5]*wx1*wx2+f[0]*Bmag[4]*wx1*wx2+0.8944271909999159*Bmag[1]*f[1]*wx1*wx2+0.1844277783908294*Bmag[8]*f[54]*dv1*wx2+0.2886751345948129*Bmag[5]*f[54]*dv1*wx2+0.2581988897471611*Bmag[7]*f[33]*dv1*wx2+0.1844277783908294*Bmag[6]*f[32]*dv1*wx2+0.2886751345948129*Bmag[2]*f[32]*dv1*wx2+0.2886751345948129*Bmag[8]*f[22]*dv1*wx2+0.1844277783908294*Bmag[4]*f[21]*dv1*wx2+0.2886751345948129*Bmag[0]*f[21]*dv1*wx2+0.2581988897471612*Bmag[3]*f[15]*dv1*wx2+0.2886751345948129*Bmag[6]*f[7]*dv1*wx2+0.2581988897471612*Bmag[1]*f[6]*dv1*wx2+0.2886751345948129*f[3]*Bmag[4]*dv1*wx2+0.1844277783908294*Bmag[8]*f[57]*dv2*wx1+0.2886751345948129*Bmag[5]*f[57]*dv2*wx1+0.2581988897471611*Bmag[7]*f[36]*dv2*wx1+0.1844277783908294*Bmag[6]*f[35]*dv2*wx1+0.2886751345948129*Bmag[2]*f[35]*dv2*wx1+0.2886751345948129*Bmag[8]*f[26]*dv2*wx1+0.1844277783908294*Bmag[4]*f[25]*dv2*wx1+0.2886751345948129*Bmag[0]*f[25]*dv2*wx1+0.2581988897471612*Bmag[3]*f[16]*dv2*wx1+0.2886751345948129*Bmag[6]*f[9]*dv2*wx1+0.2581988897471612*Bmag[1]*f[8]*dv2*wx1+0.2886751345948129*Bmag[4]*f[4]*dv2*wx1+0.05323971374999499*Bmag[8]*f[66]*dv1*dv2+0.08333333333333333*Bmag[5]*f[66]*dv1*dv2+0.07453559924999298*Bmag[7]*f[51]*dv1*dv2+0.05323971374999499*Bmag[6]*f[50]*dv1*dv2+0.08333333333333336*Bmag[2]*f[50]*dv1*dv2+0.08333333333333333*Bmag[8]*f[38]*dv1*dv2+0.05323971374999499*Bmag[4]*f[37]*dv1*dv2+0.08333333333333333*Bmag[0]*f[37]*dv1*dv2+0.07453559924999298*Bmag[3]*f[31]*dv1*dv2+0.08333333333333336*Bmag[6]*f[18]*dv1*dv2+0.07453559924999298*Bmag[1]*f[17]*dv1*dv2+0.08333333333333333*Bmag[4]*f[10]*dv1*dv2))/m_; 
  out[5] += (volFact*(0.6388765649999399*Bmag[8]*f[44]*wx1*wx2+Bmag[4]*f[44]*wx1*wx2+0.6388765649999399*Bmag[7]*f[20]*wx1*wx2+1.0*Bmag[1]*f[20]*wx1*wx2+0.8944271909999159*Bmag[6]*f[19]*wx1*wx2+0.6388765649999399*Bmag[5]*f[12]*wx1*wx2+Bmag[0]*f[12]*wx1*wx2+Bmag[8]*f[11]*wx1*wx2+1.0*f[1]*Bmag[7]*wx1*wx2+0.8944271909999159*Bmag[3]*f[5]*wx1*wx2+f[0]*Bmag[5]*wx1*wx2+0.8944271909999159*Bmag[2]*f[2]*wx1*wx2+0.1844277783908294*Bmag[8]*f[54]*dv1*wx2+0.2886751345948129*Bmag[4]*f[54]*dv1*wx2+0.1844277783908294*Bmag[7]*f[33]*dv1*wx2+0.2886751345948129*Bmag[1]*f[33]*dv1*wx2+0.2581988897471611*Bmag[6]*f[32]*dv1*wx2+0.1844277783908294*Bmag[5]*f[22]*dv1*wx2+0.2886751345948129*Bmag[0]*f[22]*dv1*wx2+0.2886751345948129*Bmag[8]*f[21]*dv1*wx2+0.2581988897471612*Bmag[3]*f[15]*dv1*wx2+0.2581988897471612*Bmag[2]*f[7]*dv1*wx2+0.2886751345948129*f[6]*Bmag[7]*dv1*wx2+0.2886751345948129*f[3]*Bmag[5]*dv1*wx2+0.1844277783908294*Bmag[8]*f[57]*dv2*wx1+0.2886751345948129*Bmag[4]*f[57]*dv2*wx1+0.1844277783908294*Bmag[7]*f[36]*dv2*wx1+0.2886751345948129*Bmag[1]*f[36]*dv2*wx1+0.2581988897471611*Bmag[6]*f[35]*dv2*wx1+0.1844277783908294*Bmag[5]*f[26]*dv2*wx1+0.2886751345948129*Bmag[0]*f[26]*dv2*wx1+0.2886751345948129*Bmag[8]*f[25]*dv2*wx1+0.2581988897471612*Bmag[3]*f[16]*dv2*wx1+0.2581988897471612*Bmag[2]*f[9]*dv2*wx1+0.2886751345948129*Bmag[7]*f[8]*dv2*wx1+0.2886751345948129*f[4]*Bmag[5]*dv2*wx1+0.05323971374999499*Bmag[8]*f[66]*dv1*dv2+0.08333333333333333*Bmag[4]*f[66]*dv1*dv2+0.05323971374999499*Bmag[7]*f[51]*dv1*dv2+0.08333333333333336*Bmag[1]*f[51]*dv1*dv2+0.07453559924999298*Bmag[6]*f[50]*dv1*dv2+0.05323971374999499*Bmag[5]*f[38]*dv1*dv2+0.08333333333333333*Bmag[0]*f[38]*dv1*dv2+0.08333333333333333*Bmag[8]*f[37]*dv1*dv2+0.07453559924999298*Bmag[3]*f[31]*dv1*dv2+0.07453559924999298*Bmag[2]*f[18]*dv1*dv2+0.08333333333333336*Bmag[7]*f[17]*dv1*dv2+0.08333333333333333*Bmag[5]*f[10]*dv1*dv2))/m_; 
  out[6] += (volFact*(0.5714285714285714*Bmag[6]*f[44]*wx1*wx2+0.8944271909999161*Bmag[2]*f[44]*wx1*wx2+0.8*Bmag[3]*f[20]*wx1*wx2+0.5714285714285714*Bmag[8]*f[19]*wx1*wx2+0.8944271909999159*Bmag[5]*f[19]*wx1*wx2+0.6388765649999399*Bmag[4]*f[19]*wx1*wx2+Bmag[0]*f[19]*wx1*wx2+0.8944271909999159*Bmag[6]*f[12]*wx1*wx2+0.6388765649999399*Bmag[6]*f[11]*wx1*wx2+1.0*Bmag[2]*f[11]*wx1*wx2+0.8944271909999161*f[2]*Bmag[8]*wx1*wx2+0.8*f[5]*Bmag[7]*wx1*wx2+f[0]*Bmag[6]*wx1*wx2+0.8944271909999161*Bmag[1]*f[5]*wx1*wx2+1.0*f[2]*Bmag[4]*wx1*wx2+0.8944271909999161*f[1]*Bmag[3]*wx1*wx2+0.1649572197684645*Bmag[6]*f[54]*dv1*wx2+0.2581988897471611*Bmag[2]*f[54]*dv1*wx2+0.2309401076758503*Bmag[3]*f[33]*dv1*wx2+0.1649572197684645*Bmag[8]*f[32]*dv1*wx2+0.2581988897471611*Bmag[5]*f[32]*dv1*wx2+0.1844277783908294*Bmag[4]*f[32]*dv1*wx2+0.2886751345948129*Bmag[0]*f[32]*dv1*wx2+0.2581988897471611*Bmag[6]*f[22]*dv1*wx2+0.1844277783908294*Bmag[6]*f[21]*dv1*wx2+0.2886751345948129*Bmag[2]*f[21]*dv1*wx2+0.2309401076758504*Bmag[7]*f[15]*dv1*wx2+0.2581988897471611*Bmag[1]*f[15]*dv1*wx2+0.2581988897471611*f[7]*Bmag[8]*dv1*wx2+0.2886751345948129*Bmag[4]*f[7]*dv1*wx2+0.2581988897471611*Bmag[3]*f[6]*dv1*wx2+0.2886751345948129*f[3]*Bmag[6]*dv1*wx2+0.1649572197684645*Bmag[6]*f[57]*dv2*wx1+0.2581988897471611*Bmag[2]*f[57]*dv2*wx1+0.2309401076758503*Bmag[3]*f[36]*dv2*wx1+0.1649572197684645*Bmag[8]*f[35]*dv2*wx1+0.2581988897471611*Bmag[5]*f[35]*dv2*wx1+0.1844277783908294*Bmag[4]*f[35]*dv2*wx1+0.2886751345948129*Bmag[0]*f[35]*dv2*wx1+0.2581988897471611*Bmag[6]*f[26]*dv2*wx1+0.1844277783908294*Bmag[6]*f[25]*dv2*wx1+0.2886751345948129*Bmag[2]*f[25]*dv2*wx1+0.2309401076758504*Bmag[7]*f[16]*dv2*wx1+0.2581988897471611*Bmag[1]*f[16]*dv2*wx1+0.2581988897471611*Bmag[8]*f[9]*dv2*wx1+0.2886751345948129*Bmag[4]*f[9]*dv2*wx1+0.2581988897471611*Bmag[3]*f[8]*dv2*wx1+0.2886751345948129*f[4]*Bmag[6]*dv2*wx1+0.04761904761904762*Bmag[6]*f[66]*dv1*dv2+0.07453559924999302*Bmag[2]*f[66]*dv1*dv2+0.06666666666666667*Bmag[3]*f[51]*dv1*dv2+0.04761904761904762*Bmag[8]*f[50]*dv1*dv2+0.07453559924999298*Bmag[5]*f[50]*dv1*dv2+0.05323971374999499*Bmag[4]*f[50]*dv1*dv2+0.08333333333333333*Bmag[0]*f[50]*dv1*dv2+0.07453559924999298*Bmag[6]*f[38]*dv1*dv2+0.05323971374999499*Bmag[6]*f[37]*dv1*dv2+0.08333333333333336*Bmag[2]*f[37]*dv1*dv2+0.06666666666666667*Bmag[7]*f[31]*dv1*dv2+0.07453559924999302*Bmag[1]*f[31]*dv1*dv2+0.07453559924999302*Bmag[8]*f[18]*dv1*dv2+0.08333333333333336*Bmag[4]*f[18]*dv1*dv2+0.07453559924999302*Bmag[3]*f[17]*dv1*dv2+0.08333333333333333*Bmag[6]*f[10]*dv1*dv2))/m_; 
  out[7] += (volFact*(0.5714285714285714*Bmag[7]*f[44]*wx1*wx2+0.8944271909999161*Bmag[1]*f[44]*wx1*wx2+0.5714285714285714*Bmag[8]*f[20]*wx1*wx2+0.6388765649999399*Bmag[5]*f[20]*wx1*wx2+0.8944271909999159*Bmag[4]*f[20]*wx1*wx2+Bmag[0]*f[20]*wx1*wx2+0.8*Bmag[3]*f[19]*wx1*wx2+0.6388765649999399*Bmag[7]*f[12]*wx1*wx2+1.0*Bmag[1]*f[12]*wx1*wx2+0.8944271909999159*Bmag[7]*f[11]*wx1*wx2+0.8944271909999161*f[1]*Bmag[8]*wx1*wx2+f[0]*Bmag[7]*wx1*wx2+0.8*f[5]*Bmag[6]*wx1*wx2+0.8944271909999161*Bmag[2]*f[5]*wx1*wx2+1.0*f[1]*Bmag[5]*wx1*wx2+0.8944271909999161*f[2]*Bmag[3]*wx1*wx2+0.1649572197684645*Bmag[7]*f[54]*dv1*wx2+0.2581988897471611*Bmag[1]*f[54]*dv1*wx2+0.1649572197684645*Bmag[8]*f[33]*dv1*wx2+0.1844277783908294*Bmag[5]*f[33]*dv1*wx2+0.2581988897471611*Bmag[4]*f[33]*dv1*wx2+0.2886751345948129*Bmag[0]*f[33]*dv1*wx2+0.2309401076758503*Bmag[3]*f[32]*dv1*wx2+0.1844277783908294*Bmag[7]*f[22]*dv1*wx2+0.2886751345948129*Bmag[1]*f[22]*dv1*wx2+0.2581988897471611*Bmag[7]*f[21]*dv1*wx2+0.2309401076758504*Bmag[6]*f[15]*dv1*wx2+0.2581988897471611*Bmag[2]*f[15]*dv1*wx2+0.2581988897471611*f[6]*Bmag[8]*dv1*wx2+0.2581988897471611*Bmag[3]*f[7]*dv1*wx2+0.2886751345948129*f[3]*Bmag[7]*dv1*wx2+0.2886751345948129*Bmag[5]*f[6]*dv1*wx2+0.1649572197684645*Bmag[7]*f[57]*dv2*wx1+0.2581988897471611*Bmag[1]*f[57]*dv2*wx1+0.1649572197684645*Bmag[8]*f[36]*dv2*wx1+0.1844277783908294*Bmag[5]*f[36]*dv2*wx1+0.2581988897471611*Bmag[4]*f[36]*dv2*wx1+0.2886751345948129*Bmag[0]*f[36]*dv2*wx1+0.2309401076758503*Bmag[3]*f[35]*dv2*wx1+0.1844277783908294*Bmag[7]*f[26]*dv2*wx1+0.2886751345948129*Bmag[1]*f[26]*dv2*wx1+0.2581988897471611*Bmag[7]*f[25]*dv2*wx1+0.2309401076758504*Bmag[6]*f[16]*dv2*wx1+0.2581988897471611*Bmag[2]*f[16]*dv2*wx1+0.2581988897471611*Bmag[3]*f[9]*dv2*wx1+0.2581988897471611*Bmag[8]*f[8]*dv2*wx1+0.2886751345948129*Bmag[5]*f[8]*dv2*wx1+0.2886751345948129*f[4]*Bmag[7]*dv2*wx1+0.04761904761904762*Bmag[7]*f[66]*dv1*dv2+0.07453559924999302*Bmag[1]*f[66]*dv1*dv2+0.04761904761904762*Bmag[8]*f[51]*dv1*dv2+0.05323971374999499*Bmag[5]*f[51]*dv1*dv2+0.07453559924999298*Bmag[4]*f[51]*dv1*dv2+0.08333333333333333*Bmag[0]*f[51]*dv1*dv2+0.06666666666666667*Bmag[3]*f[50]*dv1*dv2+0.05323971374999499*Bmag[7]*f[38]*dv1*dv2+0.08333333333333336*Bmag[1]*f[38]*dv1*dv2+0.07453559924999298*Bmag[7]*f[37]*dv1*dv2+0.06666666666666667*Bmag[6]*f[31]*dv1*dv2+0.07453559924999302*Bmag[2]*f[31]*dv1*dv2+0.07453559924999302*Bmag[3]*f[18]*dv1*dv2+0.07453559924999302*Bmag[8]*f[17]*dv1*dv2+0.08333333333333336*Bmag[5]*f[17]*dv1*dv2+0.08333333333333333*Bmag[7]*f[10]*dv1*dv2))/m_; 
  out[8] += (volFact*(0.4081632653061225*Bmag[8]*f[44]*wx1*wx2+0.6388765649999399*Bmag[5]*f[44]*wx1*wx2+0.6388765649999399*Bmag[4]*f[44]*wx1*wx2+Bmag[0]*f[44]*wx1*wx2+0.5714285714285714*Bmag[7]*f[20]*wx1*wx2+0.8944271909999161*Bmag[1]*f[20]*wx1*wx2+0.5714285714285714*Bmag[6]*f[19]*wx1*wx2+0.8944271909999161*Bmag[2]*f[19]*wx1*wx2+0.6388765649999399*Bmag[8]*f[12]*wx1*wx2+Bmag[4]*f[12]*wx1*wx2+0.6388765649999399*Bmag[8]*f[11]*wx1*wx2+Bmag[5]*f[11]*wx1*wx2+f[0]*Bmag[8]*wx1*wx2+0.8944271909999161*f[1]*Bmag[7]*wx1*wx2+0.8944271909999161*f[2]*Bmag[6]*wx1*wx2+0.8*Bmag[3]*f[5]*wx1*wx2+0.1178265855489032*Bmag[8]*f[54]*dv1*wx2+0.1844277783908294*Bmag[5]*f[54]*dv1*wx2+0.1844277783908294*Bmag[4]*f[54]*dv1*wx2+0.2886751345948129*Bmag[0]*f[54]*dv1*wx2+0.1649572197684645*Bmag[7]*f[33]*dv1*wx2+0.2581988897471612*Bmag[1]*f[33]*dv1*wx2+0.1649572197684645*Bmag[6]*f[32]*dv1*wx2+0.2581988897471612*Bmag[2]*f[32]*dv1*wx2+0.1844277783908294*Bmag[8]*f[22]*dv1*wx2+0.2886751345948129*Bmag[4]*f[22]*dv1*wx2+0.1844277783908294*Bmag[8]*f[21]*dv1*wx2+0.2886751345948129*Bmag[5]*f[21]*dv1*wx2+0.2309401076758504*Bmag[3]*f[15]*dv1*wx2+0.2886751345948129*f[3]*Bmag[8]*dv1*wx2+0.2581988897471611*Bmag[6]*f[7]*dv1*wx2+0.2581988897471611*f[6]*Bmag[7]*dv1*wx2+0.1178265855489032*Bmag[8]*f[57]*dv2*wx1+0.1844277783908294*Bmag[5]*f[57]*dv2*wx1+0.1844277783908294*Bmag[4]*f[57]*dv2*wx1+0.2886751345948129*Bmag[0]*f[57]*dv2*wx1+0.1649572197684645*Bmag[7]*f[36]*dv2*wx1+0.2581988897471612*Bmag[1]*f[36]*dv2*wx1+0.1649572197684645*Bmag[6]*f[35]*dv2*wx1+0.2581988897471612*Bmag[2]*f[35]*dv2*wx1+0.1844277783908294*Bmag[8]*f[26]*dv2*wx1+0.2886751345948129*Bmag[4]*f[26]*dv2*wx1+0.1844277783908294*Bmag[8]*f[25]*dv2*wx1+0.2886751345948129*Bmag[5]*f[25]*dv2*wx1+0.2309401076758504*Bmag[3]*f[16]*dv2*wx1+0.2581988897471611*Bmag[6]*f[9]*dv2*wx1+0.2581988897471611*Bmag[7]*f[8]*dv2*wx1+0.2886751345948129*f[4]*Bmag[8]*dv2*wx1+0.03401360544217687*Bmag[8]*f[66]*dv1*dv2+0.05323971374999499*Bmag[5]*f[66]*dv1*dv2+0.05323971374999499*Bmag[4]*f[66]*dv1*dv2+0.08333333333333333*Bmag[0]*f[66]*dv1*dv2+0.04761904761904762*Bmag[7]*f[51]*dv1*dv2+0.07453559924999302*Bmag[1]*f[51]*dv1*dv2+0.04761904761904762*Bmag[6]*f[50]*dv1*dv2+0.07453559924999302*Bmag[2]*f[50]*dv1*dv2+0.05323971374999499*Bmag[8]*f[38]*dv1*dv2+0.08333333333333333*Bmag[4]*f[38]*dv1*dv2+0.05323971374999499*Bmag[8]*f[37]*dv1*dv2+0.08333333333333333*Bmag[5]*f[37]*dv1*dv2+0.06666666666666667*Bmag[3]*f[31]*dv1*dv2+0.07453559924999302*Bmag[6]*f[18]*dv1*dv2+0.07453559924999302*Bmag[7]*f[17]*dv1*dv2+0.08333333333333333*Bmag[8]*f[10]*dv1*dv2))/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
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
__host__ __device__ void GkMomentCalc2x2vTensor_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
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
  outM0[8] += 2.0*f[44]*volFact; 
  outM1[0] += 0.3333333333333333*volFact*(6.0*f[0]*wx1+1.732050807568877*f[3]*dv1); 
  outM1[1] += 0.3333333333333333*volFact*(6.0*f[1]*wx1+1.732050807568877*f[6]*dv1); 
  outM1[2] += 0.3333333333333333*volFact*(6.0*f[2]*wx1+1.732050807568877*f[7]*dv1); 
  outM1[3] += 0.3333333333333333*volFact*(6.0*f[5]*wx1+1.732050807568877*f[15]*dv1); 
  outM1[4] += 0.06666666666666667*volFact*(30.0*f[11]*wx1+8.660254037844387*f[21]*dv1); 
  outM1[5] += 0.06666666666666667*volFact*(30.0*f[12]*wx1+8.660254037844387*f[22]*dv1); 
  outM1[6] += 0.06666666666666667*volFact*(30.0*f[19]*wx1+8.660254037844387*f[32]*dv1); 
  outM1[7] += 0.06666666666666667*volFact*(30.0*f[20]*wx1+8.660254037844387*f[33]*dv1); 
  outM1[8] += 0.3333333333333333*volFact*(6.0*f[44]*wx1+1.732050807568877*f[54]*dv1); 
  outM2[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  outM2[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[23]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  outM2[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[24]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  outM2[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[15]*dv1*wx1+0.149071198499986*f[34]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
  outM2[4] += volFact*(2.0*f[11]*wx1_sq+1.154700538379251*f[21]*dv1*wx1+0.149071198499986*f[45]*dv1_sq+0.1666666666666667*f[11]*dv1_sq); 
  outM2[5] += volFact*(2.0*f[12]*wx1_sq+1.154700538379251*f[22]*dv1*wx1+0.149071198499986*f[46]*dv1_sq+0.1666666666666667*f[12]*dv1_sq); 
  outM2[6] += volFact*(2.0*f[19]*wx1_sq+1.154700538379251*f[32]*dv1*wx1+0.149071198499986*f[55]*dv1_sq+0.1666666666666667*f[19]*dv1_sq); 
  outM2[7] += volFact*(2.0*f[20]*wx1_sq+1.154700538379251*f[33]*dv1*wx1+0.149071198499986*f[56]*dv1_sq+0.1666666666666667*f[20]*dv1_sq); 
  outM2[8] += volFact*(2.0*f[44]*wx1_sq+1.154700538379252*f[54]*dv1*wx1+0.149071198499986*f[72]*dv1_sq+0.1666666666666667*f[44]*dv1_sq); 
  double tmp[9]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[16]*dv2; 
  tmp[4] = 2.0*f[11]*wx2+0.5773502691896257*f[25]*dv2; 
  tmp[5] = 2.0*f[12]*wx2+0.5773502691896257*f[26]*dv2; 
  tmp[6] = 2.0*f[19]*wx2+0.5773502691896257*f[35]*dv2; 
  tmp[7] = 2.0*f[20]*wx2+0.5773502691896257*f[36]*dv2; 
  tmp[8] = 2.0*f[44]*wx2+0.5773502691896258*f[57]*dv2; 
  outM2[0] += (2.0*(0.5*Bmag[8]*tmp[8]+0.5*Bmag[7]*tmp[7]+0.5*Bmag[6]*tmp[6]+0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.447213595499958*Bmag[7]*tmp[8]+0.447213595499958*tmp[7]*Bmag[8]+0.5000000000000001*Bmag[5]*tmp[7]+0.5000000000000001*tmp[5]*Bmag[7]+0.447213595499958*Bmag[3]*tmp[6]+0.447213595499958*tmp[3]*Bmag[6]+0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  outM2[2] += (2.0*(0.447213595499958*Bmag[6]*tmp[8]+0.447213595499958*tmp[6]*Bmag[8]+0.447213595499958*Bmag[3]*tmp[7]+0.447213595499958*tmp[3]*Bmag[7]+0.5000000000000001*Bmag[4]*tmp[6]+0.5000000000000001*tmp[4]*Bmag[6]+0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  outM2[3] += (2.0*(0.4*Bmag[3]*tmp[8]+0.4*tmp[3]*Bmag[8]+0.4*Bmag[6]*tmp[7]+0.447213595499958*Bmag[2]*tmp[7]+0.4*tmp[6]*Bmag[7]+0.447213595499958*tmp[2]*Bmag[7]+0.447213595499958*Bmag[1]*tmp[6]+0.447213595499958*tmp[1]*Bmag[6]+0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  outM2[4] += (2.0*(0.31943828249997*Bmag[8]*tmp[8]+0.5*Bmag[5]*tmp[8]+0.5*tmp[5]*Bmag[8]+0.4472135954999579*Bmag[7]*tmp[7]+0.31943828249997*Bmag[6]*tmp[6]+0.5000000000000001*Bmag[2]*tmp[6]+0.5000000000000001*tmp[2]*Bmag[6]+0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  outM2[5] += (2.0*(0.31943828249997*Bmag[8]*tmp[8]+0.5*Bmag[4]*tmp[8]+0.5*tmp[4]*Bmag[8]+0.31943828249997*Bmag[7]*tmp[7]+0.5000000000000001*Bmag[1]*tmp[7]+0.5000000000000001*tmp[1]*Bmag[7]+0.4472135954999579*Bmag[6]*tmp[6]+0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
  outM2[6] += (2.0*(0.2857142857142857*Bmag[6]*tmp[8]+0.447213595499958*Bmag[2]*tmp[8]+0.2857142857142857*tmp[6]*Bmag[8]+0.447213595499958*tmp[2]*Bmag[8]+0.4*Bmag[3]*tmp[7]+0.4*tmp[3]*Bmag[7]+0.4472135954999579*Bmag[5]*tmp[6]+0.31943828249997*Bmag[4]*tmp[6]+0.5*Bmag[0]*tmp[6]+0.4472135954999579*tmp[5]*Bmag[6]+0.31943828249997*tmp[4]*Bmag[6]+0.5*tmp[0]*Bmag[6]+0.5000000000000001*Bmag[2]*tmp[4]+0.5000000000000001*tmp[2]*Bmag[4]+0.447213595499958*Bmag[1]*tmp[3]+0.447213595499958*tmp[1]*Bmag[3])*volFact)/m_; 
  outM2[7] += (2.0*(0.2857142857142857*Bmag[7]*tmp[8]+0.447213595499958*Bmag[1]*tmp[8]+0.2857142857142857*tmp[7]*Bmag[8]+0.447213595499958*tmp[1]*Bmag[8]+0.31943828249997*Bmag[5]*tmp[7]+0.4472135954999579*Bmag[4]*tmp[7]+0.5*Bmag[0]*tmp[7]+0.31943828249997*tmp[5]*Bmag[7]+0.4472135954999579*tmp[4]*Bmag[7]+0.5*tmp[0]*Bmag[7]+0.4*Bmag[3]*tmp[6]+0.4*tmp[3]*Bmag[6]+0.5000000000000001*Bmag[1]*tmp[5]+0.5000000000000001*tmp[1]*Bmag[5]+0.447213595499958*Bmag[2]*tmp[3]+0.447213595499958*tmp[2]*Bmag[3])*volFact)/m_; 
  outM2[8] += (2.0*(0.2040816326530612*Bmag[8]*tmp[8]+0.31943828249997*Bmag[5]*tmp[8]+0.31943828249997*Bmag[4]*tmp[8]+0.5*Bmag[0]*tmp[8]+0.31943828249997*tmp[5]*Bmag[8]+0.31943828249997*tmp[4]*Bmag[8]+0.5*tmp[0]*Bmag[8]+0.2857142857142857*Bmag[7]*tmp[7]+0.447213595499958*Bmag[1]*tmp[7]+0.447213595499958*tmp[1]*Bmag[7]+0.2857142857142857*Bmag[6]*tmp[6]+0.447213595499958*Bmag[2]*tmp[6]+0.447213595499958*tmp[2]*Bmag[6]+0.5*Bmag[4]*tmp[5]+0.5*tmp[4]*Bmag[5]+0.4*Bmag[3]*tmp[3])*volFact)/m_; 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M0_step1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
__host__ __device__ void GkMomentCalc2x2vTensor_M0_step1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
  out[20] += 1.414213562373095*f[44]*volFact; 
  out[21] += 1.414213562373095*f[47]*volFact; 
  out[22] += 1.414213562373095*f[48]*volFact; 
  out[23] += 1.414213562373095*f[57]*volFact; 
  out[24] += 1.414213562373095*f[60]*volFact; 
  out[25] += 1.414213562373095*f[61]*volFact; 
  out[26] += 1.414213562373095*f[73]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M0_step2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]/2; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[4]*volFact; 
} 
__host__ __device__ void GkMomentCalc2x2vTensor_M0_step2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
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
  out[8] += 2.828427124746191*f[20]*volFact; 
} 
