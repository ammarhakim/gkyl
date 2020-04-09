#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment1x1vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 2.0*fld[0]*volFac*wv1+0.5773502691896258*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x1vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 2.0*fld[0]*volFac*wv1+0.5773502691896258*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x1vSer_v1_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 2.0*fld[0]*volFac*wv1+0.5773502691896258*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x1vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.0*fld[0]*volFac*wv1R2+1.154700538379252*fld[2]*dv1*volFac*wv1+0.1666666666666667*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x1vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.0*fld[0]*volFac*wv1R2+1.154700538379252*fld[2]*dv1*volFac*wv1+0.149071198499986*fld[5]*dv1R2*volFac+0.1666666666666667*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x1vSer_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.0*fld[0]*volFac*wv1R2+1.154700538379252*fld[2]*dv1*volFac*wv1+0.149071198499986*fld[5]*dv1R2*volFac+0.1666666666666667*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x1vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 2.0*fld[0]*volFac*wv1+0.5773502691896258*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x1vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 2.0*fld[0]*volFac*wv1+0.5773502691896258*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x1vSer_vi_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 2.0*fld[0]*volFac*wv1+0.5773502691896258*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x1vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.0*fld[0]*volFac*wv1R2+1.154700538379252*fld[2]*dv1*volFac*wv1+0.1666666666666667*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x1vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.0*fld[0]*volFac*wv1R2+1.154700538379252*fld[2]*dv1*volFac*wv1+0.149071198499986*fld[5]*dv1R2*volFac+0.1666666666666667*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x1vSer_vSq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.0*fld[0]*volFac*wv1R2+1.154700538379252*fld[2]*dv1*volFac*wv1+0.149071198499986*fld[5]*dv1R2*volFac+0.1666666666666667*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x1vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.0*fld[0]*volFac; 
  out[1] += 2.0*fld[0]*volFac*wv1+0.5773502691896258*fld[2]*dv1*volFac; 
  out[2] += 2.0*fld[0]*volFac*wv1R2+1.154700538379252*fld[2]*dv1*volFac*wv1+0.1666666666666667*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x1vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.0*fld[0]*volFac; 
  out[1] += 2.0*fld[0]*volFac*wv1+0.5773502691896258*fld[2]*dv1*volFac; 
  out[2] += 2.0*fld[0]*volFac*wv1R2+1.154700538379252*fld[2]*dv1*volFac*wv1+0.149071198499986*fld[5]*dv1R2*volFac+0.1666666666666667*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x1vSer_intM_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.0*fld[0]*volFac; 
  out[1] += 2.0*fld[0]*volFac*wv1+0.5773502691896258*fld[2]*dv1*volFac; 
  out[2] += 2.0*fld[0]*volFac*wv1R2+1.154700538379252*fld[2]*dv1*volFac*wv1+0.149071198499986*fld[5]*dv1R2*volFac+0.1666666666666667*fld[0]*dv1R2*volFac; 

} 
