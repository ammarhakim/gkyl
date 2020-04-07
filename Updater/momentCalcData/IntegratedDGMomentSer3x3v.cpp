#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment3x3vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv1 = w[3]; 
  const double dv1 = dx[3]; 
 

  out[0] += 8.0*fld[0]*volFac*wv1+2.309401076758503*fld[4]*dv1*volFac; 

} 
void IntDGMoment3x3vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv1 = w[3]; 
  const double dv1 = dx[3]; 
 

  out[0] += 8.0*fld[0]*volFac*wv1+2.309401076758503*fld[4]*dv1*volFac; 

} 
void IntDGMoment3x3vSer_v2_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv2 = w[4]; 
  const double dv2 = dx[4]; 
 

  out[0] += 8.0*fld[0]*volFac*wv2+2.309401076758503*fld[5]*dv2*volFac; 

} 
void IntDGMoment3x3vSer_v2_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv2 = w[4]; 
  const double dv2 = dx[4]; 
 

  out[0] += 8.0*fld[0]*volFac*wv2+2.309401076758503*fld[5]*dv2*volFac; 

} 
void IntDGMoment3x3vSer_v3_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv3 = w[5]; 
  const double dv3 = dx[5]; 
 

  out[0] += 8.0*fld[0]*volFac*wv3+2.309401076758503*fld[6]*dv3*volFac; 

} 
void IntDGMoment3x3vSer_v3_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv3 = w[5]; 
  const double dv3 = dx[5]; 
 

  out[0] += 8.0*fld[0]*volFac*wv3+2.309401076758503*fld[6]*dv3*volFac; 

} 
void IntDGMoment3x3vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv1 = w[3]; 
  const double dv1 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 8.0*fld[0]*volFac*wv1R2+4.618802153517007*fld[4]*dv1*volFac*wv1+0.6666666666666666*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment3x3vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv1 = w[3]; 
  const double dv1 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 8.0*fld[0]*volFac*wv1R2+4.618802153517007*fld[4]*dv1*volFac*wv1+0.5962847939999438*fld[25]*dv1R2*volFac+0.6666666666666666*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment3x3vSer_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv2 = w[4]; 
  const double dv2 = dx[4]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 8.0*fld[0]*volFac*wv2R2+4.618802153517007*fld[5]*dv2*volFac*wv2+0.6666666666666666*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment3x3vSer_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv2 = w[4]; 
  const double dv2 = dx[4]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 8.0*fld[0]*volFac*wv2R2+4.618802153517007*fld[5]*dv2*volFac*wv2+0.5962847939999438*fld[26]*dv2R2*volFac+0.6666666666666666*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment3x3vSer_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv3 = w[5]; 
  const double dv3 = dx[5]; 
 
  const double wv3R2 = std::pow(wv3,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 8.0*fld[0]*volFac*wv3R2+4.618802153517007*fld[6]*dv3*volFac*wv3+0.6666666666666666*fld[0]*dv3R2*volFac; 

} 
void IntDGMoment3x3vSer_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv3 = w[5]; 
  const double dv3 = dx[5]; 
 
  const double wv3R2 = std::pow(wv3,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 8.0*fld[0]*volFac*wv3R2+4.618802153517007*fld[6]*dv3*volFac*wv3+0.5962847939999438*fld[27]*dv3R2*volFac+0.6666666666666666*fld[0]*dv3R2*volFac; 

} 
void IntDGMoment3x3vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv1 = w[3]; 
  const double dv1 = dx[3]; 
  const double wv2 = w[4]; 
  const double dv2 = dx[4]; 
  const double wv3 = w[5]; 
  const double dv3 = dx[5]; 
 

  out[0] += 8.0*fld[0]*volFac*wv1+2.309401076758503*fld[4]*dv1*volFac; 
  out[1] += 8.0*fld[0]*volFac*wv2+2.309401076758503*fld[5]*dv2*volFac; 
  out[2] += 8.0*fld[0]*volFac*wv3+2.309401076758503*fld[6]*dv3*volFac; 

} 
void IntDGMoment3x3vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv1 = w[3]; 
  const double dv1 = dx[3]; 
  const double wv2 = w[4]; 
  const double dv2 = dx[4]; 
  const double wv3 = w[5]; 
  const double dv3 = dx[5]; 
 

  out[0] += 8.0*fld[0]*volFac*wv1+2.309401076758503*fld[4]*dv1*volFac; 
  out[1] += 8.0*fld[0]*volFac*wv2+2.309401076758503*fld[5]*dv2*volFac; 
  out[2] += 8.0*fld[0]*volFac*wv3+2.309401076758503*fld[6]*dv3*volFac; 

} 
void IntDGMoment3x3vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv1 = w[3]; 
  const double dv1 = dx[3]; 
  const double wv2 = w[4]; 
  const double dv2 = dx[4]; 
  const double wv3 = w[5]; 
  const double dv3 = dx[5]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 8.0*fld[0]*volFac*wv3R2+4.618802153517007*fld[6]*dv3*volFac*wv3+8.0*fld[0]*volFac*wv2R2+4.618802153517007*fld[5]*dv2*volFac*wv2+8.0*fld[0]*volFac*wv1R2+4.618802153517007*fld[4]*dv1*volFac*wv1+0.6666666666666666*fld[0]*dv3R2*volFac+0.6666666666666666*fld[0]*dv2R2*volFac+0.6666666666666666*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment3x3vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv1 = w[3]; 
  const double dv1 = dx[3]; 
  const double wv2 = w[4]; 
  const double dv2 = dx[4]; 
  const double wv3 = w[5]; 
  const double dv3 = dx[5]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 8.0*fld[0]*volFac*wv3R2+4.618802153517007*fld[6]*dv3*volFac*wv3+8.0*fld[0]*volFac*wv2R2+4.618802153517007*fld[5]*dv2*volFac*wv2+8.0*fld[0]*volFac*wv1R2+4.618802153517007*fld[4]*dv1*volFac*wv1+0.5962847939999438*fld[27]*dv3R2*volFac+0.6666666666666666*fld[0]*dv3R2*volFac+0.5962847939999438*fld[26]*dv2R2*volFac+0.6666666666666666*fld[0]*dv2R2*volFac+0.5962847939999438*fld[25]*dv1R2*volFac+0.6666666666666666*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment3x3vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv1 = w[3]; 
  const double dv1 = dx[3]; 
  const double wv2 = w[4]; 
  const double dv2 = dx[4]; 
  const double wv3 = w[5]; 
  const double dv3 = dx[5]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 8.0*fld[0]*volFac; 
  out[1] += 8.0*fld[0]*volFac*wv1+2.309401076758503*fld[4]*dv1*volFac; 
  out[2] += 8.0*fld[0]*volFac*wv2+2.309401076758503*fld[5]*dv2*volFac; 
  out[3] += 8.0*fld[0]*volFac*wv3+2.309401076758503*fld[6]*dv3*volFac; 
  out[4] += 8.0*fld[0]*volFac*wv3R2+4.618802153517007*fld[6]*dv3*volFac*wv3+8.0*fld[0]*volFac*wv2R2+4.618802153517007*fld[5]*dv2*volFac*wv2+8.0*fld[0]*volFac*wv1R2+4.618802153517007*fld[4]*dv1*volFac*wv1+0.6666666666666666*fld[0]*dv3R2*volFac+0.6666666666666666*fld[0]*dv2R2*volFac+0.6666666666666666*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment3x3vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wv1 = w[3]; 
  const double dv1 = dx[3]; 
  const double wv2 = w[4]; 
  const double dv2 = dx[4]; 
  const double wv3 = w[5]; 
  const double dv3 = dx[5]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 8.0*fld[0]*volFac; 
  out[1] += 8.0*fld[0]*volFac*wv1+2.309401076758503*fld[4]*dv1*volFac; 
  out[2] += 8.0*fld[0]*volFac*wv2+2.309401076758503*fld[5]*dv2*volFac; 
  out[3] += 8.0*fld[0]*volFac*wv3+2.309401076758503*fld[6]*dv3*volFac; 
  out[4] += 8.0*fld[0]*volFac*wv3R2+4.618802153517007*fld[6]*dv3*volFac*wv3+8.0*fld[0]*volFac*wv2R2+4.618802153517007*fld[5]*dv2*volFac*wv2+8.0*fld[0]*volFac*wv1R2+4.618802153517007*fld[4]*dv1*volFac*wv1+0.5962847939999438*fld[27]*dv3R2*volFac+0.6666666666666666*fld[0]*dv3R2*volFac+0.5962847939999438*fld[26]*dv2R2*volFac+0.6666666666666666*fld[0]*dv2R2*volFac+0.5962847939999438*fld[25]*dv1R2*volFac+0.6666666666666666*fld[0]*dv1R2*volFac; 

} 
