#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment1x3vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x3vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x3vSer_v1_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x3vSer_v2_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 

  out[0] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[3]*dv2*volFac; 

} 
void IntDGMoment1x3vSer_v2_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 

  out[0] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[3]*dv2*volFac; 

} 
void IntDGMoment1x3vSer_v2_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 

  out[0] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[3]*dv2*volFac; 

} 
void IntDGMoment1x3vSer_v3_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv3+1.154700538379252*fld[4]*dv3*volFac; 

} 
void IntDGMoment1x3vSer_v3_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv3+1.154700538379252*fld[4]*dv3*volFac; 

} 
void IntDGMoment1x3vSer_v3_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv3+1.154700538379252*fld[4]*dv3*volFac; 

} 
void IntDGMoment1x3vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[2]*dv1*volFac*wv1+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x3vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[2]*dv1*volFac*wv1+0.2981423969999719*fld[12]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x3vSer_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[2]*dv1*volFac*wv1+0.2981423969999719*fld[12]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x3vSer_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[3]*dv2*volFac*wv2+0.3333333333333333*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment1x3vSer_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[3]*dv2*volFac*wv2+0.2981423969999719*fld[13]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment1x3vSer_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[3]*dv2*volFac*wv2+0.2981423969999719*fld[13]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment1x3vSer_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 
  const double wv3R2 = std::pow(wv3,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 4.0*fld[0]*volFac*wv3R2+2.309401076758503*fld[4]*dv3*volFac*wv3+0.3333333333333333*fld[0]*dv3R2*volFac; 

} 
void IntDGMoment1x3vSer_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 
  const double wv3R2 = std::pow(wv3,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 4.0*fld[0]*volFac*wv3R2+2.309401076758503*fld[4]*dv3*volFac*wv3+0.2981423969999719*fld[14]*dv3R2*volFac+0.3333333333333333*fld[0]*dv3R2*volFac; 

} 
void IntDGMoment1x3vSer_v3Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 
  const double wv3R2 = std::pow(wv3,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 4.0*fld[0]*volFac*wv3R2+2.309401076758503*fld[4]*dv3*volFac*wv3+0.2981423969999719*fld[14]*dv3R2*volFac+0.3333333333333333*fld[0]*dv3R2*volFac; 

} 
void IntDGMoment1x3vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[2]*dv1*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[3]*dv2*volFac; 
  out[2] += 4.0*fld[0]*volFac*wv3+1.154700538379252*fld[4]*dv3*volFac; 

} 
void IntDGMoment1x3vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[2]*dv1*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[3]*dv2*volFac; 
  out[2] += 4.0*fld[0]*volFac*wv3+1.154700538379252*fld[4]*dv3*volFac; 

} 
void IntDGMoment1x3vSer_vi_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[2]*dv1*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[3]*dv2*volFac; 
  out[2] += 4.0*fld[0]*volFac*wv3+1.154700538379252*fld[4]*dv3*volFac; 

} 
void IntDGMoment1x3vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 4.0*fld[0]*volFac*wv3R2+2.309401076758503*fld[4]*dv3*volFac*wv3+4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[3]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[2]*dv1*volFac*wv1+0.3333333333333333*fld[0]*dv3R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x3vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 4.0*fld[0]*volFac*wv3R2+2.309401076758503*fld[4]*dv3*volFac*wv3+4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[3]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[2]*dv1*volFac*wv1+0.2981423969999719*fld[14]*dv3R2*volFac+0.3333333333333333*fld[0]*dv3R2*volFac+0.2981423969999719*fld[13]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac+0.2981423969999719*fld[12]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x3vSer_vSq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 4.0*fld[0]*volFac*wv3R2+2.309401076758503*fld[4]*dv3*volFac*wv3+4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[3]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[2]*dv1*volFac*wv1+0.2981423969999719*fld[14]*dv3R2*volFac+0.3333333333333333*fld[0]*dv3R2*volFac+0.2981423969999719*fld[13]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac+0.2981423969999719*fld[12]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x3vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 4.0*fld[0]*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[2]*dv1*volFac; 
  out[2] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[3]*dv2*volFac; 
  out[3] += 4.0*fld[0]*volFac*wv3+1.154700538379252*fld[4]*dv3*volFac; 
  out[4] += 4.0*fld[0]*volFac*wv3R2+2.309401076758503*fld[4]*dv3*volFac*wv3+4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[3]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[2]*dv1*volFac*wv1+0.3333333333333333*fld[0]*dv3R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x3vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 4.0*fld[0]*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[2]*dv1*volFac; 
  out[2] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[3]*dv2*volFac; 
  out[3] += 4.0*fld[0]*volFac*wv3+1.154700538379252*fld[4]*dv3*volFac; 
  out[4] += 4.0*fld[0]*volFac*wv3R2+2.309401076758503*fld[4]*dv3*volFac*wv3+4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[3]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[2]*dv1*volFac*wv1+0.2981423969999719*fld[14]*dv3R2*volFac+0.3333333333333333*fld[0]*dv3R2*volFac+0.2981423969999719*fld[13]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac+0.2981423969999719*fld[12]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x3vSer_intM_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
  const double wv3 = w[3]; 
  const double dv3 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 4.0*fld[0]*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[2]*dv1*volFac; 
  out[2] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[3]*dv2*volFac; 
  out[3] += 4.0*fld[0]*volFac*wv3+1.154700538379252*fld[4]*dv3*volFac; 
  out[4] += 4.0*fld[0]*volFac*wv3R2+2.309401076758503*fld[4]*dv3*volFac*wv3+4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[3]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[2]*dv1*volFac*wv1+0.2981423969999719*fld[14]*dv3R2*volFac+0.3333333333333333*fld[0]*dv3R2*volFac+0.2981423969999719*fld[13]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac+0.2981423969999719*fld[12]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
