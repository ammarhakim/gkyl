#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment2x2vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[3]*dv1*volFac; 

} 
void IntDGMoment2x2vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[3]*dv1*volFac; 

} 
void IntDGMoment2x2vSer_v1_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[3]*dv1*volFac; 

} 
void IntDGMoment2x2vSer_v2_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[4]*dv2*volFac; 

} 
void IntDGMoment2x2vSer_v2_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[4]*dv2*volFac; 

} 
void IntDGMoment2x2vSer_v2_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[4]*dv2*volFac; 

} 
void IntDGMoment2x2vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[3]*dv1*volFac*wv1+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x2vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[3]*dv1*volFac*wv1+0.2981423969999719*fld[13]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x2vSer_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[3]*dv1*volFac*wv1+0.2981423969999719*fld[13]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x2vSer_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[4]*dv2*volFac*wv2+0.3333333333333333*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment2x2vSer_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[4]*dv2*volFac*wv2+0.2981423969999719*fld[14]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment2x2vSer_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[4]*dv2*volFac*wv2+0.2981423969999719*fld[14]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment2x2vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[3]*dv1*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[4]*dv2*volFac; 

} 
void IntDGMoment2x2vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[3]*dv1*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[4]*dv2*volFac; 

} 
void IntDGMoment2x2vSer_vi_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[3]*dv1*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[4]*dv2*volFac; 

} 
void IntDGMoment2x2vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[4]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[3]*dv1*volFac*wv1+0.3333333333333333*fld[0]*dv2R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x2vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[4]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[3]*dv1*volFac*wv1+0.2981423969999719*fld[14]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac+0.2981423969999719*fld[13]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x2vSer_vSq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[4]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[3]*dv1*volFac*wv1+0.2981423969999719*fld[14]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac+0.2981423969999719*fld[13]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x2vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[3]*dv1*volFac; 
  out[2] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[4]*dv2*volFac; 
  out[3] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[4]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[3]*dv1*volFac*wv1+0.3333333333333333*fld[0]*dv2R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x2vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[3]*dv1*volFac; 
  out[2] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[4]*dv2*volFac; 
  out[3] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[4]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[3]*dv1*volFac*wv1+0.2981423969999719*fld[14]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac+0.2981423969999719*fld[13]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x2vSer_intM_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 4.0*fld[0]*volFac; 
  out[1] += 4.0*fld[0]*volFac*wv1+1.154700538379252*fld[3]*dv1*volFac; 
  out[2] += 4.0*fld[0]*volFac*wv2+1.154700538379252*fld[4]*dv2*volFac; 
  out[3] += 4.0*fld[0]*volFac*wv2R2+2.309401076758503*fld[4]*dv2*volFac*wv2+4.0*fld[0]*volFac*wv1R2+2.309401076758503*fld[3]*dv1*volFac*wv1+0.2981423969999719*fld[14]*dv2R2*volFac+0.3333333333333333*fld[0]*dv2R2*volFac+0.2981423969999719*fld[13]*dv1R2*volFac+0.3333333333333333*fld[0]*dv1R2*volFac; 

} 
