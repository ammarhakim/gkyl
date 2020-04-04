#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment2xSer_one_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
 
  out[0] += 2.0*fld[0]*volFac; 

} 
void IntDGMoment2xSer_one_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
 
  out[0] += 2.0*fld[0]*volFac; 

} 
void IntDGMoment2xSer_one_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
 
  out[0] += 2.0*fld[0]*volFac; 

} 
void IntDGMoment2xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 2.0*fld[0]*volFac*wx1+0.5773502691896258*fld[1]*dx1*volFac; 

} 
void IntDGMoment2xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 2.0*fld[0]*volFac*wx1+0.5773502691896258*fld[1]*dx1*volFac; 

} 
void IntDGMoment2xSer_x1_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 2.0*fld[0]*volFac*wx1+0.5773502691896258*fld[1]*dx1*volFac; 

} 
void IntDGMoment2xSer_x2_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 2.0*fld[0]*volFac*wx2+0.5773502691896258*fld[2]*dx2*volFac; 

} 
void IntDGMoment2xSer_x2_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 2.0*fld[0]*volFac*wx2+0.5773502691896258*fld[2]*dx2*volFac; 

} 
void IntDGMoment2xSer_x2_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 2.0*fld[0]*volFac*wx2+0.5773502691896258*fld[2]*dx2*volFac; 

} 
void IntDGMoment2xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 2.0*fld[0]*volFac*wx1R2+1.154700538379252*fld[1]*dx1*volFac*wx1+0.1666666666666667*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment2xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 2.0*fld[0]*volFac*wx1R2+1.154700538379252*fld[1]*dx1*volFac*wx1+0.149071198499986*fld[4]*dx1R2*volFac+0.1666666666666667*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment2xSer_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 2.0*fld[0]*volFac*wx1R2+1.154700538379252*fld[1]*dx1*volFac*wx1+0.149071198499986*fld[4]*dx1R2*volFac+0.1666666666666667*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment2xSer_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 2.0*fld[0]*volFac*wx2R2+1.154700538379252*fld[2]*dx2*volFac*wx2+0.1666666666666667*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment2xSer_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 2.0*fld[0]*volFac*wx2R2+1.154700538379252*fld[2]*dx2*volFac*wx2+0.149071198499986*fld[5]*dx2R2*volFac+0.1666666666666667*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment2xSer_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 2.0*fld[0]*volFac*wx2R2+1.154700538379252*fld[2]*dx2*volFac*wx2+0.149071198499986*fld[5]*dx2R2*volFac+0.1666666666666667*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment2xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 2.0*fld[0]*volFac*wx2R2+1.154700538379252*fld[2]*dx2*volFac*wx2+2.0*fld[0]*volFac*wx1R2+1.154700538379252*fld[1]*dx1*volFac*wx1+0.1666666666666667*fld[0]*dx2R2*volFac+0.1666666666666667*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment2xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 2.0*fld[0]*volFac*wx2R2+1.154700538379252*fld[2]*dx2*volFac*wx2+2.0*fld[0]*volFac*wx1R2+1.154700538379252*fld[1]*dx1*volFac*wx1+0.149071198499986*fld[5]*dx2R2*volFac+0.1666666666666667*fld[0]*dx2R2*volFac+0.149071198499986*fld[4]*dx1R2*volFac+0.1666666666666667*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment2xSer_xSq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*0.25; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 2.0*fld[0]*volFac*wx2R2+1.154700538379252*fld[2]*dx2*volFac*wx2+2.0*fld[0]*volFac*wx1R2+1.154700538379252*fld[1]*dx1*volFac*wx1+0.149071198499986*fld[5]*dx2R2*volFac+0.1666666666666667*fld[0]*dx2R2*volFac+0.149071198499986*fld[4]*dx1R2*volFac+0.1666666666666667*fld[0]*dx1R2*volFac; 

} 
