#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment6xSer_one_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
 
  out[0] += 8.0*fld[0]*volFac; 

} 
void IntDGMoment6xSer_one_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
 
  out[0] += 8.0*fld[0]*volFac; 

} 
void IntDGMoment6xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 8.0*fld[0]*volFac*wx1+2.309401076758503*fld[1]*dx1*volFac; 

} 
void IntDGMoment6xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 8.0*fld[0]*volFac*wx1+2.309401076758503*fld[1]*dx1*volFac; 

} 
void IntDGMoment6xSer_x2_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 8.0*fld[0]*volFac*wx2+2.309401076758503*fld[2]*dx2*volFac; 

} 
void IntDGMoment6xSer_x2_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 8.0*fld[0]*volFac*wx2+2.309401076758503*fld[2]*dx2*volFac; 

} 
void IntDGMoment6xSer_x3_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 8.0*fld[0]*volFac*wx3+2.309401076758503*fld[3]*dx3*volFac; 

} 
void IntDGMoment6xSer_x3_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 8.0*fld[0]*volFac*wx3+2.309401076758503*fld[3]*dx3*volFac; 

} 
void IntDGMoment6xSer_x4_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 

  out[0] += 8.0*fld[0]*volFac*wx4+2.309401076758503*fld[4]*dx4*volFac; 

} 
void IntDGMoment6xSer_x4_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 

  out[0] += 8.0*fld[0]*volFac*wx4+2.309401076758503*fld[4]*dx4*volFac; 

} 
void IntDGMoment6xSer_x5_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx5 = w[4]; 
  const double dx5 = dx[4]; 
 

  out[0] += 8.0*fld[0]*volFac*wx5+2.309401076758503*fld[5]*dx5*volFac; 

} 
void IntDGMoment6xSer_x5_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx5 = w[4]; 
  const double dx5 = dx[4]; 
 

  out[0] += 8.0*fld[0]*volFac*wx5+2.309401076758503*fld[5]*dx5*volFac; 

} 
void IntDGMoment6xSer_x6_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx6 = w[5]; 
  const double dx6 = dx[5]; 
 

  out[0] += 8.0*fld[0]*volFac*wx6+2.309401076758503*fld[6]*dx6*volFac; 

} 
void IntDGMoment6xSer_x6_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx6 = w[5]; 
  const double dx6 = dx[5]; 
 

  out[0] += 8.0*fld[0]*volFac*wx6+2.309401076758503*fld[6]*dx6*volFac; 

} 
void IntDGMoment6xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 8.0*fld[0]*volFac*wx1R2+4.618802153517007*fld[1]*dx1*volFac*wx1+0.6666666666666666*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment6xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 8.0*fld[0]*volFac*wx1R2+4.618802153517007*fld[1]*dx1*volFac*wx1+0.5962847939999438*fld[22]*dx1R2*volFac+0.6666666666666666*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment6xSer_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 8.0*fld[0]*volFac*wx2R2+4.618802153517007*fld[2]*dx2*volFac*wx2+0.6666666666666666*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment6xSer_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 8.0*fld[0]*volFac*wx2R2+4.618802153517007*fld[2]*dx2*volFac*wx2+0.5962847939999438*fld[23]*dx2R2*volFac+0.6666666666666666*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment6xSer_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx3R2 = std::pow(wx3,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 8.0*fld[0]*volFac*wx3R2+4.618802153517007*fld[3]*dx3*volFac*wx3+0.6666666666666666*fld[0]*dx3R2*volFac; 

} 
void IntDGMoment6xSer_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx3R2 = std::pow(wx3,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 8.0*fld[0]*volFac*wx3R2+4.618802153517007*fld[3]*dx3*volFac*wx3+0.5962847939999438*fld[24]*dx3R2*volFac+0.6666666666666666*fld[0]*dx3R2*volFac; 

} 
void IntDGMoment6xSer_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 
  const double wx4R2 = std::pow(wx4,2);
  const double dx4R2 = std::pow(dx4,2);

  out[0] += 8.0*fld[0]*volFac*wx4R2+4.618802153517007*fld[4]*dx4*volFac*wx4+0.6666666666666666*fld[0]*dx4R2*volFac; 

} 
void IntDGMoment6xSer_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 
  const double wx4R2 = std::pow(wx4,2);
  const double dx4R2 = std::pow(dx4,2);

  out[0] += 8.0*fld[0]*volFac*wx4R2+4.618802153517007*fld[4]*dx4*volFac*wx4+0.5962847939999438*fld[25]*dx4R2*volFac+0.6666666666666666*fld[0]*dx4R2*volFac; 

} 
void IntDGMoment6xSer_x5Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx5 = w[4]; 
  const double dx5 = dx[4]; 
 
  const double wx5R2 = std::pow(wx5,2);
  const double dx5R2 = std::pow(dx5,2);

  out[0] += 8.0*fld[0]*volFac*wx5R2+4.618802153517007*fld[5]*dx5*volFac*wx5+0.6666666666666666*fld[0]*dx5R2*volFac; 

} 
void IntDGMoment6xSer_x5Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx5 = w[4]; 
  const double dx5 = dx[4]; 
 
  const double wx5R2 = std::pow(wx5,2);
  const double dx5R2 = std::pow(dx5,2);

  out[0] += 8.0*fld[0]*volFac*wx5R2+4.618802153517007*fld[5]*dx5*volFac*wx5+0.5962847939999438*fld[26]*dx5R2*volFac+0.6666666666666666*fld[0]*dx5R2*volFac; 

} 
void IntDGMoment6xSer_x6Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx6 = w[5]; 
  const double dx6 = dx[5]; 
 
  const double wx6R2 = std::pow(wx6,2);
  const double dx6R2 = std::pow(dx6,2);

  out[0] += 8.0*fld[0]*volFac*wx6R2+4.618802153517007*fld[6]*dx6*volFac*wx6+0.6666666666666666*fld[0]*dx6R2*volFac; 

} 
void IntDGMoment6xSer_x6Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx6 = w[5]; 
  const double dx6 = dx[5]; 
 
  const double wx6R2 = std::pow(wx6,2);
  const double dx6R2 = std::pow(dx6,2);

  out[0] += 8.0*fld[0]*volFac*wx6R2+4.618802153517007*fld[6]*dx6*volFac*wx6+0.5962847939999438*fld[27]*dx6R2*volFac+0.6666666666666666*fld[0]*dx6R2*volFac; 

} 
void IntDGMoment6xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
  const double wx5 = w[4]; 
  const double dx5 = dx[4]; 
  const double wx6 = w[5]; 
  const double dx6 = dx[5]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double wx3R2 = std::pow(wx3,2);
  const double wx4R2 = std::pow(wx4,2);
  const double wx5R2 = std::pow(wx5,2);
  const double wx6R2 = std::pow(wx6,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);
  const double dx3R2 = std::pow(dx3,2);
  const double dx4R2 = std::pow(dx4,2);
  const double dx5R2 = std::pow(dx5,2);
  const double dx6R2 = std::pow(dx6,2);

  out[0] += 8.0*fld[0]*volFac*wx6R2+4.618802153517007*fld[6]*dx6*volFac*wx6+8.0*fld[0]*volFac*wx5R2+4.618802153517007*fld[5]*dx5*volFac*wx5+8.0*fld[0]*volFac*wx4R2+4.618802153517007*fld[4]*dx4*volFac*wx4+8.0*fld[0]*volFac*wx3R2+4.618802153517007*fld[3]*dx3*volFac*wx3+8.0*fld[0]*volFac*wx2R2+4.618802153517007*fld[2]*dx2*volFac*wx2+8.0*fld[0]*volFac*wx1R2+4.618802153517007*fld[1]*dx1*volFac*wx1+0.6666666666666666*fld[0]*dx6R2*volFac+0.6666666666666666*fld[0]*dx5R2*volFac+0.6666666666666666*fld[0]*dx4R2*volFac+0.6666666666666666*fld[0]*dx3R2*volFac+0.6666666666666666*fld[0]*dx2R2*volFac+0.6666666666666666*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment6xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*dx[5]*0.015625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
  const double wx5 = w[4]; 
  const double dx5 = dx[4]; 
  const double wx6 = w[5]; 
  const double dx6 = dx[5]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double wx3R2 = std::pow(wx3,2);
  const double wx4R2 = std::pow(wx4,2);
  const double wx5R2 = std::pow(wx5,2);
  const double wx6R2 = std::pow(wx6,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);
  const double dx3R2 = std::pow(dx3,2);
  const double dx4R2 = std::pow(dx4,2);
  const double dx5R2 = std::pow(dx5,2);
  const double dx6R2 = std::pow(dx6,2);

  out[0] += 8.0*fld[0]*volFac*wx6R2+4.618802153517007*fld[6]*dx6*volFac*wx6+8.0*fld[0]*volFac*wx5R2+4.618802153517007*fld[5]*dx5*volFac*wx5+8.0*fld[0]*volFac*wx4R2+4.618802153517007*fld[4]*dx4*volFac*wx4+8.0*fld[0]*volFac*wx3R2+4.618802153517007*fld[3]*dx3*volFac*wx3+8.0*fld[0]*volFac*wx2R2+4.618802153517007*fld[2]*dx2*volFac*wx2+8.0*fld[0]*volFac*wx1R2+4.618802153517007*fld[1]*dx1*volFac*wx1+0.5962847939999438*fld[27]*dx6R2*volFac+0.6666666666666666*fld[0]*dx6R2*volFac+0.5962847939999438*fld[26]*dx5R2*volFac+0.6666666666666666*fld[0]*dx5R2*volFac+0.5962847939999438*fld[25]*dx4R2*volFac+0.6666666666666666*fld[0]*dx4R2*volFac+0.5962847939999438*fld[24]*dx3R2*volFac+0.6666666666666666*fld[0]*dx3R2*volFac+0.5962847939999438*fld[23]*dx2R2*volFac+0.6666666666666666*fld[0]*dx2R2*volFac+0.5962847939999438*fld[22]*dx1R2*volFac+0.6666666666666666*fld[0]*dx1R2*volFac; 

} 
