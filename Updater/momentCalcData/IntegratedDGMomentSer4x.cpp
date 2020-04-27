#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment4xSer_one_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
 
  out[0] += 4.0*fld[0]*volFac; 

} 
void IntDGMoment4xSer_one_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
 
  out[0] += 4.0*fld[0]*volFac; 

} 
void IntDGMoment4xSer_one_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
 
  out[0] += 4.0*fld[0]*volFac; 

} 
void IntDGMoment4xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 4.0*fld[0]*volFac*wx1+1.154700538379252*fld[1]*dx1*volFac; 

} 
void IntDGMoment4xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 4.0*fld[0]*volFac*wx1+1.154700538379252*fld[1]*dx1*volFac; 

} 
void IntDGMoment4xSer_x1_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 4.0*fld[0]*volFac*wx1+1.154700538379252*fld[1]*dx1*volFac; 

} 
void IntDGMoment4xSer_x2_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 4.0*fld[0]*volFac*wx2+1.154700538379252*fld[2]*dx2*volFac; 

} 
void IntDGMoment4xSer_x2_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 4.0*fld[0]*volFac*wx2+1.154700538379252*fld[2]*dx2*volFac; 

} 
void IntDGMoment4xSer_x2_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 4.0*fld[0]*volFac*wx2+1.154700538379252*fld[2]*dx2*volFac; 

} 
void IntDGMoment4xSer_x3_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 4.0*fld[0]*volFac*wx3+1.154700538379252*fld[3]*dx3*volFac; 

} 
void IntDGMoment4xSer_x3_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 4.0*fld[0]*volFac*wx3+1.154700538379252*fld[3]*dx3*volFac; 

} 
void IntDGMoment4xSer_x3_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 4.0*fld[0]*volFac*wx3+1.154700538379252*fld[3]*dx3*volFac; 

} 
void IntDGMoment4xSer_x4_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wx4+1.154700538379252*fld[4]*dx4*volFac; 

} 
void IntDGMoment4xSer_x4_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wx4+1.154700538379252*fld[4]*dx4*volFac; 

} 
void IntDGMoment4xSer_x4_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wx4+1.154700538379252*fld[4]*dx4*volFac; 

} 
void IntDGMoment4xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 4.0*fld[0]*volFac*wx1R2+2.309401076758503*fld[1]*dx1*volFac*wx1+0.3333333333333333*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment4xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 4.0*fld[0]*volFac*wx1R2+2.309401076758503*fld[1]*dx1*volFac*wx1+0.2981423969999719*fld[11]*dx1R2*volFac+0.3333333333333333*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment4xSer_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 4.0*fld[0]*volFac*wx1R2+2.309401076758503*fld[1]*dx1*volFac*wx1+0.2981423969999719*fld[11]*dx1R2*volFac+0.3333333333333333*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment4xSer_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 4.0*fld[0]*volFac*wx2R2+2.309401076758503*fld[2]*dx2*volFac*wx2+0.3333333333333333*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment4xSer_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 4.0*fld[0]*volFac*wx2R2+2.309401076758503*fld[2]*dx2*volFac*wx2+0.2981423969999719*fld[12]*dx2R2*volFac+0.3333333333333333*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment4xSer_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 4.0*fld[0]*volFac*wx2R2+2.309401076758503*fld[2]*dx2*volFac*wx2+0.2981423969999719*fld[12]*dx2R2*volFac+0.3333333333333333*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment4xSer_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx3R2 = std::pow(wx3,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 4.0*fld[0]*volFac*wx3R2+2.309401076758503*fld[3]*dx3*volFac*wx3+0.3333333333333333*fld[0]*dx3R2*volFac; 

} 
void IntDGMoment4xSer_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx3R2 = std::pow(wx3,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 4.0*fld[0]*volFac*wx3R2+2.309401076758503*fld[3]*dx3*volFac*wx3+0.2981423969999719*fld[13]*dx3R2*volFac+0.3333333333333333*fld[0]*dx3R2*volFac; 

} 
void IntDGMoment4xSer_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx3R2 = std::pow(wx3,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 4.0*fld[0]*volFac*wx3R2+2.309401076758503*fld[3]*dx3*volFac*wx3+0.2981423969999719*fld[13]*dx3R2*volFac+0.3333333333333333*fld[0]*dx3R2*volFac; 

} 
void IntDGMoment4xSer_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 
  const double wx4R2 = std::pow(wx4,2);
  const double dx4R2 = std::pow(dx4,2);

  out[0] += 4.0*fld[0]*volFac*wx4R2+2.309401076758503*fld[4]*dx4*volFac*wx4+0.3333333333333333*fld[0]*dx4R2*volFac; 

} 
void IntDGMoment4xSer_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 
  const double wx4R2 = std::pow(wx4,2);
  const double dx4R2 = std::pow(dx4,2);

  out[0] += 4.0*fld[0]*volFac*wx4R2+2.309401076758503*fld[4]*dx4*volFac*wx4+0.2981423969999719*fld[14]*dx4R2*volFac+0.3333333333333333*fld[0]*dx4R2*volFac; 

} 
void IntDGMoment4xSer_x4Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 
  const double wx4R2 = std::pow(wx4,2);
  const double dx4R2 = std::pow(dx4,2);

  out[0] += 4.0*fld[0]*volFac*wx4R2+2.309401076758503*fld[4]*dx4*volFac*wx4+0.2981423969999719*fld[14]*dx4R2*volFac+0.3333333333333333*fld[0]*dx4R2*volFac; 

} 
void IntDGMoment4xSer_xi_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wx1+1.154700538379252*fld[1]*dx1*volFac; 
  out[1] += 4.0*fld[0]*volFac*wx2+1.154700538379252*fld[2]*dx2*volFac; 
  out[2] += 4.0*fld[0]*volFac*wx3+1.154700538379252*fld[3]*dx3*volFac; 
  out[3] += 4.0*fld[0]*volFac*wx4+1.154700538379252*fld[4]*dx4*volFac; 

} 
void IntDGMoment4xSer_xi_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wx1+1.154700538379252*fld[1]*dx1*volFac; 
  out[1] += 4.0*fld[0]*volFac*wx2+1.154700538379252*fld[2]*dx2*volFac; 
  out[2] += 4.0*fld[0]*volFac*wx3+1.154700538379252*fld[3]*dx3*volFac; 
  out[3] += 4.0*fld[0]*volFac*wx4+1.154700538379252*fld[4]*dx4*volFac; 

} 
void IntDGMoment4xSer_xi_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 

  out[0] += 4.0*fld[0]*volFac*wx1+1.154700538379252*fld[1]*dx1*volFac; 
  out[1] += 4.0*fld[0]*volFac*wx2+1.154700538379252*fld[2]*dx2*volFac; 
  out[2] += 4.0*fld[0]*volFac*wx3+1.154700538379252*fld[3]*dx3*volFac; 
  out[3] += 4.0*fld[0]*volFac*wx4+1.154700538379252*fld[4]*dx4*volFac; 

} 
void IntDGMoment4xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double wx3R2 = std::pow(wx3,2);
  const double wx4R2 = std::pow(wx4,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);
  const double dx3R2 = std::pow(dx3,2);
  const double dx4R2 = std::pow(dx4,2);

  out[0] += 4.0*fld[0]*volFac*wx4R2+2.309401076758503*fld[4]*dx4*volFac*wx4+4.0*fld[0]*volFac*wx3R2+2.309401076758503*fld[3]*dx3*volFac*wx3+4.0*fld[0]*volFac*wx2R2+2.309401076758503*fld[2]*dx2*volFac*wx2+4.0*fld[0]*volFac*wx1R2+2.309401076758503*fld[1]*dx1*volFac*wx1+0.3333333333333333*fld[0]*dx4R2*volFac+0.3333333333333333*fld[0]*dx3R2*volFac+0.3333333333333333*fld[0]*dx2R2*volFac+0.3333333333333333*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment4xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double wx3R2 = std::pow(wx3,2);
  const double wx4R2 = std::pow(wx4,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);
  const double dx3R2 = std::pow(dx3,2);
  const double dx4R2 = std::pow(dx4,2);

  out[0] += 4.0*fld[0]*volFac*wx4R2+2.309401076758503*fld[4]*dx4*volFac*wx4+4.0*fld[0]*volFac*wx3R2+2.309401076758503*fld[3]*dx3*volFac*wx3+4.0*fld[0]*volFac*wx2R2+2.309401076758503*fld[2]*dx2*volFac*wx2+4.0*fld[0]*volFac*wx1R2+2.309401076758503*fld[1]*dx1*volFac*wx1+0.2981423969999719*fld[14]*dx4R2*volFac+0.3333333333333333*fld[0]*dx4R2*volFac+0.2981423969999719*fld[13]*dx3R2*volFac+0.3333333333333333*fld[0]*dx3R2*volFac+0.2981423969999719*fld[12]*dx2R2*volFac+0.3333333333333333*fld[0]*dx2R2*volFac+0.2981423969999719*fld[11]*dx1R2*volFac+0.3333333333333333*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment4xSer_xSq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*0.0625; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double wx3R2 = std::pow(wx3,2);
  const double wx4R2 = std::pow(wx4,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);
  const double dx3R2 = std::pow(dx3,2);
  const double dx4R2 = std::pow(dx4,2);

  out[0] += 4.0*fld[0]*volFac*wx4R2+2.309401076758503*fld[4]*dx4*volFac*wx4+4.0*fld[0]*volFac*wx3R2+2.309401076758503*fld[3]*dx3*volFac*wx3+4.0*fld[0]*volFac*wx2R2+2.309401076758503*fld[2]*dx2*volFac*wx2+4.0*fld[0]*volFac*wx1R2+2.309401076758503*fld[1]*dx1*volFac*wx1+0.2981423969999719*fld[14]*dx4R2*volFac+0.3333333333333333*fld[0]*dx4R2*volFac+0.2981423969999719*fld[13]*dx3R2*volFac+0.3333333333333333*fld[0]*dx3R2*volFac+0.2981423969999719*fld[12]*dx2R2*volFac+0.3333333333333333*fld[0]*dx2R2*volFac+0.2981423969999719*fld[11]*dx1R2*volFac+0.3333333333333333*fld[0]*dx1R2*volFac; 

} 
