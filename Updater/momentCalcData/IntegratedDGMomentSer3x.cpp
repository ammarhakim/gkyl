#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment3xSer_one_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
 
  out[0] += 2.828427124746191*fld[0]*volFac; 

} 
void IntDGMoment3xSer_one_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
 
  out[0] += 2.828427124746191*fld[0]*volFac; 

} 
void IntDGMoment3xSer_one_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
 
  out[0] += 2.828427124746191*fld[0]*volFac; 

} 
void IntDGMoment3xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx1+0.8164965809277261*fld[1]*dx1*volFac; 

} 
void IntDGMoment3xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx1+0.8164965809277261*fld[1]*dx1*volFac; 

} 
void IntDGMoment3xSer_x1_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx1+0.8164965809277261*fld[1]*dx1*volFac; 

} 
void IntDGMoment3xSer_x2_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx2+0.8164965809277261*fld[2]*dx2*volFac; 

} 
void IntDGMoment3xSer_x2_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx2+0.8164965809277261*fld[2]*dx2*volFac; 

} 
void IntDGMoment3xSer_x2_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx2+0.8164965809277261*fld[2]*dx2*volFac; 

} 
void IntDGMoment3xSer_x3_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx3+0.8164965809277261*fld[3]*dx3*volFac; 

} 
void IntDGMoment3xSer_x3_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx3+0.8164965809277261*fld[3]*dx3*volFac; 

} 
void IntDGMoment3xSer_x3_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx3+0.8164965809277261*fld[3]*dx3*volFac; 

} 
void IntDGMoment3xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx1R2+1.632993161855453*fld[1]*dx1*volFac*wx1+0.2357022603955158*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment3xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx1R2+1.632993161855453*fld[1]*dx1*volFac*wx1+0.210818510677892*fld[7]*dx1R2*volFac+0.2357022603955158*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment3xSer_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx1R2+1.632993161855453*fld[1]*dx1*volFac*wx1+0.210818510677892*fld[7]*dx1R2*volFac+0.2357022603955158*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment3xSer_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx2R2+1.632993161855453*fld[2]*dx2*volFac*wx2+0.2357022603955158*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment3xSer_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx2R2+1.632993161855453*fld[2]*dx2*volFac*wx2+0.210818510677892*fld[8]*dx2R2*volFac+0.2357022603955158*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment3xSer_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx2R2+1.632993161855453*fld[2]*dx2*volFac*wx2+0.210818510677892*fld[8]*dx2R2*volFac+0.2357022603955158*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment3xSer_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx3R2 = std::pow(wx3,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx3R2+1.632993161855453*fld[3]*dx3*volFac*wx3+0.2357022603955158*fld[0]*dx3R2*volFac; 

} 
void IntDGMoment3xSer_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx3R2 = std::pow(wx3,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx3R2+1.632993161855453*fld[3]*dx3*volFac*wx3+0.210818510677892*fld[9]*dx3R2*volFac+0.2357022603955158*fld[0]*dx3R2*volFac; 

} 
void IntDGMoment3xSer_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx3R2 = std::pow(wx3,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx3R2+1.632993161855453*fld[3]*dx3*volFac*wx3+0.210818510677892*fld[9]*dx3R2*volFac+0.2357022603955158*fld[0]*dx3R2*volFac; 

} 
void IntDGMoment3xSer_xi_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx1+0.8164965809277261*fld[1]*dx1*volFac; 
  out[1] += 2.828427124746191*fld[0]*volFac*wx2+0.8164965809277261*fld[2]*dx2*volFac; 
  out[2] += 2.828427124746191*fld[0]*volFac*wx3+0.8164965809277261*fld[3]*dx3*volFac; 

} 
void IntDGMoment3xSer_xi_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx1+0.8164965809277261*fld[1]*dx1*volFac; 
  out[1] += 2.828427124746191*fld[0]*volFac*wx2+0.8164965809277261*fld[2]*dx2*volFac; 
  out[2] += 2.828427124746191*fld[0]*volFac*wx3+0.8164965809277261*fld[3]*dx3*volFac; 

} 
void IntDGMoment3xSer_xi_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wx1+0.8164965809277261*fld[1]*dx1*volFac; 
  out[1] += 2.828427124746191*fld[0]*volFac*wx2+0.8164965809277261*fld[2]*dx2*volFac; 
  out[2] += 2.828427124746191*fld[0]*volFac*wx3+0.8164965809277261*fld[3]*dx3*volFac; 

} 
void IntDGMoment3xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double wx3R2 = std::pow(wx3,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx3R2+1.632993161855453*fld[3]*dx3*volFac*wx3+2.828427124746191*fld[0]*volFac*wx2R2+1.632993161855453*fld[2]*dx2*volFac*wx2+2.828427124746191*fld[0]*volFac*wx1R2+1.632993161855453*fld[1]*dx1*volFac*wx1+0.2357022603955158*fld[0]*dx3R2*volFac+0.2357022603955158*fld[0]*dx2R2*volFac+0.2357022603955158*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment3xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double wx3R2 = std::pow(wx3,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx3R2+1.632993161855453*fld[3]*dx3*volFac*wx3+2.828427124746191*fld[0]*volFac*wx2R2+1.632993161855453*fld[2]*dx2*volFac*wx2+2.828427124746191*fld[0]*volFac*wx1R2+1.632993161855453*fld[1]*dx1*volFac*wx1+0.210818510677892*fld[9]*dx3R2*volFac+0.2357022603955158*fld[0]*dx3R2*volFac+0.210818510677892*fld[8]*dx2R2*volFac+0.2357022603955158*fld[0]*dx2R2*volFac+0.210818510677892*fld[7]*dx1R2*volFac+0.2357022603955158*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment3xSer_xSq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double wx3R2 = std::pow(wx3,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wx3R2+1.632993161855453*fld[3]*dx3*volFac*wx3+2.828427124746191*fld[0]*volFac*wx2R2+1.632993161855453*fld[2]*dx2*volFac*wx2+2.828427124746191*fld[0]*volFac*wx1R2+1.632993161855453*fld[1]*dx1*volFac*wx1+0.210818510677892*fld[9]*dx3R2*volFac+0.2357022603955158*fld[0]*dx3R2*volFac+0.210818510677892*fld[8]*dx2R2*volFac+0.2357022603955158*fld[0]*dx2R2*volFac+0.210818510677892*fld[7]*dx1R2*volFac+0.2357022603955158*fld[0]*dx1R2*volFac; 

} 
