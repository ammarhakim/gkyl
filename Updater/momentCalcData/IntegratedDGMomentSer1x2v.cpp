#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment1x2vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wv1+0.8164965809277261*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x2vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wv1+0.8164965809277261*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x2vSer_v1_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wv1+0.8164965809277261*fld[2]*dv1*volFac; 

} 
void IntDGMoment1x2vSer_v2_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wv2+0.8164965809277261*fld[3]*dv2*volFac; 

} 
void IntDGMoment1x2vSer_v2_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wv2+0.8164965809277261*fld[3]*dv2*volFac; 

} 
void IntDGMoment1x2vSer_v2_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wv2+0.8164965809277261*fld[3]*dv2*volFac; 

} 
void IntDGMoment1x2vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wv1R2+1.632993161855453*fld[2]*dv1*volFac*wv1+0.2357022603955158*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x2vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wv1R2+1.632993161855453*fld[2]*dv1*volFac*wv1+0.210818510677892*fld[8]*dv1R2*volFac+0.2357022603955158*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x2vSer_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wv1R2+1.632993161855453*fld[2]*dv1*volFac*wv1+0.210818510677892*fld[8]*dv1R2*volFac+0.2357022603955158*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x2vSer_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wv2R2+1.632993161855453*fld[3]*dv2*volFac*wv2+0.2357022603955158*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment1x2vSer_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wv2R2+1.632993161855453*fld[3]*dv2*volFac*wv2+0.210818510677892*fld[9]*dv2R2*volFac+0.2357022603955158*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment1x2vSer_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wv2R2+1.632993161855453*fld[3]*dv2*volFac*wv2+0.210818510677892*fld[9]*dv2R2*volFac+0.2357022603955158*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment1x2vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wv1+0.8164965809277261*fld[2]*dv1*volFac; 
  out[1] += 2.828427124746191*fld[0]*volFac*wv2+0.8164965809277261*fld[3]*dv2*volFac; 

} 
void IntDGMoment1x2vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wv1+0.8164965809277261*fld[2]*dv1*volFac; 
  out[1] += 2.828427124746191*fld[0]*volFac*wv2+0.8164965809277261*fld[3]*dv2*volFac; 

} 
void IntDGMoment1x2vSer_vi_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 

  out[0] += 2.828427124746191*fld[0]*volFac*wv1+0.8164965809277261*fld[2]*dv1*volFac; 
  out[1] += 2.828427124746191*fld[0]*volFac*wv2+0.8164965809277261*fld[3]*dv2*volFac; 

} 
void IntDGMoment1x2vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wv2R2+1.632993161855453*fld[3]*dv2*volFac*wv2+2.828427124746191*fld[0]*volFac*wv1R2+1.632993161855453*fld[2]*dv1*volFac*wv1+0.2357022603955158*fld[0]*dv2R2*volFac+0.2357022603955158*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x2vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wv2R2+1.632993161855453*fld[3]*dv2*volFac*wv2+2.828427124746191*fld[0]*volFac*wv1R2+1.632993161855453*fld[2]*dv1*volFac*wv1+0.210818510677892*fld[9]*dv2R2*volFac+0.2357022603955158*fld[0]*dv2R2*volFac+0.210818510677892*fld[8]*dv1R2*volFac+0.2357022603955158*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x2vSer_vSq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 2.828427124746191*fld[0]*volFac*wv2R2+1.632993161855453*fld[3]*dv2*volFac*wv2+2.828427124746191*fld[0]*volFac*wv1R2+1.632993161855453*fld[2]*dv1*volFac*wv1+0.210818510677892*fld[9]*dv2R2*volFac+0.2357022603955158*fld[0]*dv2R2*volFac+0.210818510677892*fld[8]*dv1R2*volFac+0.2357022603955158*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x2vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 2.828427124746191*fld[0]*volFac; 
  out[1] += 2.828427124746191*fld[0]*volFac*wv1+0.8164965809277261*fld[2]*dv1*volFac; 
  out[2] += 2.828427124746191*fld[0]*volFac*wv2+0.8164965809277261*fld[3]*dv2*volFac; 
  out[3] += 2.828427124746191*fld[0]*volFac*wv2R2+1.632993161855453*fld[3]*dv2*volFac*wv2+2.828427124746191*fld[0]*volFac*wv1R2+1.632993161855453*fld[2]*dv1*volFac*wv1+0.2357022603955158*fld[0]*dv2R2*volFac+0.2357022603955158*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x2vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 2.828427124746191*fld[0]*volFac; 
  out[1] += 2.828427124746191*fld[0]*volFac*wv1+0.8164965809277261*fld[2]*dv1*volFac; 
  out[2] += 2.828427124746191*fld[0]*volFac*wv2+0.8164965809277261*fld[3]*dv2*volFac; 
  out[3] += 2.828427124746191*fld[0]*volFac*wv2R2+1.632993161855453*fld[3]*dv2*volFac*wv2+2.828427124746191*fld[0]*volFac*wv1R2+1.632993161855453*fld[2]*dv1*volFac*wv1+0.210818510677892*fld[9]*dv2R2*volFac+0.2357022603955158*fld[0]*dv2R2*volFac+0.210818510677892*fld[8]*dv1R2*volFac+0.2357022603955158*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment1x2vSer_intM_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*0.125; 
  const double wv1 = w[1]; 
  const double dv1 = dx[1]; 
  const double wv2 = w[2]; 
  const double dv2 = dx[2]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 2.828427124746191*fld[0]*volFac; 
  out[1] += 2.828427124746191*fld[0]*volFac*wv1+0.8164965809277261*fld[2]*dv1*volFac; 
  out[2] += 2.828427124746191*fld[0]*volFac*wv2+0.8164965809277261*fld[3]*dv2*volFac; 
  out[3] += 2.828427124746191*fld[0]*volFac*wv2R2+1.632993161855453*fld[3]*dv2*volFac*wv2+2.828427124746191*fld[0]*volFac*wv1R2+1.632993161855453*fld[2]*dv1*volFac*wv1+0.210818510677892*fld[9]*dv2R2*volFac+0.2357022603955158*fld[0]*dv2R2*volFac+0.210818510677892*fld[8]*dv1R2*volFac+0.2357022603955158*fld[0]*dv1R2*volFac; 

} 
