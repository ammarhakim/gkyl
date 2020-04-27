#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment1xSer_one_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*0.5; 
 
  out[0] += 1.414213562373095*fld[0]*volFac; 

} 
void IntDGMoment1xSer_one_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*0.5; 
 
  out[0] += 1.414213562373095*fld[0]*volFac; 

} 
void IntDGMoment1xSer_one_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*0.5; 
 
  out[0] += 1.414213562373095*fld[0]*volFac; 

} 
void IntDGMoment1xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 1.414213562373095*fld[0]*volFac*wx1+0.408248290463863*fld[1]*dx1*volFac; 

} 
void IntDGMoment1xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 1.414213562373095*fld[0]*volFac*wx1+0.408248290463863*fld[1]*dx1*volFac; 

} 
void IntDGMoment1xSer_x1_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 1.414213562373095*fld[0]*volFac*wx1+0.408248290463863*fld[1]*dx1*volFac; 

} 
void IntDGMoment1xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 1.414213562373095*fld[0]*volFac*wx1R2+0.8164965809277261*fld[1]*dx1*volFac*wx1+0.1178511301977579*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment1xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 1.414213562373095*fld[0]*volFac*wx1R2+0.8164965809277261*fld[1]*dx1*volFac*wx1+0.105409255338946*fld[2]*dx1R2*volFac+0.1178511301977579*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment1xSer_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 1.414213562373095*fld[0]*volFac*wx1R2+0.8164965809277261*fld[1]*dx1*volFac*wx1+0.105409255338946*fld[2]*dx1R2*volFac+0.1178511301977579*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment1xSer_xi_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 1.414213562373095*fld[0]*volFac*wx1+0.408248290463863*fld[1]*dx1*volFac; 

} 
void IntDGMoment1xSer_xi_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 1.414213562373095*fld[0]*volFac*wx1+0.408248290463863*fld[1]*dx1*volFac; 

} 
void IntDGMoment1xSer_xi_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 1.414213562373095*fld[0]*volFac*wx1+0.408248290463863*fld[1]*dx1*volFac; 

} 
void IntDGMoment1xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 1.414213562373095*fld[0]*volFac*wx1R2+0.8164965809277261*fld[1]*dx1*volFac*wx1+0.1178511301977579*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment1xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 1.414213562373095*fld[0]*volFac*wx1R2+0.8164965809277261*fld[1]*dx1*volFac*wx1+0.105409255338946*fld[2]*dx1R2*volFac+0.1178511301977579*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment1xSer_xSq_P3(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*0.5; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 1.414213562373095*fld[0]*volFac*wx1R2+0.8164965809277261*fld[1]*dx1*volFac*wx1+0.105409255338946*fld[2]*dx1R2*volFac+0.1178511301977579*fld[0]*dx1R2*volFac; 

} 
