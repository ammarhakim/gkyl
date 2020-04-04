#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment5xSer_one_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
 
  out[0] += 5.656854249492382*fld[0]*volFac; 

} 
void IntDGMoment5xSer_one_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
 
  out[0] += 5.656854249492382*fld[0]*volFac; 

} 
void IntDGMoment5xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wx1+1.632993161855453*fld[1]*dx1*volFac; 

} 
void IntDGMoment5xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wx1+1.632993161855453*fld[1]*dx1*volFac; 

} 
void IntDGMoment5xSer_x2_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wx2+1.632993161855453*fld[2]*dx2*volFac; 

} 
void IntDGMoment5xSer_x2_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wx2+1.632993161855453*fld[2]*dx2*volFac; 

} 
void IntDGMoment5xSer_x3_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wx3+1.632993161855453*fld[3]*dx3*volFac; 

} 
void IntDGMoment5xSer_x3_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wx3+1.632993161855453*fld[3]*dx3*volFac; 

} 
void IntDGMoment5xSer_x4_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wx4+1.632993161855453*fld[4]*dx4*volFac; 

} 
void IntDGMoment5xSer_x4_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wx4+1.632993161855453*fld[4]*dx4*volFac; 

} 
void IntDGMoment5xSer_x5_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx5 = w[4]; 
  const double dx5 = dx[4]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wx5+1.632993161855453*fld[5]*dx5*volFac; 

} 
void IntDGMoment5xSer_x5_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx5 = w[4]; 
  const double dx5 = dx[4]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wx5+1.632993161855453*fld[5]*dx5*volFac; 

} 
void IntDGMoment5xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx1R2+3.265986323710906*fld[1]*dx1*volFac*wx1+0.4714045207910317*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment5xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx1 = w[0]; 
  const double dx1 = dx[0]; 
 
  const double wx1R2 = std::pow(wx1,2);
  const double dx1R2 = std::pow(dx1,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx1R2+3.265986323710906*fld[1]*dx1*volFac*wx1+0.421637021355784*fld[16]*dx1R2*volFac+0.4714045207910317*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment5xSer_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx2R2+3.265986323710906*fld[2]*dx2*volFac*wx2+0.4714045207910317*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment5xSer_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx2 = w[1]; 
  const double dx2 = dx[1]; 
 
  const double wx2R2 = std::pow(wx2,2);
  const double dx2R2 = std::pow(dx2,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx2R2+3.265986323710906*fld[2]*dx2*volFac*wx2+0.421637021355784*fld[17]*dx2R2*volFac+0.4714045207910317*fld[0]*dx2R2*volFac; 

} 
void IntDGMoment5xSer_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx3R2 = std::pow(wx3,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx3R2+3.265986323710906*fld[3]*dx3*volFac*wx3+0.4714045207910317*fld[0]*dx3R2*volFac; 

} 
void IntDGMoment5xSer_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx3 = w[2]; 
  const double dx3 = dx[2]; 
 
  const double wx3R2 = std::pow(wx3,2);
  const double dx3R2 = std::pow(dx3,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx3R2+3.265986323710906*fld[3]*dx3*volFac*wx3+0.421637021355784*fld[18]*dx3R2*volFac+0.4714045207910317*fld[0]*dx3R2*volFac; 

} 
void IntDGMoment5xSer_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 
  const double wx4R2 = std::pow(wx4,2);
  const double dx4R2 = std::pow(dx4,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx4R2+3.265986323710906*fld[4]*dx4*volFac*wx4+0.4714045207910317*fld[0]*dx4R2*volFac; 

} 
void IntDGMoment5xSer_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx4 = w[3]; 
  const double dx4 = dx[3]; 
 
  const double wx4R2 = std::pow(wx4,2);
  const double dx4R2 = std::pow(dx4,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx4R2+3.265986323710906*fld[4]*dx4*volFac*wx4+0.421637021355784*fld[19]*dx4R2*volFac+0.4714045207910317*fld[0]*dx4R2*volFac; 

} 
void IntDGMoment5xSer_x5Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx5 = w[4]; 
  const double dx5 = dx[4]; 
 
  const double wx5R2 = std::pow(wx5,2);
  const double dx5R2 = std::pow(dx5,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx5R2+3.265986323710906*fld[5]*dx5*volFac*wx5+0.4714045207910317*fld[0]*dx5R2*volFac; 

} 
void IntDGMoment5xSer_x5Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wx5 = w[4]; 
  const double dx5 = dx[4]; 
 
  const double wx5R2 = std::pow(wx5,2);
  const double dx5R2 = std::pow(dx5,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx5R2+3.265986323710906*fld[5]*dx5*volFac*wx5+0.421637021355784*fld[20]*dx5R2*volFac+0.4714045207910317*fld[0]*dx5R2*volFac; 

} 
void IntDGMoment5xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
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
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double wx3R2 = std::pow(wx3,2);
  const double wx4R2 = std::pow(wx4,2);
  const double wx5R2 = std::pow(wx5,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);
  const double dx3R2 = std::pow(dx3,2);
  const double dx4R2 = std::pow(dx4,2);
  const double dx5R2 = std::pow(dx5,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx5R2+3.265986323710906*fld[5]*dx5*volFac*wx5+5.656854249492382*fld[0]*volFac*wx4R2+3.265986323710906*fld[4]*dx4*volFac*wx4+5.656854249492382*fld[0]*volFac*wx3R2+3.265986323710906*fld[3]*dx3*volFac*wx3+5.656854249492382*fld[0]*volFac*wx2R2+3.265986323710906*fld[2]*dx2*volFac*wx2+5.656854249492382*fld[0]*volFac*wx1R2+3.265986323710906*fld[1]*dx1*volFac*wx1+0.4714045207910317*fld[0]*dx5R2*volFac+0.4714045207910317*fld[0]*dx4R2*volFac+0.4714045207910317*fld[0]*dx3R2*volFac+0.4714045207910317*fld[0]*dx2R2*volFac+0.4714045207910317*fld[0]*dx1R2*volFac; 

} 
void IntDGMoment5xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
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
 
  const double wx1R2 = std::pow(wx1,2);
  const double wx2R2 = std::pow(wx2,2);
  const double wx3R2 = std::pow(wx3,2);
  const double wx4R2 = std::pow(wx4,2);
  const double wx5R2 = std::pow(wx5,2);
  const double dx1R2 = std::pow(dx1,2);
  const double dx2R2 = std::pow(dx2,2);
  const double dx3R2 = std::pow(dx3,2);
  const double dx4R2 = std::pow(dx4,2);
  const double dx5R2 = std::pow(dx5,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wx5R2+3.265986323710906*fld[5]*dx5*volFac*wx5+5.656854249492382*fld[0]*volFac*wx4R2+3.265986323710906*fld[4]*dx4*volFac*wx4+5.656854249492382*fld[0]*volFac*wx3R2+3.265986323710906*fld[3]*dx3*volFac*wx3+5.656854249492382*fld[0]*volFac*wx2R2+3.265986323710906*fld[2]*dx2*volFac*wx2+5.656854249492382*fld[0]*volFac*wx1R2+3.265986323710906*fld[1]*dx1*volFac*wx1+0.421637021355784*fld[20]*dx5R2*volFac+0.4714045207910317*fld[0]*dx5R2*volFac+0.421637021355784*fld[19]*dx4R2*volFac+0.4714045207910317*fld[0]*dx4R2*volFac+0.421637021355784*fld[18]*dx3R2*volFac+0.4714045207910317*fld[0]*dx3R2*volFac+0.421637021355784*fld[17]*dx2R2*volFac+0.4714045207910317*fld[0]*dx2R2*volFac+0.421637021355784*fld[16]*dx1R2*volFac+0.4714045207910317*fld[0]*dx1R2*volFac; 

} 
