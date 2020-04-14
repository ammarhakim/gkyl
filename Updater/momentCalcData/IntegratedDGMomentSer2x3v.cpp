#include <IntegratedDGMomentModDecl.h> 
 
void IntDGMoment2x3vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wv1+1.632993161855453*fld[3]*dv1*volFac; 

} 
void IntDGMoment2x3vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wv1+1.632993161855453*fld[3]*dv1*volFac; 

} 
void IntDGMoment2x3vSer_v2_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wv2+1.632993161855453*fld[4]*dv2*volFac; 

} 
void IntDGMoment2x3vSer_v2_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wv2+1.632993161855453*fld[4]*dv2*volFac; 

} 
void IntDGMoment2x3vSer_v3_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv3 = w[4]; 
  const double dv3 = dx[4]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wv3+1.632993161855453*fld[5]*dv3*volFac; 

} 
void IntDGMoment2x3vSer_v3_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv3 = w[4]; 
  const double dv3 = dx[4]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wv3+1.632993161855453*fld[5]*dv3*volFac; 

} 
void IntDGMoment2x3vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wv1R2+3.265986323710906*fld[3]*dv1*volFac*wv1+0.4714045207910317*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x3vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double dv1R2 = std::pow(dv1,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wv1R2+3.265986323710906*fld[3]*dv1*volFac*wv1+0.421637021355784*fld[18]*dv1R2*volFac+0.4714045207910317*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x3vSer_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wv2R2+3.265986323710906*fld[4]*dv2*volFac*wv2+0.4714045207910317*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment2x3vSer_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
 
  const double wv2R2 = std::pow(wv2,2);
  const double dv2R2 = std::pow(dv2,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wv2R2+3.265986323710906*fld[4]*dv2*volFac*wv2+0.421637021355784*fld[19]*dv2R2*volFac+0.4714045207910317*fld[0]*dv2R2*volFac; 

} 
void IntDGMoment2x3vSer_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv3 = w[4]; 
  const double dv3 = dx[4]; 
 
  const double wv3R2 = std::pow(wv3,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wv3R2+3.265986323710906*fld[5]*dv3*volFac*wv3+0.4714045207910317*fld[0]*dv3R2*volFac; 

} 
void IntDGMoment2x3vSer_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv3 = w[4]; 
  const double dv3 = dx[4]; 
 
  const double wv3R2 = std::pow(wv3,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wv3R2+3.265986323710906*fld[5]*dv3*volFac*wv3+0.421637021355784*fld[20]*dv3R2*volFac+0.4714045207910317*fld[0]*dv3R2*volFac; 

} 
void IntDGMoment2x3vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
  const double wv3 = w[4]; 
  const double dv3 = dx[4]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wv1+1.632993161855453*fld[3]*dv1*volFac; 
  out[1] += 5.656854249492382*fld[0]*volFac*wv2+1.632993161855453*fld[4]*dv2*volFac; 
  out[2] += 5.656854249492382*fld[0]*volFac*wv3+1.632993161855453*fld[5]*dv3*volFac; 

} 
void IntDGMoment2x3vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
  const double wv3 = w[4]; 
  const double dv3 = dx[4]; 
 

  out[0] += 5.656854249492382*fld[0]*volFac*wv1+1.632993161855453*fld[3]*dv1*volFac; 
  out[1] += 5.656854249492382*fld[0]*volFac*wv2+1.632993161855453*fld[4]*dv2*volFac; 
  out[2] += 5.656854249492382*fld[0]*volFac*wv3+1.632993161855453*fld[5]*dv3*volFac; 

} 
void IntDGMoment2x3vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
  const double wv3 = w[4]; 
  const double dv3 = dx[4]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wv3R2+3.265986323710906*fld[5]*dv3*volFac*wv3+5.656854249492382*fld[0]*volFac*wv2R2+3.265986323710906*fld[4]*dv2*volFac*wv2+5.656854249492382*fld[0]*volFac*wv1R2+3.265986323710906*fld[3]*dv1*volFac*wv1+0.4714045207910317*fld[0]*dv3R2*volFac+0.4714045207910317*fld[0]*dv2R2*volFac+0.4714045207910317*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x3vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
  const double wv3 = w[4]; 
  const double dv3 = dx[4]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 5.656854249492382*fld[0]*volFac*wv3R2+3.265986323710906*fld[5]*dv3*volFac*wv3+5.656854249492382*fld[0]*volFac*wv2R2+3.265986323710906*fld[4]*dv2*volFac*wv2+5.656854249492382*fld[0]*volFac*wv1R2+3.265986323710906*fld[3]*dv1*volFac*wv1+0.421637021355784*fld[20]*dv3R2*volFac+0.4714045207910317*fld[0]*dv3R2*volFac+0.421637021355784*fld[19]*dv2R2*volFac+0.4714045207910317*fld[0]*dv2R2*volFac+0.421637021355784*fld[18]*dv1R2*volFac+0.4714045207910317*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x3vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
  const double wv3 = w[4]; 
  const double dv3 = dx[4]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 5.656854249492382*fld[0]*volFac; 
  out[1] += 5.656854249492382*fld[0]*volFac*wv1+1.632993161855453*fld[3]*dv1*volFac; 
  out[2] += 5.656854249492382*fld[0]*volFac*wv2+1.632993161855453*fld[4]*dv2*volFac; 
  out[3] += 5.656854249492382*fld[0]*volFac*wv3+1.632993161855453*fld[5]*dv3*volFac; 
  out[4] += 5.656854249492382*fld[0]*volFac*wv3R2+3.265986323710906*fld[5]*dv3*volFac*wv3+5.656854249492382*fld[0]*volFac*wv2R2+3.265986323710906*fld[4]*dv2*volFac*wv2+5.656854249492382*fld[0]*volFac*wv1R2+3.265986323710906*fld[3]*dv1*volFac*wv1+0.4714045207910317*fld[0]*dv3R2*volFac+0.4714045207910317*fld[0]*dv2R2*volFac+0.4714045207910317*fld[0]*dv1R2*volFac; 

} 
void IntDGMoment2x3vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out) 
{ 
  const double volFac = dx[0]*dx[1]*dx[2]*dx[3]*dx[4]*0.03125; 
  const double wv1 = w[2]; 
  const double dv1 = dx[2]; 
  const double wv2 = w[3]; 
  const double dv2 = dx[3]; 
  const double wv3 = w[4]; 
  const double dv3 = dx[4]; 
 
  const double wv1R2 = std::pow(wv1,2);
  const double wv2R2 = std::pow(wv2,2);
  const double wv3R2 = std::pow(wv3,2);
  const double dv1R2 = std::pow(dv1,2);
  const double dv2R2 = std::pow(dv2,2);
  const double dv3R2 = std::pow(dv3,2);

  out[0] += 5.656854249492382*fld[0]*volFac; 
  out[1] += 5.656854249492382*fld[0]*volFac*wv1+1.632993161855453*fld[3]*dv1*volFac; 
  out[2] += 5.656854249492382*fld[0]*volFac*wv2+1.632993161855453*fld[4]*dv2*volFac; 
  out[3] += 5.656854249492382*fld[0]*volFac*wv3+1.632993161855453*fld[5]*dv3*volFac; 
  out[4] += 5.656854249492382*fld[0]*volFac*wv3R2+3.265986323710906*fld[5]*dv3*volFac*wv3+5.656854249492382*fld[0]*volFac*wv2R2+3.265986323710906*fld[4]*dv2*volFac*wv2+5.656854249492382*fld[0]*volFac*wv1R2+3.265986323710906*fld[3]*dv1*volFac*wv1+0.421637021355784*fld[20]*dv3R2*volFac+0.4714045207910317*fld[0]*dv3R2*volFac+0.421637021355784*fld[19]*dv2R2*volFac+0.4714045207910317*fld[0]*dv2R2*volFac+0.421637021355784*fld[18]*dv1R2*volFac+0.4714045207910317*fld[0]*dv1R2*volFac; 

} 
