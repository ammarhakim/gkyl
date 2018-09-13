#ifndef DIST_FUNC_MOMENT_CALC_MOD_DELC_H 
#define DIST_FUNC_MOMENT_CALC_MOD_DELC_H 
extern "C" { 
void MomentCalc1x1vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc1x1vSer_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc1x1vSer_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vSer_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vSer_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vSer_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vSer_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vSer_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void GkMomentCalc1x1vSer_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x1vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc1x1vSer_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vSer_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vSer_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vSer_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vSer_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vSer_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void IntMomentCalc1x1vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x1vSer_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x2vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc1x2vSer_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc1x2vSer_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void GkMomentCalc1x2vSer_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x2vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc1x2vSer_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vSer_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void IntMomentCalc1x2vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x2vSer_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc1x3vSer_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x3vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x3vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x3vSer_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x2vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc2x2vSer_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc2x2vSer_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M0_step1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M0_step2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void GkMomentCalc2x2vSer_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x2vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc2x2vSer_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M0_step1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M0_step2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vSer_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void IntMomentCalc2x2vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x2vSer_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc2x3vSer_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x3vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x3vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x3vSer_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc3x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc3x3vSer_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc3x3vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc3x3vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc3x3vSer_P3(const double *w, const double *dxv, const double *f, double *out); 


void GkMomentCalc3x2vSer_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M0_step1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M0_step2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void GkMomentCalc3x2vSer_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 

void GkMomentCalc3x2vSer_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M0_step1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M0_step2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vSer_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 

void MomentCalc1x1vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc1x1vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc1x1vMax_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void GkMomentCalc1x1vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x1vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc1x1vMax_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x1vMax_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void IntMomentCalc1x1vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x1vMax_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x2vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc1x2vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc1x2vMax_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void GkMomentCalc1x2vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x2vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc1x2vMax_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc1x2vMax_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void IntMomentCalc1x2vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x2vMax_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x3vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc1x3vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x3vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x3vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x3vMax_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x2vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc2x2vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc2x2vMax_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void GkMomentCalc2x2vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x2vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void GkMomentCalc2x2vMax_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc2x2vMax_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void IntMomentCalc2x2vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x2vMax_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x3vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc2x3vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x3vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x3vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x3vMax_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc3x3vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void MomentCalc3x3vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc3x3vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc3x3vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc3x3vMax_P3(const double *w, const double *dxv, const double *f, double *out); 



 
void GkMomentCalc3x2vMax_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 
void GkMomentCalc3x2vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 

void GkMomentCalc3x2vMax_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out); 
void GkMomentCalc3x2vMax_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3); 

void MomentCalc1x1vTensor_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x1vTensor_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 

void MomentCalc1x1vTensor_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x1vTensor_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vTensor_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vTensor_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x1vTensor_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x2vTensor_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x2vTensor_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 

void MomentCalc1x2vTensor_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x2vTensor_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vTensor_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vTensor_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x2vTensor_P3(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x3vTensor_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vTensor_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vTensor_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vTensor_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vTensor_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vTensor_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x3vTensor_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 

void MomentCalc1x3vTensor_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vTensor_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vTensor_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vTensor_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vTensor_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vTensor_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc1x3vTensor_P2(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x2vTensor_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vTensor_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vTensor_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vTensor_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vTensor_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vTensor_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x2vTensor_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 

void MomentCalc2x2vTensor_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vTensor_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vTensor_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vTensor_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vTensor_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vTensor_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2); 
void IntMomentCalc2x2vTensor_P2(const double *w, const double *dxv, const double *f, double *out);


} 
#endif 
