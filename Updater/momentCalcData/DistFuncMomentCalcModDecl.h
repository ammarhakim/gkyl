#ifndef DIST_FUNC_MOMENT_CALC_MOD_DELC_H 
#define DIST_FUNC_MOMENT_CALC_MOD_DELC_H 
extern "C" { 
void MomentCalc1x1vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x2vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x2vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc3x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 



 
void MomentCalc1x1vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x2vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x3vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x2vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x3vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc3x3vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 



 
} 
#endif 
