#ifndef CHARGEEXCHANGE_MOD_DELC_H 
#define CHARGEEXCHANGE_MOD_DELC_H 
#include <cmath> 

#include <algorithm> 

extern "C" { 
void SigmaCXcellAvSer1x1v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvSer1x2v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvSer1x3v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvSer2x2v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvSer2x3v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvSer3x3v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvSer1x1v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvSer1x2v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvSer1x3v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvSer2x2v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvSer2x3v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvSer3x3v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvSer1x1v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvSer1x2v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvSer1x3v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvSer2x2v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvSer2x3v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvSer3x3v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvMax1x1v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvMax1x2v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvMax1x3v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvMax2x2v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvMax2x3v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvMax3x3v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvMax1x1v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvMax1x2v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvMax1x3v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvMax2x2v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvMax2x3v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvMax3x3v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvMax1x1v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvMax1x2v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvMax1x3v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvMax2x2v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 

void SigmaCXcellAvMax2x3v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 


void SigmaCXcellAvMax3x3v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX); 



 
} 
#endif 
