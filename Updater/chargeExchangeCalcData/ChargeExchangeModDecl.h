#ifndef CHARGEEXCHANGE_MOD_DELC_H 
#define CHARGEEXCHANGE_MOD_DELC_H 
#include <cmath> 

#include <algorithm> 

extern "C" { 
double SigmaCXcellAvSer1x1v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 

double SigmaCXcellAvSer1x2v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 

double SigmaCXcellAvSer1x3v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 


double SigmaCXcellAvSer2x2v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 

double SigmaCXcellAvSer2x3v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 


double SigmaCXcellAvSer3x3v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 


double SigmaCXcellAvSer1x1v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 

double SigmaCXcellAvSer1x2v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 

double SigmaCXcellAvSer1x3v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 


double SigmaCXcellAvSer2x2v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 

double SigmaCXcellAvSer2x3v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 


double SigmaCXcellAvSer3x3v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 


double SigmaCXcellAvSer1x1v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 

double SigmaCXcellAvSer1x2v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 

double SigmaCXcellAvSer1x3v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 


double SigmaCXcellAvSer2x2v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 

double SigmaCXcellAvSer2x3v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 


double SigmaCXcellAvSer3x3v_P3(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX); 



 
} 
#endif 
