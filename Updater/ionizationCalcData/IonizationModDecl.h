#ifndef IONIZATION_MOD_DECL_H 
#define IONIZATION_MOD_DECL_H 
#include <cmath> 

#include <algorithm> 

extern "C" { 
double VoronovReactRateCellAv1xSer_P1(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp1xSer_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv2xSer_P1(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp2xSer_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv3xSer_P1(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp3xSer_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 


double VoronovReactRateCellAv1xSer_P2(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp1xSer_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv2xSer_P2(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp2xSer_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv3xSer_P2(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp3xSer_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 


double VoronovReactRateCellAv1xSer_P3(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp1xSer_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv2xSer_P3(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp2xSer_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv3xSer_P3(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp3xSer_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 


double VoronovReactRateCellAv1xMax_P1(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp1xMax_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv2xMax_P1(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp2xMax_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv3xMax_P1(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp3xMax_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 


double VoronovReactRateCellAv1xMax_P2(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp1xMax_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv2xMax_P2(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp2xMax_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv3xMax_P2(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp3xMax_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 


double VoronovReactRateCellAv1xMax_P3(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp1xMax_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv2xMax_P3(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp2xMax_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 

double VoronovReactRateCellAv3xMax_P3(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz); 

void IonizationTemp3xMax_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz); 



 
} 
#endif 
