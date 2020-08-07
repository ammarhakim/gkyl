#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void SigmaCXcellAvMax2x2v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uIon[6]:        ion fluid velocity. 
  // uNeut[6]:       neutral fluid velocity. 
  // vtSqIon[3]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[3]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uNeut[3],2)-0.5*uIon[3]*uNeut[3]+0.25*pow(uIon[3],2)+0.25*pow(uNeut[0],2)-0.5*uIon[0]*uNeut[0]+0.25*pow(uIon[0],2); 
  vSigmaCX[0] = 1.128379167095513*a*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)-1.128379167095513*b*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)*log(0.5641895835477563*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)); 
 
} 
void SigmaCXcellAvMax2x2v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uIon[12]:        ion fluid velocity. 
  // uNeut[12]:       neutral fluid velocity. 
  // vtSqIon[6]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[6]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uNeut[6],2)-0.5*uIon[6]*uNeut[6]+0.25*pow(uIon[6],2)+0.25*pow(uNeut[0],2)-0.5*uIon[0]*uNeut[0]+0.25*pow(uIon[0],2); 
  vSigmaCX[0] = 1.128379167095513*a*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)-1.128379167095513*b*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)*log(0.5641895835477563*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)); 
 
} 
