#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void SigmaCXcellAvMax3x3v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0[4]:         neutral particle density. 
  // uIon[12]:        ion fluid velocity. 
  // uNeut[12]:       neutral fluid velocity. 
  // vtSqIon[4]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[4]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double m0NeutAv = 0.3535533905932738*m0[0]; 
  double vtSqIonAv = 0.3535533905932738*vtSqIon[0]; 
  double vtSqNeutAv = 0.3535533905932738*vtSqNeut[0]; 
  double vINSqAv = 0.125*pow(uNeut[8],2)-0.25*uIon[8]*uNeut[8]+0.125*pow(uIon[8],2)+0.125*pow(uNeut[4],2)-0.25*uIon[4]*uNeut[4]+0.125*pow(uIon[4],2)+0.125*pow(uNeut[0],2)-0.25*uIon[0]*uNeut[0]+0.125*pow(uIon[0],2); 
 
  vSigmaCX[0] = 2.828427124746191*a*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)-2.828427124746191*b*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)*log(sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)); 
 
  if (m0NeutAv <= 0) { 
    vSigmaCX[0] = 0.0;
  }
} 
void SigmaCXcellAvMax3x3v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0[10]:         neutral particle density. 
  // uIon[30]:        ion fluid velocity. 
  // uNeut[30]:       neutral fluid velocity. 
  // vtSqIon[10]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[10]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double m0NeutAv = 0.3535533905932738*m0[0]; 
  double vtSqIonAv = 0.3535533905932738*vtSqIon[0]; 
  double vtSqNeutAv = 0.3535533905932738*vtSqNeut[0]; 
  double vINSqAv = 0.125*pow(uNeut[20],2)-0.25*uIon[20]*uNeut[20]+0.125*pow(uIon[20],2)+0.125*pow(uNeut[10],2)-0.25*uIon[10]*uNeut[10]+0.125*pow(uIon[10],2)+0.125*pow(uNeut[0],2)-0.25*uIon[0]*uNeut[0]+0.125*pow(uIon[0],2); 
 
  vSigmaCX[0] = 2.828427124746191*a*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)-2.828427124746191*b*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)*log(sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)); 
 
  if (m0NeutAv <= 0) { 
    vSigmaCX[0] = 0.0;
  }
} 
