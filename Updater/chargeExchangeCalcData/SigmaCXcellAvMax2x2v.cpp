#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void SigmaCXcellAvMax2x2v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0[3]:         neutral particle density. 
  // uIon[6]:        ion fluid velocity. 
  // uNeut[6]:       neutral fluid velocity. 
  // vtSqIon[3]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[3]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double m0NeutAv = 0.5*m0[0]; 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uNeut[3],2)-0.5*uIon[3]*uNeut[3]+0.25*pow(uIon[3],2)+0.25*pow(uNeut[0],2)-0.5*uIon[0]*uNeut[0]+0.25*pow(uIon[0],2); 
 
  vSigmaCX[0] = 2.0*a*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)-2.0*b*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)*log(sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)); 
 
  if (m0NeutAv <= 0) { 
    vSigmaCX[0] = 0.0;
  }
} 
void SigmaCXcellAvMax2x2v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0[6]:         neutral particle density. 
  // uIon[12]:        ion fluid velocity. 
  // uNeut[12]:       neutral fluid velocity. 
  // vtSqIon[6]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[6]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double m0NeutAv = 0.5*m0[0]; 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uNeut[6],2)-0.5*uIon[6]*uNeut[6]+0.25*pow(uIon[6],2)+0.25*pow(uNeut[0],2)-0.5*uIon[0]*uNeut[0]+0.25*pow(uIon[0],2); 
 
  vSigmaCX[0] = 2.0*a*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)-2.0*b*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)*log(sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)); 
 
  if (m0NeutAv <= 0) { 
    vSigmaCX[0] = 0.0;
  }
} 
