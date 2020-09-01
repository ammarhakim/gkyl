#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void SigmaCXcellAvSer2x2v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0[4]:         neutral particle density. 
  // uIon[8]:        ion fluid velocity. 
  // uNeut[8]:       neutral fluid velocity. 
  // vtSqIon[4]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[4]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double m0NeutAv = 0.5*m0[0]; 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uNeut[4],2)-0.5*uIon[4]*uNeut[4]+0.25*pow(uIon[4],2)+0.25*pow(uNeut[0],2)-0.5*uIon[0]*uNeut[0]+0.25*pow(uIon[0],2); 
 
  vSigmaCX[0] = 2.0*a*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)-2.0*b*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)*log(sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)); 
 
  if (m0NeutAv <= 0) { 
    vSigmaCX[0] = 0.0;
  }
} 
void SigmaCXcellAvSer2x2v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0[8]:         neutral particle density. 
  // uIon[16]:        ion fluid velocity. 
  // uNeut[16]:       neutral fluid velocity. 
  // vtSqIon[8]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[8]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double m0NeutAv = 0.5*m0[0]; 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uNeut[8],2)-0.5*uIon[8]*uNeut[8]+0.25*pow(uIon[8],2)+0.25*pow(uNeut[0],2)-0.5*uIon[0]*uNeut[0]+0.25*pow(uIon[0],2); 
 
  vSigmaCX[0] = 2.0*a*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)-2.0*b*sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)*log(sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)); 
 
  if (m0NeutAv <= 0) { 
    vSigmaCX[0] = 0.0;
  }
} 
