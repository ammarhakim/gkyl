#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void VmSigmaCXcellAvMax2x3v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uIon[9]:        ion fluid velocity. 
  // uNeut[9]:       neutral fluid velocity. 
  // vtSqIon[3]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[3]:    neutral squared thermal speed, sqrt(T/m). 
  // sigmaCX:          cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uNeut[6],2)-0.5*uIon[6]*uNeut[6]+0.25*pow(uIon[6],2)+0.25*pow(uNeut[3],2)-0.5*uIon[3]*uNeut[3]+0.25*pow(uIon[3],2)+0.25*pow(uNeut[0],2)-0.5*uIon[0]*uNeut[0]+0.25*pow(uIon[0],2); 
  sigmaCX[0] = 2.0*a-2.0*b*log(sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)); 
 
} 
void VmSigmaCXcellAvMax2x3v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uIon[18]:        ion fluid velocity. 
  // uNeut[18]:       neutral fluid velocity. 
  // vtSqIon[6]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[6]:    neutral squared thermal speed, sqrt(T/m). 
  // sigmaCX:          cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uNeut[12],2)-0.5*uIon[12]*uNeut[12]+0.25*pow(uIon[12],2)+0.25*pow(uNeut[6],2)-0.5*uIon[6]*uNeut[6]+0.25*pow(uIon[6],2)+0.25*pow(uNeut[0],2)-0.5*uIon[0]*uNeut[0]+0.25*pow(uIon[0],2); 
  sigmaCX[0] = 2.0*a-2.0*b*log(sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)); 
 
} 
void VmSigmaCXcellAvMax2x3v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uIon[30]:        ion fluid velocity. 
  // uNeut[30]:       neutral fluid velocity. 
  // vtSqIon[10]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[10]:    neutral squared thermal speed, sqrt(T/m). 
  // sigmaCX:          cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uNeut[20],2)-0.5*uIon[20]*uNeut[20]+0.25*pow(uIon[20],2)+0.25*pow(uNeut[10],2)-0.5*uIon[10]*uNeut[10]+0.25*pow(uIon[10],2)+0.25*pow(uNeut[0],2)-0.5*uIon[0]*uNeut[0]+0.25*pow(uIon[0],2); 
  sigmaCX[0] = 2.0*a-2.0*b*log(sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv)); 
 
} 
