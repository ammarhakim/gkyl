#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void GkSigmaCXcellAvMax2x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uParIon[6]:    ion fluid velocity. 
  // uParNeut[6]:    neutral fluid velocity. 
  // vtSqIon[3]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[3]:    neutral squared thermal speed, sqrt(T/m). 
  // sigmaCX:          cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uParNeut[0],2)-0.5*uParIon[0]*uParNeut[0]+0.25*pow(uParIon[0],2); 
  sigmaCX[0] = 2.0*a-2.0*b*log(0.5641895835477563*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)); 
 
} 
void GkSigmaCXcellAvMax2x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uParIon[12]:    ion fluid velocity. 
  // uParNeut[12]:    neutral fluid velocity. 
  // vtSqIon[6]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[6]:    neutral squared thermal speed, sqrt(T/m). 
  // sigmaCX:          cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uParNeut[0],2)-0.5*uParIon[0]*uParNeut[0]+0.25*pow(uParIon[0],2); 
  sigmaCX[0] = 2.0*a-2.0*b*log(0.5641895835477563*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)); 
 
} 
void GkSigmaCXcellAvMax2x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uParIon[20]:    ion fluid velocity. 
  // uParNeut[20]:    neutral fluid velocity. 
  // vtSqIon[10]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[10]:    neutral squared thermal speed, sqrt(T/m). 
  // sigmaCX:          cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.5*vtSqIon[0]; 
  double vtSqNeutAv = 0.5*vtSqNeut[0]; 
  double vINSqAv = 0.25*pow(uParNeut[0],2)-0.5*uParIon[0]*uParNeut[0]+0.25*pow(uParIon[0],2); 
  sigmaCX[0] = 2.0*a-2.0*b*log(0.5641895835477563*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)); 
 
} 
