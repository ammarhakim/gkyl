#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void SigmaCXcellAvSer3x3v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uIon[24]:        ion fluid velocity. 
  // uNeut[24]:       neutral fluid velocity. 
  // vtSqIon[8]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[8]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.3535533905932738*vtSqIon[0]; 
  double vtSqNeutAv = 0.3535533905932738*vtSqNeut[0]; 
  double vINSqAv = 0.125*pow(uNeut[16],2)-0.25*uIon[16]*uNeut[16]+0.125*pow(uIon[16],2)+0.125*pow(uNeut[8],2)-0.25*uIon[8]*uNeut[8]+0.125*pow(uIon[8],2)+0.125*pow(uNeut[0],2)-0.25*uIon[0]*uNeut[0]+0.125*pow(uIon[0],2); 
  vSigmaCX[0] = 1.595769121605731*a*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)-1.595769121605731*b*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)*log(0.5641895835477563*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)); 
 
} 
void SigmaCXcellAvSer3x3v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uIon[60]:        ion fluid velocity. 
  // uNeut[60]:       neutral fluid velocity. 
  // vtSqIon[20]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[20]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.3535533905932738*vtSqIon[0]; 
  double vtSqNeutAv = 0.3535533905932738*vtSqNeut[0]; 
  double vINSqAv = 0.125*pow(uNeut[40],2)-0.25*uIon[40]*uNeut[40]+0.125*pow(uIon[40],2)+0.125*pow(uNeut[20],2)-0.25*uIon[20]*uNeut[20]+0.125*pow(uIon[20],2)+0.125*pow(uNeut[0],2)-0.25*uIon[0]*uNeut[0]+0.125*pow(uIon[0],2); 
  vSigmaCX[0] = 1.595769121605731*a*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)-1.595769121605731*b*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)*log(0.5641895835477563*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)); 
 
} 
