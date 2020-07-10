#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void GkSigmaCXcellAvSer1x1v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uParIon[2]:    ion fluid velocity. 
  // uParNeut[2]:   neutral fluid velocity. 
  // vtSqIon[2]:    ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[2]:   neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:       cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.7071067811865476*vtSqIon[0]; 
  double vtSqNeutAv = 0.7071067811865476*vtSqNeut[0]; 
  double vINSqAv = 0.5*pow(uParNeut[0],2)-1.0*uParIon[0]*uParNeut[0]+0.5*pow(uParIon[0],2); 
  vSigmaCX[0] = 0.7978845608028654*a*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)-0.7978845608028654*b*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)*log(0.5641895835477563*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)); 
 
} 
void GkSigmaCXcellAvSer1x1v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uParIon[3]:    ion fluid velocity. 
  // uParNeut[3]:   neutral fluid velocity. 
  // vtSqIon[3]:    ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[3]:   neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:       cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.7071067811865476*vtSqIon[0]; 
  double vtSqNeutAv = 0.7071067811865476*vtSqNeut[0]; 
  double vINSqAv = 0.5*pow(uParNeut[0],2)-1.0*uParIon[0]*uParNeut[0]+0.5*pow(uParIon[0],2); 
  vSigmaCX[0] = 0.7978845608028654*a*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)-0.7978845608028654*b*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)*log(0.5641895835477563*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)); 
 
} 
void GkSigmaCXcellAvSer1x1v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // uParIon[4]:    ion fluid velocity. 
  // uParNeut[4]:   neutral fluid velocity. 
  // vtSqIon[4]:    ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[4]:   neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:       cell ave cross section fitting eqn. 
 
  double vtSqIonAv = 0.7071067811865476*vtSqIon[0]; 
  double vtSqNeutAv = 0.7071067811865476*vtSqNeut[0]; 
  double vINSqAv = 0.5*pow(uParNeut[0],2)-1.0*uParIon[0]*uParNeut[0]+0.5*pow(uParIon[0],2); 
  vSigmaCX[0] = 0.7978845608028654*a*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)-0.7978845608028654*b*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)*log(0.5641895835477563*sqrt(4.0*vtSqNeutAv+4.0*vtSqIonAv+3.141592653589793*vINSqAv)); 
 
} 
