#include <ChargeExchangeModDecl.h> 
#include <math.h> 
double SigmaCXcellAvSer2x2v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX) 
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
  if ((vtSqIonAv > 0.) && (vtSqIonAv < vtSqIonMin)) vtSqIonAv = vtSqIonMin;
  if ((vtSqNeutAv > 0.) && (vtSqNeutAv < vtSqNeutMin)) vtSqNeutAv = vtSqNeutMin;
  
  if (m0NeutAv <= 0 || vtSqNeutAv <= 0 || vtSqIonAv <= 0) { 
    vSigmaCX[0] = 0.0;
    return 0.0; 
  } else {
  double vINSqAv = 0.25*pow(uNeut[4],2)-0.5*uIon[4]*uNeut[4]+0.25*pow(uIon[4],2)+0.25*pow(uNeut[0],2)-0.5*uIon[0]*uNeut[0]+0.25*pow(uIon[0],2); 
 
  double Vcx = sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv);
  vSigmaCX[0] = 2.0*Vcx*a-2.0*Vcx*log(Vcx)*b; 
 
    return 0.08333333333333333*m0[0]*vSigmaCX[0]; 
  }
} 
double SigmaCXcellAvSer2x2v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX) 
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
  if ((vtSqIonAv > 0.) && (vtSqIonAv < vtSqIonMin)) vtSqIonAv = vtSqIonMin;
  if ((vtSqNeutAv > 0.) && (vtSqNeutAv < vtSqNeutMin)) vtSqNeutAv = vtSqNeutMin;
  
  if (m0NeutAv <= 0 || vtSqNeutAv <= 0 || vtSqIonAv <= 0) { 
    vSigmaCX[0] = 0.0;
    return 0.0; 
  } else {
  double vINSqAv = 0.25*pow(uNeut[8],2)-0.5*uIon[8]*uNeut[8]+0.25*pow(uIon[8],2)+0.25*pow(uNeut[0],2)-0.5*uIon[0]*uNeut[0]+0.25*pow(uIon[0],2); 
 
  double Vcx = sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv);
  vSigmaCX[0] = 2.0*Vcx*a-2.0*Vcx*log(Vcx)*b; 
 
    return 0.05*m0[0]*vSigmaCX[0]; 
  }
} 
