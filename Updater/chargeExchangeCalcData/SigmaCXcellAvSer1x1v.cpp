#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void SigmaCXcellAvSer1x1v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0[2]:         neutral particle density. 
  // uIon[2]:        ion fluid velocity. 
  // uNeut[2]:       neutral fluid velocity. 
  // vtSqIon[2]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[2]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double m0NeutAv = 0.7071067811865476*m0[0]; 
  double vtSqIonAv = 0.7071067811865476*vtSqIon[0]; 
  double vtSqNeutAv = 0.7071067811865476*vtSqNeut[0]; 
  if ((vtSqIonAv > 0.) && (vtSqIonAv < vtSqIonMin)) vtSqIonAv = vtSqIonMin;
  if ((vtSqNeutAv > 0.) && (vtSqNeutAv < vtSqNeutMin)) vtSqNeutAv = vtSqNeutMin;
  
  if (m0NeutAv <= 0 || vtSqNeutAv <= 0 || vtSqIonAv <= 0) { 
    vSigmaCX[0] = 0.0;
  } else {
  double vINSqAv = 0.5*pow(uNeut[0],2)-1.0*uIon[0]*uNeut[0]+0.5*pow(uIon[0],2); 
 
  double Vcx = sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv);
  vSigmaCX[0] = 1.414213562373095*Vcx*a-1.414213562373095*Vcx*log(Vcx)*b; 
 
  }
} 
void SigmaCXcellAvSer1x1v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0[3]:         neutral particle density. 
  // uIon[3]:        ion fluid velocity. 
  // uNeut[3]:       neutral fluid velocity. 
  // vtSqIon[3]:     ion squared thermal speed, sqrt(T/m). 
  // vtSqNeut[3]:    neutral squared thermal speed, sqrt(T/m). 
  // vSigmaCX:          cell ave cross section fitting eqn. 
 
  double m0NeutAv = 0.7071067811865476*m0[0]; 
  double vtSqIonAv = 0.7071067811865476*vtSqIon[0]; 
  double vtSqNeutAv = 0.7071067811865476*vtSqNeut[0]; 
  if ((vtSqIonAv > 0.) && (vtSqIonAv < vtSqIonMin)) vtSqIonAv = vtSqIonMin;
  if ((vtSqNeutAv > 0.) && (vtSqNeutAv < vtSqNeutMin)) vtSqNeutAv = vtSqNeutMin;
  
  if (m0NeutAv <= 0 || vtSqNeutAv <= 0 || vtSqIonAv <= 0) { 
    vSigmaCX[0] = 0.0;
  } else {
  double vINSqAv = 0.5*pow(uNeut[0],2)-1.0*uIon[0]*uNeut[0]+0.5*pow(uIon[0],2); 
 
  double Vcx = sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv);
  vSigmaCX[0] = 1.414213562373095*Vcx*a-1.414213562373095*Vcx*log(Vcx)*b; 
 
  }
} 
