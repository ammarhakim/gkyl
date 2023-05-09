#include <ChargeExchangeModDecl.h> 
#include <math.h> 
double SigmaCXcellAvSer1x3v_P1(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0[2]:         neutral particle density. 
  // uIon[6]:        ion fluid velocity. 
  // uNeut[6]:       neutral fluid velocity. 
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
    return 0.0; 
  } else {
  double vINSqAv = 0.5*pow(uNeut[4],2)-1.0*uIon[4]*uNeut[4]+0.5*pow(uIon[4],2)+0.5*pow(uNeut[2],2)-1.0*uIon[2]*uNeut[2]+0.5*pow(uIon[2],2)+0.5*pow(uNeut[0],2)-1.0*uIon[0]*uNeut[0]+0.5*pow(uIon[0],2); 
 
  double Vcx = sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv);
  vSigmaCX[0] = 1.414213562373095*Vcx*a-1.414213562373095*Vcx*log(Vcx)*b; 
 
    return 0.1666666666666667*m0[0]*vSigmaCX[0]; 
  }
} 
double SigmaCXcellAvSer1x3v_P2(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX) 
{ 
  // a               constant in fitting function. 
  // b               constant in fitting function. 
  // m0[3]:         neutral particle density. 
  // uIon[9]:        ion fluid velocity. 
  // uNeut[9]:       neutral fluid velocity. 
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
    return 0.0; 
  } else {
  double vINSqAv = 0.5*pow(uNeut[6],2)-1.0*uIon[6]*uNeut[6]+0.5*pow(uIon[6],2)+0.5*pow(uNeut[3],2)-1.0*uIon[3]*uNeut[3]+0.5*pow(uIon[3],2)+0.5*pow(uNeut[0],2)-1.0*uIon[0]*uNeut[0]+0.5*pow(uIon[0],2); 
 
  double Vcx = sqrt(1.273239544735163*vtSqNeutAv+1.273239544735163*vtSqIonAv+vINSqAv);
  vSigmaCX[0] = 1.414213562373095*Vcx*a-1.414213562373095*Vcx*log(Vcx)*b; 
 
    return 0.1*m0[0]*vSigmaCX[0]; 
  }
} 
