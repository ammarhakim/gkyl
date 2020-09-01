#include <IonizationModDecl.h> 
#include <math.h> 
double VoronovReactRateCellAv2xSer_P1(const double elemCharge, const double m_, const double *m0, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	  Voronov ionization energy. 
  // A:   	  Voronov constant. 
  // K:   	  Voronov constant. 
  // P:   	  Voronov constant. 
  // X:   	  Voronov constant. 
  // m_:          mass of electron. 
  // m0[4]:      neutral density. 
  // vtSq[4]:    elc squared thermal speed, sqrt(T/m). 
  // coefIz[4]:  ionization reaction rate. 
 
  double m0NeutAv = 0.5*m0[0]; 
  double vtSq0 = 0.5*vtSq[0]; 
  double T0 = (0.5*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0) { 
    coefIz[0] = 0.0;
  }
  return 0.08333333333333333*coefIz[0]*m0[0]; 
} 
 
double VoronovReactRateCellAv2xSer_P2(const double elemCharge, const double m_, const double *m0, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	  Voronov ionization energy. 
  // A:   	  Voronov constant. 
  // K:   	  Voronov constant. 
  // P:   	  Voronov constant. 
  // X:   	  Voronov constant. 
  // m_:          mass of electron. 
  // m0[8]:      neutral density. 
  // vtSq[8]:    elc squared thermal speed, sqrt(T/m). 
  // coefIz[8]:  ionization reaction rate. 
 
  double m0NeutAv = 0.5*m0[0]; 
  double vtSq0 = 0.5*vtSq[0]; 
  double T0 = (0.5*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0) { 
    coefIz[0] = 0.0;
  }
  return 0.05*coefIz[0]*m0[0]; 
} 
 
double VoronovReactRateCellAv2xSer_P3(const double elemCharge, const double m_, const double *m0, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	  Voronov ionization energy. 
  // A:   	  Voronov constant. 
  // K:   	  Voronov constant. 
  // P:   	  Voronov constant. 
  // X:   	  Voronov constant. 
  // m_:          mass of electron. 
  // m0[12]:      neutral density. 
  // vtSq[12]:    elc squared thermal speed, sqrt(T/m). 
  // coefIz[12]:  ionization reaction rate. 
 
  double m0NeutAv = 0.5*m0[0]; 
  double vtSq0 = 0.5*vtSq[0]; 
  double T0 = (0.5*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0) { 
    coefIz[0] = 0.0;
  }
  return 0.03571428571428571*coefIz[0]*m0[0]; 
} 
 
