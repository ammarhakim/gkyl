#include <IonizationModDecl.h> 
#include <math.h> 
double VoronovReactRateCellAv2xSer_P1(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	     Voronov ionization energy. 
  // A:   	     Voronov constant. 
  // K:   	     Voronov constant. 
  // P:   	     Voronov constant. 
  // X:   	     Voronov constant. 
  // m_:             mass of electron. 
  // m0[4]:         neutral density. 
  // vtSqElc[4]:    electron squared thermal speed, sqrt(T/m). 
  // vtSqNeut[4]:   neutral squared thermal speed, sqrt(T/m). 
  // coefIz[4]:     ionization reaction rate. 
 
  double m0NeutAv = 0.5*m0[0]; 
  double vtSqElc0 = 0.5*vtSqElc[0]; 
  double vtSqNeut0 = 0.5*vtSqNeut[0]; 
  double T0 = (0.5*vtSqElc[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0 || vtSqNeut0 <= 0 || vtSqElc0 <= 0) { 
    coefIz[0] = 0.0;
    return 0.0;
  }
  else { 
    coefIz[0] = (A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
    return 0.08333333333333333*coefIz[0]*m0[0]; 
  }
 
} 
 
double VoronovReactRateCellAv2xSer_P2(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	     Voronov ionization energy. 
  // A:   	     Voronov constant. 
  // K:   	     Voronov constant. 
  // P:   	     Voronov constant. 
  // X:   	     Voronov constant. 
  // m_:             mass of electron. 
  // m0[8]:         neutral density. 
  // vtSqElc[8]:    electron squared thermal speed, sqrt(T/m). 
  // vtSqNeut[8]:   neutral squared thermal speed, sqrt(T/m). 
  // coefIz[8]:     ionization reaction rate. 
 
  double m0NeutAv = 0.5*m0[0]; 
  double vtSqElc0 = 0.5*vtSqElc[0]; 
  double vtSqNeut0 = 0.5*vtSqNeut[0]; 
  double T0 = (0.5*vtSqElc[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0 || vtSqNeut0 <= 0 || vtSqElc0 <= 0) { 
    coefIz[0] = 0.0;
    return 0.0;
  }
  else { 
    coefIz[0] = (A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
    return 0.05*coefIz[0]*m0[0]; 
  }
 
} 
 
double VoronovReactRateCellAv2xSer_P3(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	     Voronov ionization energy. 
  // A:   	     Voronov constant. 
  // K:   	     Voronov constant. 
  // P:   	     Voronov constant. 
  // X:   	     Voronov constant. 
  // m_:             mass of electron. 
  // m0[12]:         neutral density. 
  // vtSqElc[12]:    electron squared thermal speed, sqrt(T/m). 
  // vtSqNeut[12]:   neutral squared thermal speed, sqrt(T/m). 
  // coefIz[12]:     ionization reaction rate. 
 
  double m0NeutAv = 0.5*m0[0]; 
  double vtSqElc0 = 0.5*vtSqElc[0]; 
  double vtSqNeut0 = 0.5*vtSqNeut[0]; 
  double T0 = (0.5*vtSqElc[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0 || vtSqNeut0 <= 0 || vtSqElc0 <= 0) { 
    coefIz[0] = 0.0;
    return 0.0;
  }
  else { 
    coefIz[0] = (A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
    return 0.03571428571428571*coefIz[0]*m0[0]; 
  }
 
} 
 
