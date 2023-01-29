#include <IonizationModDecl.h> 
#include <math.h> 
double VoronovReactRateCellAv3xSer_P1(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, double vtSqNeutMin, const double *vtSqElc, double vtSqElcMin, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
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
 
  double m0NeutAv = 0.3535533905932738*m0[0]; 
  double vtSqElc0 = 0.3535533905932738*vtSqElc[0]; 
  double vtSqNeut0 = 0.3535533905932738*vtSqNeut[0]; 
  if ((vtSqNeut0 > 0.) && (vtSqNeut0 < vtSqNeutMin)) vtSqNeut0 = vtSqNeutMin;
  if ((vtSqElc0 > 0.) && (vtSqElc0 < vtSqElcMin)) vtSqElc0 = vtSqElcMin;
 
  double T0 = (m_*vtSqElc0)/elemCharge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0. || vtSqNeut0 <= 0. || vtSqElc0 <= 0.) { 
    coefIz[0] = 0.0;
    return 0.0; 
  } else {
    coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
    return 0.04166666666666666*coefIz[0]*m0[0]; 
  };
} 
 
double VoronovReactRateCellAv3xSer_P2(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, double vtSqNeutMin, const double *vtSqElc, double vtSqElcMin, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	     Voronov ionization energy. 
  // A:   	     Voronov constant. 
  // K:   	     Voronov constant. 
  // P:   	     Voronov constant. 
  // X:   	     Voronov constant. 
  // m_:             mass of electron. 
  // m0[20]:         neutral density. 
  // vtSqElc[20]:    electron squared thermal speed, sqrt(T/m). 
  // vtSqNeut[20]:   neutral squared thermal speed, sqrt(T/m). 
  // coefIz[20]:     ionization reaction rate. 
 
  double m0NeutAv = 0.3535533905932738*m0[0]; 
  double vtSqElc0 = 0.3535533905932738*vtSqElc[0]; 
  double vtSqNeut0 = 0.3535533905932738*vtSqNeut[0]; 
  if ((vtSqNeut0 > 0.) && (vtSqNeut0 < vtSqNeutMin)) vtSqNeut0 = vtSqNeutMin;
  if ((vtSqElc0 > 0.) && (vtSqElc0 < vtSqElcMin)) vtSqElc0 = vtSqElcMin;
 
  double T0 = (m_*vtSqElc0)/elemCharge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0. || vtSqNeut0 <= 0. || vtSqElc0 <= 0.) { 
    coefIz[0] = 0.0;
    return 0.0; 
  } else {
    coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
    return 0.025*coefIz[0]*m0[0]; 
  };
} 
 
double VoronovReactRateCellAv3xSer_P3(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, double vtSqNeutMin, const double *vtSqElc, double vtSqElcMin, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	     Voronov ionization energy. 
  // A:   	     Voronov constant. 
  // K:   	     Voronov constant. 
  // P:   	     Voronov constant. 
  // X:   	     Voronov constant. 
  // m_:             mass of electron. 
  // m0[32]:         neutral density. 
  // vtSqElc[32]:    electron squared thermal speed, sqrt(T/m). 
  // vtSqNeut[32]:   neutral squared thermal speed, sqrt(T/m). 
  // coefIz[32]:     ionization reaction rate. 
 
  double m0NeutAv = 0.3535533905932738*m0[0]; 
  double vtSqElc0 = 0.3535533905932738*vtSqElc[0]; 
  double vtSqNeut0 = 0.3535533905932738*vtSqNeut[0]; 
  if ((vtSqNeut0 > 0.) && (vtSqNeut0 < vtSqNeutMin)) vtSqNeut0 = vtSqNeutMin;
  if ((vtSqElc0 > 0.) && (vtSqElc0 < vtSqElcMin)) vtSqElc0 = vtSqElcMin;
 
  double T0 = (m_*vtSqElc0)/elemCharge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0. || vtSqNeut0 <= 0. || vtSqElc0 <= 0.) { 
    coefIz[0] = 0.0;
    return 0.0; 
  } else {
    coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
    return 0.01785714285714286*coefIz[0]*m0[0]; 
  };
} 
 
