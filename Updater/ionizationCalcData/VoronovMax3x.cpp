#include <IonizationModDecl.h> 
#include <math.h> 
double VoronovReactRateCellAv3xMax_P1(const double elemCharge, const double m_, const double *m0, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
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
 
  double m0NeutAv = 0.3535533905932738*m0[0]; 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  double T0 = (0.3535533905932738*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0) { 
    coefIz[0] = 0.0;
  }
  return 0.04166666666666666*coefIz[0]*m0[0]; 
} 
 
double VoronovReactRateCellAv3xMax_P2(const double elemCharge, const double m_, const double *m0, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	  Voronov ionization energy. 
  // A:   	  Voronov constant. 
  // K:   	  Voronov constant. 
  // P:   	  Voronov constant. 
  // X:   	  Voronov constant. 
  // m_:          mass of electron. 
  // m0[10]:      neutral density. 
  // vtSq[10]:    elc squared thermal speed, sqrt(T/m). 
  // coefIz[10]:  ionization reaction rate. 
 
  double m0NeutAv = 0.3535533905932738*m0[0]; 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  double T0 = (0.3535533905932738*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0) { 
    coefIz[0] = 0.0;
  }
  return 0.025*coefIz[0]*m0[0]; 
} 
 
double VoronovReactRateCellAv3xMax_P3(const double elemCharge, const double m_, const double *m0, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	  Voronov ionization energy. 
  // A:   	  Voronov constant. 
  // K:   	  Voronov constant. 
  // P:   	  Voronov constant. 
  // X:   	  Voronov constant. 
  // m_:          mass of electron. 
  // m0[20]:      neutral density. 
  // vtSq[20]:    elc squared thermal speed, sqrt(T/m). 
  // coefIz[20]:  ionization reaction rate. 
 
  double m0NeutAv = 0.3535533905932738*m0[0]; 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  double T0 = (0.3535533905932738*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0) { 
    coefIz[0] = 0.0;
  }
  return 0.01785714285714286*coefIz[0]*m0[0]; 
} 
 
