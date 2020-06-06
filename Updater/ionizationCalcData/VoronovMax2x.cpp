#include <IonizationModDecl.h> 
#include <math.h> 
void VoronovReactRateCellAv2xMax_P1(const double elemCharge, const double m_, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	  Voronov ionization energy. 
  // A:   	  Voronov constant. 
  // K:   	  Voronov constant. 
  // P:   	  Voronov constant. 
  // X:   	  Voronov constant. 
  // m_:          mass of electron. 
  // vtSq[3]:    squared thermal speed, sqrt(T/m) 
  // coefIz[3]:  ionization reaction rate. 
 
  double vtSq0 = 0.5*vtSq[0]; 
  double T0 = (0.5*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
 
  if (U > 3.0/2.0) { 
    coefIz[0] = 0.0;
  }
} 
 
void VoronovReactRateCellAv2xMax_P2(const double elemCharge, const double m_, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	  Voronov ionization energy. 
  // A:   	  Voronov constant. 
  // K:   	  Voronov constant. 
  // P:   	  Voronov constant. 
  // X:   	  Voronov constant. 
  // m_:          mass of electron. 
  // vtSq[6]:    squared thermal speed, sqrt(T/m) 
  // coefIz[6]:  ionization reaction rate. 
 
  double vtSq0 = 0.5*vtSq[0]; 
  double T0 = (0.5*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
 
  if (U > 3.0/2.0) { 
    coefIz[0] = 0.0;
  }
} 
 
void VoronovReactRateCellAv2xMax_P3(const double elemCharge, const double m_, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	  Voronov ionization energy. 
  // A:   	  Voronov constant. 
  // K:   	  Voronov constant. 
  // P:   	  Voronov constant. 
  // X:   	  Voronov constant. 
  // m_:          mass of electron. 
  // vtSq[10]:    squared thermal speed, sqrt(T/m) 
  // coefIz[10]:  ionization reaction rate. 
 
  double vtSq0 = 0.5*vtSq[0]; 
  double T0 = (0.5*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (A*P*pow(U,K+1/2))/(500000.0*X*exp(U)+500000.0*U*exp(U))+(A*pow(U,K))/(500000.0*X*exp(U)+500000.0*U*exp(U)); 
 
  if (U > 3.0/2.0) { 
    coefIz[0] = 0.0;
  }
} 
 
