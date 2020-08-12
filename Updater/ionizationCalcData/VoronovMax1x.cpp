#include <IonizationModDecl.h> 
#include <math.h> 
double VoronovReactRateCellAv1xMax_P1(const double elemCharge, const double m_, const double *m0, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	  Voronov ionization energy. 
  // A:   	  Voronov constant. 
  // K:   	  Voronov constant. 
  // P:   	  Voronov constant. 
  // X:   	  Voronov constant. 
  // m_:          mass of electron. 
  // m0[2]:      neutral density. 
  // vtSq[2]:    elc squared thermal speed, sqrt(T/m). 
  // coefIz[2]:  ionization reaction rate. 
 
  double vtSq0 = 0.7071067811865476*vtSq[0]; 
  double T0 = (0.7071067811865476*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(1000000.0*X*exp(U)+1000000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(1000000.0*X*exp(U)+1000000.0*U*exp(U)); 
 
  if (U > 3.0/2.0) { 
    coefIz[0] = 0.0;
  }
  return 0.1666666666666667*coefIz[0]*m0[0]; 
} 
 
double VoronovReactRateCellAv1xMax_P2(const double elemCharge, const double m_, const double *m0, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	  Voronov ionization energy. 
  // A:   	  Voronov constant. 
  // K:   	  Voronov constant. 
  // P:   	  Voronov constant. 
  // X:   	  Voronov constant. 
  // m_:          mass of electron. 
  // m0[3]:      neutral density. 
  // vtSq[3]:    elc squared thermal speed, sqrt(T/m). 
  // coefIz[3]:  ionization reaction rate. 
 
  double vtSq0 = 0.7071067811865476*vtSq[0]; 
  double T0 = (0.7071067811865476*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(1000000.0*X*exp(U)+1000000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(1000000.0*X*exp(U)+1000000.0*U*exp(U)); 
 
  if (U > 3.0/2.0) { 
    coefIz[0] = 0.0;
  }
  return 0.1*coefIz[0]*m0[0]; 
} 
 
double VoronovReactRateCellAv1xMax_P3(const double elemCharge, const double m_, const double *m0, const double *vtSq, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
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
 
  double vtSq0 = 0.7071067811865476*vtSq[0]; 
  double T0 = (0.7071067811865476*vtSq[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(1000000.0*X*exp(U)+1000000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(1000000.0*X*exp(U)+1000000.0*U*exp(U)); 
 
  if (U > 3.0/2.0) { 
    coefIz[0] = 0.0;
  }
  return 0.07142857142857142*coefIz[0]*m0[0]; 
} 
 
