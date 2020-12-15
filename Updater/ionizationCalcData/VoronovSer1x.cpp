#include <IonizationModDecl.h> 
#include <math.h> 
double VoronovReactRateCellAv1xSer_P1(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	     Voronov ionization energy. 
  // A:   	     Voronov constant. 
  // K:   	     Voronov constant. 
  // P:   	     Voronov constant. 
  // X:   	     Voronov constant. 
  // m_:             mass of electron. 
  // m0[2]:         neutral density. 
  // vtSqElc[2]:    electron squared thermal speed, sqrt(T/m). 
  // vtSqNeut[2]:   neutral squared thermal speed, sqrt(T/m). 
  // coefIz[2]:     ionization reaction rate. 
 
  double m0NeutAv = 0.7071067811865476*m0[0]; 
  double vtSqElc0 = 0.7071067811865476*vtSqElc[0]; 
  double vtSqNeut0 = 0.7071067811865476*vtSqNeut[0]; 
  double T0 = (0.7071067811865476*vtSqElc[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0 || vtSqNeut0 <= 0 || vtSqElc0 <= 0) { 
    coefIz[0] = 0.0;
    return 0.0;
  }
  else { 
    coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(1000000.0*X*exp(U)+1000000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(1000000.0*X*exp(U)+1000000.0*U*exp(U)); 
    return 0.1666666666666667*coefIz[0]*m0[0]; 
  }
 
} 
 
double VoronovReactRateCellAv1xSer_P2(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // E:   	     Voronov ionization energy. 
  // A:   	     Voronov constant. 
  // K:   	     Voronov constant. 
  // P:   	     Voronov constant. 
  // X:   	     Voronov constant. 
  // m_:             mass of electron. 
  // m0[3]:         neutral density. 
  // vtSqElc[3]:    electron squared thermal speed, sqrt(T/m). 
  // vtSqNeut[3]:   neutral squared thermal speed, sqrt(T/m). 
  // coefIz[3]:     ionization reaction rate. 
 
  double m0NeutAv = 0.7071067811865476*m0[0]; 
  double vtSqElc0 = 0.7071067811865476*vtSqElc[0]; 
  double vtSqNeut0 = 0.7071067811865476*vtSqNeut[0]; 
  double T0 = (0.7071067811865476*vtSqElc[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0 || vtSqNeut0 <= 0 || vtSqElc0 <= 0) { 
    coefIz[0] = 0.0;
    return 0.0;
  }
  else { 
    coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(1000000.0*X*exp(U)+1000000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(1000000.0*X*exp(U)+1000000.0*U*exp(U)); 
    return 0.1*coefIz[0]*m0[0]; 
  }
 
} 
 
double VoronovReactRateCellAv1xSer_P3(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz) 
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
 
  double m0NeutAv = 0.7071067811865476*m0[0]; 
  double vtSqElc0 = 0.7071067811865476*vtSqElc[0]; 
  double vtSqNeut0 = 0.7071067811865476*vtSqNeut[0]; 
  double T0 = (0.7071067811865476*vtSqElc[0]*m_)/elemCharge; 
  double U = E/T0; 
 
  if (U >= 3.0/2.0 || m0NeutAv <= 0 || vtSqNeut0 <= 0 || vtSqElc0 <= 0) { 
    coefIz[0] = 0.0;
    return 0.0;
  }
  else { 
    coefIz[0] = (1.414213562373095*A*P*pow(U,K+1/2))/(1000000.0*X*exp(U)+1000000.0*U*exp(U))+(1.414213562373095*A*pow(U,K))/(1000000.0*X*exp(U)+1000000.0*U*exp(U)); 
    return 0.07142857142857142*coefIz[0]*m0[0]; 
  }
 
} 
 
