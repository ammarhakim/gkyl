#include <IonizationModDecl.h> 
#include <math.h> 
void IonizationTemp2xMax_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[3]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[3]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.5*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.6666666666666666*E*elemCharge)/m_; 
  }
 
  else { 
     vtSqIz[0] = 2.0e-10; 
  }
 
} 
 
void IonizationTemp2xMax_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[6]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[6]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.5*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.6666666666666666*E*elemCharge)/m_; 
  }
 
  else { 
     vtSqIz[0] = 2.0e-10; 
  }
 
} 
 
void IonizationTemp2xMax_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[10]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[10]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.5*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.6666666666666666*E*elemCharge)/m_; 
  }
 
  else { 
     vtSqIz[0] = 2.0e-10; 
  }
 
} 
 
