#include <IonizationModDecl.h> 
#include <math.h> 
void IonizationTemp2xSer_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[4]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[4]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.5*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.6666666666666666*E*elemCharge)/m_; 
  }
 
  else { 
     vtSqIz[0] = 2.0e-10; 
  }
 
} 
 
void IonizationTemp2xSer_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[8]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[8]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.5*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.6666666666666666*E*elemCharge)/m_; 
  }
 
  else { 
     vtSqIz[0] = 2.0e-10; 
  }
 
} 
 
void IonizationTemp2xSer_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[12]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[12]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.5*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.6666666666666666*E*elemCharge)/m_; 
  }
 
  else { 
     vtSqIz[0] = 2.0e-10; 
  }
 
} 
 
