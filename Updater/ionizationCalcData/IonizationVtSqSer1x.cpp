#include <IonizationModDecl.h> 
#include <math.h> 
void IonizationTemp1xSer_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[2]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[2]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.7071067811865476*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.4714045207910317*E*elemCharge)/m_; 
  }
 
  else { 
     vtSqIz[0] = 1.414213562373093e-10; 
  }
 
} 
 
void IonizationTemp1xSer_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[3]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[3]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.7071067811865476*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.4714045207910317*E*elemCharge)/m_; 
  }
 
  else { 
     vtSqIz[0] = 1.414213562373093e-10; 
  }
 
} 
 
void IonizationTemp1xSer_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[4]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[4]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.7071067811865476*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.4714045207910317*E*elemCharge)/m_; 
  }
 
  else { 
     vtSqIz[0] = 1.414213562373093e-10; 
  }
 
} 
 
