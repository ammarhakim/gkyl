#include <IonizationModDecl.h> 
#include <math.h> 
void IonizationTemp3xMax_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[4]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[4]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-0.9428090415820636*E*elemCharge; 
  }
 
  else { 
     vtSqIz[0] = 2.828427124746187e-10; 
  }
 
} 
 
void IonizationTemp3xMax_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[10]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[10]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-0.9428090415820636*E*elemCharge; 
  }
 
  else { 
     vtSqIz[0] = 2.828427124746187e-10; 
  }
 
} 
 
void IonizationTemp3xMax_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[20]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[20]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-0.9428090415820636*E*elemCharge; 
  }
 
  else { 
     vtSqIz[0] = 2.828427124746187e-10; 
  }
 
} 
 
