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
     vtSqIz[1] = 0.5*vtSq[1]; 
     vtSqIz[2] = 0.5*vtSq[2]; 
  }
 
  else { 
     vtSqIz[0] = 1.0e-10*vtSq[0]; 
     vtSqIz[1] = 1.0e-10*vtSq[1]; 
     vtSqIz[2] = 1.0e-10*vtSq[2]; 
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
     vtSqIz[1] = 0.5*vtSq[1]; 
     vtSqIz[2] = 0.5*vtSq[2]; 
     vtSqIz[3] = 0.5*vtSq[3]; 
     vtSqIz[4] = 0.5*vtSq[4]; 
     vtSqIz[5] = 0.5*vtSq[5]; 
  }
 
  else { 
     vtSqIz[0] = 1.0e-10*vtSq[0]; 
     vtSqIz[1] = 1.0e-10*vtSq[1]; 
     vtSqIz[2] = 1.0e-10*vtSq[2]; 
     vtSqIz[3] = 1.0e-10*vtSq[3]; 
     vtSqIz[4] = 1.0e-10*vtSq[4]; 
     vtSqIz[5] = 1.0e-10*vtSq[5]; 
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
     vtSqIz[1] = 0.5*vtSq[1]; 
     vtSqIz[2] = 0.5*vtSq[2]; 
     vtSqIz[3] = 0.5*vtSq[3]; 
     vtSqIz[4] = 0.5*vtSq[4]; 
     vtSqIz[5] = 0.5*vtSq[5]; 
     vtSqIz[6] = 0.5*vtSq[6]; 
     vtSqIz[7] = 0.5*vtSq[7]; 
     vtSqIz[8] = 0.5*vtSq[8]; 
     vtSqIz[9] = 0.5*vtSq[9]; 
  }
 
  else { 
     vtSqIz[0] = 1.0e-10*vtSq[0]; 
     vtSqIz[1] = 1.0e-10*vtSq[1]; 
     vtSqIz[2] = 1.0e-10*vtSq[2]; 
     vtSqIz[3] = 1.0e-10*vtSq[3]; 
     vtSqIz[4] = 1.0e-10*vtSq[4]; 
     vtSqIz[5] = 1.0e-10*vtSq[5]; 
     vtSqIz[6] = 1.0e-10*vtSq[6]; 
     vtSqIz[7] = 1.0e-10*vtSq[7]; 
     vtSqIz[8] = 1.0e-10*vtSq[8]; 
     vtSqIz[9] = 1.0e-10*vtSq[9]; 
  }
 
} 
 
