#include <IonizationModDecl.h> 
#include <math.h> 
void IonizationTemp3xSer_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[8]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[8]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.9428090415820636*E*elemCharge)/m_; 
     vtSqIz[1] = 0.5*vtSq[1]; 
     vtSqIz[2] = 0.5*vtSq[2]; 
     vtSqIz[3] = 0.5*vtSq[3]; 
     vtSqIz[4] = 0.5*vtSq[4]; 
     vtSqIz[5] = 0.5*vtSq[5]; 
     vtSqIz[6] = 0.5*vtSq[6]; 
     vtSqIz[7] = 0.5*vtSq[7]; 
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
  }
 
} 
 
void IonizationTemp3xSer_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[20]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[20]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.9428090415820636*E*elemCharge)/m_; 
     vtSqIz[1] = 0.5*vtSq[1]; 
     vtSqIz[2] = 0.5*vtSq[2]; 
     vtSqIz[3] = 0.5*vtSq[3]; 
     vtSqIz[4] = 0.5*vtSq[4]; 
     vtSqIz[5] = 0.5*vtSq[5]; 
     vtSqIz[6] = 0.5*vtSq[6]; 
     vtSqIz[7] = 0.5*vtSq[7]; 
     vtSqIz[8] = 0.5*vtSq[8]; 
     vtSqIz[9] = 0.5*vtSq[9]; 
     vtSqIz[10] = 0.5*vtSq[10]; 
     vtSqIz[11] = 0.5*vtSq[11]; 
     vtSqIz[12] = 0.5*vtSq[12]; 
     vtSqIz[13] = 0.5*vtSq[13]; 
     vtSqIz[14] = 0.5*vtSq[14]; 
     vtSqIz[15] = 0.5*vtSq[15]; 
     vtSqIz[16] = 0.5*vtSq[16]; 
     vtSqIz[17] = 0.5*vtSq[17]; 
     vtSqIz[18] = 0.5*vtSq[18]; 
     vtSqIz[19] = 0.5*vtSq[19]; 
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
     vtSqIz[10] = 1.0e-10*vtSq[10]; 
     vtSqIz[11] = 1.0e-10*vtSq[11]; 
     vtSqIz[12] = 1.0e-10*vtSq[12]; 
     vtSqIz[13] = 1.0e-10*vtSq[13]; 
     vtSqIz[14] = 1.0e-10*vtSq[14]; 
     vtSqIz[15] = 1.0e-10*vtSq[15]; 
     vtSqIz[16] = 1.0e-10*vtSq[16]; 
     vtSqIz[17] = 1.0e-10*vtSq[17]; 
     vtSqIz[18] = 1.0e-10*vtSq[18]; 
     vtSqIz[19] = 1.0e-10*vtSq[19]; 
  }
 
} 
 
void IonizationTemp3xSer_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[32]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[32]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.9428090415820636*E*elemCharge)/m_; 
     vtSqIz[1] = 0.5*vtSq[1]; 
     vtSqIz[2] = 0.5*vtSq[2]; 
     vtSqIz[3] = 0.5*vtSq[3]; 
     vtSqIz[4] = 0.5*vtSq[4]; 
     vtSqIz[5] = 0.5*vtSq[5]; 
     vtSqIz[6] = 0.5*vtSq[6]; 
     vtSqIz[7] = 0.5*vtSq[7]; 
     vtSqIz[8] = 0.5*vtSq[8]; 
     vtSqIz[9] = 0.5*vtSq[9]; 
     vtSqIz[10] = 0.5*vtSq[10]; 
     vtSqIz[11] = 0.5*vtSq[11]; 
     vtSqIz[12] = 0.5*vtSq[12]; 
     vtSqIz[13] = 0.5*vtSq[13]; 
     vtSqIz[14] = 0.5*vtSq[14]; 
     vtSqIz[15] = 0.5*vtSq[15]; 
     vtSqIz[16] = 0.5*vtSq[16]; 
     vtSqIz[17] = 0.5*vtSq[17]; 
     vtSqIz[18] = 0.5*vtSq[18]; 
     vtSqIz[19] = 0.5*vtSq[19]; 
     vtSqIz[20] = 0.5*vtSq[20]; 
     vtSqIz[21] = 0.5*vtSq[21]; 
     vtSqIz[22] = 0.5*vtSq[22]; 
     vtSqIz[23] = 0.5*vtSq[23]; 
     vtSqIz[24] = 0.5*vtSq[24]; 
     vtSqIz[25] = 0.5*vtSq[25]; 
     vtSqIz[26] = 0.5*vtSq[26]; 
     vtSqIz[27] = 0.5*vtSq[27]; 
     vtSqIz[28] = 0.5*vtSq[28]; 
     vtSqIz[29] = 0.5*vtSq[29]; 
     vtSqIz[30] = 0.5*vtSq[30]; 
     vtSqIz[31] = 0.5*vtSq[31]; 
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
     vtSqIz[10] = 1.0e-10*vtSq[10]; 
     vtSqIz[11] = 1.0e-10*vtSq[11]; 
     vtSqIz[12] = 1.0e-10*vtSq[12]; 
     vtSqIz[13] = 1.0e-10*vtSq[13]; 
     vtSqIz[14] = 1.0e-10*vtSq[14]; 
     vtSqIz[15] = 1.0e-10*vtSq[15]; 
     vtSqIz[16] = 1.0e-10*vtSq[16]; 
     vtSqIz[17] = 1.0e-10*vtSq[17]; 
     vtSqIz[18] = 1.0e-10*vtSq[18]; 
     vtSqIz[19] = 1.0e-10*vtSq[19]; 
     vtSqIz[20] = 1.0e-10*vtSq[20]; 
     vtSqIz[21] = 1.0e-10*vtSq[21]; 
     vtSqIz[22] = 1.0e-10*vtSq[22]; 
     vtSqIz[23] = 1.0e-10*vtSq[23]; 
     vtSqIz[24] = 1.0e-10*vtSq[24]; 
     vtSqIz[25] = 1.0e-10*vtSq[25]; 
     vtSqIz[26] = 1.0e-10*vtSq[26]; 
     vtSqIz[27] = 1.0e-10*vtSq[27]; 
     vtSqIz[28] = 1.0e-10*vtSq[28]; 
     vtSqIz[29] = 1.0e-10*vtSq[29]; 
     vtSqIz[30] = 1.0e-10*vtSq[30]; 
     vtSqIz[31] = 1.0e-10*vtSq[31]; 
  }
 
} 
 
