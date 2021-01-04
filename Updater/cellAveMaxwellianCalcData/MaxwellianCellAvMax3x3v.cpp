#include "MaxwellianCellAvModDecl.h" 
#include <math.h> 
void MaxwellianCellAvMax3x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[6]:      cell-center coordinates. 
  // m0[4]:     particle density. 
  // u[12]:      fluid velocity. 
  // vtSq[4]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.3535533905932738*m0[0]; 
  double vtSqAv = 0.3535533905932738*vtSq[0]; 
  double vSqAv = 0.125*pow(u[8],2.0)-0.7071067811865475*w[5]*u[8]+pow(w[5],2.0)+pow(w[4],2.0)-0.7071067811865475*u[4]*w[4]+0.125*pow(u[4],2.0)+pow(w[3],2.0)-0.7071067811865475*u[0]*w[3]+0.125*pow(u[0],2.0); 
 
  fMax[0] = (0.5079490874739284*m0Av)/(pow(vtSqAv,1.5)*exp(vSqAv/(2.0*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void MaxwellianCellAvMax3x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[6]:      cell-center coordinates. 
  // m0[10]:     particle density. 
  // u[30]:      fluid velocity. 
  // vtSq[10]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.3535533905932738*m0[0]; 
  double vtSqAv = 0.3535533905932738*vtSq[0]; 
  double vSqAv = 0.125*pow(u[20],2.0)-0.7071067811865475*w[5]*u[20]+0.125*pow(u[10],2.0)-0.7071067811865475*w[4]*u[10]+pow(w[5],2.0)+pow(w[4],2.0)+pow(w[3],2.0)-0.7071067811865475*u[0]*w[3]+0.125*pow(u[0],2.0); 
 
  fMax[0] = (0.5079490874739284*m0Av)/(pow(vtSqAv,1.5)*exp(vSqAv/(2.0*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
