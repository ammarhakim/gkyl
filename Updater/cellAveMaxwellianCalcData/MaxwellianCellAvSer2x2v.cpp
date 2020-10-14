#include "MaxwellianCellAvModDecl.h" 
#include <math.h> 
void MaxwellianCellAvSer2x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[4]:      cell-center coordinates. 
  // m0[4]:     particle density. 
  // u[8]:      fluid velocity. 
  // vtSq[4]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.5*m0[0]; 
  double vtSqAv = 0.5*vtSq[0]; 
  double vSqAv = 0.25*pow(u[4],2.0)-1.0*w[3]*u[4]+pow(w[3],2.0)+pow(w[2],2.0)-1.0*u[0]*w[2]+0.25*pow(u[0],2.0); 
 
  fMax[0] = (0.6366197723675814*m0Av)/(exp(vSqAv/(2.0*vtSqAv))*abs(vtSqAv)); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void MaxwellianCellAvSer2x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[4]:      cell-center coordinates. 
  // m0[8]:     particle density. 
  // u[16]:      fluid velocity. 
  // vtSq[8]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.5*m0[0]; 
  double vtSqAv = 0.5*vtSq[0]; 
  double vSqAv = 0.25*pow(u[8],2.0)-1.0*w[3]*u[8]+pow(w[3],2.0)+pow(w[2],2.0)-1.0*u[0]*w[2]+0.25*pow(u[0],2.0); 
 
  fMax[0] = (0.6366197723675814*m0Av)/(exp(vSqAv/(2.0*vtSqAv))*abs(vtSqAv)); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
