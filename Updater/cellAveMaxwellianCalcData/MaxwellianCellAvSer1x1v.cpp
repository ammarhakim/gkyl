#include "MaxwellianCellAvModDecl.h" 
#include <math.h> 
void MaxwellianCellAvSer1x1v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[2]:      cell-center coordinates. 
  // m0[2]:     particle density. 
  // u[2]:      fluid velocity. 
  // vtSq[2]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.7071067811865476*m0[0]; 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double vSqAv = pow(w[1],2)-1.414213562373095*u[0]*w[1]+0.5*pow(u[0],2); 
 
  fMax[0] = (0.564189583547756*m0Av)/(sqrt(vtSqAv)*exp(vSqAv/(2*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void MaxwellianCellAvSer1x1v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[2]:      cell-center coordinates. 
  // m0[3]:     particle density. 
  // u[3]:      fluid velocity. 
  // vtSq[3]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.7071067811865476*m0[0]; 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double vSqAv = pow(w[1],2)-1.414213562373095*u[0]*w[1]+0.5*pow(u[0],2); 
 
  fMax[0] = (0.564189583547756*m0Av)/(sqrt(vtSqAv)*exp(vSqAv/(2*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
