#include "MaxwellianCellAvModDecl.h" 
#include <math.h> 
void MaxwellianCellAvSer1x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[3]:      cell-center coordinates. 
  // m0[2]:     particle density. 
  // u[4]:      fluid velocity. 
  // vtSq[2]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.7071067811865476*m0[0]; 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double vSqAv = pow(w[2],2)-1.414213562373095*u[2]*w[2]+0.5*pow(u[2],2)+pow(w[1],2)-1.414213562373095*u[0]*w[1]+0.5*pow(u[0],2); 
 
  fMax[0] = (0.2250790790392765*m0Av)/(exp(vSqAv/(2*vtSqAv))*abs(vtSqAv)); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void MaxwellianCellAvSer1x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[3]:      cell-center coordinates. 
  // m0[3]:     particle density. 
  // u[6]:      fluid velocity. 
  // vtSq[3]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.7071067811865476*m0[0]; 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double vSqAv = 0.5*pow(u[3],2)-1.414213562373095*w[2]*u[3]+pow(w[2],2)+pow(w[1],2)-1.414213562373095*u[0]*w[1]+0.5*pow(u[0],2); 
 
  fMax[0] = (0.2250790790392765*m0Av)/(exp(vSqAv/(2*vtSqAv))*abs(vtSqAv)); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
