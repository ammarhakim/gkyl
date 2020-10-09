#include "MaxwellianCellAvModDecl.h" 
#include <math.h> 
void MaxwellianCellAvMax1x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[4]:      cell-center coordinates. 
  // m0[2]:     particle density. 
  // u[6]:      fluid velocity. 
  // vtSq[2]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.7071067811865476*m0[0]; 
  double vtSqAv = 0.7071067811865476*m0[0]; 
  double vSqAv = 0.5*pow(u[4],2)-1.414213562373095*w[3]*u[4]+pow(w[3],2)+pow(w[2],2)-1.414213562373095*u[2]*w[2]+0.5*pow(u[2],2)+pow(w[1],2)-1.414213562373095*u[0]*w[1]+0.5*pow(u[0],2); 
 
  fMax[0] = (0.08979356106258324*m0Av)/(pow(vtSqAv,3/2)*exp(vSqAv/(2*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void MaxwellianCellAvMax1x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[4]:      cell-center coordinates. 
  // m0[3]:     particle density. 
  // u[9]:      fluid velocity. 
  // vtSq[3]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.7071067811865476*m0[0]; 
  double vtSqAv = 0.7071067811865476*m0[0]; 
  double vSqAv = 0.5*pow(u[6],2)-1.414213562373095*w[3]*u[6]+pow(w[3],2)+0.5*pow(u[3],2)-1.414213562373095*w[2]*u[3]+pow(w[2],2)+pow(w[1],2)-1.414213562373095*u[0]*w[1]+0.5*pow(u[0],2); 
 
  fMax[0] = (0.08979356106258324*m0Av)/(pow(vtSqAv,3/2)*exp(vSqAv/(2*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
