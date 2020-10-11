#include "MaxwellianCellAvModDecl.h" 
#include <math.h> 
void MaxwellianCellAvSer3x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[6]:      cell-center coordinates. 
  // m0[8]:     particle density. 
  // u[24]:      fluid velocity. 
  // vtSq[8]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.3535533905932738*m0[0]; 
  double vtSqAv = 0.3535533905932738*vtSq[0]; 
  double vSqAv = 0.125*pow(u[16],2)-0.7071067811865475*w[5]*u[16]+0.125*pow(u[8],2)-0.7071067811865475*w[4]*u[8]+pow(w[5],2)+pow(w[4],2)+pow(w[3],2)-0.7071067811865475*u[0]*w[3]+0.125*pow(u[0],2); 
 
  fMax[0] = (0.1795871221251665*m0Av)/(pow(vtSqAv,3/2)*exp(vSqAv/(2*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void MaxwellianCellAvSer3x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[6]:      cell-center coordinates. 
  // m0[20]:     particle density. 
  // u[60]:      fluid velocity. 
  // vtSq[20]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.3535533905932738*m0[0]; 
  double vtSqAv = 0.3535533905932738*vtSq[0]; 
  double vSqAv = 0.125*pow(u[40],2)-0.7071067811865475*w[5]*u[40]+0.125*pow(u[20],2)-0.7071067811865475*w[4]*u[20]+pow(w[5],2)+pow(w[4],2)+pow(w[3],2)-0.7071067811865475*u[0]*w[3]+0.125*pow(u[0],2); 
 
  fMax[0] = (0.1795871221251665*m0Av)/(pow(vtSqAv,3/2)*exp(vSqAv/(2*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
