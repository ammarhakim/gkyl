#include "MaxwellianCellAvModDecl.h" 
#include <math.h> 
void MaxwellianCellAvSer2x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[5]:      cell-center coordinates. 
  // m0[4]:     particle density. 
  // u[12]:      fluid velocity. 
  // vtSq[4]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.5*m0[0]; 
  double vtSqAv = 0.5*m0[0]; 
  double vSqAv = 0.25*pow(u[8],2)-1.0*w[4]*u[8]+pow(w[4],2)+0.25*pow(u[4],2)-1.0*w[3]*u[4]+pow(w[3],2)+pow(w[2],2)-1.0*u[0]*w[2]+0.25*pow(u[0],2); 
 
  fMax[0] = (0.1269872718684819*m0Av)/(pow(vtSqAv,3/2)*exp(vSqAv/(2*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void MaxwellianCellAvSer2x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[5]:      cell-center coordinates. 
  // m0[8]:     particle density. 
  // u[24]:      fluid velocity. 
  // vtSq[8]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.5*m0[0]; 
  double vtSqAv = 0.5*m0[0]; 
  double vSqAv = 0.25*pow(u[16],2)-1.0*w[4]*u[16]+0.25*pow(u[8],2)-1.0*w[3]*u[8]+pow(w[4],2)+pow(w[3],2)+pow(w[2],2)-1.0*u[0]*w[2]+0.25*pow(u[0],2); 
 
  fMax[0] = (0.1269872718684819*m0Av)/(pow(vtSqAv,3/2)*exp(vSqAv/(2*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
