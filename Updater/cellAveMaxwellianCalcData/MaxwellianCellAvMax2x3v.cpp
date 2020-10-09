#include "MaxwellianCellAvModDecl.h" 
#include <math.h> 
void MaxwellianCellAvMax2x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[5]:      cell-center coordinates. 
  // m0[3]:     particle density. 
  // u[9]:      fluid velocity. 
  // vtSq[3]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.5*m0[0]; 
  double vtSqAv = 0.5*m0[0]; 
  double vSqAv = 0.25*pow(u[6],2)-1.0*w[4]*u[6]+pow(w[4],2)+pow(w[3],2)-1.0*u[3]*w[3]+0.25*pow(u[3],2)+pow(w[2],2)-1.0*u[0]*w[2]+0.25*pow(u[0],2); 
 
  fMax[0] = (0.1269872718684819*m0Av)/(pow(vtSqAv,3/2)*exp(vSqAv/(2*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void MaxwellianCellAvMax2x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax) 
{ 
  // w[5]:      cell-center coordinates. 
  // m0[6]:     particle density. 
  // u[18]:      fluid velocity. 
  // vtSq[6]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.5*m0[0]; 
  double vtSqAv = 0.5*m0[0]; 
  double vSqAv = 0.25*pow(u[12],2)-1.0*w[4]*u[12]+0.25*pow(u[6],2)-1.0*w[3]*u[6]+pow(w[4],2)+pow(w[3],2)+pow(w[2],2)-1.0*u[0]*w[2]+0.25*pow(u[0],2); 
 
  fMax[0] = (0.1269872718684819*m0Av)/(pow(vtSqAv,3/2)*exp(vSqAv/(2*vtSqAv))); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
