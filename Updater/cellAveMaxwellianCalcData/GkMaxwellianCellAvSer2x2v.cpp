#include "MaxwellianCellAvModDecl.h" 
#include <math.h> 
void GkMaxwellianCellAvSer2x2v_P1(const double m_, const double *w, const double *m0, const double *uPar, const double *vtSq, double *bmag, double *fMax) 
{ 
  // w[4]:      cell-center coordinates. 
  // m0[4]:     particle density. 
  // uPar[8]:   fluid velocity. 
  // vtSq[4]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.5*m0[0]; 
  double uParAv = 0.5*uPar[0]; 
  double vtSqAv = 0.5*vtSq[0]; 
  double bmagAv = 0.5*bmag[0]; 
  double vSqAv = pow(uParAv,2)-2.0*w[2]*uParAv+pow(w[2],2); 
 
  fMax[0] = (0.3183098861837907*bmagAv*m0Av*exp((-vSqAv/(2*vtSqAv))-(w[3]*bmagAv)/(m_*vtSqAv)))/abs(vtSqAv); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void GkMaxwellianCellAvSer2x2v_P2(const double m_, const double *w, const double *m0, const double *uPar, const double *vtSq, double *bmag, double *fMax) 
{ 
  // w[4]:      cell-center coordinates. 
  // m0[8]:     particle density. 
  // uPar[16]:   fluid velocity. 
  // vtSq[8]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.5*m0[0]; 
  double uParAv = 0.5*uPar[0]; 
  double vtSqAv = 0.5*vtSq[0]; 
  double bmagAv = 0.5*bmag[0]; 
  double vSqAv = pow(uParAv,2)-2.0*w[2]*uParAv+pow(w[2],2); 
 
  fMax[0] = (0.3183098861837907*bmagAv*m0Av*exp((-vSqAv/(2*vtSqAv))-(w[3]*bmagAv)/(m_*vtSqAv)))/abs(vtSqAv); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
