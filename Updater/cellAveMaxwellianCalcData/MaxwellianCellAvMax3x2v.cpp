#include <MaxwellianCellAvModDecl.h> 
#include <math.h> 
void MaxwellianCellAvMax3x2v_P1(const double m_, const double *w, const double *m0, const double *uPar, const double *vtSq, double *bmag, double *fMax) 
{ 
  // w[5]:      cell-center coordinates. 
  // m0[4]:     particle density. 
  // uPar[8]:   fluid velocity. 
  // vtSq[4]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.3535533905932738*m0[0]; 
  double uParAv = 0.3535533905932738*uPar[0]; 
  double vtSqAv = 0.3535533905932738*vtSq[0]; 
  double bmagAv = 0.3535533905932738*bmag[0]; 
  double vSqAv = pow(uParAv,2)-2.0*w[3]*uParAv+pow(w[3],2); 
 
  fMax[0] = (0.4501581580785531*bmagAv*m0Av*exp((-vSqAv/(2*vtSqAv))-(w[4]*bmagAv)/(m_*vtSqAv)))/abs(vtSqAv); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void MaxwellianCellAvMax3x2v_P2(const double m_, const double *w, const double *m0, const double *uPar, const double *vtSq, double *bmag, double *fMax) 
{ 
  // w[5]:      cell-center coordinates. 
  // m0[10]:     particle density. 
  // uPar[20]:   fluid velocity. 
  // vtSq[10]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.3535533905932738*m0[0]; 
  double uParAv = 0.3535533905932738*uPar[0]; 
  double vtSqAv = 0.3535533905932738*vtSq[0]; 
  double bmagAv = 0.3535533905932738*bmag[0]; 
  double vSqAv = pow(uParAv,2)-2.0*w[3]*uParAv+pow(w[3],2); 
 
  fMax[0] = (0.4501581580785531*bmagAv*m0Av*exp((-vSqAv/(2*vtSqAv))-(w[4]*bmagAv)/(m_*vtSqAv)))/abs(vtSqAv); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
