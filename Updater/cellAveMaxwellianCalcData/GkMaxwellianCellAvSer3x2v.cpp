#include <MaxwellianCellAvModDecl.h> 
#include <math.h> 
void GkMaxwellianCellAvSer3x2v_P1(const double m_, const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax) 
{ 
  // w[5]:      cell-center coordinates. 
  // m0[8]:     particle density. 
  // uPar[16]:   fluid velocity. 
  // vtSq[8]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.3535533905932738*m0[0]; 
  double uParAv = 0.3535533905932738*uPar[0]; 
  double vtSqAv = 0.3535533905932738*vtSq[0]; 
  double bmagAv = 0.3535533905932738*bmag[0]; 
  double vSqAv = 0.5*pow(uParAv,2.0)-1.0*w[3]*uParAv+(w[4]*bmagAv)/m_+0.5*pow(w[3],2.0); 
 
  fMax[0] = (0.3591742442503332*bmagAv*m0Av)/(pow(vtSqAv,1.5)*exp(vSqAv/vtSqAv)); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void GkMaxwellianCellAvSer3x2v_P2(const double m_, const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax) 
{ 
  // w[5]:      cell-center coordinates. 
  // m0[20]:     particle density. 
  // uPar[40]:   fluid velocity. 
  // vtSq[20]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.3535533905932738*m0[0]; 
  double uParAv = 0.3535533905932738*uPar[0]; 
  double vtSqAv = 0.3535533905932738*vtSq[0]; 
  double bmagAv = 0.3535533905932738*bmag[0]; 
  double vSqAv = 0.5*pow(uParAv,2.0)-1.0*w[3]*uParAv+(w[4]*bmagAv)/m_+0.5*pow(w[3],2.0); 
 
  fMax[0] = (0.3591742442503332*bmagAv*m0Av)/(pow(vtSqAv,1.5)*exp(vSqAv/vtSqAv)); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
