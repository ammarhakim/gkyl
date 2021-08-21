#include <MaxwellianCellAvModDecl.h> 
#include <math.h> 
void GkMaxwellianCellAvSer1x2v_P1(const double m_, const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax) 
{ 
  // w[3]:      cell-center coordinates. 
  // m0[2]:     particle density. 
  // uPar[4]:   fluid velocity. 
  // vtSq[2]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.7071067811865476*m0[0]; 
  double uParAv = 0.7071067811865476*uPar[0]; 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double bmagAv = 0.7071067811865476*bmag[0]; 
  double vSqAv = 0.5*pow(uParAv,2.0)-1.0*w[1]*uParAv+(w[2]*bmagAv)/m_+0.5*pow(w[1],2.0); 
 
  fMax[0] = (0.1795871221251666*bmagAv*m0Av)/(pow(vtSqAv,1.5)*exp(vSqAv/vtSqAv)); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void GkMaxwellianCellAvSer1x2v_P2(const double m_, const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax) 
{ 
  // w[3]:      cell-center coordinates. 
  // m0[3]:     particle density. 
  // uPar[6]:   fluid velocity. 
  // vtSq[3]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.7071067811865476*m0[0]; 
  double uParAv = 0.7071067811865476*uPar[0]; 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double bmagAv = 0.7071067811865476*bmag[0]; 
  double vSqAv = 0.5*pow(uParAv,2.0)-1.0*w[1]*uParAv+(w[2]*bmagAv)/m_+0.5*pow(w[1],2.0); 
 
  fMax[0] = (0.1795871221251666*bmagAv*m0Av)/(pow(vtSqAv,1.5)*exp(vSqAv/vtSqAv)); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
