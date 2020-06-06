#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void GkProdCXcellAvSer1x1v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[2]:      cell-center coordinates. 
  // m0[1]:     density. 
  // uPar[1]:    velocity. 
  // vtSq[1]:   squared thermal speed, sqrt(T/m). 
  // fOther:     distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double xSqAv = pow(w[1],2)/vtSqAv-(1.414213562373095*uPar[0]*w[1])/vtSqAv+(0.5*pow(uPar[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.7071067811865475*fOther[1]*m0[1]*vrelCX+0.7071067811865475*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.7071067811865475*fOther[0]*m0[1]*vrelCX+0.7071067811865475*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.7071067811865475*m0[1]*fOther[3]*vrelCX+0.7071067811865475*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.7071067811865475*m0[0]*fOther[3]*vrelCX+0.7071067811865475*m0[1]*fOther[2]*vrelCX; 
 
} 
void GkProdCXcellAvSer1x1v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[2]:      cell-center coordinates. 
  // m0[1]:     density. 
  // uPar[1]:    velocity. 
  // vtSq[1]:   squared thermal speed, sqrt(T/m). 
  // fOther:     distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double xSqAv = pow(w[1],2)/vtSqAv-(1.414213562373095*uPar[0]*w[1])/vtSqAv+(0.5*pow(uPar[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.7071067811865475*m0[2]*fOther[4]*vrelCX+0.7071067811865475*fOther[1]*m0[1]*vrelCX+0.7071067811865475*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.6324555320336759*m0[1]*fOther[4]*vrelCX+0.6324555320336759*fOther[1]*m0[2]*vrelCX+0.7071067811865475*fOther[0]*m0[1]*vrelCX+0.7071067811865475*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.7071067811865475*m0[2]*fOther[6]*vrelCX+0.7071067811865475*m0[1]*fOther[3]*vrelCX+0.7071067811865475*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.632455532033676*m0[1]*fOther[6]*vrelCX+0.6324555320336759*m0[2]*fOther[3]*vrelCX+0.7071067811865475*m0[0]*fOther[3]*vrelCX+0.7071067811865475*m0[1]*fOther[2]*vrelCX; 
  prodCX[4] = 0.4517539514526256*m0[2]*fOther[4]*vrelCX+0.7071067811865475*m0[0]*fOther[4]*vrelCX+0.7071067811865475*fOther[0]*m0[2]*vrelCX+0.6324555320336759*fOther[1]*m0[1]*vrelCX; 
  prodCX[5] = 0.7071067811865475*m0[1]*fOther[7]*vrelCX+0.7071067811865475*m0[0]*fOther[5]*vrelCX; 
  prodCX[6] = 0.4517539514526256*m0[2]*fOther[6]*vrelCX+0.7071067811865475*m0[0]*fOther[6]*vrelCX+0.632455532033676*m0[1]*fOther[3]*vrelCX+0.7071067811865475*fOther[2]*m0[2]*vrelCX; 
  prodCX[7] = 0.6324555320336759*m0[2]*fOther[7]*vrelCX+0.7071067811865475*m0[0]*fOther[7]*vrelCX+0.7071067811865475*m0[1]*fOther[5]*vrelCX; 
 
} 
void GkProdCXcellAvSer1x1v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[2]:      cell-center coordinates. 
  // m0[1]:     density. 
  // uPar[1]:    velocity. 
  // vtSq[1]:   squared thermal speed, sqrt(T/m). 
  // fOther:     distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double xSqAv = pow(w[1],2)/vtSqAv-(1.414213562373095*uPar[0]*w[1])/vtSqAv+(0.5*pow(uPar[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.7071067811865475*m0[3]*fOther[8]*vrelCX+0.7071067811865475*m0[2]*fOther[4]*vrelCX+0.7071067811865475*fOther[1]*m0[1]*vrelCX+0.7071067811865475*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.6210590034081186*m0[2]*fOther[8]*vrelCX+0.6210590034081186*m0[3]*fOther[4]*vrelCX+0.6324555320336759*m0[1]*fOther[4]*vrelCX+0.6324555320336759*fOther[1]*m0[2]*vrelCX+0.7071067811865475*fOther[0]*m0[1]*vrelCX+0.7071067811865475*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.7071067811865474*m0[3]*fOther[10]*vrelCX+0.7071067811865475*m0[2]*fOther[6]*vrelCX+0.7071067811865475*m0[1]*fOther[3]*vrelCX+0.7071067811865475*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.6210590034081187*m0[2]*fOther[10]*vrelCX+0.6210590034081187*m0[3]*fOther[6]*vrelCX+0.632455532033676*m0[1]*fOther[6]*vrelCX+0.6324555320336759*m0[2]*fOther[3]*vrelCX+0.7071067811865475*m0[0]*fOther[3]*vrelCX+0.7071067811865475*m0[1]*fOther[2]*vrelCX; 
  prodCX[4] = 0.421637021355784*m0[3]*fOther[8]*vrelCX+0.6210590034081186*m0[1]*fOther[8]*vrelCX+0.4517539514526256*m0[2]*fOther[4]*vrelCX+0.7071067811865475*m0[0]*fOther[4]*vrelCX+0.6210590034081186*fOther[1]*m0[3]*vrelCX+0.7071067811865475*fOther[0]*m0[2]*vrelCX+0.6324555320336759*fOther[1]*m0[1]*vrelCX; 
  prodCX[5] = 0.7071067811865475*m0[1]*fOther[7]*vrelCX+0.7071067811865475*m0[0]*fOther[5]*vrelCX; 
  prodCX[6] = 0.4216370213557839*m0[3]*fOther[10]*vrelCX+0.6210590034081187*m0[1]*fOther[10]*vrelCX+0.4517539514526256*m0[2]*fOther[6]*vrelCX+0.7071067811865475*m0[0]*fOther[6]*vrelCX+0.6210590034081187*fOther[3]*m0[3]*vrelCX+0.632455532033676*m0[1]*fOther[3]*vrelCX+0.7071067811865475*fOther[2]*m0[2]*vrelCX; 
  prodCX[7] = 0.6324555320336759*m0[2]*fOther[7]*vrelCX+0.7071067811865475*m0[0]*fOther[7]*vrelCX+0.7071067811865475*m0[1]*fOther[5]*vrelCX; 
  prodCX[8] = 0.421637021355784*m0[2]*fOther[8]*vrelCX+0.7071067811865475*m0[0]*fOther[8]*vrelCX+0.421637021355784*m0[3]*fOther[4]*vrelCX+0.6210590034081186*m0[1]*fOther[4]*vrelCX+0.7071067811865475*fOther[0]*m0[3]*vrelCX+0.6210590034081186*fOther[1]*m0[2]*vrelCX; 
  prodCX[9] = 0.7071067811865474*m0[1]*fOther[11]*vrelCX+0.7071067811865475*m0[0]*fOther[9]*vrelCX; 
  prodCX[10] = 0.421637021355784*m0[2]*fOther[10]*vrelCX+0.7071067811865475*m0[0]*fOther[10]*vrelCX+0.4216370213557839*m0[3]*fOther[6]*vrelCX+0.6210590034081187*m0[1]*fOther[6]*vrelCX+0.7071067811865474*fOther[2]*m0[3]*vrelCX+0.6210590034081187*m0[2]*fOther[3]*vrelCX; 
  prodCX[11] = 0.6324555320336759*m0[2]*fOther[11]*vrelCX+0.7071067811865475*m0[0]*fOther[11]*vrelCX+0.7071067811865474*m0[1]*fOther[9]*vrelCX; 
 
} 
