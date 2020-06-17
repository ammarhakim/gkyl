#ifndef CHARGEEXCHANGE_MOD_DELC_H 
#define CHARGEEXCHANGE_MOD_DELC_H 
#include <cmath> 

#include <algorithm> 

extern "C" { 
void GkProdCXcellAvSer1x1v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x1v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer1x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvSer2x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer2x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvSer3x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer3x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvSer1x1v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x1v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer1x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvSer2x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer2x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvSer3x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer3x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvSer1x1v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x1v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer1x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvSer2x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer2x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvSer3x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer3x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvMax1x1v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x1v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax1x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvMax2x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax2x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvMax3x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax3x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvMax1x1v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x1v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax1x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvMax2x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax2x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvMax3x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax3x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvMax1x1v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x1v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax1x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvMax2x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax2x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void GkProdCXcellAvMax3x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax3x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 



 
} 
#endif 
