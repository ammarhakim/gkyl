#ifndef RELATIVEVELOCITY_MOD_DELC_H 
#define RELATIVEVELOCITY_MOD_DELC_H 
#include <cmath> 

#include <algorithm> 

extern "C" { 
void GkProdCXcellAvSer1x1v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvSer1x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer1x1v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer1x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer1x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer2x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer2x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer2x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer3x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer3x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer1x1v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvSer1x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer1x1v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer1x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer1x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer2x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer2x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer2x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer3x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer3x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer1x1v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvSer1x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer1x1v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer1x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer1x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer2x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer2x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer2x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer3x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvSer3x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax1x1v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvMax1x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax1x1v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax1x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax1x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax2x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax2x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax2x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax3x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax3x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax1x1v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvMax1x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax1x1v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax1x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax1x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax2x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax2x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax2x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax3x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax3x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax1x1v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvMax1x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax1x1v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax1x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax1x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax2x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax2x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax2x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax3x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void VmProdCXcellAvMax3x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 



 
} 
#endif 
