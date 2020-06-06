-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in calculating charge exchange vRel and reaction rate.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

void VmProdCXcellAvSer1x1v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer1x1v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvSer1x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer1x2v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvSer1x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer1x3v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer1x1v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x1v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer1x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvSer2x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer2x2v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvSer2x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer2x3v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer2x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer2x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvSer3x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer3x3v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer3x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer3x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvSer1x1v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer1x1v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvSer1x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer1x2v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvSer1x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer1x3v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer1x1v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x1v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer1x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvSer2x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer2x2v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvSer2x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer2x3v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer2x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer2x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvSer3x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer3x3v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer3x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer3x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvSer1x1v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer1x1v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvSer1x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer1x2v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvSer1x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer1x3v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer1x1v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x1v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer1x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer1x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvSer2x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer2x2v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvSer2x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer2x3v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer2x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer2x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvSer3x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvSer3x3v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvSer3x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvSer3x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvMax1x1v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax1x1v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvMax1x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax1x2v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvMax1x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax1x3v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax1x1v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x1v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax1x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvMax2x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax2x2v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvMax2x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax2x3v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax2x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax2x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvMax3x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax3x3v_P1(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax3x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax3x2v_P1(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvMax1x1v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax1x1v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvMax1x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax1x2v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvMax1x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax1x3v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax1x1v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x1v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax1x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvMax2x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax2x2v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvMax2x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax2x3v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax2x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax2x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvMax3x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax3x3v_P2(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax3x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax3x2v_P2(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvMax1x1v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax1x1v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvMax1x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax1x2v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvMax1x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax1x3v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax1x1v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x1v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax1x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax1x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvMax2x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax2x2v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void VmProdCXcellAvMax2x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax2x3v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax2x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax2x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


void VmProdCXcellAvMax3x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX); 

void VmSigmaCXcellAvMax3x3v_P3(const double a, const double b, const double *uIon, const double *uNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 

void GkProdCXcellAvMax3x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkSigmaCXcellAvMax3x2v_P3(const double a, const double b, const double *uParIon, const double *uParNeut, const double *vtSqIon, const double *vtSqNeut, double *sigmaCX); 


]]
