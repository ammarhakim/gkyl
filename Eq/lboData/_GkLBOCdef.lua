-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the GK Lenard Bernstein solver.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

double GkLBOconstNuVol1x1vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf1x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf1x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 

double GkLBOconstNuVol1x1vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf1x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf1x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 


double GkLBOconstNuVol1x2vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf1x2vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf1x2vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuSurf1x2vSer_Mu_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf1x2vSer_Mu_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 

double GkLBOconstNuVol1x2vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf1x2vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf1x2vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuSurf1x2vSer_Mu_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf1x2vSer_Mu_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 


double GkLBOconstNuVol2x1vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf2x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf2x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 

double GkLBOconstNuVol2x1vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf2x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf2x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 


double GkLBOconstNuVol2x2vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf2x2vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf2x2vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuSurf2x2vSer_Mu_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf2x2vSer_Mu_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 

double GkLBOconstNuVol2x2vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf2x2vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf2x2vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuSurf2x2vSer_Mu_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf2x2vSer_Mu_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 


double GkLBOconstNuVol3x1vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf3x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf3x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 

double GkLBOconstNuVol3x1vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf3x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf3x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 


double GkLBOconstNuVol3x2vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf3x2vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf3x2vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuSurf3x2vSer_Mu_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf3x2vSer_Mu_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 

double GkLBOconstNuVol3x2vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out); 
double GkLBOconstNuSurf3x2vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf3x2vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuSurf3x2vSer_Mu_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 
double GkLBOconstNuBoundarySurf3x2vSer_Mu_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr); 


]]

