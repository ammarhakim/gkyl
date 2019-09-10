-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the ConstDiffusion solver
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

double ConstDiffusionVol1xMaxP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity1xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol1xMaxP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol1xMaxP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 


double ConstDiffusionVol2xMaxP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity2xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf2xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity2xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol2xMaxP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf2xMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol2xMaxP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf2xMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 


double ConstDiffusionVol3xMaxP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xMax_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xMax_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol3xMaxP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xMax_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol3xMaxP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xMax_Z_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Z_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 



 
double ConstDiffusionVol1xSerP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity1xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol1xSerP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol1xSerP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol1xSerP4(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 


double ConstDiffusionVol2xSerP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity2xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf2xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity2xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol2xSerP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf2xSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol2xSerP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf2xSer_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol2xSerP4(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf2xSer_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 


double ConstDiffusionVol3xSerP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xSer_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol3xSerP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol3xSerP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Z_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Z_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol3xSerP4(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Z_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Z_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 



]] 
