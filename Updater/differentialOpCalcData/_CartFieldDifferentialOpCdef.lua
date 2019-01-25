-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used by differential operators.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

void CartFieldDifferentialOpDxxRecoveryVol1xSer_P1(const double *w, const double *dx, const int *idx, const double *f, double *out); 
void CartFieldDifferentialOpDxxRecoverySurf1xSer_X_P1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *fl, const double *fr, double *outl, double *outr); 
void CartFieldDifferentialOpDxxRecoveryBoundarySurf1xSer_X_P1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *fl, const double *fr, double *outl, double *outr); 



void CartFieldDifferentialOpDxxRecoveryVol1xMax_P1(const double *w, const double *dx, const int *idx, const double *f, double *out); 
void CartFieldDifferentialOpDxxRecoverySurf1xMax_X_P1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *fl, const double *fr, double *outl, double *outr); 
void CartFieldDifferentialOpDxxRecoveryBoundarySurf1xMax_X_P1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *fl, const double *fr, double *outl, double *outr); 




]] 
