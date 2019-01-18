-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used by differential operators.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

void CartFieldDifferentialOpDxxVol1xSer_P1(const double *w, const double *dx, const int *idx, const double *f, double *out); 
void CartFieldDifferentialOpDxxSurf1xSer_X_P1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *fl, const double *fr, double *outl, double *outr); 
void CartFieldDifferentialOpDxxBoundarySurf1xSer_X_P1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *fl, const double *fr, double *outl, double *outr); 



void CartFieldDifferentialOpDxxVol1xMax_P1(const double *w, const double *dx, const int *idx, const double *f, double *out); 
void CartFieldDifferentialOpDxxSurf1xMax_X_P1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *fl, const double *fr, double *outl, double *outr); 
void CartFieldDifferentialOpDxxBoundarySurf1xMax_X_P1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *fl, const double *fr, double *outl, double *outr); 




]] 
