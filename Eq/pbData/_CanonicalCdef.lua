-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the canonical poisson bracket eq
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[
double CanonicalVol1x1vMaxP1(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf1x1vMax_X_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf1x1vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 

double CanonicalVol1x1vMaxP2(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf1x1vMax_X_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf1x1vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 


double CanonicalVol2x2vMaxP1(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf2x2vMax_X_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vMax_Y_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vMax_VY_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 

double CanonicalVol2x2vMaxP2(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf2x2vMax_X_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vMax_Y_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vMax_VY_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 



 
double CanonicalVol1x1vSerP1(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf1x1vSer_X_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf1x1vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 

double CanonicalVol1x1vSerP2(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf1x1vSer_X_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf1x1vSer_VX_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 


double CanonicalVol2x2vSerP1(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf2x2vSer_X_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_Y_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_VY_P1(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 

double CanonicalVol2x2vSerP2(const double *w, const double *dxv, const double *H, const double *f, double *out); 
double CanonicalSurf2x2vSer_X_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_Y_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_VX_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 
double CanonicalSurf2x2vSer_VY_P2(const double *w, const double *dxv, const double amax, const double *H, const double *fl, const double *fr, double *outl, double *outr); 

]]
