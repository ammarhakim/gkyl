-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the gyrokinetic eq
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[
double GyrokineticVol2x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *f, double *out); 
double GyrokineticSurf2x2vSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double GyrokineticSurf2x2vSer_Y_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double GyrokineticSurf2x2vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 


double GyrokineticVol3x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *f, double *out); 
double GyrokineticSurf3x2vSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double GyrokineticSurf3x2vSer_Y_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double GyrokineticSurf3x2vSer_Z_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double GyrokineticSurf3x2vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 

]]
