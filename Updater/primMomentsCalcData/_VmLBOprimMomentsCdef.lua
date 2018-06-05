-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the binary operator calc.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

void VmLBOconstNuPrimMoments1x1vSer_P1(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, double *u, double *vtSq); 

]]
