-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the primitive moment calc
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[
void GkPrimMomentCalc1x1vSer_Upar_P1(const double *m0, const double *m1, const double *m2, double *out); 

]]
