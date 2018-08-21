-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in calculating lagreange fix.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

void lagrangeFixSer1x1v1p(double *dm0, double *dm1, double *dm2, double *f, const double vc, const double L, const double Nv); 

]]

