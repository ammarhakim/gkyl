-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for C kernels used by the CartFieldInterpolate updater.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

ffi.cdef [[

void CartFieldInterpProlong1xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);
void CartFieldInterpRestrict1xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldF, double *fldC);

void CartFieldInterpProlong2xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);
void CartFieldInterpRestrict2xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldF, double *fldC);


]]
