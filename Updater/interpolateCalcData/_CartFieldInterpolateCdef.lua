-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for C kernels used by the CartFieldInterpolate updater.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

ffi.cdef [[
 
void CartFieldInterp1xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);

void CartFieldInterp2xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);


void CartFieldInterp1xSer_P2(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);

void CartFieldInterp2xSer_P2(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);

]]
