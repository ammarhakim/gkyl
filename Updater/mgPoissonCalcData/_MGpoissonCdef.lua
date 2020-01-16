-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for C kernels used by the multigrid Poisson solver.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

ffi.cdef [[

void MGpoissonProlong1xSer_P1(const double *fldC, double *fldFL, double *fldFR);
void MGpoissonRestrict1xSer_P1(const double *fldFL, const double *fldFR, double *fldC);

void MGpoissonJacobi1xSer_P1(double **dx, const double *rho, double **phi);
void MGpoissonDampedJacobi1xSer_P1(const double omega, double **dx, const double *rho, double **phi);

void MGpoissonProlong2xSer_P1(const double *fldC, double *fldFL, double *fldFR);
void MGpoissonRestrict2xSer_P1(const double *fldFL, const double *fldFR, double *fldC);

void MGpoissonJacobi2xSer_P1(double **dx, const double *rho, double **phi);
void MGpoissonDampedJacobi2xSer_P1(const double omega, double **dx, const double *rho, double **phi);

]]
