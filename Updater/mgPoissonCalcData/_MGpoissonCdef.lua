-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for C kernels used by the multigrid Poisson solver.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

ffi.cdef [[
 
void MGpoissonProlong1xSer_P1(const double *fldC, double **fldF);
void MGpoissonRestrict1xSer_P1(double **fldF, double *fldC);

void MGpoissonDampedGaussSeidel1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedGaussSeidel1xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedGaussSeidel1xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi1xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi1xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);

void MGpoissonResidue1xSer_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue1xSer_LxRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue1xSer_UxRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);

void MGpoissonProlong2xSer_P1(const double *fldC, double **fldF);
void MGpoissonRestrict2xSer_P1(double **fldF, double *fldC);

void MGpoissonDampedGaussSeidel2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi2xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi2xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi2xSer_LyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi2xSer_UyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi2xSer_LxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi2xSer_LxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi2xSer_UxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDampedJacobi2xSer_UxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);

void MGpoissonResidue2xSer_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxRobinLyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxRobinUyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxRobinLyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxRobinUyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);

]]
