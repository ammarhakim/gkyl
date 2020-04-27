#ifndef MG_POISSON_MOD_DECL_H 
#define MG_POISSON_MOD_DECL_H 
 
#include <cmath>
 
extern "C" { 
 
void MGpoissonDGProlong1xSer_P1(const double *fldC, double **fldF);
void MGpoissonDGRestrict1xSer_P1(double **fldF, double *fldC);

void MGpoissonDGDampedGaussSeidel1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel1xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel1xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi1xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi1xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);

void MGpoissonDGResidue1xSer_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue1xSer_LxRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue1xSer_UxRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);

void MGpoissonFEM_DGtoFEM_1xSer_P1(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_1xSer_Lx_P1(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_1xSer_Ux_P1(const double *dgFld, double **femOut);

void MGpoissonFEMDampedGaussSeidel1xSer_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel1xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel1xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);

void MGpoissonFEMProlong1xSer_P1(const double *fldC, double **fldF);
void MGpoissonFEMProlong1xSer_LxRobin_P1(const double *fldC, double **fldF);
void MGpoissonFEMProlong1xSer_UxRobin_P1(const double *fldC, double **fldF);
void MGpoissonFEMRestrict1xSer_P1(double **fldF, double *fldC);
void MGpoissonFEMRestrict1xSer_LxRobin_P1(double **fldF, double *fldC);
void MGpoissonFEMRestrict1xSer_UxRobin_P1(double **fldF, double *fldC);

void MGpoissonDGProlong2xSer_P1(const double *fldC, double **fldF);
void MGpoissonDGRestrict2xSer_P1(double **fldF, double *fldC);

void MGpoissonDGDampedGaussSeidel2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_LyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_UyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_LxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_LxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_UxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_UxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_LyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_UyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_LxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_LxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_UxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_UxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);

void MGpoissonDGResidue2xSer_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_LxRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_UxRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_LyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_UyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_LxRobinLyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_LxRobinUyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_UxRobinLyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_UxRobinUyRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);

void MGpoissonFEM_DGtoFEM_2xSer_P1(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_Lx_P1(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_Ux_P1(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_Ly_P1(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_Uy_P1(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_LxLy_P1(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_LxUy_P1(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_UxLy_P1(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_UxUy_P1(const double *dgFld, double **femOut);

void MGpoissonFEMDampedGaussSeidel2xSer_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_LyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_UyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_LxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_LxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_UxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_UxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);

void MGpoissonFEMProlong2xSer_P1(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_LxRobin_P1(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_UxRobin_P1(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_LyRobin_P1(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_UyRobin_P1(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_LxRobinLyRobin_P1(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_LxRobinUyRobin_P1(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_UxRobinLyRobin_P1(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_UxRobinUyRobin_P1(const double *fldC, double **fldF);
void MGpoissonFEMRestrict2xSer_P1(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_LxRobin_P1(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_UxRobin_P1(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_LyRobin_P1(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_UyRobin_P1(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_LxRobinLyRobin_P1(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_LxRobinUyRobin_P1(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_UxRobinLyRobin_P1(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_UxRobinUyRobin_P1(double **fldF, double *fldC);


void MGpoissonDGProlong1xSer_P2(const double *fldC, double **fldF);
void MGpoissonDGRestrict1xSer_P2(double **fldF, double *fldC);

void MGpoissonDGDampedGaussSeidel1xSer_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi1xSer_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel1xSer_LxRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel1xSer_UxRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi1xSer_LxRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi1xSer_UxRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);

void MGpoissonDGResidue1xSer_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue1xSer_LxRobin_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue1xSer_UxRobin_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);

void MGpoissonFEM_DGtoFEM_1xSer_P2(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_1xSer_Lx_P2(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_1xSer_Ux_P2(const double *dgFld, double **femOut);

void MGpoissonFEMDampedGaussSeidel1xSer_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel1xSer_LxRobin_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel1xSer_UxRobin_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);

void MGpoissonFEMProlong1xSer_P2(const double *fldC, double **fldF);
void MGpoissonFEMProlong1xSer_LxRobin_P2(const double *fldC, double **fldF);
void MGpoissonFEMProlong1xSer_UxRobin_P2(const double *fldC, double **fldF);
void MGpoissonFEMRestrict1xSer_P2(double **fldF, double *fldC);
void MGpoissonFEMRestrict1xSer_LxRobin_P2(double **fldF, double *fldC);
void MGpoissonFEMRestrict1xSer_UxRobin_P2(double **fldF, double *fldC);

void MGpoissonDGProlong2xSer_P2(const double *fldC, double **fldF);
void MGpoissonDGRestrict2xSer_P2(double **fldF, double *fldC);

void MGpoissonDGDampedGaussSeidel2xSer_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_LxRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_UxRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_LyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_UyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_LxRobinLyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_LxRobinUyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_UxRobinLyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedGaussSeidel2xSer_UxRobinUyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_LxRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_UxRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_LyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_UyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_LxRobinLyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_LxRobinUyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_UxRobinLyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);
void MGpoissonDGDampedJacobi2xSer_UxRobinUyRobin_P2(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi);

void MGpoissonDGResidue2xSer_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_LxRobin_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_UxRobin_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_LyRobin_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_UyRobin_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_LxRobinLyRobin_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_LxRobinUyRobin_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_UxRobinLyRobin_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonDGResidue2xSer_UxRobinUyRobin_P2(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);

void MGpoissonFEM_DGtoFEM_2xSer_P2(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_Lx_P2(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_Ux_P2(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_Ly_P2(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_Uy_P2(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_LxLy_P2(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_LxUy_P2(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_UxLy_P2(const double *dgFld, double **femOut);
void MGpoissonFEM_DGtoFEM_2xSer_UxUy_P2(const double *dgFld, double **femOut);

void MGpoissonFEMDampedGaussSeidel2xSer_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_LxRobin_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_UxRobin_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_LyRobin_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_UyRobin_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_LxRobinLyRobin_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_LxRobinUyRobin_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_UxRobinLyRobin_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);
void MGpoissonFEMDampedGaussSeidel2xSer_UxRobinUyRobin_P2(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi);

void MGpoissonFEMProlong2xSer_P2(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_LxRobin_P2(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_UxRobin_P2(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_LyRobin_P2(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_UyRobin_P2(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_LxRobinLyRobin_P2(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_LxRobinUyRobin_P2(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_UxRobinLyRobin_P2(const double *fldC, double **fldF);
void MGpoissonFEMProlong2xSer_UxRobinUyRobin_P2(const double *fldC, double **fldF);
void MGpoissonFEMRestrict2xSer_P2(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_LxRobin_P2(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_UxRobin_P2(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_LyRobin_P2(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_UyRobin_P2(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_LxRobinLyRobin_P2(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_LxRobinUyRobin_P2(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_UxRobinLyRobin_P2(double **fldF, double *fldC);
void MGpoissonFEMRestrict2xSer_UxRobinUyRobin_P2(double **fldF, double *fldC);


} 
#endif 
