#ifndef MG_POISSON_MOD_DECL_H 
#define MG_POISSON_MOD_DECL_H 
 
extern "C" { 
 
void MGpoissonProlong1xSer_P1(const double *fldC, double **fldF);
void MGpoissonRestrict1xSer_P1(double **fldF, double *fldC);

void MGpoissonGaussSeidel1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel1xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel1xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel1xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel1xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel1xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel1xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel1xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel1xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);

void MGpoissonResidue1xSer_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue1xSer_LxDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue1xSer_LxNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue1xSer_UxDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue1xSer_UxNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);

void MGpoissonProlong2xSer_P1(const double *fldC, double **fldF);
void MGpoissonRestrict2xSer_P1(double **fldF, double *fldC);

void MGpoissonGaussSeidel2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_LxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonGaussSeidel2xSer_UxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_LxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);
void MGpoissonDampedGaussSeidel2xSer_UxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi);

void MGpoissonResidue2xSer_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxDirichletLyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxDirichletLyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxNeumannLyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxNeumannLyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxDirichletUyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxDirichletUyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxNeumannUyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_LxNeumannUyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxDirichletLyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxDirichletLyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxNeumannLyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxNeumannLyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxDirichletUyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxDirichletUyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxNeumannUyDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);
void MGpoissonResidue2xSer_UxNeumannUyNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut);


} 
#endif 
