#ifndef MG_POISSON_MOD_DECL_H 
#define MG_POISSON_MOD_DECL_H 
 
extern "C" { 
 
void MGpoissonProlong1xSer_P1(const double *fldC, double **fldF);
void MGpoissonRestrict1xSer_P1(double **fldF, double *fldC);

void MGpoissonJacobi1xSer_P1(const double omega, double **dx, const double *rho, double **phi);
void MGpoissonDampedJacobi1xSer_P1(const double omega, double **dx, const double *rho, double **phi);

void MGpoissonResidue1xSer_P1(double **dx, const double *rho, double **phi, double *resOut);

void MGpoissonProlong2xSer_P1(const double *fldC, double **fldF);
void MGpoissonRestrict2xSer_P1(double **fldF, double *fldC);

void MGpoissonJacobi2xSer_P1(const double omega, double **dx, const double *rho, double **phi);
void MGpoissonDampedJacobi2xSer_P1(const double omega, double **dx, const double *rho, double **phi);

void MGpoissonResidue2xSer_P1(double **dx, const double *rho, double **phi, double *resOut);


} 
#endif 
