#ifndef MG_POISSON_MOD_DECL_H 
#define MG_POISSON_MOD_DECL_H 
 
extern "C" { 
 
void MGpoissonProlong1xSer_P1(const double *fldC, double *fldFL, double *fldFR);
void MGpoissonRestrict1xSer_P1(const double *fldFL, const double *fldFR, double *fldC);

void MGpoissonJacobi1xSer_P1(double **dx, const double *rho, double **phi);
void MGpoissonDampedJacobi1xSer_P1(const double omega, double **dx, const double *rho, double **phi);

void MGpoissonProlong2xSer_P1(const double *fldC, double *fldFL, double *fldFR);
void MGpoissonRestrict2xSer_P1(const double *fldFL, const double *fldFR, double *fldC);

void MGpoissonJacobi2xSer_P1(double **dx, const double *rho, double **phi);
void MGpoissonDampedJacobi2xSer_P1(const double omega, double **dx, const double *rho, double **phi);


} 
#endif 
