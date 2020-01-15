#ifndef MG_POISSON_MOD_DECL_H 
#define MG_POISSON_MOD_DECL_H 
 
extern "C" { 
 
void MGpoissonProlong1xSer_P1(const double *fldC, double *fldFL, double *fldFR);
void MGpoissonRestrict1xSer_P1(const double *fldFL, const double *fldFR, double *fldC);


} 
#endif 
