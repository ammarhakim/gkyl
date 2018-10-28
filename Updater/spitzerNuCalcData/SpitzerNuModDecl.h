#ifndef SPITZERNU_MOD_DELC_H 
#define SPITZERNU_MOD_DELC_H 
#include <../../Lib/gkyl_ipow.h> 

extern "C" { 
double gkyl_ipow(double base, int exp);

void SpitzerNuCellAv1xMax_P1(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv1xMax_P2(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv1xMax_P3(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAv2xMax_P1(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv2xMax_P2(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv2xMax_P3(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAv3xMax_P1(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv3xMax_P2(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv3xMax_P3(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 



 
void SpitzerNuCellAv1xSer_P1(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv1xSer_P2(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv1xSer_P3(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAv2xSer_P1(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv2xSer_P2(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv2xSer_P3(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAv3xSer_P1(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv3xSer_P2(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAv3xSer_P3(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu); 



 
} 
#endif 
