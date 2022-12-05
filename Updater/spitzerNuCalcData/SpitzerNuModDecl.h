#ifndef SPITZERNU_MOD_DELC_H 
#define SPITZERNU_MOD_DELC_H 
#include <cmath> 

#include <algorithm> 

#include <../../Lib/gkyl_ipow.h> 

extern "C" { 
double gkyl_ipow(double base, int exp);

void SpitzerNuCellAvScale1xSer_P1(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvScale2xSer_P1(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvScale3xSer_P1(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 


void SpitzerNuCellAvScale1xSer_P2(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvScale2xSer_P2(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvScale3xSer_P2(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 


void SpitzerNuCellAvScale1xSer_P3(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvScale2xSer_P3(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvScale3xSer_P3(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 


void SpitzerNuCellAvBuild1xSer_P1(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvBuild2xSer_P1(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvBuild3xSer_P1(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 


void SpitzerNuCellAvBuild1xSer_P2(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvBuild2xSer_P2(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvBuild3xSer_P2(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 


void SpitzerNuCellAvBuild1xSer_P3(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvBuild2xSer_P3(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 

void SpitzerNuCellAvBuild3xSer_P3(double elemCharge, double eps0, double hBar, double nuFrac, double qA, double massA, const double *m0A, const double *vtSqA, double vtSqMinA, double qB, double massB, const double *m0B, const double *vtSqB, double vtSqMinB, double normNu, const double *Bmag, double *nu); 



 
} 
#endif 
