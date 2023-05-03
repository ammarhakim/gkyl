#ifndef SKINGHOSTAVG_MOD_DECL_H
#define SKINGHOSTAVG_MOD_DECL_H

#include <cmath>

extern "C" { 

void skin_ghost_avg_lower_3x_Ser_p1(double *skinField, const double *ghostField); 
void skin_ghost_avg_upper_3x_Ser_p1(double *skinField, const double *ghostField); 



 
}

#endif 
