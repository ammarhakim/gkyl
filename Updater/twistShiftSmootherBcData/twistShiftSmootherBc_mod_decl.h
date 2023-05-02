#ifndef TWISTSHIFTSMOOTHERBC_MOD_DECL_H
#define TWISTSHIFTSMOOTHERBC_MOD_DECL_H

#include <cmath>

extern "C" { 

void twist_shift_smoother_bc_mod_lower_3x_Ser_p1(double *skinField, const double *ghostField); 
void twist_shift_smoother_bc_mod_upper_3x_Ser_p1(double *skinField, const double *ghostField); 



 
}

#endif 
