#ifndef CART_FIELD_INTERPOLATE_MOD_DECL_H 
#define CART_FIELD_INTERPOLATE_MOD_DECL_H 

#include <cmath>
#include <algorithm>
 
extern "C" { 
 
void CartFieldInterpProlong1xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);
void CartFieldInterpRestrict1xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldF, double *fldC);

void CartFieldInterpProlong2xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);
void CartFieldInterpRestrict2xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldF, double *fldC);


} 
#endif 
