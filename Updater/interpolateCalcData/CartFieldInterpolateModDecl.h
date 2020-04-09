#ifndef CART_FIELD_INTERPOLATE_MOD_DECL_H 
#define CART_FIELD_INTERPOLATE_MOD_DECL_H 
 
#include <cmath>
#include <algorithm>
 
extern "C" { 
 
void CartFieldInterp1xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);

void CartFieldInterp2xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);


void CartFieldInterp1xSer_P2(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);

void CartFieldInterp2xSer_P2(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF);


} 
#endif 
