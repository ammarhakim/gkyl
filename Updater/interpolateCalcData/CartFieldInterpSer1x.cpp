#include <CartFieldInterpolateModDecl.h> 
 
void CartFieldInterpProlong1xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF) 
{ 
  // wC:   cell center of coarse-grid cell.
  // wF:   cell center of fine-grid cell.
  // dxC:  cell length of coarse-grid cell.
  // dxF:  cell length of fine-grid cell.
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double eLo[1];
  double eUp[1];
  eLo[0] = std::max(-1.0,1.0-(2.0*(wF[0]-1.0*wC[0]+0.5*dxF[0]+0.5*dxC[0]))/dxF[0]);
  eUp[0] = std::min( 1.0,(2.0*((-1.0*wF[0])+wC[0]+0.5*dxF[0]+0.5*dxC[0]))/dxF[0]-1.0);

  double eLo0R2 = std::pow(eLo[0],2);
  double eLo0R3 = std::pow(eLo[0],3);
  double eUp0R2 = std::pow(eUp[0],2);
  double eUp0R3 = std::pow(eUp[0],3);

  fldF[0] += (0.25*(((6.928203230275509*eUp[0]-6.928203230275509*eLo[0])*wF[0]+(6.928203230275509*eLo[0]-6.928203230275509*eUp[0])*wC[0]+1.732050807568877*dxF[0]*eUp0R2-1.732050807568877*dxF[0]*eLo0R2)*fldC[1]+(2.0*dxC[0]*eUp[0]-2.0*dxC[0]*eLo[0])*fldC[0]))/dxC[0]; 
  fldF[1] += (0.25*(((6.0*eUp0R2-6.0*eLo0R2)*wF[0]+(6.0*eLo0R2-6.0*eUp0R2)*wC[0]+2.0*dxF[0]*eUp0R3-2.0*dxF[0]*eLo0R3)*fldC[1]+(1.732050807568877*dxC[0]*eUp0R2-1.732050807568877*dxC[0]*eLo0R2)*fldC[0]))/dxC[0]; 

}

