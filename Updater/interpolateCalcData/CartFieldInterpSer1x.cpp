#include <CartFieldInterpolateModDecl.h> 
 
void CartFieldInterp1xSer_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF) 
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

void CartFieldInterp1xSer_P2(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF) 
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

  double wC0R2 = std::pow(wC[0],2);
  double wF0R2 = std::pow(wF[0],2);
  double dxC0R2 = std::pow(dxC[0],2);
  double dxF0R2 = std::pow(dxF[0],2);
  double eLo0R2 = std::pow(eLo[0],2);
  double eLo0R3 = std::pow(eLo[0],3);
  double eLo0R4 = std::pow(eLo[0],4);
  double eLo0R5 = std::pow(eLo[0],5);
  double eUp0R2 = std::pow(eUp[0],2);
  double eUp0R3 = std::pow(eUp[0],3);
  double eUp0R4 = std::pow(eUp[0],4);
  double eUp0R5 = std::pow(eUp[0],5);

  fldF[0] += (0.25*(((26.83281572999748*eUp[0]-26.83281572999748*eLo[0])*wF0R2+((53.66563145999496*eLo[0]-53.66563145999496*eUp[0])*wC[0]+13.41640786499874*dxF[0]*eUp0R2-13.41640786499874*dxF[0]*eLo0R2)*wF[0]+(26.83281572999748*eUp[0]-26.83281572999748*eLo[0])*wC0R2+(13.41640786499874*dxF[0]*eLo0R2-13.41640786499874*dxF[0]*eUp0R2)*wC[0]+2.23606797749979*dxF0R2*eUp0R3-2.23606797749979*dxC0R2*eUp[0]-2.23606797749979*dxF0R2*eLo0R3+2.23606797749979*dxC0R2*eLo[0])*fldC[2]+((6.928203230275509*dxC[0]*eUp[0]-6.928203230275509*dxC[0]*eLo[0])*wF[0]+(6.928203230275509*dxC[0]*eLo[0]-6.928203230275509*dxC[0]*eUp[0])*wC[0]+1.732050807568877*dxC[0]*dxF[0]*eUp0R2-1.732050807568877*dxC[0]*dxF[0]*eLo0R2)*fldC[1]+(2.0*dxC0R2*eUp[0]-2.0*dxC0R2*eLo[0])*fldC[0]))/dxC0R2; 
  fldF[1] += (0.0625*(((92.951600308978*eUp0R2-92.951600308978*eLo0R2)*wF0R2+((185.903200617956*eLo0R2-185.903200617956*eUp0R2)*wC[0]+61.96773353931867*dxF[0]*eUp0R3-61.96773353931867*dxF[0]*eLo0R3)*wF[0]+(92.951600308978*eUp0R2-92.951600308978*eLo0R2)*wC0R2+(61.96773353931867*dxF[0]*eLo0R3-61.96773353931867*dxF[0]*eUp0R3)*wC[0]+11.61895003862225*dxF0R2*eUp0R4-7.745966692414834*dxC0R2*eUp0R2-11.61895003862225*dxF0R2*eLo0R4+7.745966692414834*dxC0R2*eLo0R2)*fldC[2]+((24.0*dxC[0]*eUp0R2-24.0*dxC[0]*eLo0R2)*wF[0]+(24.0*dxC[0]*eLo0R2-24.0*dxC[0]*eUp0R2)*wC[0]+8.0*dxC[0]*dxF[0]*eUp0R3-8.0*dxC[0]*dxF[0]*eLo0R3)*fldC[1]+(6.928203230275509*dxC0R2*eUp0R2-6.928203230275509*dxC0R2*eLo0R2)*fldC[0]))/dxC0R2; 
  fldF[2] += (0.0625*(((120.0*eUp0R3-120.0*eUp[0]-120.0*eLo0R3+120.0*eLo[0])*wF0R2+(((-240.0*eUp0R3)+240.0*eUp[0]+240.0*eLo0R3-240.0*eLo[0])*wC[0]+90.0*dxF[0]*eUp0R4-60.0*dxF[0]*eUp0R2-90.0*dxF[0]*eLo0R4+60.0*dxF[0]*eLo0R2)*wF[0]+(120.0*eUp0R3-120.0*eUp[0]-120.0*eLo0R3+120.0*eLo[0])*wC0R2+((-90.0*dxF[0]*eUp0R4)+60.0*dxF[0]*eUp0R2+90.0*dxF[0]*eLo0R4-60.0*dxF[0]*eLo0R2)*wC[0]+18.0*dxF0R2*eUp0R5+((-10.0*dxF0R2)-10.0*dxC0R2)*eUp0R3+10.0*dxC0R2*eUp[0]-18.0*dxF0R2*eLo0R5+(10.0*dxF0R2+10.0*dxC0R2)*eLo0R3-10.0*dxC0R2*eLo[0])*fldC[2]+((30.98386676965934*dxC[0]*eUp0R3-30.98386676965934*dxC[0]*eUp[0]-30.98386676965934*dxC[0]*eLo0R3+30.98386676965934*dxC[0]*eLo[0])*wF[0]+((-30.98386676965934*dxC[0]*eUp0R3)+30.98386676965934*dxC[0]*eUp[0]+30.98386676965934*dxC[0]*eLo0R3-30.98386676965934*dxC[0]*eLo[0])*wC[0]+11.61895003862225*dxC[0]*dxF[0]*eUp0R4-7.745966692414834*dxC[0]*dxF[0]*eUp0R2-11.61895003862225*dxC[0]*dxF[0]*eLo0R4+7.745966692414834*dxC[0]*dxF[0]*eLo0R2)*fldC[1]+(8.94427190999916*dxC0R2*eUp0R3-8.94427190999916*dxC0R2*eUp[0]-8.94427190999916*dxC0R2*eLo0R3+8.94427190999916*dxC0R2*eLo[0])*fldC[0]))/dxC0R2; 

}

void CartFieldInterp1xSer_X_P1(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF) 
{ 
  // wC:   cell center of coarse-grid cell.
  // wF:   cell center of fine-grid cell.
  // dxC:  cell length of coarse-grid cell.
  // dxF:  cell length of fine-grid cell.
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double eLo[1];
  double eUp[1];
  eLo[0] = std::max(-1.0,-1.0);
  eUp[0] = std::min( 1.0,1.0);

  double eLo0R2 = std::pow(eLo[0],2);
  double eLo0R3 = std::pow(eLo[0],3);
  double eUp0R2 = std::pow(eUp[0],2);
  double eUp0R3 = std::pow(eUp[0],3);

  fldF[0] += 0.25*((1.732050807568877*eUp0R2-1.732050807568877*eLo0R2)*fldC[1]+(2.0*eUp[0]-2.0*eLo[0])*fldC[0]); 
  fldF[1] += 0.25*((2.0*eUp0R3-2.0*eLo0R3)*fldC[1]+(1.732050807568877*eUp0R2-1.732050807568877*eLo0R2)*fldC[0]); 

}

void CartFieldInterp1xSer_X_P2(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF) 
{ 
  // wC:   cell center of coarse-grid cell.
  // wF:   cell center of fine-grid cell.
  // dxC:  cell length of coarse-grid cell.
  // dxF:  cell length of fine-grid cell.
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double eLo[1];
  double eUp[1];
  eLo[0] = std::max(-1.0,-1.0);
  eUp[0] = std::min( 1.0,1.0);

  double eLo0R2 = std::pow(eLo[0],2);
  double eLo0R3 = std::pow(eLo[0],3);
  double eLo0R4 = std::pow(eLo[0],4);
  double eLo0R5 = std::pow(eLo[0],5);
  double eUp0R2 = std::pow(eUp[0],2);
  double eUp0R3 = std::pow(eUp[0],3);
  double eUp0R4 = std::pow(eUp[0],4);
  double eUp0R5 = std::pow(eUp[0],5);

  fldF[0] += 0.25*((2.23606797749979*eUp0R3-2.23606797749979*eUp[0]-2.23606797749979*eLo0R3+2.23606797749979*eLo[0])*fldC[2]+(1.732050807568877*eUp0R2-1.732050807568877*eLo0R2)*fldC[1]+(2.0*eUp[0]-2.0*eLo[0])*fldC[0]); 
  fldF[1] += 0.0625*((11.61895003862225*eUp0R4-7.745966692414834*eUp0R2-11.61895003862225*eLo0R4+7.745966692414834*eLo0R2)*fldC[2]+(8.0*eUp0R3-8.0*eLo0R3)*fldC[1]+(6.928203230275509*eUp0R2-6.928203230275509*eLo0R2)*fldC[0]); 
  fldF[2] += 0.0625*((18.0*eUp0R5-20.0*eUp0R3+10.0*eUp[0]-18.0*eLo0R5+20.0*eLo0R3-10.0*eLo[0])*fldC[2]+(11.61895003862225*eUp0R4-7.745966692414834*eUp0R2-11.61895003862225*eLo0R4+7.745966692414834*eLo0R2)*fldC[1]+(8.94427190999916*eUp0R3-8.94427190999916*eUp[0]-8.94427190999916*eLo0R3+8.94427190999916*eLo[0])*fldC[0]); 

}

