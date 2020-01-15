#include <MGpoissonModDecl.h> 
 
void MGpoissonProlong1xSer_P1(const double *fldC, double *fldFL, double *fldFR) 
{ 
  // fldC:  coarse-grid field.
  // fldFL: fine-grid field (left).
  // fldFR: fine-grid field (right).
  fldFL[0] = fldC[0]-0.8660254037844386*fldC[1]; 
  fldFL[1] = 0.5*fldC[1]; 

  fldFR[0] = 0.8660254037844386*fldC[1]+fldC[0]; 
  fldFR[1] = 0.5*fldC[1]; 

}

void MGpoissonRestrict1xSer_P1(const double *fldFL, const double *fldFR, double *fldC) 
{ 
  // fldFL: fine-grid field (left).
  // fldFR: fine-grid field (right).
  // fldC:  coarse-grid field.
  fldC[0] = 0.5*fldFR[0]+0.5*fldFL[0]; 
  fldC[1] = 0.25*fldFR[1]+0.25*fldFL[1]+0.4330127018922193*fldFR[0]-0.4330127018922193*fldFL[0]; 

}

