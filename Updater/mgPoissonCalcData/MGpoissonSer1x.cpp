#include <MGpoissonModDecl.h> 
 
void MGpoissonProlong1xSer_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF1 = fldF[0];
  double *fldF2 = fldF[1];

  fldF1[0] = fldC[0]-0.8660254037844386*fldC[1]; 
  fldF1[1] = 0.5*fldC[1]; 

  fldF2[0] = 0.8660254037844386*fldC[1]+fldC[0]; 
  fldF2[1] = 0.5*fldC[1]; 

}

void MGpoissonRestrict1xSer_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in stencils pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF1 = fldF[0];
  double *fldF2 = fldF[1];

  fldC[0] = 0.5*fldF2[0]+0.5*fldF1[0]; 
  fldC[1] = 0.25*fldF2[1]+0.25*fldF1[1]+0.4330127018922193*fldF2[0]-0.4330127018922193*fldF1[0]; 

}

void MGpoissonGaussSeidel1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double *dxUx = dx[1]; 
  double *dxLx = dx[2]; 

  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxLx[1]; 
  double rdxUx[1]; 
  double rdxLxSq[1]; 
  double rdxUxSq[1]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxLx[0]   = 1.0/dxLx[0]; 
  rdxUx[0]   = 1.0/dxUx[0]; 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(1.0*(((17.32050807568877*rdxUx[0]-17.32050807568877*rdxLx[0])*rho[1]+((-46.0*rdxUx[0])-46.0*rdxLx[0])*rho[0])*volFac+(69.28203230275508*rdxUxSq[0]+129.9038105676658*rdxLx[0]*rdxUx[0])*phiUx[1]+((-129.9038105676658*rdxLx[0]*rdxUx[0])-69.28203230275508*rdxLxSq[0])*phiLx[1]-66.0*phiUx[0]*rdxUxSq[0]+((-141.0*phiUx[0])-141.0*phiLx[0])*rdxLx[0]*rdxUx[0]-66.0*phiLx[0]*rdxLxSq[0]))/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 
  phiC[1] = (((18.0*rdxUx[0]+18.0*rdxLx[0])*rho[1]+(45.03332099679081*rdxLx[0]-45.03332099679081*rdxUx[0])*rho[0])*volFac+(66.0*rdxUxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(66.0*rdxLxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiLx[1]-62.35382907247956*phiUx[0]*rdxUxSq[0]+(140.296115413079*phiUx[0]-140.296115413079*phiLx[0])*rdxLx[0]*rdxUx[0]+62.35382907247956*phiLx[0]*rdxLxSq[0])/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 

}

void MGpoissonGaussSeidel1xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (6.873217490352689e-4*((3.0*rho[1]+86.60254037844386*rho[0])*volFac-720.0*rdxCp2Sq[0]*phiUx[1]+(623.5382907247956*phiUx[0]+1175.755076535926*bcVals[0])*rdxCp2Sq[0]))/rdxCp2Sq[0]; 
  phiC[1] = (0.001030982623552903*((5.196152422706631*rho[1]+10.0*rho[0])*volFac-277.1281292110203*rdxCp2Sq[0]*phiUx[1]+(240.0*phiUx[0]-339.411254969543*bcVals[0])*rdxCp2Sq[0]))/rdxCp2Sq[0]; 

}

void MGpoissonGaussSeidel1xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(0.004811252243246883*((9.0*rho[1]-51.96152422706631*rho[0])*volFac+240.0*rdxCp2Sq[0]*phiUx[1]+(195.9591794226543*bcVals[0]-207.8460969082653*phiUx[0])*rdxCp2Sq[0]))/rdxCp2Sq[0]; 
  phiC[1] = (0.01443375672974065*((1.732050807568877*rho[1]-5.0*rho[0])*volFac+28.28427124746191*bcVals[0]*rdxCp2Sq[0]))/rdxCp2Sq[0]; 

}

void MGpoissonGaussSeidel1xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(6.873217490352689e-4*((3.0*rho[1]-86.60254037844386*rho[0])*volFac-720.0*rdxCp2Sq[0]*phiLx[1]-1175.755076535926*rdxCp2Sq[0]*bcVals[1]-623.5382907247956*phiLx[0]*rdxCp2Sq[0]))/rdxCp2Sq[0]; 
  phiC[1] = (0.001030982623552903*((5.196152422706631*rho[1]-10.0*rho[0])*volFac-277.1281292110203*rdxCp2Sq[0]*phiLx[1]+339.411254969543*rdxCp2Sq[0]*bcVals[1]-240.0*phiLx[0]*rdxCp2Sq[0]))/rdxCp2Sq[0]; 

}

void MGpoissonGaussSeidel1xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (0.004811252243246883*((9.0*rho[1]+51.96152422706631*rho[0])*volFac+240.0*rdxCp2Sq[0]*phiLx[1]+195.9591794226543*rdxCp2Sq[0]*bcVals[1]+207.8460969082653*phiLx[0]*rdxCp2Sq[0]))/rdxCp2Sq[0]; 
  phiC[1] = (0.01443375672974065*((1.732050807568877*rho[1]+5.0*rho[0])*volFac+28.28427124746191*rdxCp2Sq[0]*bcVals[1]))/rdxCp2Sq[0]; 

}

void MGpoissonDampedGaussSeidel1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double *dxUx = dx[1]; 
  double *dxLx = dx[2]; 

  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxLx[1]; 
  double rdxUx[1]; 
  double rdxLxSq[1]; 
  double rdxUxSq[1]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxLx[0]   = 1.0/dxLx[0]; 
  rdxUx[0]   = 1.0/dxUx[0]; 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(1.0*(((17.32050807568877*rdxUx[0]-17.32050807568877*rdxLx[0])*rho[1]+((-46.0*rdxUx[0])-46.0*rdxLx[0])*rho[0])*omega*volFac+((69.28203230275508*rdxUxSq[0]+129.9038105676658*rdxLx[0]*rdxUx[0])*phiUx[1]+((-129.9038105676658*rdxLx[0]*rdxUx[0])-69.28203230275508*rdxLxSq[0])*phiLx[1]+(6.0*phiC[0]-66.0*phiUx[0])*rdxUxSq[0]+((-141.0*phiUx[0])-141.0*phiLx[0]+402.0*phiC[0])*rdxLx[0]*rdxUx[0]+(6.0*phiC[0]-66.0*phiLx[0])*rdxLxSq[0])*omega-6.0*phiC[0]*rdxUxSq[0]-402.0*phiC[0]*rdxLx[0]*rdxUx[0]-6.0*phiC[0]*rdxLxSq[0]))/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 
  phiC[1] = (((18.0*rdxUx[0]+18.0*rdxLx[0])*rho[1]+(45.03332099679081*rdxLx[0]-45.03332099679081*rdxUx[0])*rho[0])*omega*volFac+((66.0*rdxUxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(66.0*rdxLxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiLx[1]+((-6.0*rdxUxSq[0])-402.0*rdxLx[0]*rdxUx[0]-6.0*rdxLxSq[0])*phiC[1]-62.35382907247956*phiUx[0]*rdxUxSq[0]+(140.296115413079*phiUx[0]-140.296115413079*phiLx[0])*rdxLx[0]*rdxUx[0]+62.35382907247956*phiLx[0]*rdxLxSq[0])*omega+(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0])*phiC[1])/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 

}

void MGpoissonDampedGaussSeidel1xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (6.873217490352689e-4*((3.0*rho[1]+86.60254037844386*rho[0])*omega*volFac+((623.5382907247956*phiUx[0]-1454.922678357857*phiC[0]+1175.755076535926*bcVals[0])*rdxCp2Sq[0]-720.0*rdxCp2Sq[0]*phiUx[1])*omega+1454.922678357857*phiC[0]*rdxCp2Sq[0]))/rdxCp2Sq[0]; 
  phiC[1] = (0.001030982623552903*((5.196152422706631*rho[1]+10.0*rho[0])*omega*volFac+((-277.1281292110203*rdxCp2Sq[0]*phiUx[1])-969.9484522385712*rdxCp2Sq[0]*phiC[1]+(240.0*phiUx[0]-339.411254969543*bcVals[0])*rdxCp2Sq[0])*omega+969.9484522385712*rdxCp2Sq[0]*phiC[1]))/rdxCp2Sq[0]; 

}

void MGpoissonDampedGaussSeidel1xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(0.004811252243246883*((9.0*rho[1]-51.96152422706631*rho[0])*omega*volFac+(240.0*rdxCp2Sq[0]*phiUx[1]+((-207.8460969082653*phiUx[0])+207.8460969082653*phiC[0]+195.9591794226543*bcVals[0])*rdxCp2Sq[0])*omega-207.8460969082653*phiC[0]*rdxCp2Sq[0]))/rdxCp2Sq[0]; 
  phiC[1] = (0.01443375672974065*((1.732050807568877*rho[1]-5.0*rho[0])*omega*volFac+(28.28427124746191*bcVals[0]*rdxCp2Sq[0]-69.28203230275508*rdxCp2Sq[0]*phiC[1])*omega+69.28203230275508*rdxCp2Sq[0]*phiC[1]))/rdxCp2Sq[0]; 

}

void MGpoissonDampedGaussSeidel1xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(6.873217490352689e-4*((3.0*rho[1]-86.60254037844386*rho[0])*omega*volFac+((-720.0*rdxCp2Sq[0]*phiLx[1])-1175.755076535926*rdxCp2Sq[0]*bcVals[1]+(1454.922678357857*phiC[0]-623.5382907247956*phiLx[0])*rdxCp2Sq[0])*omega-1454.922678357857*phiC[0]*rdxCp2Sq[0]))/rdxCp2Sq[0]; 
  phiC[1] = (0.001030982623552903*((5.196152422706631*rho[1]-10.0*rho[0])*omega*volFac+((-277.1281292110203*rdxCp2Sq[0]*phiLx[1])-969.9484522385712*rdxCp2Sq[0]*phiC[1]+339.411254969543*rdxCp2Sq[0]*bcVals[1]-240.0*phiLx[0]*rdxCp2Sq[0])*omega+969.9484522385712*rdxCp2Sq[0]*phiC[1]))/rdxCp2Sq[0]; 

}

void MGpoissonDampedGaussSeidel1xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phi) 
{ 
  // omega:  relaxation parameter.
  // dx:     cell lengths of cells pointed to by the stencil.
  // bcVals: values to impose as BCs.
  // rho:    right-side source in the current cell.
  // phi:    iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (0.004811252243246883*((9.0*rho[1]+51.96152422706631*rho[0])*omega*volFac+(240.0*rdxCp2Sq[0]*phiLx[1]+195.9591794226543*rdxCp2Sq[0]*bcVals[1]+(207.8460969082653*phiLx[0]-207.8460969082653*phiC[0])*rdxCp2Sq[0])*omega+207.8460969082653*phiC[0]*rdxCp2Sq[0]))/rdxCp2Sq[0]; 
  phiC[1] = (0.01443375672974065*((1.732050807568877*rho[1]+5.0*rho[0])*omega*volFac+(28.28427124746191*rdxCp2Sq[0]*bcVals[1]-69.28203230275508*rdxCp2Sq[0]*phiC[1])*omega+69.28203230275508*rdxCp2Sq[0]*phiC[1]))/rdxCp2Sq[0]; 

}

void MGpoissonResidue1xSer_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 
  double *dxUx = dx[1]; 
  double *dxLx = dx[2]; 

  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxLx[1]; 
  double rdxUx[1]; 
  double rdxLxSq[1]; 
  double rdxUxSq[1]; 
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxLx[0]   = 1.0/dxLx[0]; 
  rdxUx[0]   = 1.0/dxUx[0]; 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.125*(4.0*rho[0]*volFac-8.660254037844386*rdxUx[0]*phiUx[1]+8.660254037844386*rdxLx[0]*phiLx[1]+(8.660254037844386*rdxLx[0]-8.660254037844386*rdxUx[0])*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0])*rdxUx[0]+(9.0*phiLx[0]-9.0*phiC[0])*rdxLx[0]); 
  resOut[1] = 0.125*(4.0*rho[1]*volFac-7.0*rdxUx[0]*phiUx[1]-7.0*rdxLx[0]*phiLx[1]+((-23.0*rdxUx[0])-23.0*rdxLx[0])*phiC[1]+(8.660254037844386*phiUx[0]-22.5166604983954*phiC[0])*rdxUx[0]+(22.5166604983954*phiC[0]-8.660254037844386*phiLx[0])*rdxLx[0]); 

}

void MGpoissonResidue1xSer_LxDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.3535533905932737*(1.414213562373095*rho[0]*volFac-9.797958971132715*rdxCp2[0]*phiUx[1]+9.797958971132715*rdxCp2[0]*phiC[1]+(8.485281374238571*phiUx[0]-25.45584412271572*phiC[0]+24.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.2041241452319315*(2.449489742783178*rho[1]*volFac-97.97958971132716*rdxCp2[0]*phiUx[1]-489.8979485566358*rdxCp2[0]*phiC[1]+(84.85281374238573*phiUx[0]+84.85281374238573*phiC[0]-240.0*bcVals[0])*rdxCp2[0]); 

}

void MGpoissonResidue1xSer_LxNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.1178511301977579*(4.242640687119286*rho[0]*volFac-39.19183588453087*rdxCp2[0]*phiUx[1]-58.7877538267963*rdxCp2[0]*phiC[1]+(33.9411254969543*phiUx[0]-33.9411254969543*phiC[0]-8.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.06804138174397717*(7.348469228349534*rho[1]*volFac-195.9591794226543*rdxCp2[0]*phiUx[1]-587.877538267963*rdxCp2[0]*phiC[1]+(169.7056274847715*phiUx[0]-169.7056274847715*phiC[0]+80.0*bcVals[0])*rdxCp2[0]); 

}

void MGpoissonResidue1xSer_UxDirichlet_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.3535533905932737*(1.414213562373095*rho[0]*volFac+9.797958971132715*rdxCp2[0]*phiLx[1]-9.797958971132715*rdxCp2[0]*phiC[1]+24.0*rdxCp2[0]*bcVals[1]+(8.485281374238571*phiLx[0]-25.45584412271572*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.2041241452319315*(2.449489742783178*rho[1]*volFac-97.97958971132716*rdxCp2[0]*phiLx[1]-489.8979485566358*rdxCp2[0]*phiC[1]+240.0*rdxCp2[0]*bcVals[1]+((-84.85281374238573*phiLx[0])-84.85281374238573*phiC[0])*rdxCp2[0]); 

}

void MGpoissonResidue1xSer_UxNeumann_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdxCp2[1]; 
  double rdxCp2Sq[1]; 
  double rdxCp2R3[1]; 
  double rdxCp2R4[1]; 
  double rdxCp2R6[1]; 
  double rdxCp2R8[1]; 
  rdxCp2[0]  = volFac/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.1178511301977579*(4.242640687119286*rho[0]*volFac+39.19183588453087*rdxCp2[0]*phiLx[1]+58.7877538267963*rdxCp2[0]*phiC[1]+8.0*rdxCp2[0]*bcVals[1]+(33.9411254969543*phiLx[0]-33.9411254969543*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.06804138174397717*(7.348469228349534*rho[1]*volFac-195.9591794226543*rdxCp2[0]*phiLx[1]-587.877538267963*rdxCp2[0]*phiC[1]+80.0*rdxCp2[0]*bcVals[1]+(169.7056274847715*phiC[0]-169.7056274847715*phiLx[0])*rdxCp2[0]); 

}

