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
  rdxLx[0]   = volFac*4.0/(dxLx[0]*dxLx[0]); 
  rdxUx[0]   = volFac*4.0/(dxUx[0]*dxUx[0]); 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(1.0*(((69.28203230275508*rdxUx[0]-69.28203230275508*rdxLx[0])*rho[1]+((-184.0*rdxUx[0])-184.0*rdxLx[0])*rho[0])*volFac+(69.28203230275508*rdxUxSq[0]+129.9038105676658*rdxLx[0]*rdxUx[0])*phiUx[1]+((-129.9038105676658*rdxLx[0]*rdxUx[0])-69.28203230275508*rdxLxSq[0])*phiLx[1]-66.0*phiUx[0]*rdxUxSq[0]+((-141.0*phiUx[0])-141.0*phiLx[0])*rdxLx[0]*rdxUx[0]-66.0*phiLx[0]*rdxLxSq[0]))/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 
  phiC[1] = (((72.0*rdxUx[0]+72.0*rdxLx[0])*rho[1]+(180.1332839871632*rdxLx[0]-180.1332839871632*rdxUx[0])*rho[0])*volFac+(66.0*rdxUxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(66.0*rdxLxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiLx[1]-62.35382907247956*phiUx[0]*rdxUxSq[0]+(140.296115413079*phiUx[0]-140.296115413079*phiLx[0])*rdxLx[0]*rdxUx[0]+62.35382907247956*phiLx[0]*rdxLxSq[0])/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (8.951166964180246e-4*((324.0*rho[1]+692.8203230275509*rho[0])*volFac-360.0*rdxCp2[0]*phiUx[1]+(441.6729559300637*phiUx[0]+955.3009996854395*bcVals[0])*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.002685350089254074*((48.49742261192856*rho[1]+40.0*rho[0])*volFac-95.26279441628824*rdxCp2[0]*phiUx[1]+(90.0*phiUx[0]-127.2792206135786*bcVals[0])*rdxCp2[0]))/rdxCp2[0]; 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(0.06415002990995844*((6.0*rho[1]-36.37306695894642*rho[0])*volFac+15.0*rdxCp2[0]*phiUx[1]+(31.84336665618132*bcVals[0]-15.58845726811989*phiUx[0])*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.1154700538379252*((3.464101615137754*rho[1]-5.0*rho[0])*volFac+7.071067811865476*bcVals[0]*rdxCp2[0]))/rdxCp2[0]; 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(8.951166964180246e-4*((324.0*rho[1]-692.8203230275509*rho[0])*volFac-360.0*rdxCp2[0]*phiLx[1]-955.3009996854395*rdxCp2[0]*bcVals[1]-441.6729559300637*phiLx[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.002685350089254074*((48.49742261192856*rho[1]-40.0*rho[0])*volFac-95.26279441628824*rdxCp2[0]*phiLx[1]+127.2792206135786*rdxCp2[0]*bcVals[1]-90.0*phiLx[0]*rdxCp2[0]))/rdxCp2[0]; 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (0.06415002990995844*((6.0*rho[1]+36.37306695894642*rho[0])*volFac+15.0*rdxCp2[0]*phiLx[1]+31.84336665618132*rdxCp2[0]*bcVals[1]+15.58845726811989*phiLx[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.1154700538379252*((3.464101615137754*rho[1]+5.0*rho[0])*volFac+7.071067811865476*rdxCp2[0]*bcVals[1]))/rdxCp2[0]; 

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
  rdxLx[0]   = volFac*4.0/(dxLx[0]*dxLx[0]); 
  rdxUx[0]   = volFac*4.0/(dxUx[0]*dxUx[0]); 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(1.0*(((69.28203230275508*rdxUx[0]-69.28203230275508*rdxLx[0])*rho[1]+((-184.0*rdxUx[0])-184.0*rdxLx[0])*rho[0])*omega*volFac+((69.28203230275508*rdxUxSq[0]+129.9038105676658*rdxLx[0]*rdxUx[0])*phiUx[1]+((-129.9038105676658*rdxLx[0]*rdxUx[0])-69.28203230275508*rdxLxSq[0])*phiLx[1]+(6.0*phiC[0]-66.0*phiUx[0])*rdxUxSq[0]+((-141.0*phiUx[0])-141.0*phiLx[0]+402.0*phiC[0])*rdxLx[0]*rdxUx[0]+(6.0*phiC[0]-66.0*phiLx[0])*rdxLxSq[0])*omega-6.0*phiC[0]*rdxUxSq[0]-402.0*phiC[0]*rdxLx[0]*rdxUx[0]-6.0*phiC[0]*rdxLxSq[0]))/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 
  phiC[1] = (((72.0*rdxUx[0]+72.0*rdxLx[0])*rho[1]+(180.1332839871632*rdxLx[0]-180.1332839871632*rdxUx[0])*rho[0])*omega*volFac+((66.0*rdxUxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(66.0*rdxLxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiLx[1]+((-6.0*rdxUxSq[0])-402.0*rdxLx[0]*rdxUx[0]-6.0*rdxLxSq[0])*phiC[1]-62.35382907247956*phiUx[0]*rdxUxSq[0]+(140.296115413079*phiUx[0]-140.296115413079*phiLx[0])*rdxLx[0]*rdxUx[0]+62.35382907247956*phiLx[0]*rdxLxSq[0])*omega+(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0])*phiC[1])/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (8.951166964180246e-4*((324.0*rho[1]+692.8203230275509*rho[0])*omega*volFac+((441.6729559300637*phiUx[0]-1117.172770881926*phiC[0]+955.3009996854395*bcVals[0])*rdxCp2[0]-360.0*rdxCp2[0]*phiUx[1])*omega+1117.172770881926*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.002685350089254074*((48.49742261192856*rho[1]+40.0*rho[0])*omega*volFac+((-95.26279441628824*rdxCp2[0]*phiUx[1])-372.3909236273086*rdxCp2[0]*phiC[1]+(90.0*phiUx[0]-127.2792206135786*bcVals[0])*rdxCp2[0])*omega+372.3909236273086*rdxCp2[0]*phiC[1]))/rdxCp2[0]; 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(0.06415002990995844*((6.0*rho[1]-36.37306695894642*rho[0])*omega*volFac+(15.0*rdxCp2[0]*phiUx[1]+((-15.58845726811989*phiUx[0])+15.58845726811989*phiC[0]+31.84336665618132*bcVals[0])*rdxCp2[0])*omega-15.58845726811989*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.1154700538379252*((3.464101615137754*rho[1]-5.0*rho[0])*omega*volFac+(7.071067811865476*bcVals[0]*rdxCp2[0]-8.660254037844386*rdxCp2[0]*phiC[1])*omega+8.660254037844386*rdxCp2[0]*phiC[1]))/rdxCp2[0]; 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(8.951166964180246e-4*((324.0*rho[1]-692.8203230275509*rho[0])*omega*volFac+((-360.0*rdxCp2[0]*phiLx[1])-955.3009996854395*rdxCp2[0]*bcVals[1]+(1117.172770881926*phiC[0]-441.6729559300637*phiLx[0])*rdxCp2[0])*omega-1117.172770881926*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.002685350089254074*((48.49742261192856*rho[1]-40.0*rho[0])*omega*volFac+((-95.26279441628824*rdxCp2[0]*phiLx[1])-372.3909236273086*rdxCp2[0]*phiC[1]+127.2792206135786*rdxCp2[0]*bcVals[1]-90.0*phiLx[0]*rdxCp2[0])*omega+372.3909236273086*rdxCp2[0]*phiC[1]))/rdxCp2[0]; 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (0.06415002990995844*((6.0*rho[1]+36.37306695894642*rho[0])*omega*volFac+(15.0*rdxCp2[0]*phiLx[1]+31.84336665618132*rdxCp2[0]*bcVals[1]+(15.58845726811989*phiLx[0]-15.58845726811989*phiC[0])*rdxCp2[0])*omega+15.58845726811989*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.1154700538379252*((3.464101615137754*rho[1]+5.0*rho[0])*omega*volFac+(7.071067811865476*rdxCp2[0]*bcVals[1]-8.660254037844386*rdxCp2[0]*phiC[1])*omega+8.660254037844386*rdxCp2[0]*phiC[1]))/rdxCp2[0]; 

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
  rdxLx[0]   = volFac*4.0/(dxLx[0]*dxLx[0]); 
  rdxUx[0]   = volFac*4.0/(dxUx[0]*dxUx[0]); 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdxUx[0]*phiUx[1]+8.660254037844386*rdxLx[0]*phiLx[1]+(8.660254037844386*rdxLx[0]-8.660254037844386*rdxUx[0])*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0])*rdxUx[0]+(9.0*phiLx[0]-9.0*phiC[0])*rdxLx[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-7.0*rdxUx[0]*phiUx[1]-7.0*rdxLx[0]*phiLx[1]+((-23.0*rdxUx[0])-23.0*rdxLx[0])*phiC[1]+(8.660254037844386*phiUx[0]-22.5166604983954*phiC[0])*rdxUx[0]+(22.5166604983954*phiC[0]-8.660254037844386*phiLx[0])*rdxLx[0]); 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.0883883476483184*(11.31370849898477*rho[0]*volFac+7.348469228349534*rdxCp2[0]*phiUx[1]+66.13622305514579*rdxCp2[0]*phiC[1]+((-4.242640687119286*phiUx[0])-29.698484809835*phiC[0]+48.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.1020620726159657*(9.797958971132715*rho[1]*volFac-24.49489742783179*rdxCp2[0]*phiUx[1]-122.4744871391589*rdxCp2[0]*phiC[1]+(21.21320343559643*phiUx[0]+21.21320343559643*phiC[0]-60.0*bcVals[0])*rdxCp2[0]); 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.0441941738241592*(22.62741699796953*rho[0]*volFac-12.24744871391589*rdxCp2[0]*phiUx[1]-12.24744871391589*rdxCp2[0]*phiC[1]+(12.72792206135786*phiUx[0]-12.72792206135786*phiC[0]-16.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.01275775907699571*(78.38367176906175*rho[1]*volFac-61.23724356957945*rdxCp2[0]*phiUx[1]-257.1964229922337*rdxCp2[0]*phiC[1]+(63.63961030678928*phiUx[0]-63.63961030678928*phiC[0]+80.0*bcVals[0])*rdxCp2[0]); 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.0883883476483184*(11.31370849898477*rho[0]*volFac-7.348469228349534*rdxCp2[0]*phiLx[1]-66.13622305514579*rdxCp2[0]*phiC[1]+48.0*rdxCp2[0]*bcVals[1]+((-4.242640687119286*phiLx[0])-29.698484809835*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.1020620726159657*(9.797958971132715*rho[1]*volFac-24.49489742783179*rdxCp2[0]*phiLx[1]-122.4744871391589*rdxCp2[0]*phiC[1]+60.0*rdxCp2[0]*bcVals[1]+((-21.21320343559643*phiLx[0])-21.21320343559643*phiC[0])*rdxCp2[0]); 

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
  rdxCp2[0]  = volFac*4.0/(dxC[0]*dxC[0]); 
  rdxCp2Sq[0]  = rdxCp2[0]*rdxCp2[0]; 
  rdxCp2R3[0]  = rdxCp2[0]*rdxCp2Sq[0]; 
  rdxCp2R4[0]  = rdxCp2Sq[0]*rdxCp2Sq[0]; 
  rdxCp2R6[0]  = rdxCp2Sq[0]*rdxCp2R4[0]; 
  rdxCp2R8[0]  = rdxCp2R4[0]*rdxCp2R4[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.0441941738241592*(22.62741699796953*rho[0]*volFac+12.24744871391589*rdxCp2[0]*phiLx[1]+12.24744871391589*rdxCp2[0]*phiC[1]+16.0*rdxCp2[0]*bcVals[1]+(12.72792206135786*phiLx[0]-12.72792206135786*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.01275775907699571*(78.38367176906175*rho[1]*volFac-61.23724356957945*rdxCp2[0]*phiLx[1]-257.1964229922337*rdxCp2[0]*phiC[1]+80.0*rdxCp2[0]*bcVals[1]+(63.63961030678928*phiC[0]-63.63961030678928*phiLx[0])*rdxCp2[0]); 

}

