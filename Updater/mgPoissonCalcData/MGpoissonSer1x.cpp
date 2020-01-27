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

  phiC[0] = (0.005498573992282151*((6.0*rho[1]+86.60254037844386*rho[0])*volFac-90.0*rdxCp2[0]*phiUx[1]+(77.94228634059945*phiUx[0]+146.9693845669907*bcVals[0])*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.01649572197684645*((5.196152422706631*rho[1]+5.0*rho[0])*volFac-17.32050807568877*rdxCp2[0]*phiUx[1]+(15.0*phiUx[0]-21.21320343559643*bcVals[0])*rdxCp2[0]))/rdxCp2[0]; 

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

  phiC[0] = -(0.03849001794597506*((18.0*rho[1]-51.96152422706631*rho[0])*volFac+30.0*rdxCp2[0]*phiUx[1]+(48.98979485566358*bcVals[0]-25.98076211353316*phiUx[0])*rdxCp2[0]))/rdxCp2[0]; 
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

  phiC[0] = -(0.005498573992282151*((6.0*rho[1]-86.60254037844386*rho[0])*volFac-90.0*rdxCp2[0]*phiLx[1]-146.9693845669907*rdxCp2[0]*bcVals[1]-77.94228634059945*phiLx[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.01649572197684645*((5.196152422706631*rho[1]-5.0*rho[0])*volFac-17.32050807568877*rdxCp2[0]*phiLx[1]+21.21320343559643*rdxCp2[0]*bcVals[1]-15.0*phiLx[0]*rdxCp2[0]))/rdxCp2[0]; 

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

  phiC[0] = (0.03849001794597506*((18.0*rho[1]+51.96152422706631*rho[0])*volFac+30.0*rdxCp2[0]*phiLx[1]+48.98979485566358*rdxCp2[0]*bcVals[1]+25.98076211353316*phiLx[0]*rdxCp2[0]))/rdxCp2[0]; 
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

  phiC[0] = (0.005498573992282151*((6.0*rho[1]+86.60254037844386*rho[0])*omega*volFac+((77.94228634059945*phiUx[0]-181.8653347947321*phiC[0]+146.9693845669907*bcVals[0])*rdxCp2[0]-90.0*rdxCp2[0]*phiUx[1])*omega+181.8653347947321*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.01649572197684645*((5.196152422706631*rho[1]+5.0*rho[0])*omega*volFac+((-17.32050807568877*rdxCp2[0]*phiUx[1])-60.6217782649107*rdxCp2[0]*phiC[1]+(15.0*phiUx[0]-21.21320343559643*bcVals[0])*rdxCp2[0])*omega+60.6217782649107*rdxCp2[0]*phiC[1]))/rdxCp2[0]; 

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

  phiC[0] = -(0.03849001794597506*((18.0*rho[1]-51.96152422706631*rho[0])*omega*volFac+(30.0*rdxCp2[0]*phiUx[1]+((-25.98076211353316*phiUx[0])+25.98076211353316*phiC[0]+48.98979485566358*bcVals[0])*rdxCp2[0])*omega-25.98076211353316*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
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

  phiC[0] = -(0.005498573992282151*((6.0*rho[1]-86.60254037844386*rho[0])*omega*volFac+((-90.0*rdxCp2[0]*phiLx[1])-146.9693845669907*rdxCp2[0]*bcVals[1]+(181.8653347947321*phiC[0]-77.94228634059945*phiLx[0])*rdxCp2[0])*omega-181.8653347947321*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.01649572197684645*((5.196152422706631*rho[1]-5.0*rho[0])*omega*volFac+((-17.32050807568877*rdxCp2[0]*phiLx[1])-60.6217782649107*rdxCp2[0]*phiC[1]+21.21320343559643*rdxCp2[0]*bcVals[1]-15.0*phiLx[0]*rdxCp2[0])*omega+60.6217782649107*rdxCp2[0]*phiC[1]))/rdxCp2[0]; 

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

  phiC[0] = (0.03849001794597506*((18.0*rho[1]+51.96152422706631*rho[0])*omega*volFac+(30.0*rdxCp2[0]*phiLx[1]+48.98979485566358*rdxCp2[0]*bcVals[1]+(25.98076211353316*phiLx[0]-25.98076211353316*phiC[0])*rdxCp2[0])*omega+25.98076211353316*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
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

  resOut[0] = 0.1767766952966368*(5.656854249492382*rho[0]*volFac-4.898979485566357*rdxCp2[0]*phiUx[1]+4.898979485566357*rdxCp2[0]*phiC[1]+(4.242640687119286*phiUx[0]-12.72792206135786*phiC[0]+12.0*bcVals[0])*rdxCp2[0]); 
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

  resOut[0] = 0.2357022603955158*(4.242640687119286*rho[0]*volFac-4.898979485566357*rdxCp2[0]*phiUx[1]-7.348469228349534*rdxCp2[0]*phiC[1]+(4.242640687119286*phiUx[0]-4.242640687119286*phiC[0]-2.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.06804138174397717*(14.69693845669907*rho[1]*volFac-24.49489742783179*rdxCp2[0]*phiUx[1]-73.48469228349535*rdxCp2[0]*phiC[1]+(21.21320343559643*phiUx[0]-21.21320343559643*phiC[0]+20.0*bcVals[0])*rdxCp2[0]); 

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

  resOut[0] = 0.1767766952966368*(5.656854249492382*rho[0]*volFac+4.898979485566357*rdxCp2[0]*phiLx[1]-4.898979485566357*rdxCp2[0]*phiC[1]+12.0*rdxCp2[0]*bcVals[1]+(4.242640687119286*phiLx[0]-12.72792206135786*phiC[0])*rdxCp2[0]); 
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

  resOut[0] = 0.2357022603955158*(4.242640687119286*rho[0]*volFac+4.898979485566357*rdxCp2[0]*phiLx[1]+7.348469228349534*rdxCp2[0]*phiC[1]+2.0*rdxCp2[0]*bcVals[1]+(4.242640687119286*phiLx[0]-4.242640687119286*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.06804138174397717*(14.69693845669907*rho[1]*volFac-24.49489742783179*rdxCp2[0]*phiLx[1]-73.48469228349535*rdxCp2[0]*phiC[1]+20.0*rdxCp2[0]*bcVals[1]+(21.21320343559643*phiC[0]-21.21320343559643*phiLx[0])*rdxCp2[0]); 

}

