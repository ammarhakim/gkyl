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

void MGpoissonJacobi1xSer_P1(double **dx, const double *rho, double **phi) 
{ 
  // dx:  cell lengths of cells pointed to by the stencil.
  // rho: right-side source in the current cell.
  // phi: iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double *dxUx = dx[1]; 
  double *dxLx = dx[2]; 

  double rdxC2[1]; 
  double rdxLx[1]; 
  double rdxUx[1]; 
  double rdxLxSq[1]; 
  double rdxUxSq[1]; 
  rdxC2[0]   = 2.0/dxC[0]; 
  rdxLx[0]   = 1.0/dxLx[0]; 
  rdxUx[0]   = 1.0/dxUx[0]; 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 

  double volFac = 1.0/rdxC2[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(1.0*(((34.64101615137754*rdxUx[0]-34.64101615137754*rdxLx[0])*rho[1]+((-92.0*rdxUx[0])-92.0*rdxLx[0])*rho[0])*volFac+(69.28203230275508*rdxUxSq[0]+129.9038105676658*rdxLx[0]*rdxUx[0])*phiUx[1]+((-129.9038105676658*rdxLx[0]*rdxUx[0])-69.28203230275508*rdxLxSq[0])*phiLx[1]-66.0*phiUx[0]*rdxUxSq[0]+((-141.0*phiUx[0])-141.0*phiLx[0])*rdxLx[0]*rdxUx[0]-66.0*phiLx[0]*rdxLxSq[0]))/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 
  phiC[1] = (((36.0*rdxUx[0]+36.0*rdxLx[0])*rho[1]+(90.06664199358161*rdxLx[0]-90.06664199358161*rdxUx[0])*rho[0])*volFac+(66.0*rdxUxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(66.0*rdxLxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiLx[1]-62.35382907247956*phiUx[0]*rdxUxSq[0]+(140.296115413079*phiUx[0]-140.296115413079*phiLx[0])*rdxLx[0]*rdxUx[0]+62.35382907247956*phiLx[0]*rdxLxSq[0])/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 

}

void MGpoissonDampedJacobi1xSer_P1(const double omega, double **dx, const double *rho, double **phi) 
{ 
  // omega:     relaxation parameter.
  // dx:  cell lengths of cells pointed to by the stencil.
  // rho: right-side source in the current cell.
  // phi: iterate cells pointed to by the stencil.

  double *dxC  = dx[0]; 
  double *dxUx = dx[1]; 
  double *dxLx = dx[2]; 

  double rdxC2[1]; 
  double rdxLx[1]; 
  double rdxUx[1]; 
  double rdxLxSq[1]; 
  double rdxUxSq[1]; 
  rdxC2[0]   = 2.0/dxC[0]; 
  rdxLx[0]   = 1.0/dxLx[0]; 
  rdxUx[0]   = 1.0/dxUx[0]; 
  rdxLxSq[0] = rdxLx[0]*rdxLx[0]; 
  rdxUxSq[0] = rdxUx[0]*rdxUx[0]; 

  double volFac = 1.0/rdxC2[0]; 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(1.0*(((34.64101615137754*rdxUx[0]-34.64101615137754*rdxLx[0])*rho[1]+((-92.0*rdxUx[0])-92.0*rdxLx[0])*rho[0])*omega*volFac+((69.28203230275508*rdxUxSq[0]+129.9038105676658*rdxLx[0]*rdxUx[0])*phiUx[1]+((-129.9038105676658*rdxLx[0]*rdxUx[0])-69.28203230275508*rdxLxSq[0])*phiLx[1]+(6.0*phiC[0]-66.0*phiUx[0])*rdxUxSq[0]+((-141.0*phiUx[0])-141.0*phiLx[0]+402.0*phiC[0])*rdxLx[0]*rdxUx[0]+(6.0*phiC[0]-66.0*phiLx[0])*rdxLxSq[0])*omega-6.0*phiC[0]*rdxUxSq[0]-402.0*phiC[0]*rdxLx[0]*rdxUx[0]-6.0*phiC[0]*rdxLxSq[0]))/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 
  phiC[1] = (((36.0*rdxUx[0]+36.0*rdxLx[0])*rho[1]+(90.06664199358161*rdxLx[0]-90.06664199358161*rdxUx[0])*rho[0])*omega*volFac+((66.0*rdxUxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(66.0*rdxLxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiLx[1]+((-6.0*rdxUxSq[0])-402.0*rdxLx[0]*rdxUx[0]-6.0*rdxLxSq[0])*phiC[1]-62.35382907247956*phiUx[0]*rdxUxSq[0]+(140.296115413079*phiUx[0]-140.296115413079*phiLx[0])*rdxLx[0]*rdxUx[0]+62.35382907247956*phiLx[0]*rdxLxSq[0])*omega+(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0])*phiC[1])/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 

}

