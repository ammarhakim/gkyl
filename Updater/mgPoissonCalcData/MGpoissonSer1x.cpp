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

  phiC[0] = -(1.0*(((34.64101615137754*rdxUx[0]-34.64101615137754*rdxLx[0])*rho[1]+((-92.0*rdxUx[0])-92.0*rdxLx[0])*rho[0])*volFac+(69.28203230275508*rdxUxSq[0]+129.9038105676658*rdxLx[0]*rdxUx[0])*phiUx[1]+((-129.9038105676658*rdxLx[0]*rdxUx[0])-69.28203230275508*rdxLxSq[0])*phiLx[1]-66.0*phiUx[0]*rdxUxSq[0]+((-141.0*phiUx[0])-141.0*phiLx[0])*rdxLx[0]*rdxUx[0]-66.0*phiLx[0]*rdxLxSq[0]))/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 
  phiC[1] = (((36.0*rdxUx[0]+36.0*rdxLx[0])*rho[1]+(90.06664199358161*rdxLx[0]-90.06664199358161*rdxUx[0])*rho[0])*volFac+(66.0*rdxUxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(66.0*rdxLxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiLx[1]-62.35382907247956*phiUx[0]*rdxUxSq[0]+(140.296115413079*phiUx[0]-140.296115413079*phiLx[0])*rdxLx[0]*rdxUx[0]+62.35382907247956*phiLx[0]*rdxLxSq[0])/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 

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

  phiC[0] = (0.001374643498070538*((3.0*rho[1]+86.60254037844386*rho[0])*volFac-360.0*rdxCp2[0]*phiUx[1]+(311.7691453623978*phiUx[0]+587.877538267963*bcVals[0])*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.002061965247105807*((5.196152422706631*rho[1]+10.0*rho[0])*volFac-138.5640646055102*rdxCp2[0]*phiUx[1]+(120.0*phiUx[0]-169.7056274847715*bcVals[0])*rdxCp2[0]))/rdxCp2[0]; 

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

  phiC[0] = -(0.009622504486493766*((9.0*rho[1]-51.96152422706631*rho[0])*volFac+120.0*rdxCp2[0]*phiUx[1]+(97.97958971132716*bcVals[0]-103.9230484541326*phiUx[0])*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.02886751345948129*((1.732050807568877*rho[1]-5.0*rho[0])*volFac+14.14213562373095*bcVals[0]*rdxCp2[0]))/rdxCp2[0]; 

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

  phiC[0] = -(0.001374643498070538*((3.0*rho[1]-86.60254037844386*rho[0])*volFac-360.0*rdxCp2[0]*phiLx[1]-587.877538267963*rdxCp2[0]*bcVals[1]-311.7691453623978*phiLx[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.002061965247105807*((5.196152422706631*rho[1]-10.0*rho[0])*volFac-138.5640646055102*rdxCp2[0]*phiLx[1]+169.7056274847715*rdxCp2[0]*bcVals[1]-120.0*phiLx[0]*rdxCp2[0]))/rdxCp2[0]; 

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

  phiC[0] = (0.009622504486493766*((9.0*rho[1]+51.96152422706631*rho[0])*volFac+120.0*rdxCp2[0]*phiLx[1]+97.97958971132716*rdxCp2[0]*bcVals[1]+103.9230484541326*phiLx[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.02886751345948129*((1.732050807568877*rho[1]+5.0*rho[0])*volFac+14.14213562373095*rdxCp2[0]*bcVals[1]))/rdxCp2[0]; 

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

  phiC[0] = -(1.0*(((34.64101615137754*rdxUx[0]-34.64101615137754*rdxLx[0])*rho[1]+((-92.0*rdxUx[0])-92.0*rdxLx[0])*rho[0])*omega*volFac+((69.28203230275508*rdxUxSq[0]+129.9038105676658*rdxLx[0]*rdxUx[0])*phiUx[1]+((-129.9038105676658*rdxLx[0]*rdxUx[0])-69.28203230275508*rdxLxSq[0])*phiLx[1]+(6.0*phiC[0]-66.0*phiUx[0])*rdxUxSq[0]+((-141.0*phiUx[0])-141.0*phiLx[0]+402.0*phiC[0])*rdxLx[0]*rdxUx[0]+(6.0*phiC[0]-66.0*phiLx[0])*rdxLxSq[0])*omega-6.0*phiC[0]*rdxUxSq[0]-402.0*phiC[0]*rdxLx[0]*rdxUx[0]-6.0*phiC[0]*rdxLxSq[0]))/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 
  phiC[1] = (((36.0*rdxUx[0]+36.0*rdxLx[0])*rho[1]+(90.06664199358161*rdxLx[0]-90.06664199358161*rdxUx[0])*rho[0])*omega*volFac+((66.0*rdxUxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiUx[1]+(66.0*rdxLxSq[0]-129.0*rdxLx[0]*rdxUx[0])*phiLx[1]+((-6.0*rdxUxSq[0])-402.0*rdxLx[0]*rdxUx[0]-6.0*rdxLxSq[0])*phiC[1]-62.35382907247956*phiUx[0]*rdxUxSq[0]+(140.296115413079*phiUx[0]-140.296115413079*phiLx[0])*rdxLx[0]*rdxUx[0]+62.35382907247956*phiLx[0]*rdxLxSq[0])*omega+(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0])*phiC[1])/(6.0*rdxUxSq[0]+402.0*rdxLx[0]*rdxUx[0]+6.0*rdxLxSq[0]); 

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

  phiC[0] = (0.001374643498070538*((3.0*rho[1]+86.60254037844386*rho[0])*omega*volFac+((311.7691453623978*phiUx[0]-727.4613391789284*phiC[0]+587.877538267963*bcVals[0])*rdxCp2[0]-360.0*rdxCp2[0]*phiUx[1])*omega+727.4613391789284*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.002061965247105807*((5.196152422706631*rho[1]+10.0*rho[0])*omega*volFac+((-138.5640646055102*rdxCp2[0]*phiUx[1])-484.9742261192856*rdxCp2[0]*phiC[1]+(120.0*phiUx[0]-169.7056274847715*bcVals[0])*rdxCp2[0])*omega+484.9742261192856*rdxCp2[0]*phiC[1]))/rdxCp2[0]; 

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

  phiC[0] = -(0.009622504486493766*((9.0*rho[1]-51.96152422706631*rho[0])*omega*volFac+(120.0*rdxCp2[0]*phiUx[1]+((-103.9230484541326*phiUx[0])+103.9230484541326*phiC[0]+97.97958971132716*bcVals[0])*rdxCp2[0])*omega-103.9230484541326*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.02886751345948129*((1.732050807568877*rho[1]-5.0*rho[0])*omega*volFac+(14.14213562373095*bcVals[0]*rdxCp2[0]-34.64101615137754*rdxCp2[0]*phiC[1])*omega+34.64101615137754*rdxCp2[0]*phiC[1]))/rdxCp2[0]; 

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

  phiC[0] = -(0.001374643498070538*((3.0*rho[1]-86.60254037844386*rho[0])*omega*volFac+((-360.0*rdxCp2[0]*phiLx[1])-587.877538267963*rdxCp2[0]*bcVals[1]+(727.4613391789284*phiC[0]-311.7691453623978*phiLx[0])*rdxCp2[0])*omega-727.4613391789284*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.002061965247105807*((5.196152422706631*rho[1]-10.0*rho[0])*omega*volFac+((-138.5640646055102*rdxCp2[0]*phiLx[1])-484.9742261192856*rdxCp2[0]*phiC[1]+169.7056274847715*rdxCp2[0]*bcVals[1]-120.0*phiLx[0]*rdxCp2[0])*omega+484.9742261192856*rdxCp2[0]*phiC[1]))/rdxCp2[0]; 

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

  phiC[0] = (0.009622504486493766*((9.0*rho[1]+51.96152422706631*rho[0])*omega*volFac+(120.0*rdxCp2[0]*phiLx[1]+97.97958971132716*rdxCp2[0]*bcVals[1]+(103.9230484541326*phiLx[0]-103.9230484541326*phiC[0])*rdxCp2[0])*omega+103.9230484541326*phiC[0]*rdxCp2[0]))/rdxCp2[0]; 
  phiC[1] = (0.02886751345948129*((1.732050807568877*rho[1]+5.0*rho[0])*omega*volFac+(14.14213562373095*rdxCp2[0]*bcVals[1]-34.64101615137754*rdxCp2[0]*phiC[1])*omega+34.64101615137754*rdxCp2[0]*phiC[1]))/rdxCp2[0]; 

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

  resOut[0] = 0.125*(8.0*rho[0]*volFac-8.660254037844386*rdxUx[0]*phiUx[1]+8.660254037844386*rdxLx[0]*phiLx[1]+(8.660254037844386*rdxLx[0]-8.660254037844386*rdxUx[0])*phiC[1]+(9.0*phiUx[0]-9.0*phiC[0])*rdxUx[0]+(9.0*phiLx[0]-9.0*phiC[0])*rdxLx[0]); 
  resOut[1] = 0.125*(8.0*rho[1]*volFac-7.0*rdxUx[0]*phiUx[1]-7.0*rdxLx[0]*phiLx[1]+((-23.0*rdxUx[0])-23.0*rdxLx[0])*phiC[1]+(8.660254037844386*phiUx[0]-22.5166604983954*phiC[0])*rdxUx[0]+(22.5166604983954*phiC[0]-8.660254037844386*phiLx[0])*rdxLx[0]); 

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

  resOut[0] = 0.7071067811865475*(1.414213562373095*rho[0]*volFac-4.898979485566357*rdxCp2[0]*phiUx[1]+4.898979485566357*rdxCp2[0]*phiC[1]+(4.242640687119286*phiUx[0]-12.72792206135786*phiC[0]+12.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.408248290463863*(2.449489742783178*rho[1]*volFac-48.98979485566358*rdxCp2[0]*phiUx[1]-244.9489742783179*rdxCp2[0]*phiC[1]+(42.42640687119286*phiUx[0]+42.42640687119286*phiC[0]-120.0*bcVals[0])*rdxCp2[0]); 

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

  resOut[0] = 0.2357022603955158*(4.242640687119286*rho[0]*volFac-19.59591794226543*rdxCp2[0]*phiUx[1]-29.39387691339815*rdxCp2[0]*phiC[1]+(16.97056274847715*phiUx[0]-16.97056274847715*phiC[0]-4.0*bcVals[0])*rdxCp2[0]); 
  resOut[1] = 0.1360827634879543*(7.348469228349534*rho[1]*volFac-97.97958971132716*rdxCp2[0]*phiUx[1]-293.9387691339815*rdxCp2[0]*phiC[1]+(84.85281374238573*phiUx[0]-84.85281374238573*phiC[0]+40.0*bcVals[0])*rdxCp2[0]); 

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

  resOut[0] = 0.7071067811865475*(1.414213562373095*rho[0]*volFac+4.898979485566357*rdxCp2[0]*phiLx[1]-4.898979485566357*rdxCp2[0]*phiC[1]+12.0*rdxCp2[0]*bcVals[1]+(4.242640687119286*phiLx[0]-12.72792206135786*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.408248290463863*(2.449489742783178*rho[1]*volFac-48.98979485566358*rdxCp2[0]*phiLx[1]-244.9489742783179*rdxCp2[0]*phiC[1]+120.0*rdxCp2[0]*bcVals[1]+((-42.42640687119286*phiLx[0])-42.42640687119286*phiC[0])*rdxCp2[0]); 

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

  resOut[0] = 0.2357022603955158*(4.242640687119286*rho[0]*volFac+19.59591794226543*rdxCp2[0]*phiLx[1]+29.39387691339815*rdxCp2[0]*phiC[1]+4.0*rdxCp2[0]*bcVals[1]+(16.97056274847715*phiLx[0]-16.97056274847715*phiC[0])*rdxCp2[0]); 
  resOut[1] = 0.1360827634879543*(7.348469228349534*rho[1]*volFac-97.97958971132716*rdxCp2[0]*phiLx[1]-293.9387691339815*rdxCp2[0]*phiC[1]+40.0*rdxCp2[0]*bcVals[1]+(84.85281374238573*phiC[0]-84.85281374238573*phiLx[0])*rdxCp2[0]); 

}

