#include <MGpoissonModDecl.h> 
 
void MGpoissonFEMProlong1xSer_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];

  fldF0[0] += 0.5*fldC[0]; 
  fldF1[0] += fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
}

void MGpoissonFEM_DGtoFEM_1xSer_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 

  femFld[0] += 0.3535533905932737*dgFld[0]-0.6123724356957944*dgFld[1]; 
  femFldUx[0] += 0.6123724356957944*dgFld[1]+0.3535533905932737*dgFld[0]; 

}

void MGpoissonFEM_DGtoFEM_1xSer_Lx_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 

  femFld[0] += 0.7071067811865475*dgFld[0]-1.224744871391589*dgFld[1]; 
  femFldUx[0] += 0.6123724356957944*dgFld[1]+0.3535533905932737*dgFld[0]; 

}

void MGpoissonFEM_DGtoFEM_1xSer_Ux_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 

  femFld[0] += 0.3535533905932737*dgFld[0]-0.6123724356957944*dgFld[1]; 
  femFld[1] += 1.224744871391589*dgFld[1]+0.7071067811865475*dgFld[0]; 

}

void MGpoissonFEMDampedGaussSeidel1xSer_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 


  phiC[0] = (0.1666666666666667*((2.0*rhoUx[0]+2.0*rhoLx[0]+8.0*rhoC[0])*omega*volFac+(3.0*phiUx[0]+3.0*phiLx[0]-6.0*phiC[0])*rdx2SqVol[0]*omega+6.0*phiC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonFEMDampedGaussSeidel1xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 


  phiC[0] = -(1.0*((2.0*bcVals[2]+(phiC[0]-1.0*phiUx[0])*bcVals[1]-2.0*bcVals[0]*phiC[0])*omega-1.0*phiC[0]*bcVals[1]+2.0*bcVals[0]*phiC[0]))/(bcVals[1]-2.0*bcVals[0]); 

}

void MGpoissonFEMDampedGaussSeidel1xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 


  phiC[0] = (((2.0*rhoC[1]+2.0*rhoLx[0]+8.0*rhoC[0])*bcVals[4]+(4.0*rhoC[1]+4.0*rhoLx[0]+16.0*rhoC[0])*bcVals[3])*omega*volFac+(6.0*rdx2SqVol[0]*bcVals[5]+(3.0*phiLx[0]-3.0*phiC[0])*rdx2SqVol[0]*bcVals[4]+(6.0*phiLx[0]-12.0*phiC[0])*rdx2SqVol[0]*bcVals[3])*omega+3.0*phiC[0]*rdx2SqVol[0]*bcVals[4]+12.0*phiC[0]*rdx2SqVol[0]*bcVals[3])/(3.0*rdx2SqVol[0]*bcVals[4]+12.0*rdx2SqVol[0]*bcVals[3]); 
  phiC[1] = ((2.0*rhoC[1]+2.0*rhoLx[0]+8.0*rhoC[0])*bcVals[4]*omega*volFac+(12.0*rdx2SqVol[0]*bcVals[5]+(3.0*phiLx[0]*rdx2SqVol[0]-3.0*rdx2SqVol[0]*phiC[1])*bcVals[4]-12.0*rdx2SqVol[0]*phiC[1]*bcVals[3])*omega+3.0*rdx2SqVol[0]*phiC[1]*bcVals[4]+12.0*rdx2SqVol[0]*phiC[1]*bcVals[3])/(3.0*rdx2SqVol[0]*bcVals[4]+12.0*rdx2SqVol[0]*bcVals[3]); 

}

