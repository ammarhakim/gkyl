#include <MGpoissonModDecl.h> 
 
void MGpoissonFEMprolong1xSer_P1(double *fldCC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];
  double *fldFLxx = fldF[2];

  fldFC[0] += 0.5*fldCC[0]; 
  fldFLx[0] += fldCC[0]; 
  fldFLxx[0] += 0.5*fldCC[0]; 
}

void MGpoissonFEMprolong1xSer_LxDirichlet_P1(double *fldCC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];

  fldFC[0] += 0.5*fldCC[0]; 
  fldFLx[0] += fldCC[0]; 
}

void MGpoissonFEMprolong1xSer_LxNeumann_P1(double *fldCC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];

  fldFC[0] += 0.5*fldCC[0]; 
  fldFLx[0] += fldCC[0]; 
}

void MGpoissonFEMprolong1xSer_LxRobin_P1(double *fldCC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];

  fldFC[0] += 0.5*fldCC[0]; 
  fldFLx[0] += fldCC[0]; 
}

void MGpoissonFEMprolong1xSer_UxDirichlet_P1(double *fldCC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];
  double *fldFLxx = fldF[2];

  fldFC[0] += 0.5*fldCC[1]+0.5*fldCC[0]; 
  fldFC[1] += fldCC[1]; 
  fldFLx[0] += fldCC[0]; 
  fldFLxx[0] += 0.5*fldCC[0]; 
}

void MGpoissonFEMprolong1xSer_UxNeumann_P1(double *fldCC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];
  double *fldFLxx = fldF[2];

  fldFC[0] += 0.5*fldCC[1]+0.5*fldCC[0]; 
  fldFC[1] += fldCC[1]; 
  fldFLx[0] += fldCC[0]; 
  fldFLxx[0] += 0.5*fldCC[0]; 
}

void MGpoissonFEMprolong1xSer_UxRobin_P1(double *fldCC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];
  double *fldFLxx = fldF[2];

  fldFC[0] += 0.5*fldCC[1]+0.5*fldCC[0]; 
  fldFC[1] += fldCC[1]; 
  fldFLx[0] += fldCC[0]; 
  fldFLxx[0] += 0.5*fldCC[0]; 
}

void MGpoissonFEMrestrict1xSer_P1(double **fldF, double *fldCC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];
  double *fldFLxx = fldF[2];

  fldCC[0] += 0.5*fldFLxx[0]+fldFLx[0]+0.5*fldFC[0]; 
}

void MGpoissonFEMrestrict1xSer_LxDirichlet_P1(double **fldF, double *fldCC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldFLx = fldF[1];

  fldCC[0] += fldFLx[0]; 
}

void MGpoissonFEMrestrict1xSer_LxNeumann_P1(double **fldF, double *fldCC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];

  fldCC[0] += fldFLx[0]+0.5*fldFC[0]; 
}

void MGpoissonFEMrestrict1xSer_LxRobin_P1(double **fldF, double *fldCC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];

  fldCC[0] += fldFLx[0]+0.5*fldFC[0]; 
}

void MGpoissonFEMrestrict1xSer_UxDirichlet_P1(double **fldF, double *fldCC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];
  double *fldFLxx = fldF[2];

  fldCC[0] += 0.5*fldFLxx[0]+fldFLx[0]+0.5*fldFC[0]; 
  fldCC[1] += fldFC[1]; 
}

void MGpoissonFEMrestrict1xSer_UxNeumann_P1(double **fldF, double *fldCC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];
  double *fldFLxx = fldF[2];

  fldCC[0] += 0.5*fldFLxx[0]+fldFLx[0]+0.5*fldFC[0]; 
  fldCC[1] += fldFC[1]+0.5*fldFC[0]; 
}

void MGpoissonFEMrestrict1xSer_UxRobin_P1(double **fldF, double *fldCC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldFC = fldF[0];
  double *fldFLx = fldF[1];
  double *fldFLxx = fldF[2];

  fldCC[0] += 0.5*fldFLxx[0]+fldFLx[0]+0.5*fldFC[0]; 
  fldCC[1] += fldFC[1]+0.5*fldFC[0]; 
}

void MGpoissonFEM_DGtoFEM_1xSer_P1(double **dgFld, double *femFld) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femFld: FEM (nodal) field coefficients.

  double *dgFldC = dgFld[0]; 
  double *dgFldLx = dgFld[1]; 

  femFld[0] = 0.6123724356957944*dgFldLx[1]-0.6123724356957944*dgFldC[1]+0.3535533905932737*dgFldLx[0]+0.3535533905932737*dgFldC[0]; 

}

void MGpoissonFEM_DGtoFEM_1xSer_LxNonPeriodic_P1(double **dgFld, double *femFld) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femFld: FEM (nodal) field coefficients.

  double *dgFldC = dgFld[0]; 

  femFld[0] = 0.7071067811865475*dgFldC[0]-1.224744871391589*dgFldC[1]; 

}

void MGpoissonFEM_DGtoFEM_1xSer_UxNonPeriodic_P1(double **dgFld, double *femFld) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femFld: FEM (nodal) field coefficients.

  double *dgFldC = dgFld[0]; 
  double *dgFldLx = dgFld[1]; 

  femFld[0] = 0.6123724356957944*dgFldLx[1]-0.6123724356957944*dgFldC[1]+0.3535533905932737*dgFldLx[0]+0.3535533905932737*dgFldC[0]; 
  femFld[1] = 1.224744871391589*dgFldC[1]+0.7071067811865475*dgFldC[0]; 

}

void MGpoissonFEM_FEMtoDG_1xSer_P1(double **femFld, double *dgFld) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femFld: FEM (nodal) field coefficients.

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 

  dgFld[0] = 0.7071067811865475*femFldUx[0]+0.7071067811865475*femFldC[0]; 
  dgFld[1] = 0.408248290463863*femFldUx[0]-0.408248290463863*femFldC[0]; 

}

void MGpoissonFEM_FEMtoDG_1xSer_LxNonPeriodic_P1(double **femFld, double *dgFld) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femFld: FEM (nodal) field coefficients.

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 

  dgFld[0] = 0.7071067811865475*femFldUx[0]+0.7071067811865475*femFldC[0]; 
  dgFld[1] = 0.408248290463863*femFldUx[0]-0.408248290463863*femFldC[0]; 

}

void MGpoissonFEM_FEMtoDG_1xSer_UxNonPeriodic_P1(double **femFld, double *dgFld) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femFld: FEM (nodal) field coefficients.

  double *femFldC = femFld[0]; 

  dgFld[0] = 0.7071067811865475*femFldC[1]+0.7071067811865475*femFldC[0]; 
  dgFld[1] = 0.408248290463863*femFldC[1]-0.408248290463863*femFldC[0]; 

}

void MGpoissonFEMproject1xSer_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 

  out[0] = 0.3333333333333333*(femFldUx[0]+femFldLx[0]+4.0*femFldC[0])*volFac; 

}

void MGpoissonFEMproject1xSer_LxNonPeriodic_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 

  out[0] = 0.3333333333333333*(femFldUx[0]+2.0*femFldC[0])*volFac; 

}

void MGpoissonFEMproject1xSer_UxNonPeriodic_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 

  out[0] = 0.3333333333333333*(femFldC[1]+femFldLx[0]+4.0*femFldC[0])*volFac; 
  out[1] = 0.3333333333333333*(2.0*femFldC[1]+femFldC[0])*volFac; 

}

void MGpoissonFEMDampedJacobi1xSer_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
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
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *phiUx = phiPrev[2]; 


  phiC[0] = (0.5*((2.0*rhoC[0]+(phiUx[0]-2.0*phiPrevC[0]+phiLx[0])*rdx2SqVol[0])*omega+2.0*phiPrevC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonFEMDampedJacobi1xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
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
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *phiUx = phiPrev[2]; 


  phiC[0] = (bcVals[2]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 

}

void MGpoissonFEMDampedJacobi1xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
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
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *phiUx = phiPrev[2]; 


  phiC[0] = -(1.0*((2.0*rdx2SqVol[0]*bcVals[2]-2.0*rhoC[0]+(phiPrevC[0]-1.0*phiUx[0])*rdx2SqVol[0])*omega-1.0*phiPrevC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonFEMDampedJacobi1xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
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
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *phiUx = phiPrev[2]; 


  phiC[0] = -(1.0*((2.0*rdx2SqVol[0]*bcVals[2]+((phiPrevC[0]-1.0*phiUx[0])*rdx2SqVol[0]-2.0*rhoC[0])*bcVals[1]-2.0*bcVals[0]*phiPrevC[0]*rdx2SqVol[0])*omega-1.0*phiPrevC[0]*rdx2SqVol[0]*bcVals[1]+2.0*bcVals[0]*phiPrevC[0]*rdx2SqVol[0]))/(rdx2SqVol[0]*bcVals[1]-2.0*bcVals[0]*rdx2SqVol[0]); 

}

void MGpoissonFEMDampedJacobi1xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
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
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *phiUx = phiPrev[2]; 


  phiC[0] = (0.5*((rdx2SqVol[0]*bcVals[5]+2.0*rhoC[0]+(phiLx[0]-2.0*phiPrevC[0])*rdx2SqVol[0])*omega+2.0*phiPrevC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (bcVals[5]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi1xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
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
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *phiUx = phiPrev[2]; 


  phiC[0] = ((2.0*rdx2SqVol[0]*bcVals[5]+2.0*rhoC[1]+2.0*rhoC[0]+(phiLx[0]-1.0*phiPrevC[0])*rdx2SqVol[0])*omega+phiPrevC[0]*rdx2SqVol[0])/rdx2SqVol[0]; 
  phiC[1] = ((4.0*rdx2SqVol[0]*bcVals[5]+4.0*rhoC[1]-1.0*rdx2SqVol[0]*phiPrevC[1]+2.0*rhoC[0]+phiLx[0]*rdx2SqVol[0])*omega+rdx2SqVol[0]*phiPrevC[1])/rdx2SqVol[0]; 

}

void MGpoissonFEMDampedJacobi1xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
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
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *phiUx = phiPrev[2]; 


  phiC[0] = ((2.0*rdx2SqVol[0]*bcVals[5]+(2.0*rhoC[1]+2.0*rhoC[0]+(phiLx[0]-1.0*phiPrevC[0])*rdx2SqVol[0])*bcVals[4]+(4.0*rhoC[0]+(2.0*phiLx[0]-4.0*phiPrevC[0])*rdx2SqVol[0])*bcVals[3])*omega+phiPrevC[0]*rdx2SqVol[0]*bcVals[4]+4.0*phiPrevC[0]*rdx2SqVol[0]*bcVals[3])/(rdx2SqVol[0]*bcVals[4]+4.0*rdx2SqVol[0]*bcVals[3]); 
  phiC[1] = ((4.0*rdx2SqVol[0]*bcVals[5]+(4.0*rhoC[1]-1.0*rdx2SqVol[0]*phiPrevC[1]+2.0*rhoC[0]+phiLx[0]*rdx2SqVol[0])*bcVals[4]-4.0*rdx2SqVol[0]*phiPrevC[1]*bcVals[3])*omega+rdx2SqVol[0]*phiPrevC[1]*bcVals[4]+4.0*rdx2SqVol[0]*phiPrevC[1]*bcVals[3])/(rdx2SqVol[0]*bcVals[4]+4.0*rdx2SqVol[0]*bcVals[3]); 

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
  double *phiUx = phi[2]; 


  phiC[0] = (0.5*((2.0*rhoC[0]+(phiUx[0]+phiLx[0]-2.0*phiC[0])*rdx2SqVol[0])*omega+2.0*phiC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonFEMDampedGaussSeidel1xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[2]; 


  phiC[0] = (bcVals[2]-1.0*phiC[0])*omega+phiC[0]; 

}

void MGpoissonFEMDampedGaussSeidel1xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[2]; 


  phiC[0] = -(1.0*((2.0*rdx2SqVol[0]*bcVals[2]-2.0*rhoC[0]+(phiC[0]-1.0*phiUx[0])*rdx2SqVol[0])*omega-1.0*phiC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 

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
  double *phiUx = phi[2]; 


  phiC[0] = -(1.0*((2.0*rdx2SqVol[0]*bcVals[2]+((phiC[0]-1.0*phiUx[0])*rdx2SqVol[0]-2.0*rhoC[0])*bcVals[1]-2.0*bcVals[0]*phiC[0]*rdx2SqVol[0])*omega-1.0*phiC[0]*rdx2SqVol[0]*bcVals[1]+2.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/(rdx2SqVol[0]*bcVals[1]-2.0*bcVals[0]*rdx2SqVol[0]); 

}

void MGpoissonFEMDampedGaussSeidel1xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[2]; 


  phiC[0] = (0.5*((rdx2SqVol[0]*bcVals[5]+2.0*rhoC[0]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0])*omega+2.0*phiC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (bcVals[5]-1.0*phiC[1])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel1xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[2]; 


  phiC[0] = ((2.0*rdx2SqVol[0]*bcVals[5]+2.0*rhoC[1]+2.0*rhoC[0]+(phiLx[0]-1.0*phiC[0])*rdx2SqVol[0])*omega+phiC[0]*rdx2SqVol[0])/rdx2SqVol[0]; 
  phiC[1] = ((4.0*rdx2SqVol[0]*bcVals[5]+4.0*rhoC[1]-1.0*rdx2SqVol[0]*phiC[1]+2.0*rhoC[0]+phiLx[0]*rdx2SqVol[0])*omega+rdx2SqVol[0]*phiC[1])/rdx2SqVol[0]; 

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
  double *phiUx = phi[2]; 


  phiC[0] = ((2.0*rdx2SqVol[0]*bcVals[5]+(2.0*rhoC[1]+2.0*rhoC[0]+(phiLx[0]-1.0*phiC[0])*rdx2SqVol[0])*bcVals[4]+(4.0*rhoC[0]+(2.0*phiLx[0]-4.0*phiC[0])*rdx2SqVol[0])*bcVals[3])*omega+phiC[0]*rdx2SqVol[0]*bcVals[4]+4.0*phiC[0]*rdx2SqVol[0]*bcVals[3])/(rdx2SqVol[0]*bcVals[4]+4.0*rdx2SqVol[0]*bcVals[3]); 
  phiC[1] = ((4.0*rdx2SqVol[0]*bcVals[5]+(4.0*rhoC[1]-1.0*rdx2SqVol[0]*phiC[1]+2.0*rhoC[0]+phiLx[0]*rdx2SqVol[0])*bcVals[4]-4.0*rdx2SqVol[0]*phiC[1]*bcVals[3])*omega+rdx2SqVol[0]*phiC[1]*bcVals[4]+4.0*rdx2SqVol[0]*phiC[1]*bcVals[3])/(rdx2SqVol[0]*bcVals[4]+4.0*rdx2SqVol[0]*bcVals[3]); 

}

void MGpoissonFEMresidual1xSer_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residual in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  resOut[0] = 0.5*(2.0*rhoC[0]+(phiUx[0]+phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 

}

void MGpoissonFEMresidual1xSer_LxDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residual in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  resOut[0] = phiC[0]-1.0*bcVals[2]; 

}

void MGpoissonFEMresidual1xSer_LxNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residual in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  resOut[0] = -0.5*(2.0*rdx2SqVol[0]*bcVals[2]-2.0*rhoC[0]+(phiC[0]-1.0*phiUx[0])*rdx2SqVol[0]); 

}

void MGpoissonFEMresidual1xSer_LxRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residual in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  resOut[0] = -(0.5*(2.0*rdx2SqVol[0]*bcVals[2]+((phiC[0]-1.0*phiUx[0])*rdx2SqVol[0]-2.0*rhoC[0])*bcVals[1]-2.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/bcVals[1]; 

}

void MGpoissonFEMresidual1xSer_UxDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residual in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  resOut[0] = 0.5*(rdx2SqVol[0]*phiC[1]+2.0*rhoC[0]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[5]; 

}

void MGpoissonFEMresidual1xSer_UxNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residual in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  resOut[0] = 0.5*(rdx2SqVol[0]*phiC[1]+2.0*rhoC[0]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.5*(2.0*rdx2SqVol[0]*bcVals[5]+2.0*rhoC[1]-1.0*rdx2SqVol[0]*phiC[1]+phiC[0]*rdx2SqVol[0]); 

}

void MGpoissonFEMresidual1xSer_UxRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residual in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  resOut[0] = 0.5*(rdx2SqVol[0]*phiC[1]+2.0*rhoC[0]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = (0.5*(2.0*rdx2SqVol[0]*bcVals[5]+(2.0*rhoC[1]-1.0*rdx2SqVol[0]*phiC[1]+phiC[0]*rdx2SqVol[0])*bcVals[4]-2.0*rdx2SqVol[0]*phiC[1]*bcVals[3]))/bcVals[4]; 

}

void MGpoissonFEML2norm1xSer_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.5*dxC[0]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldUx0R2 = std::pow(femFldUx[0],2);

  normOut[0] += 0.3333333333333333*(2.0*femFldUx0R2+2.0*femFldC[0]*femFldUx[0]+2.0*femFldC0R2)*volFac; 
}

void MGpoissonFEML2norm1xSer_UxNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.5*dxC[0]; 

  double *femFldC = femFld[0]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldC1R2 = std::pow(femFldC[1],2);

  normOut[0] += 0.3333333333333333*(2.0*femFldC1R2+2.0*femFldC[0]*femFldC[1]+2.0*femFldC0R2)*volFac; 
}

void MGpoissonFEMM0norm1xSer_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.5*dxC[0]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 


  normOut[0] += (femFldUx[0]+femFldC[0])*volFac; 
}

void MGpoissonFEMM0norm1xSer_UxNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.5*dxC[0]; 

  double *femFldC = femFld[0]; 


  normOut[0] += (femFldC[1]+femFldC[0])*volFac; 
}

void MGpoissonFEMaccuConst1xSer_P1(const double constIn, double *femFld) 
{ 
  // constIn: constant to accumulate.
  // femFld:  FEM field to accumulate.

  femFld[0] += constIn; 
}

void MGpoissonFEMaccuConst1xSer_UxNonPeriodic_P1(const double constIn, double *femFld) 
{ 
  // constIn: constant to accumulate.
  // femFld:  FEM field to accumulate.

  femFld[0] += constIn; 
  femFld[1] += constIn; 
}

