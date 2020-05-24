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

void MGpoissonFEMProlong1xSer_LxDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];

  fldF0[0] += 0.5*fldC[0]; 
  fldF1[0] += fldC[0]; 
}

void MGpoissonFEMProlong1xSer_LxNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];

  fldF0[0] += 0.5*fldC[0]; 
  fldF1[0] += fldC[0]; 
}

void MGpoissonFEMProlong1xSer_LxRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];

  fldF0[0] += 0.5*fldC[0]; 
  fldF1[0] += fldC[0]; 
}

void MGpoissonFEMProlong1xSer_UxDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];

  fldF0[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF1[0] += fldC[0]; 
  fldF0[1] += fldC[1]; 
  fldF2[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong1xSer_UxNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];

  fldF0[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF1[0] += fldC[0]; 
  fldF0[1] += fldC[1]; 
  fldF2[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong1xSer_UxRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];

  fldF0[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF1[0] += fldC[0]; 
  fldF0[1] += fldC[1]; 
  fldF2[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMRestrict1xSer_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];

  fldC[0] += 0.5*fldF2[0]+fldF1[0]+0.5*fldF0[0]; 
}

void MGpoissonFEMRestrict1xSer_LxDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF1 = fldF[1];

  fldC[0] += fldF1[0]; 
}

void MGpoissonFEMRestrict1xSer_LxNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];

  fldC[0] += fldF1[0]+0.5*fldF0[0]; 
}

void MGpoissonFEMRestrict1xSer_LxRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];

  fldC[0] += fldF1[0]+0.5*fldF0[0]; 
}

void MGpoissonFEMRestrict1xSer_UxDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];

  fldC[0] += 0.5*fldF2[0]+fldF1[0]+0.5*fldF0[0]; 
  fldC[1] += fldF0[1]; 
}

void MGpoissonFEMRestrict1xSer_UxNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];

  fldC[0] += 0.5*fldF2[0]+fldF1[0]+0.5*fldF0[0]; 
  fldC[1] += fldF0[1]+0.5*fldF0[0]; 
}

void MGpoissonFEMRestrict1xSer_UxRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];

  fldC[0] += 0.5*fldF2[0]+fldF1[0]+0.5*fldF0[0]; 
  fldC[1] += fldF0[1]+0.5*fldF0[0]; 
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

void MGpoissonFEM_DGtoFEM_1xSer_LxNonPeriodic_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 

  femFld[0] += 0.7071067811865475*dgFld[0]-1.224744871391589*dgFld[1]; 
  femFldUx[0] += 0.6123724356957944*dgFld[1]+0.3535533905932737*dgFld[0]; 

}

void MGpoissonFEM_DGtoFEM_1xSer_UxNonPeriodic_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 

  femFld[0] += 0.3535533905932737*dgFld[0]-0.6123724356957944*dgFld[1]; 
  femFld[1] += 1.224744871391589*dgFld[1]+0.7071067811865475*dgFld[0]; 

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

void MGpoissonFEMResidue1xSer_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

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

void MGpoissonFEMResidue1xSer_LxDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

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

void MGpoissonFEMResidue1xSer_LxNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

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

void MGpoissonFEMResidue1xSer_LxRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

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

void MGpoissonFEMResidue1xSer_UxDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

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

void MGpoissonFEMResidue1xSer_UxNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

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

void MGpoissonFEMResidue1xSer_UxRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

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

void MGpoissonFEML2norm1xSer_P1(double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldUx0R2 = std::pow(femFldUx[0],2);

  normOut[0] += 0.3333333333333333*(2.0*femFldUx0R2+2.0*femFldC[0]*femFldUx[0]+2.0*femFldC0R2); 
}

void MGpoissonFEML2norm1xSer_UxNonPeriodic_P1(double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double *femFldC = femFld[0]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldC1R2 = std::pow(femFldC[1],2);

  normOut[0] += 0.3333333333333333*(2.0*femFldC1R2+2.0*femFldC[0]*femFldC[1]+2.0*femFldC0R2); 
}

