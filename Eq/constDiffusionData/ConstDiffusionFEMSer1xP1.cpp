void ConstDiffusionFEM1xSer_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  diffOut[0] = 0.5*(phiUx[0]+phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]; 

}

void ConstDiffusionFEM1xSer_LxDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  diffOut[0] = phiC[0]-1.0*bcVals[2]; 

}

void ConstDiffusionFEM1xSer_LxNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  diffOut[0] = -0.5*(2.0*rdx2SqVol[0]*bcVals[2]+(phiC[0]-1.0*phiUx[0])*rdx2SqVol[0]); 

}

void ConstDiffusionFEM1xSer_LxRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  diffOut[0] = -(0.5*(2.0*rdx2SqVol[0]*bcVals[2]+(phiC[0]-1.0*phiUx[0])*rdx2SqVol[0]*bcVals[1]-2.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/bcVals[1]; 

}

void ConstDiffusionFEM1xSer_UxDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  diffOut[0] = 0.5*(rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = phiC[1]-1.0*bcVals[5]; 

}

void ConstDiffusionFEM1xSer_UxNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  diffOut[0] = 0.5*(rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = 0.5*(2.0*rdx2SqVol[0]*bcVals[5]-1.0*rdx2SqVol[0]*phiC[1]+phiC[0]*rdx2SqVol[0]); 

}

void ConstDiffusionFEM1xSer_UxRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 


  diffOut[0] = 0.5*(rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = (0.5*(2.0*rdx2SqVol[0]*bcVals[5]+(phiC[0]*rdx2SqVol[0]-1.0*rdx2SqVol[0]*phiC[1])*bcVals[4]-2.0*rdx2SqVol[0]*phiC[1]*bcVals[3]))/bcVals[4]; 

}

