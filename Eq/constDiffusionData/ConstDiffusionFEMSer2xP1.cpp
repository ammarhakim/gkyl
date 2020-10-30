void ConstDiffusionFEM2xSer_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((4.0*phiUy[0]+phiUxUy[0]+phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+((-2.0*phiUy[0])+phiUxUy[0]+phiUxLy[0]+4.0*phiUx[0]-2.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 

}

void ConstDiffusionFEM2xSer_LxDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[2]; 

}

void ConstDiffusionFEM2xSer_LxNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]-1.0*phiUxLy[0]+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*rdx2SqVol[1]+(phiUy[0]-1.0*phiUxUy[0]-1.0*phiUxLy[0]-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 

}

void ConstDiffusionFEM2xSer_LxRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]-1.0*phiUxLy[0]+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*bcVals[1]*rdx2SqVol[1]+(phiUy[0]-1.0*phiUxUy[0]-1.0*phiUxLy[0]-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiC[0])*rdx2SqVol[0]))/bcVals[1]; 

}

void ConstDiffusionFEM2xSer_UxDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*bcVals[5]+(phiLy[1]-2.0*phiC[1]+4.0*phiUy[0]+4.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiUy[0])-2.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = phiC[1]-1.0*bcVals[5]; 

}

void ConstDiffusionFEM2xSer_UxNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((phiUy[1]+phiLy[1]-2.0*phiC[1]+4.0*phiUy[0]+4.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+rdx2SqVol[0]*phiLy[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiUy[0])-2.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+(2.0*phiUy[1]+2.0*phiLy[1]-4.0*phiC[1]+phiUy[0]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiUy[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 

}

void ConstDiffusionFEM2xSer_UxRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((phiUy[1]+phiLy[1]-2.0*phiC[1]+4.0*phiUy[0]+4.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+rdx2SqVol[0]*phiLy[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiUy[0])-2.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+((2.0*phiUy[1]+2.0*phiLy[1]-4.0*phiC[1]+phiUy[0]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiUy[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[4]+((-2.0*rdx2SqVol[0]*phiUy[1])-4.0*rdx2SqVol[0]*phiC[1])*bcVals[3]))/bcVals[4]; 

}

void ConstDiffusionFEM2xSer_LyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[8]; 

}

void ConstDiffusionFEM2xSer_LyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+((-4.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+(2.0*phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0]); 

}

void ConstDiffusionFEM2xSer_LyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+(((-4.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+(2.0*phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 

}

void ConstDiffusionFEM2xSer_UyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*bcVals[11]+(phiLx[1]+4.0*phiC[1]+phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLx[1]-2.0*rdx2SqVol[0]*phiC[1]+(phiUxLy[0]+4.0*phiUx[0]-2.0*phiLy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = phiC[1]-1.0*bcVals[11]; 

}

void ConstDiffusionFEM2xSer_UyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((phiUx[1]+phiLx[1]+4.0*phiC[1]+phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUx[1]+rdx2SqVol[0]*phiLx[1]-2.0*rdx2SqVol[0]*phiC[1]+(phiUxLy[0]+4.0*phiUx[0]-2.0*phiLy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+((-1.0*phiUx[1])-1.0*phiLx[1]-4.0*phiC[1]+phiUx[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiUx[1]+2.0*rdx2SqVol[0]*phiLx[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiUx[0]+phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 

}

void ConstDiffusionFEM2xSer_UyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((phiUx[1]+phiLx[1]+4.0*phiC[1]+phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUx[1]+rdx2SqVol[0]*phiLx[1]-2.0*rdx2SqVol[0]*phiC[1]+(phiUxLy[0]+4.0*phiUx[0]-2.0*phiLy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+(((-1.0*phiUx[1])-1.0*phiLx[1]-4.0*phiC[1]+phiUx[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiUx[1]+2.0*rdx2SqVol[0]*phiLx[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiUx[0]+phiLx[0]-2.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((-2.0*phiUx[1])-4.0*phiC[1])*rdx2SqVol[1]*bcVals[9]))/bcVals[10]; 

}

void ConstDiffusionFEM2xSer_LxDirichletLyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -(1.0*(dxC[1]*bcVals[8]+dxC[0]*bcVals[2]-1.0*phiC[0]*dxC[1]-1.0*dxC[0]*phiC[0]))/(dxC[1]+dxC[0]); 

}

void ConstDiffusionFEM2xSer_LxDirichletLyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[2]; 

}

void ConstDiffusionFEM2xSer_LxDirichletLyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[2]; 

}

void ConstDiffusionFEM2xSer_LxNeumannLyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[8]; 

}

void ConstDiffusionFEM2xSer_LxNeumannLyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*rdx2SqVol[1]+(phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0]); 

}

void ConstDiffusionFEM2xSer_LxNeumannLyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*rdx2SqVol[1]+(phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 

}

void ConstDiffusionFEM2xSer_LxRobinLyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[8]; 

}

void ConstDiffusionFEM2xSer_LxRobinLyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -(0.1666666666666667*(6.0*bcVals[1]*rdx2SqVol[1]*bcVals[8]+6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*bcVals[1]*rdx2SqVol[1]+(phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0]*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiC[0])*rdx2SqVol[0]))/bcVals[1]; 

}

void ConstDiffusionFEM2xSer_LxRobinLyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -(0.1666666666666667*(6.0*bcVals[1]*rdx2SqVol[1]*bcVals[8]+(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*bcVals[1]*rdx2SqVol[1]+(phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0]*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiC[0])*bcVals[1]*rdx2SqVol[1]*bcVals[6]))/(bcVals[1]*bcVals[7]); 

}

void ConstDiffusionFEM2xSer_LxDirichletUyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[2]; 
  diffOut[1] = -(1.0*(dxC[1]*bcVals[11]+dxC[0]*bcVals[2]+((-1.0*dxC[1])-1.0*dxC[0])*phiC[1]))/(dxC[1]+dxC[0]); 

}

void ConstDiffusionFEM2xSer_LxDirichletUyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[2]; 
  diffOut[1] = phiC[1]-1.0*bcVals[2]; 

}

void ConstDiffusionFEM2xSer_LxDirichletUyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[2]; 
  diffOut[1] = phiC[1]-1.0*bcVals[2]; 

}

void ConstDiffusionFEM2xSer_LxNeumannUyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+(2.0*phiC[1]+phiUxLy[0]-2.0*phiUx[0]+2.0*phiLy[0]-4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiC[1]+(phiUxLy[0]+4.0*phiUx[0]-1.0*phiLy[0]-4.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = phiC[1]-1.0*bcVals[11]; 

}

void ConstDiffusionFEM2xSer_LxNeumannUyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-1.0*phiUx[1])-2.0*phiC[1]-1.0*phiUxLy[0]+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUx[1]+rdx2SqVol[0]*phiC[1]+((-1.0*phiUxLy[0])-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+((-1.0*phiUx[1])-2.0*phiC[1]+phiUx[0]+2.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiUx[1]-2.0*rdx2SqVol[0]*phiC[1]+(phiUx[0]-1.0*phiC[0])*rdx2SqVol[0]); 

}

void ConstDiffusionFEM2xSer_LxNeumannUyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-1.0*phiUx[1])-2.0*phiC[1]-1.0*phiUxLy[0]+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUx[1]+rdx2SqVol[0]*phiC[1]+((-1.0*phiUxLy[0])-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+((-6.0*rdx2SqVol[0]*bcVals[2])+((-1.0*phiUx[1])-2.0*phiC[1]+phiUx[0]+2.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiUx[1]-2.0*rdx2SqVol[0]*phiC[1]+(phiUx[0]-1.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((-2.0*phiUx[1])-4.0*phiC[1])*rdx2SqVol[1]*bcVals[9]))/bcVals[10]; 

}

void ConstDiffusionFEM2xSer_LxRobinUyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = (0.1666666666666667*((bcVals[1]*rdx2SqVol[1]+rdx2SqVol[0]*bcVals[1])*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+(2.0*bcVals[1]*phiC[1]+(phiUxLy[0]-2.0*phiUx[0]+2.0*phiLy[0]-4.0*phiC[0])*bcVals[1])*rdx2SqVol[1]+(2.0*bcVals[0]*rdx2SqVol[0]-1.0*rdx2SqVol[0]*bcVals[1])*phiC[1]+(phiUxLy[0]+4.0*phiUx[0]-1.0*phiLy[0]-4.0*phiC[0])*rdx2SqVol[0]*bcVals[1]+4.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/bcVals[1]; 
  diffOut[1] = phiC[1]-1.0*bcVals[11]; 

}

void ConstDiffusionFEM2xSer_LxRobinUyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-1.0*bcVals[1]*phiUx[1])-2.0*bcVals[1]*phiC[1]+((-1.0*phiUxLy[0])+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*bcVals[1])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*bcVals[1]*phiUx[1]+(rdx2SqVol[0]*bcVals[1]-2.0*bcVals[0]*rdx2SqVol[0])*phiC[1]+((-1.0*phiUxLy[0])-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/bcVals[1]; 
  diffOut[1] = (0.1666666666666667*(6.0*bcVals[1]*rdx2SqVol[1]*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+((-1.0*bcVals[1]*phiUx[1])-2.0*bcVals[1]*phiC[1]+(phiUx[0]+2.0*phiC[0])*bcVals[1])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*bcVals[1]*phiUx[1]+(4.0*bcVals[0]*rdx2SqVol[0]-2.0*rdx2SqVol[0]*bcVals[1])*phiC[1]+(phiUx[0]-1.0*phiC[0])*rdx2SqVol[0]*bcVals[1]+2.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/bcVals[1]; 

}

void ConstDiffusionFEM2xSer_LxRobinUyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-1.0*bcVals[1]*phiUx[1])-2.0*bcVals[1]*phiC[1]+((-1.0*phiUxLy[0])+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*bcVals[1])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*bcVals[1]*phiUx[1]+(rdx2SqVol[0]*bcVals[1]-2.0*bcVals[0]*rdx2SqVol[0])*phiC[1]+((-1.0*phiUxLy[0])-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/bcVals[1]; 
  diffOut[1] = (0.1666666666666667*(6.0*bcVals[1]*rdx2SqVol[1]*bcVals[11]+((-6.0*rdx2SqVol[0]*bcVals[2])+((-1.0*bcVals[1]*phiUx[1])-2.0*bcVals[1]*phiC[1]+(phiUx[0]+2.0*phiC[0])*bcVals[1])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*bcVals[1]*phiUx[1]+(4.0*bcVals[0]*rdx2SqVol[0]-2.0*rdx2SqVol[0]*bcVals[1])*phiC[1]+(phiUx[0]-1.0*phiC[0])*rdx2SqVol[0]*bcVals[1]+2.0*bcVals[0]*phiC[0]*rdx2SqVol[0])*bcVals[10]+((-2.0*bcVals[1]*phiUx[1])-4.0*bcVals[1]*phiC[1])*rdx2SqVol[1]*bcVals[9]))/(bcVals[1]*bcVals[10]); 

}

void ConstDiffusionFEM2xSer_UxDirichletLyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[8]; 
  diffOut[1] = -(1.0*(dxC[1]*bcVals[8]+dxC[0]*bcVals[5]+((-1.0*dxC[1])-1.0*dxC[0])*phiC[1]))/(dxC[1]+dxC[0]); 

}

void ConstDiffusionFEM2xSer_UxDirichletLyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+((-1.0*rdx2SqVol[1])-1.0*rdx2SqVol[0])*bcVals[5]+(phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-2.0*rdx2SqVol[0]*phiC[1]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = phiC[1]-1.0*bcVals[5]; 

}

void ConstDiffusionFEM2xSer_UxDirichletLyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+(((-1.0*rdx2SqVol[1])-1.0*rdx2SqVol[0])*bcVals[5]+(phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-2.0*rdx2SqVol[0]*phiC[1]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiC[1])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 
  diffOut[1] = phiC[1]-1.0*bcVals[5]; 

}

void ConstDiffusionFEM2xSer_UxNeumannLyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[8]; 
  diffOut[1] = phiC[1]-1.0*bcVals[8]; 

}

void ConstDiffusionFEM2xSer_UxNeumannLyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+((-1.0*phiUy[1])+phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-2.0*rdx2SqVol[0]*phiC[1]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]-6.0*rdx2SqVol[0]*bcVals[5]+((-2.0*phiUy[1])+2.0*phiC[1]-1.0*phiUy[0]+phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+2.0*rdx2SqVol[0]*phiC[1]+((-1.0*phiUy[0])-2.0*phiC[0])*rdx2SqVol[0]); 

}

void ConstDiffusionFEM2xSer_UxNeumannLyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+(((-1.0*phiUy[1])+phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-2.0*rdx2SqVol[0]*phiC[1]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiC[1])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 
  diffOut[1] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+((-6.0*rdx2SqVol[0]*bcVals[5])+((-2.0*phiUy[1])+2.0*phiC[1]-1.0*phiUy[0]+phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+2.0*rdx2SqVol[0]*phiC[1]+((-1.0*phiUy[0])-2.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-4.0*phiC[1])-2.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 

}

void ConstDiffusionFEM2xSer_UxRobinLyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = phiC[0]-1.0*bcVals[8]; 
  diffOut[1] = phiC[1]-1.0*bcVals[8]; 

}

void ConstDiffusionFEM2xSer_UxRobinLyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+((-1.0*phiUy[1])+phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-2.0*rdx2SqVol[0]*phiC[1]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[4]*bcVals[8]-6.0*rdx2SqVol[0]*bcVals[5]+(((-2.0*phiUy[1])+2.0*phiC[1]-1.0*phiUy[0]+phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+2.0*rdx2SqVol[0]*phiC[1]+((-1.0*phiUy[0])-2.0*phiC[0])*rdx2SqVol[0])*bcVals[4]+(2.0*rdx2SqVol[0]*phiUy[1]+4.0*rdx2SqVol[0]*phiC[1])*bcVals[3]))/bcVals[4]; 

}

void ConstDiffusionFEM2xSer_UxRobinLyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+(((-1.0*phiUy[1])+phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-2.0*rdx2SqVol[0]*phiC[1]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiC[1])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 
  diffOut[1] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[4]*bcVals[8]+((-6.0*rdx2SqVol[0]*bcVals[5])+(((-2.0*phiUy[1])+2.0*phiC[1]-1.0*phiUy[0]+phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+2.0*rdx2SqVol[0]*phiC[1]+((-1.0*phiUy[0])-2.0*phiC[0])*rdx2SqVol[0])*bcVals[4]+(2.0*rdx2SqVol[0]*phiUy[1]+4.0*rdx2SqVol[0]*phiC[1])*bcVals[3])*bcVals[7]+((-4.0*phiC[1])-2.0*phiC[0])*rdx2SqVol[1]*bcVals[4]*bcVals[6]))/(bcVals[4]*bcVals[7]); 

}

void ConstDiffusionFEM2xSer_UxDirichletUyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = phiC[1]-1.0*bcVals[5]; 
  diffOut[2] = phiC[2]-1.0*bcVals[11]; 
  diffOut[3] = -(1.0*(dxC[1]*bcVals[11]+dxC[0]*bcVals[5]+((-1.0*dxC[1])-1.0*dxC[0])*phiC[3]))/(dxC[1]+dxC[0]); 

}

void ConstDiffusionFEM2xSer_UxDirichletUyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = phiC[1]-1.0*bcVals[5]; 
  diffOut[2] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[3] = phiC[3]-1.0*bcVals[5]; 

}

void ConstDiffusionFEM2xSer_UxDirichletUyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = phiC[1]-1.0*bcVals[5]; 
  diffOut[2] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+((2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((-2.0*rdx2SqVol[1]*phiC[3])-4.0*rdx2SqVol[1]*phiC[2])*bcVals[9]))/bcVals[10]; 
  diffOut[3] = phiC[3]-1.0*bcVals[5]; 

}

void ConstDiffusionFEM2xSer_UxNeumannUyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+(2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[2] = phiC[2]-1.0*bcVals[11]; 
  diffOut[3] = phiC[3]-1.0*bcVals[11]; 

}

void ConstDiffusionFEM2xSer_UxNeumannUyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+(2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[2] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[3] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+6.0*rdx2SqVol[0]*bcVals[5]+((-2.0*rdx2SqVol[1])-2.0*rdx2SqVol[0])*phiC[3]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[2]+(2.0*phiC[1]+phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiC[1]+phiC[0]*rdx2SqVol[0]); 

}

void ConstDiffusionFEM2xSer_UxNeumannUyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+(2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[2] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+((2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((-2.0*rdx2SqVol[1]*phiC[3])-4.0*rdx2SqVol[1]*phiC[2])*bcVals[9]))/bcVals[10]; 
  diffOut[3] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+(6.0*rdx2SqVol[0]*bcVals[5]+((-2.0*rdx2SqVol[1])-2.0*rdx2SqVol[0])*phiC[3]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[2]+(2.0*phiC[1]+phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiC[1]+phiC[0]*rdx2SqVol[0])*bcVals[10]+((-4.0*rdx2SqVol[1]*phiC[3])-2.0*rdx2SqVol[1]*phiC[2])*bcVals[9]))/bcVals[10]; 

}

void ConstDiffusionFEM2xSer_UxRobinUyDirichlet_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+((2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[4]-2.0*rdx2SqVol[0]*bcVals[3]*phiC[3]-4.0*rdx2SqVol[0]*phiC[1]*bcVals[3]))/bcVals[4]; 
  diffOut[2] = phiC[2]-1.0*bcVals[11]; 
  diffOut[3] = phiC[3]-1.0*bcVals[11]; 

}

void ConstDiffusionFEM2xSer_UxRobinUyNeumann_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+((2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[4]-2.0*rdx2SqVol[0]*bcVals[3]*phiC[3]-4.0*rdx2SqVol[0]*phiC[1]*bcVals[3]))/bcVals[4]; 
  diffOut[2] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[3] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[4]*bcVals[11]+6.0*rdx2SqVol[0]*bcVals[5]+(((-2.0*rdx2SqVol[1])-2.0*rdx2SqVol[0])*phiC[3]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[2]+(2.0*phiC[1]+phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiC[1]+phiC[0]*rdx2SqVol[0])*bcVals[4]-4.0*rdx2SqVol[0]*bcVals[3]*phiC[3]-2.0*rdx2SqVol[0]*phiC[1]*bcVals[3]))/bcVals[4]; 

}

void ConstDiffusionFEM2xSer_UxRobinUyRobin_P1(double **dx, const double *bcVals, double **phi, double *diffOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // phi:     iterate cells pointed to by the stencil.
  // diffOut: diffusion term.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  diffOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  diffOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+((2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[4]-2.0*rdx2SqVol[0]*bcVals[3]*phiC[3]-4.0*rdx2SqVol[0]*phiC[1]*bcVals[3]))/bcVals[4]; 
  diffOut[2] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+((2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((-2.0*rdx2SqVol[1]*phiC[3])-4.0*rdx2SqVol[1]*phiC[2])*bcVals[9]))/bcVals[10]; 
  diffOut[3] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[4]*bcVals[11]+(6.0*rdx2SqVol[0]*bcVals[5]+(((-2.0*rdx2SqVol[1])-2.0*rdx2SqVol[0])*phiC[3]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[2]+(2.0*phiC[1]+phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiC[1]+phiC[0]*rdx2SqVol[0])*bcVals[4]-4.0*rdx2SqVol[0]*bcVals[3]*phiC[3]-2.0*rdx2SqVol[0]*phiC[1]*bcVals[3])*bcVals[10]+((-4.0*rdx2SqVol[1]*phiC[3])-2.0*rdx2SqVol[1]*phiC[2])*bcVals[4]*bcVals[9]))/(bcVals[4]*bcVals[10]); 

}

