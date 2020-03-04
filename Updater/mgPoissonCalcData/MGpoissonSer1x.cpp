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

void MGpoissonDampedGaussSeidel1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 

  phiC[0] = (0.05555555555555555*(16.0*rho[0]*omega*volFac+((-8.660254037844386*rdx2SqVol[0]*phiUx[1])+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdx2SqVol[0])*omega+18.0*phiC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.02173913043478261*(16.0*rho[1]*omega*volFac+((-7.0*rdx2SqVol[0]*phiUx[1])-7.0*rdx2SqVol[0]*phiLx[1]-46.0*rdx2SqVol[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0])*omega+46.0*rdx2SqVol[0]*phiC[1]))/rdx2SqVol[0]; 

}

void MGpoissonDampedGaussSeidel1xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double bcValsSq[6]; 
  bcValsSq[0] = bcVals[0]*bcVals[0];
  bcValsSq[1] = bcVals[1]*bcVals[1];
  bcValsSq[2] = bcVals[2]*bcVals[2];
  bcValsSq[3] = bcVals[3]*bcVals[3];
  bcValsSq[4] = bcVals[4]*bcVals[4];
  bcValsSq[5] = bcVals[5]*bcVals[5];

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 

  phiC[0] = -(1.0*(((220.454076850486*bcVals[1]+303.7367281051142*bcVals[0])*rho[1]-1196.424673767639*rho[0]*bcVals[1]+741.0479066835021*bcVals[0]*rho[0])*omega*volFac+(1056.0*rdx2SqVol[0]*bcVals[2]+(514.3928459844674*rdx2SqVol[0]*bcVals[1]-440.9081537009721*bcVals[0]*rdx2SqVol[0])*phiUx[1]+(521.8448045156721*phiC[0]-521.8448045156721*phiUx[0])*rdx2SqVol[0]*bcVals[1]+(500.6316010800758*bcVals[0]*phiUx[0]-1247.33636201307*bcVals[0]*phiC[0])*rdx2SqVol[0])*omega-521.8448045156721*phiC[0]*rdx2SqVol[0]*bcVals[1]+1247.33636201307*bcVals[0]*phiC[0]*rdx2SqVol[0]))/(521.8448045156721*rdx2SqVol[0]*bcVals[1]-1247.33636201307*bcVals[0]*rdx2SqVol[0]); 
  phiC[1] = (((229.1025971044415*bcVals[1]-220.6173157302029*bcVals[0])*rho[1]-279.2418306772823*rho[0]*bcVals[1]-166.5653025092562*bcVals[0]*rho[0])*omega*volFac+(415.6921938165305*rdx2SqVol[0]*bcVals[2]+(12.72792206135786*rdx2SqVol[0]*bcVals[1]+280.0142853498729*bcVals[0]*rdx2SqVol[0])*phiUx[1]+(1247.33636201307*bcVals[0]*rdx2SqVol[0]-521.8448045156721*rdx2SqVol[0]*bcVals[1])*phiC[1]-293.9387691339815*bcVals[0]*phiUx[0]*rdx2SqVol[0])*omega+(521.8448045156721*rdx2SqVol[0]*bcVals[1]-1247.33636201307*bcVals[0]*rdx2SqVol[0])*phiC[1])/(521.8448045156721*rdx2SqVol[0]*bcVals[1]-1247.33636201307*bcVals[0]*rdx2SqVol[0]); 

}

void MGpoissonDampedGaussSeidel1xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double bcValsSq[6]; 
  bcValsSq[0] = bcVals[0]*bcVals[0];
  bcValsSq[1] = bcVals[1]*bcVals[1];
  bcValsSq[2] = bcVals[2]*bcVals[2];
  bcValsSq[3] = bcVals[3]*bcVals[3];
  bcValsSq[4] = bcVals[4]*bcVals[4];
  bcValsSq[5] = bcVals[5]*bcVals[5];

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 

  phiC[0] = (((220.454076850486*rho[1]+1196.424673767639*rho[0])*bcVals[4]+(741.0479066835021*rho[0]-303.7367281051142*rho[1])*bcVals[3])*omega*volFac+(1056.0*rdx2SqVol[0]*bcVals[5]+(514.3928459844674*rdx2SqVol[0]*phiLx[1]+(521.8448045156721*phiLx[0]-521.8448045156721*phiC[0])*rdx2SqVol[0])*bcVals[4]+(440.9081537009721*rdx2SqVol[0]*phiLx[1]+(500.6316010800758*phiLx[0]-1247.33636201307*phiC[0])*rdx2SqVol[0])*bcVals[3])*omega+521.8448045156721*phiC[0]*rdx2SqVol[0]*bcVals[4]+1247.33636201307*phiC[0]*rdx2SqVol[0]*bcVals[3])/(521.8448045156721*rdx2SqVol[0]*bcVals[4]+1247.33636201307*rdx2SqVol[0]*bcVals[3]); 
  phiC[1] = (((229.1025971044415*rho[1]+279.2418306772823*rho[0])*bcVals[4]+(220.6173157302029*rho[1]-166.5653025092562*rho[0])*bcVals[3])*omega*volFac+(415.6921938165305*rdx2SqVol[0]*bcVals[5]+(12.72792206135786*rdx2SqVol[0]*phiLx[1]-521.8448045156721*rdx2SqVol[0]*phiC[1])*bcVals[4]+((-280.0142853498729*rdx2SqVol[0]*phiLx[1])-1247.33636201307*rdx2SqVol[0]*phiC[1]-293.9387691339815*phiLx[0]*rdx2SqVol[0])*bcVals[3])*omega+521.8448045156721*rdx2SqVol[0]*phiC[1]*bcVals[4]+1247.33636201307*rdx2SqVol[0]*phiC[1]*bcVals[3])/(521.8448045156721*rdx2SqVol[0]*bcVals[4]+1247.33636201307*rdx2SqVol[0]*bcVals[3]); 

}

void MGpoissonDampedJacobi1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiLx = phiPrev[1]; 
  double *phiUx = phiPrev[2]; 

  phiC[0] = (0.05555555555555555*(16.0*rho[0]*omega*volFac+((-8.660254037844386*rdx2SqVol[0]*phiUx[1])+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]-18.0*phiPrevC[0]+9.0*phiLx[0])*rdx2SqVol[0])*omega+18.0*phiPrevC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.02173913043478261*(16.0*rho[1]*omega*volFac+((-7.0*rdx2SqVol[0]*phiUx[1])-46.0*rdx2SqVol[0]*phiPrevC[1]-7.0*rdx2SqVol[0]*phiLx[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0])*omega+46.0*rdx2SqVol[0]*phiPrevC[1]))/rdx2SqVol[0]; 

}

void MGpoissonDampedJacobi1xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double bcValsSq[6]; 
  bcValsSq[0] = bcVals[0]*bcVals[0];
  bcValsSq[1] = bcVals[1]*bcVals[1];
  bcValsSq[2] = bcVals[2]*bcVals[2];
  bcValsSq[3] = bcVals[3]*bcVals[3];
  bcValsSq[4] = bcVals[4]*bcVals[4];
  bcValsSq[5] = bcVals[5]*bcVals[5];

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiLx = phiPrev[1]; 
  double *phiUx = phiPrev[2]; 

  phiC[0] = -(1.0*(((220.454076850486*bcVals[1]+303.7367281051142*bcVals[0])*rho[1]-1196.424673767639*rho[0]*bcVals[1]+741.0479066835021*bcVals[0]*rho[0])*omega*volFac+(1056.0*rdx2SqVol[0]*bcVals[2]+(514.3928459844674*rdx2SqVol[0]*bcVals[1]-440.9081537009721*bcVals[0]*rdx2SqVol[0])*phiUx[1]+(521.8448045156721*phiPrevC[0]-521.8448045156721*phiUx[0])*rdx2SqVol[0]*bcVals[1]+(500.6316010800758*bcVals[0]*phiUx[0]-1247.33636201307*bcVals[0]*phiPrevC[0])*rdx2SqVol[0])*omega-521.8448045156721*phiPrevC[0]*rdx2SqVol[0]*bcVals[1]+1247.33636201307*bcVals[0]*phiPrevC[0]*rdx2SqVol[0]))/(521.8448045156721*rdx2SqVol[0]*bcVals[1]-1247.33636201307*bcVals[0]*rdx2SqVol[0]); 
  phiC[1] = (((229.1025971044415*bcVals[1]-220.6173157302029*bcVals[0])*rho[1]-279.2418306772823*rho[0]*bcVals[1]-166.5653025092562*bcVals[0]*rho[0])*omega*volFac+(415.6921938165305*rdx2SqVol[0]*bcVals[2]+(12.72792206135786*rdx2SqVol[0]*bcVals[1]+280.0142853498729*bcVals[0]*rdx2SqVol[0])*phiUx[1]+(1247.33636201307*bcVals[0]*rdx2SqVol[0]-521.8448045156721*rdx2SqVol[0]*bcVals[1])*phiPrevC[1]-293.9387691339815*bcVals[0]*phiUx[0]*rdx2SqVol[0])*omega+(521.8448045156721*rdx2SqVol[0]*bcVals[1]-1247.33636201307*bcVals[0]*rdx2SqVol[0])*phiPrevC[1])/(521.8448045156721*rdx2SqVol[0]*bcVals[1]-1247.33636201307*bcVals[0]*rdx2SqVol[0]); 

}

void MGpoissonDampedJacobi1xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double bcValsSq[6]; 
  bcValsSq[0] = bcVals[0]*bcVals[0];
  bcValsSq[1] = bcVals[1]*bcVals[1];
  bcValsSq[2] = bcVals[2]*bcVals[2];
  bcValsSq[3] = bcVals[3]*bcVals[3];
  bcValsSq[4] = bcVals[4]*bcVals[4];
  bcValsSq[5] = bcVals[5]*bcVals[5];

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *phiLx = phiPrev[1]; 
  double *phiUx = phiPrev[2]; 

  phiC[0] = (((220.454076850486*rho[1]+1196.424673767639*rho[0])*bcVals[4]+(741.0479066835021*rho[0]-303.7367281051142*rho[1])*bcVals[3])*omega*volFac+(1056.0*rdx2SqVol[0]*bcVals[5]+(514.3928459844674*rdx2SqVol[0]*phiLx[1]+(521.8448045156721*phiLx[0]-521.8448045156721*phiPrevC[0])*rdx2SqVol[0])*bcVals[4]+(440.9081537009721*rdx2SqVol[0]*phiLx[1]+(500.6316010800758*phiLx[0]-1247.33636201307*phiPrevC[0])*rdx2SqVol[0])*bcVals[3])*omega+521.8448045156721*phiPrevC[0]*rdx2SqVol[0]*bcVals[4]+1247.33636201307*phiPrevC[0]*rdx2SqVol[0]*bcVals[3])/(521.8448045156721*rdx2SqVol[0]*bcVals[4]+1247.33636201307*rdx2SqVol[0]*bcVals[3]); 
  phiC[1] = (((229.1025971044415*rho[1]+279.2418306772823*rho[0])*bcVals[4]+(220.6173157302029*rho[1]-166.5653025092562*rho[0])*bcVals[3])*omega*volFac+(415.6921938165305*rdx2SqVol[0]*bcVals[5]+(12.72792206135786*rdx2SqVol[0]*phiLx[1]-521.8448045156721*rdx2SqVol[0]*phiPrevC[1])*bcVals[4]+((-1247.33636201307*rdx2SqVol[0]*phiPrevC[1])-280.0142853498729*rdx2SqVol[0]*phiLx[1]-293.9387691339815*phiLx[0]*rdx2SqVol[0])*bcVals[3])*omega+521.8448045156721*rdx2SqVol[0]*phiPrevC[1]*bcVals[4]+1247.33636201307*rdx2SqVol[0]*phiPrevC[1]*bcVals[3])/(521.8448045156721*rdx2SqVol[0]*bcVals[4]+1247.33636201307*rdx2SqVol[0]*bcVals[3]); 

}

void MGpoissonResidue1xSer_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]-46.0*rdx2SqVol[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0]); 

}

void MGpoissonResidue1xSer_LxRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 

  resOut[0] = ((203.6467529817258*rho[0]*bcVals[1]-45.25483399593907*bcVals[0]*rho[0])*volFac-144.0*rdx2SqVol[0]*bcVals[2]+(4.898979485566357*bcVals[0]*rdx2SqVol[0]-110.227038425243*rdx2SqVol[0]*bcVals[1])*phiUx[1]+((-110.227038425243*rdx2SqVol[0]*bcVals[1])-151.8683640525571*bcVals[0]*rdx2SqVol[0])*phiC[1]+(114.5512985522207*phiUx[0]-114.5512985522207*phiC[0])*rdx2SqVol[0]*bcVals[1]+(110.3086578651014*bcVals[0]*phiC[0]-8.485281374238571*bcVals[0]*phiUx[0])*rdx2SqVol[0])/(203.6467529817258*bcVals[1]-45.25483399593907*bcVals[0]); 
  resOut[1] = ((144.0*bcVals[1]-32.0*bcVals[0])*rho[1]*volFac+137.171425595858*rdx2SqVol[0]*bcVals[2]+(38.0*bcVals[0]*rdx2SqVol[0]-87.0*rdx2SqVol[0]*bcVals[1])*phiUx[1]+(262.0*bcVals[0]*rdx2SqVol[0]-423.0*rdx2SqVol[0]*bcVals[1])*phiC[1]+(98.726896031426*phiUx[0]-98.726896031426*phiC[0])*rdx2SqVol[0]*bcVals[1]+((-38.1051177665153*bcVals[0]*phiUx[0])-58.88972745734183*bcVals[0]*phiC[0])*rdx2SqVol[0])/(144.0*bcVals[1]-32.0*bcVals[0]); 

}

void MGpoissonResidue1xSer_UxRobin_P1(double **dx, const double *bcVals, const double *rho, double **phi, double *resOut) 
{ 
  // dx:     cell lengths of cells pointed to by the stencil.
  // rho:    right-side source in the current cell.
  // bcVals: values to impose as BCs.
  // phi:    iterate cells pointed to by the stencil.
  // resOut: residue in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.5*dxC[0]; 

  double rdx2SqVol[1]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiLx = phi[1]; 
  double *phiUx = phi[2]; 

  resOut[0] = ((203.6467529817258*rho[0]*bcVals[4]+45.25483399593907*rho[0]*bcVals[3])*volFac+144.0*rdx2SqVol[0]*bcVals[5]+(110.227038425243*rdx2SqVol[0]*phiLx[1]+110.227038425243*rdx2SqVol[0]*phiC[1]+(114.5512985522207*phiLx[0]-114.5512985522207*phiC[0])*rdx2SqVol[0])*bcVals[4]+(4.898979485566357*rdx2SqVol[0]*phiLx[1]-151.8683640525571*rdx2SqVol[0]*phiC[1]+(8.485281374238571*phiLx[0]-110.3086578651014*phiC[0])*rdx2SqVol[0])*bcVals[3])/(203.6467529817258*bcVals[4]+45.25483399593907*bcVals[3]); 
  resOut[1] = ((144.0*rho[1]*bcVals[4]+32.0*rho[1]*bcVals[3])*volFac+137.171425595858*rdx2SqVol[0]*bcVals[5]+((-87.0*rdx2SqVol[0]*phiLx[1])-423.0*rdx2SqVol[0]*phiC[1]+(98.726896031426*phiC[0]-98.726896031426*phiLx[0])*rdx2SqVol[0])*bcVals[4]+((-38.0*rdx2SqVol[0]*phiLx[1])-262.0*rdx2SqVol[0]*phiC[1]+((-38.1051177665153*phiLx[0])-58.88972745734183*phiC[0])*rdx2SqVol[0])*bcVals[3])/(144.0*bcVals[4]+32.0*bcVals[3]); 

}

