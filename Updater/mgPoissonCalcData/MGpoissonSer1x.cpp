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

void MGpoissonGaussSeidel1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (0.05555555555555555*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0])*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.02173913043478261*(16.0*rho[1]*volFac-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonGaussSeidel1xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (0.001603416737384461*((151.8683640525571*rho[1]+370.523953341751*rho[0])*volFac-220.454076850486*rdx2SqVol[0]*phiUx[1]+(250.3158005400378*phiUx[0]+528.0*bcVals[0])*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.001603416737384461*((110.3086578651014*rho[1]+83.28265125462806*rho[0])*volFac-140.0071426749364*rdx2SqVol[0]*phiUx[1]+(146.9693845669907*phiUx[0]-207.8460969082653*bcVals[0])*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonGaussSeidel1xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(0.003319091792389131*((127.2792206135786*rho[1]-690.7561074648562*rho[0])*volFac+296.98484809835*rdx2SqVol[0]*phiUx[1]+(609.6818842642448*bcVals[0]-301.2872383623309*phiUx[0])*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.009957275377167391*((44.0908153700972*rho[1]-53.74011537017763*rho[0])*volFac+2.449489742783178*rdx2SqVol[0]*phiUx[1]+80.0*bcVals[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonGaussSeidel1xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(0.001603416737384461*((151.8683640525571*rho[1]-370.523953341751*rho[0])*volFac-220.454076850486*rdx2SqVol[0]*phiLx[1]-528.0*rdx2SqVol[0]*bcVals[1]-250.3158005400378*phiLx[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.001603416737384461*((110.3086578651014*rho[1]-83.28265125462806*rho[0])*volFac-140.0071426749364*rdx2SqVol[0]*phiLx[1]+207.8460969082653*rdx2SqVol[0]*bcVals[1]-146.9693845669907*phiLx[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonGaussSeidel1xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (0.003319091792389131*((127.2792206135786*rho[1]+690.7561074648562*rho[0])*volFac+296.98484809835*rdx2SqVol[0]*phiLx[1]+609.6818842642448*rdx2SqVol[0]*bcVals[1]+301.2872383623309*phiLx[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.009957275377167391*((44.0908153700972*rho[1]+53.74011537017763*rho[0])*volFac+2.449489742783178*rdx2SqVol[0]*phiLx[1]+80.0*rdx2SqVol[0]*bcVals[1]))/rdx2SqVol[0]; 

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
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (0.05555555555555555*(16.0*rho[0]*omega*volFac+((-8.660254037844386*rdx2SqVol[0]*phiUx[1])+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdx2SqVol[0])*omega+18.0*phiC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.02173913043478261*(16.0*rho[1]*omega*volFac+((-7.0*rdx2SqVol[0]*phiUx[1])-7.0*rdx2SqVol[0]*phiLx[1]-46.0*rdx2SqVol[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0])*omega+46.0*rdx2SqVol[0]*phiC[1]))/rdx2SqVol[0]; 

}

void MGpoissonDampedGaussSeidel1xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (0.001603416737384461*((151.8683640525571*rho[1]+370.523953341751*rho[0])*omega*volFac+((250.3158005400378*phiUx[0]-623.668181006535*phiC[0]+528.0*bcVals[0])*rdx2SqVol[0]-220.454076850486*rdx2SqVol[0]*phiUx[1])*omega+623.668181006535*phiC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.001603416737384461*((110.3086578651014*rho[1]+83.28265125462806*rho[0])*omega*volFac+((-140.0071426749364*rdx2SqVol[0]*phiUx[1])-623.668181006535*rdx2SqVol[0]*phiC[1]+(146.9693845669907*phiUx[0]-207.8460969082653*bcVals[0])*rdx2SqVol[0])*omega+623.668181006535*rdx2SqVol[0]*phiC[1]))/rdx2SqVol[0]; 

}

void MGpoissonDampedGaussSeidel1xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(0.003319091792389131*((127.2792206135786*rho[1]-690.7561074648562*rho[0])*omega*volFac+(296.98484809835*rdx2SqVol[0]*phiUx[1]+((-301.2872383623309*phiUx[0])+301.2872383623309*phiC[0]+609.6818842642448*bcVals[0])*rdx2SqVol[0])*omega-301.2872383623309*phiC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.009957275377167391*((44.0908153700972*rho[1]-53.74011537017763*rho[0])*omega*volFac+(2.449489742783178*rdx2SqVol[0]*phiUx[1]-100.4290794541103*rdx2SqVol[0]*phiC[1]+80.0*bcVals[0]*rdx2SqVol[0])*omega+100.4290794541103*rdx2SqVol[0]*phiC[1]))/rdx2SqVol[0]; 

}

void MGpoissonDampedGaussSeidel1xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = -(0.001603416737384461*((151.8683640525571*rho[1]-370.523953341751*rho[0])*omega*volFac+((-220.454076850486*rdx2SqVol[0]*phiLx[1])-528.0*rdx2SqVol[0]*bcVals[1]+(623.668181006535*phiC[0]-250.3158005400378*phiLx[0])*rdx2SqVol[0])*omega-623.668181006535*phiC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.001603416737384461*((110.3086578651014*rho[1]-83.28265125462806*rho[0])*omega*volFac+((-140.0071426749364*rdx2SqVol[0]*phiLx[1])-623.668181006535*rdx2SqVol[0]*phiC[1]+207.8460969082653*rdx2SqVol[0]*bcVals[1]-146.9693845669907*phiLx[0]*rdx2SqVol[0])*omega+623.668181006535*rdx2SqVol[0]*phiC[1]))/rdx2SqVol[0]; 

}

void MGpoissonDampedGaussSeidel1xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  phiC[0] = (0.003319091792389131*((127.2792206135786*rho[1]+690.7561074648562*rho[0])*omega*volFac+(296.98484809835*rdx2SqVol[0]*phiLx[1]+609.6818842642448*rdx2SqVol[0]*bcVals[1]+(301.2872383623309*phiLx[0]-301.2872383623309*phiC[0])*rdx2SqVol[0])*omega+301.2872383623309*phiC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.009957275377167391*((44.0908153700972*rho[1]+53.74011537017763*rho[0])*omega*volFac+(2.449489742783178*rdx2SqVol[0]*phiLx[1]-100.4290794541103*rdx2SqVol[0]*phiC[1]+80.0*rdx2SqVol[0]*bcVals[1])*omega+100.4290794541103*rdx2SqVol[0]*phiC[1]))/rdx2SqVol[0]; 

}

void MGpoissonJacobi1xSer_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 

  phiC[0] = (0.05555555555555555*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0])*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.02173913043478261*(16.0*rho[1]*volFac-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonJacobi1xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 

  phiC[0] = (0.001603416737384461*((151.8683640525571*rho[1]+370.523953341751*rho[0])*volFac-220.454076850486*rdx2SqVol[0]*phiUx[1]+(250.3158005400378*phiUx[0]+528.0*bcVals[0])*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.001603416737384461*((110.3086578651014*rho[1]+83.28265125462806*rho[0])*volFac-140.0071426749364*rdx2SqVol[0]*phiUx[1]+(146.9693845669907*phiUx[0]-207.8460969082653*bcVals[0])*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonJacobi1xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 

  phiC[0] = -(0.003319091792389131*((127.2792206135786*rho[1]-690.7561074648562*rho[0])*volFac+296.98484809835*rdx2SqVol[0]*phiUx[1]+(609.6818842642448*bcVals[0]-301.2872383623309*phiUx[0])*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.009957275377167391*((44.0908153700972*rho[1]-53.74011537017763*rho[0])*volFac+2.449489742783178*rdx2SqVol[0]*phiUx[1]+80.0*bcVals[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonJacobi1xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 

  phiC[0] = -(0.001603416737384461*((151.8683640525571*rho[1]-370.523953341751*rho[0])*volFac-220.454076850486*rdx2SqVol[0]*phiLx[1]-528.0*rdx2SqVol[0]*bcVals[1]-250.3158005400378*phiLx[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.001603416737384461*((110.3086578651014*rho[1]-83.28265125462806*rho[0])*volFac-140.0071426749364*rdx2SqVol[0]*phiLx[1]+207.8460969082653*rdx2SqVol[0]*bcVals[1]-146.9693845669907*phiLx[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 

}

void MGpoissonJacobi1xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 

  phiC[0] = (0.003319091792389131*((127.2792206135786*rho[1]+690.7561074648562*rho[0])*volFac+296.98484809835*rdx2SqVol[0]*phiLx[1]+609.6818842642448*rdx2SqVol[0]*bcVals[1]+301.2872383623309*phiLx[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.009957275377167391*((44.0908153700972*rho[1]+53.74011537017763*rho[0])*volFac+2.449489742783178*rdx2SqVol[0]*phiLx[1]+80.0*rdx2SqVol[0]*bcVals[1]))/rdx2SqVol[0]; 

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
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 

  phiC[0] = (0.05555555555555555*(16.0*rho[0]*omega*volFac+((-8.660254037844386*rdx2SqVol[0]*phiUx[1])+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]-18.0*phiPrevC[0]+9.0*phiLx[0])*rdx2SqVol[0])*omega+18.0*phiPrevC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.02173913043478261*(16.0*rho[1]*omega*volFac+((-7.0*rdx2SqVol[0]*phiUx[1])-46.0*rdx2SqVol[0]*phiPrevC[1]-7.0*rdx2SqVol[0]*phiLx[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0])*omega+46.0*rdx2SqVol[0]*phiPrevC[1]))/rdx2SqVol[0]; 

}

void MGpoissonDampedJacobi1xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 

  phiC[0] = (0.001603416737384461*((151.8683640525571*rho[1]+370.523953341751*rho[0])*omega*volFac+((250.3158005400378*phiUx[0]-623.668181006535*phiPrevC[0]+528.0*bcVals[0])*rdx2SqVol[0]-220.454076850486*rdx2SqVol[0]*phiUx[1])*omega+623.668181006535*phiPrevC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.001603416737384461*((110.3086578651014*rho[1]+83.28265125462806*rho[0])*omega*volFac+((-140.0071426749364*rdx2SqVol[0]*phiUx[1])-623.668181006535*rdx2SqVol[0]*phiPrevC[1]+(146.9693845669907*phiUx[0]-207.8460969082653*bcVals[0])*rdx2SqVol[0])*omega+623.668181006535*rdx2SqVol[0]*phiPrevC[1]))/rdx2SqVol[0]; 

}

void MGpoissonDampedJacobi1xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 

  phiC[0] = -(0.003319091792389131*((127.2792206135786*rho[1]-690.7561074648562*rho[0])*omega*volFac+(296.98484809835*rdx2SqVol[0]*phiUx[1]+((-301.2872383623309*phiUx[0])+301.2872383623309*phiPrevC[0]+609.6818842642448*bcVals[0])*rdx2SqVol[0])*omega-301.2872383623309*phiPrevC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.009957275377167391*((44.0908153700972*rho[1]-53.74011537017763*rho[0])*omega*volFac+(2.449489742783178*rdx2SqVol[0]*phiUx[1]-100.4290794541103*rdx2SqVol[0]*phiPrevC[1]+80.0*bcVals[0]*rdx2SqVol[0])*omega+100.4290794541103*rdx2SqVol[0]*phiPrevC[1]))/rdx2SqVol[0]; 

}

void MGpoissonDampedJacobi1xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 

  phiC[0] = -(0.001603416737384461*((151.8683640525571*rho[1]-370.523953341751*rho[0])*omega*volFac+((-220.454076850486*rdx2SqVol[0]*phiLx[1])-528.0*rdx2SqVol[0]*bcVals[1]+(623.668181006535*phiPrevC[0]-250.3158005400378*phiLx[0])*rdx2SqVol[0])*omega-623.668181006535*phiPrevC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.001603416737384461*((110.3086578651014*rho[1]-83.28265125462806*rho[0])*omega*volFac+((-623.668181006535*rdx2SqVol[0]*phiPrevC[1])-140.0071426749364*rdx2SqVol[0]*phiLx[1]+207.8460969082653*rdx2SqVol[0]*bcVals[1]-146.9693845669907*phiLx[0]*rdx2SqVol[0])*omega+623.668181006535*rdx2SqVol[0]*phiPrevC[1]))/rdx2SqVol[0]; 

}

void MGpoissonDampedJacobi1xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, const double *rho, double **phiPrev, double **phi) 
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
  double *phiUx = phiPrev[1]; 
  double *phiLx = phiPrev[2]; 

  phiC[0] = (0.003319091792389131*((127.2792206135786*rho[1]+690.7561074648562*rho[0])*omega*volFac+(296.98484809835*rdx2SqVol[0]*phiLx[1]+609.6818842642448*rdx2SqVol[0]*bcVals[1]+(301.2872383623309*phiLx[0]-301.2872383623309*phiPrevC[0])*rdx2SqVol[0])*omega+301.2872383623309*phiPrevC[0]*rdx2SqVol[0]))/rdx2SqVol[0]; 
  phiC[1] = (0.009957275377167391*((44.0908153700972*rho[1]+53.74011537017763*rho[0])*omega*volFac+((-100.4290794541103*rdx2SqVol[0]*phiPrevC[1])+2.449489742783178*rdx2SqVol[0]*phiLx[1]+80.0*rdx2SqVol[0]*bcVals[1])*omega+100.4290794541103*rdx2SqVol[0]*phiPrevC[1]))/rdx2SqVol[0]; 

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
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.0625*(16.0*rho[0]*volFac-8.660254037844386*rdx2SqVol[0]*phiUx[1]+8.660254037844386*rdx2SqVol[0]*phiLx[1]+(9.0*phiUx[0]+9.0*phiLx[0]-18.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-7.0*rdx2SqVol[0]*phiUx[1]-7.0*rdx2SqVol[0]*phiLx[1]-46.0*rdx2SqVol[0]*phiC[1]+(8.660254037844386*phiUx[0]-8.660254037844386*phiLx[0])*rdx2SqVol[0]); 

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

  double rdx2SqVol[1]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.0441941738241592*(22.62741699796953*rho[0]*volFac-2.449489742783178*rdx2SqVol[0]*phiUx[1]+75.93418202627852*rdx2SqVol[0]*phiC[1]+(4.242640687119286*phiUx[0]-55.15432893255071*phiC[0]+72.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-19.0*rdx2SqVol[0]*phiUx[1]-131.0*rdx2SqVol[0]*phiC[1]+(19.05255888325765*phiUx[0]+29.44486372867091*phiC[0]-68.585712797929*bcVals[0])*rdx2SqVol[0]); 

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

  double rdx2SqVol[1]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.0441941738241592*(22.62741699796953*rho[0]*volFac-12.24744871391589*rdx2SqVol[0]*phiUx[1]-12.24744871391589*rdx2SqVol[0]*phiC[1]+(12.72792206135786*phiUx[0]-12.72792206135786*phiC[0]-16.0*bcVals[0])*rdx2SqVol[0]); 
  resOut[1] = 0.01202813060811721*(83.1384387633061*rho[1]*volFac-50.22947341949744*rdx2SqVol[0]*phiUx[1]-244.2191638672117*rdx2SqVol[0]*phiC[1]+(57.0*phiUx[0]-57.0*phiC[0]+79.19595949289337*bcVals[0])*rdx2SqVol[0]); 

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

  double rdx2SqVol[1]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.0441941738241592*(22.62741699796953*rho[0]*volFac+2.449489742783178*rdx2SqVol[0]*phiLx[1]-75.93418202627852*rdx2SqVol[0]*phiC[1]+72.0*rdx2SqVol[0]*bcVals[1]+(4.242640687119286*phiLx[0]-55.15432893255071*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.0625*(16.0*rho[1]*volFac-19.0*rdx2SqVol[0]*phiLx[1]-131.0*rdx2SqVol[0]*phiC[1]+68.585712797929*rdx2SqVol[0]*bcVals[1]+((-19.05255888325765*phiLx[0])-29.44486372867091*phiC[0])*rdx2SqVol[0]); 

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

  double rdx2SqVol[1]; 
  rdx2SqVol[0]   = volFac*4.0/(dxC[0]*dxC[0]); 

  double *phiC = phi[0]; 
  double *phiUx = phi[1]; 
  double *phiLx = phi[2]; 

  resOut[0] = 0.0441941738241592*(22.62741699796953*rho[0]*volFac+12.24744871391589*rdx2SqVol[0]*phiLx[1]+12.24744871391589*rdx2SqVol[0]*phiC[1]+16.0*rdx2SqVol[0]*bcVals[1]+(12.72792206135786*phiLx[0]-12.72792206135786*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.01202813060811721*(83.1384387633061*rho[1]*volFac-50.22947341949744*rdx2SqVol[0]*phiLx[1]-244.2191638672117*rdx2SqVol[0]*phiC[1]+79.19595949289337*rdx2SqVol[0]*bcVals[1]+(57.0*phiC[0]-57.0*phiLx[0])*rdx2SqVol[0]); 

}

